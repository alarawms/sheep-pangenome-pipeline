include { DOWNLOAD_GENOME } from '../../modules/local/download_genome'
include { VALIDATE_GENOME } from '../../modules/local/validate_genome'

workflow INPUT_CHECK {
    take:
    samplesheet // channel: samplesheet read in from --input

    main:
    // Parse the samplesheet
    parsed_samplesheet = samplesheet
        .splitCsv(header: true, sep: ',')
        .map { create_channels_from_samplesheet(it) }

    // Separate samples with accession numbers vs local files
    parsed_samplesheet
        .branch {
            download: it[2] != null && it[2] != '' // has accession
            local: it[1] != null && it[1] != ''   // has local file path
        }
        .set { samples_branched }

    // Download genomes from NCBI
    DOWNLOAD_GENOME(
        samples_branched.download.map { meta, fasta, accession ->
            [meta, accession]
        }
    )

    // Prepare local files
    local_genomes = samples_branched.local.map { meta, fasta, accession ->
        [meta, file(fasta)]
    }

    // Combine downloaded and local genomes
    all_genomes = DOWNLOAD_GENOME.out.fasta.mix(local_genomes)

    // Validate all genomes
    VALIDATE_GENOME(all_genomes)

    emit:
    genomes    = all_genomes
    validation = VALIDATE_GENOME.out.validation
    statistics = VALIDATE_GENOME.out.stats
    downloads  = DOWNLOAD_GENOME.out.log
    versions   = DOWNLOAD_GENOME.out.versions.mix(VALIDATE_GENOME.out.versions)
}

// Function to create channels from samplesheet
def create_channels_from_samplesheet(LinkedHashMap row) {
    // Check for required columns
    if (!row.containsKey('sample')) {
        error("Samplesheet missing required 'sample' column")
    }

    def meta = [:]
    meta.id = row.sample

    // Add optional metadata
    if (row.containsKey('breed')) meta.breed = row.breed
    if (row.containsKey('population')) meta.population = row.population
    if (row.containsKey('geographic_origin')) meta.geographic_origin = row.geographic_origin
    if (row.containsKey('quality_score')) meta.quality_score = row.quality_score

    // Determine input type
    def fasta_file = row.containsKey('fasta') ? row.fasta : null
    def accession = row.containsKey('accession') ? row.accession : null

    // Validate that we have either fasta or accession
    if (!fasta_file && !accession) {
        error("Sample ${meta.id} must have either 'fasta' path or 'accession' number")
    }
    if (fasta_file && accession) {
        error("Sample ${meta.id} cannot have both 'fasta' and 'accession' - choose one")
    }

    return [meta, fasta_file, accession]
}