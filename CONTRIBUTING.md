# Contributing to Sheep Pangenome Pipeline

Thank you for your interest in contributing to the Sheep Pangenome Pipeline! 🐑🧬

## Development Philosophy

This project follows a **staged implementation approach**:
- Each stage is fully implemented, tested, and validated before proceeding
- Scientific rigor and reproducibility are paramount
- Comprehensive documentation and testing are required for all stages

## Current Status

- ✅ **Stage 1**: Data Acquisition & Preparation (Complete)
- 🔄 **Stage 2**: Genome Preprocessing & Indexing (Next)
- ⏳ **Stages 3-7**: Planned future development

## How to Contribute

### 1. Issues and Bug Reports

- Check existing issues before creating new ones
- Use the provided issue templates
- Include detailed information:
  - Pipeline stage and version
  - Error messages and logs
  - System information (cluster, resources)
  - Reproducible example

### 2. Feature Requests

- Discuss new features in GitHub Discussions first
- Consider impact on staged implementation approach
- Ensure alignment with scientific best practices

### 3. Code Contributions

#### Prerequisites
- Nextflow ≥23.04.0
- Understanding of nf-core best practices
- Experience with SLURM and HPC systems
- Knowledge of genomics workflows

#### Development Workflow

1. **Fork the repository**
2. **Create a feature branch**: `git checkout -b feature/your-feature-name`
3. **Make changes following our standards**:
   - Follow existing code style and structure
   - Add comprehensive tests
   - Update documentation
   - Maintain KAUST Ibex compatibility
4. **Test thoroughly**:
   ```bash
   python3 scripts/stage1_test.py
   nextflow run . -profile test
   ```
5. **Commit with descriptive messages**
6. **Push and create a Pull Request**

#### Code Standards

- **Nextflow DSL2** syntax
- **nf-core** best practices and conventions
- **Comprehensive documentation** for all modules
- **Test coverage** for new functionality
- **SLURM optimization** for HPC environments
- **Scientific validation** with appropriate quality gates

### 4. Documentation

- Update relevant documentation for any changes
- Follow existing documentation style
- Include usage examples and troubleshooting guides
- Ensure KAUST Ibex-specific instructions are accurate

### 5. Testing

All contributions must include appropriate tests:

- **Unit tests**: For individual modules
- **Integration tests**: For workflow components
- **Validation tests**: For scientific correctness
- **Performance tests**: For resource optimization

Run the complete test suite:
```bash
# Stage 1 comprehensive tests
python3 scripts/stage1_test.py

# Nextflow syntax validation
nextflow inspect main.nf

# Test profile execution
nextflow run . --input assets/sheep_genomes_catalog.csv -profile test
```

## Stage Development Guidelines

### For Stage 2+ Development

When contributing to future stages:

1. **Follow staged approach**: Complete implementation before moving forward
2. **Maintain input-output contracts**: Ensure proper data flow between stages
3. **Add comprehensive validation**: Scientific and technical quality gates
4. **Document thoroughly**: Usage guides and troubleshooting
5. **Test extensively**: Multiple scales and configurations

### Stage-Specific Considerations

- **Stage 2**: Genome preprocessing and standardization
- **Stage 3**: PGGB pangenome construction (high-memory optimization)
- **Stage 4**: Graph analysis and visualization
- **Stage 5**: Variant calling and genotyping
- **Stage 6**: Population genomics analysis
- **Stage 7**: Publication-ready reporting

## Review Process

1. **Automated checks**: Tests, syntax validation, documentation
2. **Code review**: Maintainer review for quality and standards
3. **Scientific review**: Validation of biological correctness
4. **Integration testing**: Compatibility with existing stages
5. **Documentation review**: Clarity and completeness

## Communication

- **GitHub Issues**: Bug reports and feature requests
- **GitHub Discussions**: General questions and ideas
- **Pull Request reviews**: Code-specific discussions

## Recognition

Contributors will be:
- Added to the AUTHORS file
- Mentioned in release notes
- Acknowledged in publications using the pipeline

## Code of Conduct

Please be respectful and professional in all interactions. This project aims to advance sheep genomics research through collaborative, open science.

## Getting Help

- Check existing documentation first
- Search GitHub Issues for similar problems
- Use GitHub Discussions for general questions
- Contact maintainers for urgent issues

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

Thank you for contributing to sheep genomics research! 🐑🧬✨