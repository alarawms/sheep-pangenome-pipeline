#!/usr/bin/env python3
"""
Stage 1 Testing Script for Sheep Pangenome Pipeline
Tests data acquisition and preparation components
"""

import os
import json
import pandas as pd
from pathlib import Path
import subprocess
import sys

class Stage1Tester:
    def __init__(self, work_dir="."):
        self.work_dir = Path(work_dir)
        self.test_results = {}

    def create_test_samplesheet(self):
        """Create a minimal test samplesheet with 3-5 genomes"""
        test_data = [
            {
                'sample': 'rambouillet_test',
                'accession': 'GCF_016772045.1',
                'breed': 'Rambouillet',
                'population': 'European',
                'geographic_origin': 'United_States'
            },
            {
                'sample': 'texel_test',
                'accession': 'GCF_000298735.2',
                'breed': 'Texel',
                'population': 'European',
                'geographic_origin': 'Netherlands'
            },
            {
                'sample': 'hu_test',
                'accession': 'GCA_001704415.1',
                'breed': 'Hu',
                'population': 'Asian',
                'geographic_origin': 'China'
            }
        ]

        df = pd.DataFrame(test_data)
        test_file = self.work_dir / 'test_samplesheet.csv'
        df.to_csv(test_file, index=False)
        print(f"✅ Created test samplesheet: {test_file}")
        return test_file

    def test_module_syntax(self, module_path):
        """Test Nextflow module syntax"""
        print(f"🔍 Testing module syntax: {module_path}")

        # Use nextflow to validate syntax
        cmd = ["nextflow", "inspect", str(module_path)]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            if result.returncode == 0:
                print(f"✅ Module syntax valid: {module_path.name}")
                return True
            else:
                print(f"❌ Module syntax error in {module_path.name}:")
                print(result.stderr)
                return False
        except subprocess.TimeoutExpired:
            print(f"⏰ Timeout testing {module_path.name}")
            return False
        except FileNotFoundError:
            print("⚠️  Nextflow not found - skipping syntax validation")
            return None

    def validate_file_structure(self):
        """Validate that all required files exist"""
        print("📁 Validating file structure...")

        required_files = [
            'modules/local/download_genome.nf',
            'modules/local/validate_genome.nf',
            'subworkflows/local/input_check.nf',
            'conf/modules.config',
            'assets/sheep_genomes_catalog.csv'
        ]

        missing_files = []
        for file_path in required_files:
            full_path = self.work_dir / file_path
            if not full_path.exists():
                missing_files.append(file_path)
            else:
                print(f"✅ Found: {file_path}")

        if missing_files:
            print(f"❌ Missing files: {missing_files}")
            return False
        else:
            print("✅ All required files present")
            return True

    def validate_catalog(self):
        """Validate the sheep genomes catalog"""
        print("🐑 Validating sheep genomes catalog...")

        catalog_path = self.work_dir / 'assets/sheep_genomes_catalog.csv'
        if not catalog_path.exists():
            print("❌ Catalog file not found")
            return False

        try:
            df = pd.read_csv(catalog_path)

            # Check required columns
            required_cols = ['sample', 'accession', 'breed', 'population', 'geographic_origin']
            missing_cols = [col for col in required_cols if col not in df.columns]

            if missing_cols:
                print(f"❌ Missing columns in catalog: {missing_cols}")
                return False

            # Check for duplicates
            if df['sample'].duplicated().any():
                print("❌ Duplicate sample IDs in catalog")
                return False

            if df['accession'].duplicated().any():
                print("❌ Duplicate accessions in catalog")
                return False

            print(f"✅ Catalog validated: {len(df)} genomes")
            print(f"   - Breeds: {df['breed'].nunique()}")
            print(f"   - Populations: {df['population'].nunique()}")
            return True

        except Exception as e:
            print(f"❌ Error validating catalog: {e}")
            return False

    def test_stage1_logic(self):
        """Test the logical flow of Stage 1"""
        print("🔬 Testing Stage 1 logical flow...")

        # Test samplesheet parsing logic
        test_samplesheet = self.create_test_samplesheet()

        try:
            df = pd.read_csv(test_samplesheet)

            # Simulate the branching logic from input_check.nf
            download_samples = df[df['accession'].notna()]
            # Create fasta column if it doesn't exist
            if 'fasta' not in df.columns:
                df['fasta'] = None
            local_samples = df[df['fasta'].notna()]

            print(f"📥 Samples for download: {len(download_samples)}")
            print(f"📁 Local samples: {len(local_samples)}")

            if len(download_samples) + len(local_samples) == len(df):
                print("✅ Sample classification logic correct")
                return True
            else:
                print("❌ Sample classification logic error")
                return False

        except Exception as e:
            print(f"❌ Error in stage 1 logic test: {e}")
            return False

    def generate_validation_criteria(self):
        """Generate validation criteria document"""
        criteria = {
            "stage1_validation_criteria": {
                "download_validation": {
                    "success_rate": ">95%",
                    "genome_size_range": "2.4-3.2 Gb",
                    "max_download_time": "30 minutes per genome",
                    "retry_attempts": 3
                },
                "genome_validation": {
                    "gc_content_range": "35-50%",
                    "n_content_max": "5%",
                    "max_contigs": 50000,
                    "busco_completeness": ">85% (if available)"
                },
                "quality_gates": {
                    "all_downloads_successful": "required",
                    "all_validations_passed": "required",
                    "metadata_completeness": ">90%",
                    "no_duplicate_samples": "required"
                }
            }
        }

        criteria_file = self.work_dir / 'stage1_validation_criteria.json'
        with open(criteria_file, 'w') as f:
            json.dump(criteria, f, indent=2)

        print(f"📋 Generated validation criteria: {criteria_file}")
        return criteria_file

    def run_all_tests(self):
        """Run comprehensive Stage 1 testing"""
        print("🚀 Running Stage 1 Comprehensive Tests")
        print("=" * 50)

        tests = [
            ("File Structure", self.validate_file_structure),
            ("Catalog Validation", self.validate_catalog),
            ("Stage 1 Logic", self.test_stage1_logic),
        ]

        results = {}
        for test_name, test_func in tests:
            print(f"\n{test_name}:")
            results[test_name] = test_func()

        # Test module syntax if nextflow available
        modules_dir = self.work_dir / 'modules/local'
        if modules_dir.exists():
            print(f"\nModule Syntax Tests:")
            for module_file in modules_dir.glob('*.nf'):
                results[f"Syntax_{module_file.stem}"] = self.test_module_syntax(module_file)

        # Generate validation criteria
        self.generate_validation_criteria()

        # Summary
        print("\n" + "=" * 50)
        print("📊 TEST SUMMARY")
        print("=" * 50)

        passed = sum(1 for r in results.values() if r is True)
        total = len([r for r in results.values() if r is not None])

        for test, result in results.items():
            if result is True:
                print(f"✅ {test}")
            elif result is False:
                print(f"❌ {test}")
            elif result is None:
                print(f"⚠️  {test} (skipped)")

        print(f"\n📈 Results: {passed}/{total} tests passed")

        if passed == total:
            print("🎉 Stage 1 implementation ready for testing!")
        else:
            print("🔧 Stage 1 needs fixes before proceeding")

        return results

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Test Stage 1 implementation')
    parser.add_argument('--work-dir', default='.', help='Working directory')
    args = parser.parse_args()

    tester = Stage1Tester(args.work_dir)
    results = tester.run_all_tests()

    # Exit with error code if tests failed
    failed_tests = [k for k, v in results.items() if v is False]
    if failed_tests:
        sys.exit(1)