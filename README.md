# Test Data Generator for SAM/BAM/CRAM and FASTQ Files

![Test Coverage](https://via.placeholder.com/800x200.png?text=42+SAM/BAM+tests+%7C+12+FASTQ+tests+%7C+3+reference+genomes)

Generates standardized test files for validating bioinformatics tools. Contains 54 pre-built scenarios covering common edge cases and format requirements.

## Key Features

```plaintext
📂 42 alignment tests        🔍 12 FASTQ variations
🔄 Reproducible outputs      🧬 3 reference genomes
📊 Metadata documentation   🐍 Python-powered
```

## Quick Start

**Prerequisites**  
- Python 3.8+  
- `pysam` library (`pip install pysam`)

**Generate All Tests**  
```bash
python -m test_data_generator.generate_tests
```

Files appear in:  
`test_data/  
├── SAM_01/  
├── SAM_02/  
...  
└── FASTQ_12/`

## Test Categories

| Category          | Tests | Example Cases                 |
|-------------------|-------|-------------------------------|
| SAM/BAM/CRAM      | 42    | Mapped pairs, QC failures     |
| FASTQ             | 12    | GZIP files, mixed read lengths|
| Special Cases     | 14    | Empty reads, circular genomes |

## Customization  
Modify `config.py` to:  
- Change output locations  
- Select specific test cases  
- Adjust reference genomes

## Why This Matters  
Properly formatted test data helps catch tool errors before real data analysis. This suite validates:  
- File format compliance  
- Edge case handling  
- Pipeline robustness

---

*"Never trust a bioinformatics tool you haven't tried to break first"* - Genomics proverb
