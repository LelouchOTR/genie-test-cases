# Test Data Generator for SAM/BAM/CRAM and FASTQ Files

Generates standardized test files for validating bioinformatics tools. Contains 54 pre-built scenarios covering common edge cases and format requirements.

```text
SAM/BAM/CRAM (42)      FASTQ (12)        Special (14)
█████████████████████  ████████████      ████████████
Mapped Pairs (18)      Paired End (6)    Empty Reads (4)
Flags/Formats (24)     Compression (2)   Edge Cases (10)
```

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

```text
SAM/BAM/CRAM (42)       FASTQ (12)
├─ Mapped Pairs        ├─ Paired End
├─ Alignment Flags     ├─ Single End  
├─ CIGAR Operations    ├─ Compression
├─ Optional Tags       └─ Read Types
└─ Format Conversions

Special Cases (14)
├─ Empty Reads
├─ Circular References
├─ Edge Case Handling
└─ Boundary Conditions
```

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
