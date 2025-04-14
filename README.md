# Test Data Generator for SAM/BAM/CRAM and FASTQ Files

Generates standardized test files for validating bioinformatics tools. Contains 54 pre-built scenarios covering common edge cases and format requirements.

```text
SAM/BAM/CRAM (42)                  FASTQ (12)       
█████████████████████              ████████████
│─ Core Alignment Tests (28)       ├─ Core FASTQ (10)
│─ Special Cases* (14)             ╰─ Compression (2)
╰─ Format Conversions (4)

*Special cases include: Empty reads, circular references, 
 edge cases, and boundary conditions
```

## Key Features

```plaintext
📂 54 total tests (42 alignments + 12 FASTQ)
🔄 Reproducible outputs      🧬 3 reference genomes
📊 Metadata documentation   🐍 Python-powered
```

## Why This Matters  
Properly formatted test data helps catch tool errors before real data analysis. This suite validates:  
- File format compliance  
- Edge case handling  
- Pipeline robustness

## Quick Start

**1. Create Conda Environment**
```bash
conda env create -f environment.yml
conda activate genie-test
```

**2. Generate Test Files**
```bash
python -m test_data_generator.generate
```

**Output Structure**
```text
test_data/
├── SAM_01/          # SAM format tests
├── SAM_02/  
...
└── FASTQ_12/        # FASTQ edge cases
```

## Test Categories

```text
SAM/BAM/CRAM (42)              FASTQ (12)
├─ Core Alignment Features     ├─ Read Types
│  ├─ Mapped Pairs             ├─ Paired/Single End  
│  ├─ Alignment Flags          ├─ Quality Scores
│  ├─ CIGAR Operations         ╰─ Compression
│  ╰─ Optional Tags
│
├─ Special Cases
│  ├─ Empty Reads
│  ├─ Circular References
│  ├─ Edge Cases
│  ╰─ Boundary Conditions
│
╰─ Format Support
   ├─ SAM/BAM/CRAM Input
   ╰─ SAM/BAM/CRAM Output
```

**Detailed test case descriptions:** [test_cases_descriptions.md](test_cases_descriptions.md)

## Customization  
Modify `config.py` to:  
- Change output locations  
- Select specific test cases  
- Adjust reference genomes

---

*"Never trust a bioinformatics tool you haven't tried to break first"* - Genomics proverb
