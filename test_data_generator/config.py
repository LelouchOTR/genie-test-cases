from pathlib import Path

# Base directory where test data will be generated
BASE_OUTPUT_DIR = Path("test_data")

# Test case configurations
TEST_CASES = [
    # === SAM/BAM Test Cases ===
    {
        "id": "SAM_01",
        "name": "Unmapped read - single end",
        "format": "sam",
        "output_subdir": "SAM_01",
        "generator_func": "sam_bam_generators.generate_sam_01",
        "params": {},
        "description": "Single-end unmapped read with basic flags",
        "notes": "Verify handling of unmapped reads",
        "output_files": ["alignment.sam"]
    },
    {
        "id": "SAM_02",
        "name": "Unmapped pair",
        "format": "sam",
        "output_subdir": "SAM_02",
        "generator_func": "sam_bam_generators.generate_sam_02",
        "params": {},
        "description": "Paired-end reads where both mates are unmapped",
        "output_files": ["alignment.sam"]
    },
    {
        "id": "SAM_03",
        "name": "Half-mapped read pair",
        "format": "sam",
        "output_subdir": "SAM_03",
        "generator_func": "sam_bam_generators.generate_sam_03",
        "params": {},
        "description": "Paired-end reads where one mate is mapped and one is unmapped",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_04",
        "name": "Mapped read single end",
        "format": "sam",
        "output_subdir": "SAM_04",
        "generator_func": "sam_bam_generators.generate_sam_04",
        "params": {},
        "description": "Single-end mapped reads with different orientations",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_05",
        "name": "Mapped read pair - same position + TLEN",
        "format": "sam",
        "output_subdir": "SAM_05",
        "generator_func": "sam_bam_generators.generate_sam_05",
        "params": {},
        "description": "Paired-end reads mapped to same position with template length",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_39",
        "name": "(BAM input) - Generate BAM",
        "format": "bam",
        "output_subdir": "SAM_39",
        "generator_func": "sam_bam_generators.generate_sam_39",
        "params": {},
        "description": "Basic BAM file generation",
        "output_files": ["alignment.bam", "simple_ref.fa"]
    },

    # === FASTQ Test Cases ===
    {
        "id": "FASTQ_01",
        "name": "Single End - constant read length",
        "format": "fastq",
        "output_subdir": "FASTQ_01",
        "generator_func": "fastq_generators.generate_fastq_01",
        "params": {},
        "description": "Single-end reads with consistent length",
        "output_files": ["reads.fq"]
    },
    {
        "id": "FASTQ_02",
        "name": "Single End - variable read length",
        "format": "fastq",
        "output_subdir": "FASTQ_02",
        "generator_func": "fastq_generators.generate_fastq_02",
        "params": {},
        "description": "Single-end reads with varying lengths",
        "output_files": ["reads.fq"]
    },
    {
        "id": "FASTQ_03",
        "name": "Paired End - both mates same length",
        "format": "fastq",
        "output_subdir": "FASTQ_03",
        "generator_func": "fastq_generators.generate_fastq_03",
        "params": {},
        "description": "Paired-end reads where both mates have same length",
        "output_files": ["reads_1.fq", "reads_2.fq"]
    },
    {
        "id": "FASTQ_12",
        "name": "Paired End - unequal file lengths",
        "format": "fastq",
        "output_subdir": "FASTQ_12",
        "generator_func": "fastq_generators.generate_fastq_12",
        "params": {},
        "description": "Paired-end files with unequal numbers of reads",
        "output_files": ["reads_1.fq", "reads_2.fq"]
    }
]
