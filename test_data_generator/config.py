"""Configuration for test data generation."""
from pathlib import Path

# Base directory where test data will be generated (relative to project root)
BASE_OUTPUT_DIR = Path("test_data")

# Test case configurations
TEST_CASES = [
    {
        "id": "01",
        "name": "Single End - constant read length", 
        "format": "fastq",
        "output_subdir": "01",
        "generator_func": "fastq_generators.generate_fastq_01",
        "params": {},
        "description": "Single-end reads with consistent length",
        "output_files": ["reads.fastq"]
    },
    {
        "id": "02",
        "name": "Single End - variable read length",
        "format": "fastq",
        "output_subdir": "02",
        "generator_func": "fastq_generators.generate_fastq_02",
        "params": {},
        "description": "Single-end reads with varying lengths",
        "output_files": ["reads.fastq"]
    },
    {
        "id": "03",
        "name": "Paired End - both mates same length",
        "format": "fastq",
        "output_subdir": "03",
        "generator_func": "fastq_generators.generate_fastq_03",
        "params": {},
        "description": "Paired-end reads where both mates have same length",
        "output_files": ["reads_1.fastq", "reads_2.fastq"]
    },
    {
        "id": "04",
        "name": "Paired End - mates different length",
        "format": "fastq", 
        "output_subdir": "04",
        "generator_func": "fastq_generators.generate_fastq_04",
        "params": {},
        "description": "Paired-end reads with different lengths per mate",
        "output_files": ["reads_1.fastq", "reads_2.fastq"]
    },
    {
        "id": "05",
        "name": "GZIP compressed input",
        "format": "fastq",
        "output_subdir": "05",
        "generator_func": "fastq_generators.generate_fastq_05",
        "params": {},
        "description": "FASTQ input compressed with gzip",
        "output_files": ["reads.fastq.gz"]
    },
    {
        "id": "06",
        "name": "GZIP compressed output",
        "format": "fastq",
        "output_subdir": "06",
        "generator_func": "fastq_generators.generate_fastq_06",
        "params": {},
        "description": "Test tool's ability to output gzipped FASTQ",
        "output_files": ["reads.fastq"]
    },
    {
        "id": "07",
        "name": "Paired End - different read names", 
        "format": "fastq",
        "output_subdir": "07",
        "generator_func": "fastq_generators.generate_fastq_07",
        "params": {},
        "description": "Paired-end files with non-matching read names",
        "output_files": ["reads_1.fastq", "reads_2.fastq"]
    },
    {
        "id": "08",
        "name": "Full quality score range",
        "format": "fastq",
        "output_subdir": "08",
        "generator_func": "fastq_generators.generate_fastq_08",
        "params": {},
        "description": "All valid Phred+33 quality scores",
        "output_files": ["reads.fastq"]
    },
    {
        "id": "09",
        "name": "ACGT nucleotide only",
        "format": "fastq",
        "output_subdir": "09",
        "generator_func": "fastq_generators.generate_fastq_09",
        "params": {},
        "description": "Reads containing only ACGT bases",
        "output_files": ["reads.fastq"]
    },
    {
        "id": "10",
        "name": "IUPAC ambiguity codes",
        "format": "fastq",
        "output_subdir": "10",
        "generator_func": "fastq_generators.generate_fastq_10",
        "params": {},
        "description": "Reads with IUPAC ambiguity characters",
        "output_files": ["reads.fastq"]
    },
    {
        "id": "11",
        "name": "Special read name characters",
        "format": "fastq",
        "output_subdir": "11",
        "generator_func": "fastq_generators.generate_fastq_11",
        "params": {},
        "description": "Read names with special characters",
        "output_files": ["reads.fastq"]
    },
    {
        "id": "12",
        "name": "Paired End - unequal file lengths",
        "format": "fastq",
        "output_subdir": "12",
        "generator_func": "fastq_generators.generate_fastq_12",
        "params": {},
        "description": "Paired-end files with unequal numbers of reads",
        "output_files": ["reads_1.fastq", "reads_2.fastq"]
    },
    # === SAM/BAM Test Cases ===
    # === SAM/BAM/CRAM Test Cases 13-54 ===
    {
        "id": "13",
        "name": "Unmapped read - single end",
        "format": "sam",
        "output_subdir": "13",
        "generator_func": "sam_bam_generators.generate_sam_13",
        "params": {"special_reference": "large_ref.fa"},
        "description": "Single-end unmapped read with basic flags",
        "notes": "Verify handling of unmapped reads",
        "output_files": ["alignment.sam"]
    },
    {
        "id": "14",
        "name": "Unmapped pair",
        "format": "sam",
        "output_subdir": "14",
        "generator_func": "sam_bam_generators.generate_sam_14",
        "params": {},
        "description": "Paired-end reads where both mates are unmapped",
        "output_files": ["alignment.sam"]
    },
    {
        "id": "15",
        "name": "Half-mapped read pair",
        "format": "sam",
        "output_subdir": "15",
        "generator_func": "sam_bam_generators.generate_sam_15",
        "params": {},
        "description": "Paired-end reads where one mate is mapped and one is unmapped",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "16",
        "name": "Mapped read single end",
        "format": "sam",
        "output_subdir": "16",
        "generator_func": "sam_bam_generators.generate_sam_16",
        "params": {},
        "description": "Single-end mapped reads with different orientations",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "17",
        "name": "Mapped read pair - same position + TLEN",
        "format": "sam",
        "output_subdir": "17",
        "generator_func": "sam_bam_generators.generate_sam_17",
        "params": {},
        "description": "Paired-end reads mapped to same position with template length",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "18",
        "name": "Mapped read pair – enclosed + TLEN",
        "format": "sam",
        "output_subdir": "18",
        "generator_func": "sam_bam_generators.generate_sam_18",
        "params": {},
        "description": "Paired reads where one read encloses the other",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "19",
        "name": "Mapped read pair – overlapping + TLEN",
        "format": "sam",
        "output_subdir": "19",
        "generator_func": "sam_bam_generators.generate_sam_19",
        "params": {},
        "description": "Paired reads with overlapping alignment",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "20",
        "name": "Mapped read pair – no overlapping + TLEN",
        "format": "sam",
        "output_subdir": "20",
        "generator_func": "sam_bam_generators.generate_sam_20",
        "params": {},
        "description": "Paired reads with no overlap between mates",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "21",
        "name": "Mapped read pair – long distance + TLEN",
        "format": "sam",
        "output_subdir": "21",
        "generator_func": "sam_bam_generators.generate_sam_21",
        "params": {"special_reference": "large_ref.fa"},
        "description": "Paired reads mapped far apart (>1M bases)",
        "output_files": ["alignment.sam", "large_ref.fa"]
    },
    {
        "id": "22",
        "name": "Mapped read pair – different reference + TLEN",
        "format": "sam",
        "output_subdir": "22",
        "generator_func": "sam_bam_generators.generate_sam_22",
        "params": {},
        "description": "Paired reads mapped to different reference sequences",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "23",
        "name": "Secondary alignment",
        "format": "sam",
        "output_subdir": "23",
        "generator_func": "sam_bam_generators.generate_sam_23",
        "params": {},
        "description": "Read with primary and secondary alignments",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "24",
        "name": "Supplementary / chimeric alignment",
        "format": "sam",
        "output_subdir": "24",
        "generator_func": "sam_bam_generators.generate_sam_24",
        "params": {},
        "description": "Chimeric read with supplementary alignments",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "25",
        "name": "Base substitution (M, =, X)",
        "format": "sam",
        "output_subdir": "25",
        "generator_func": "sam_bam_generators.generate_sam_13",
        "params": {},
        "description": "Read with M, =, X in CIGAR string",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "26",
        "name": "Base insertion",
        "format": "sam",
        "output_subdir": "26",
        "generator_func": "sam_bam_generators.generate_sam_14",
        "params": {},
        "description": "Read with I in CIGAR string",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "27",
        "name": "Base deletion",
        "format": "sam",
        "output_subdir": "27",
        "generator_func": "sam_bam_generators.generate_sam_15",
        "params": {},
        "description": "Read with D in CIGAR string",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "28",
        "name": "Softclips",
        "format": "sam",
        "output_subdir": "28",
        "generator_func": "sam_bam_generators.generate_sam_16",
        "params": {},
        "description": "Read with S in CIGAR string",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "29",
        "name": "Padding (P)",
        "format": "sam",
        "output_subdir": "29",
        "generator_func": "sam_bam_generators.generate_sam_17",
        "params": {},
        "description": "Read with P in CIGAR string",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "30",
        "name": "Hardclips",
        "format": "sam",
        "output_subdir": "30",
        "generator_func": "sam_bam_generators.generate_sam_18",
        "params": {},
        "description": "Read with H in CIGAR string",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "31",
        "name": "PCR duplicate flag",
        "format": "sam",
        "output_subdir": "31",
        "generator_func": "sam_bam_generators.generate_sam_19",
        "params": {},
        "description": "Reads with different PCR duplicate flag values",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "32",
        "name": "Paired end – different flags per mate",
        "format": "sam",
        "output_subdir": "32",
        "generator_func": "sam_bam_generators.generate_sam_20",
        "params": {},
        "description": "Pairs with different flag combinations per mate",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "33",
        "name": "next read … flags – unmapped",
        "format": "sam",
        "output_subdir": "33",
        "generator_func": "sam_bam_generators.generate_sam_21",
        "params": {},
        "description": "Unmapped read with mate flags set",
        "output_files": ["alignment.sam"]
    },
    {
        "id": "34",
        "name": "next read … flags – half mapped",
        "format": "sam",
        "output_subdir": "34",
        "generator_func": "sam_bam_generators.generate_sam_22",
        "params": {},
        "description": "Mapped read with unmapped mate flags",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "35",
        "name": "next read … flags – short distance",
        "format": "sam",
        "output_subdir": "35",
        "generator_func": "sam_bam_generators.generate_sam_23",
        "params": {},
        "description": "Properly paired reads with nearby mates",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "36",
        "name": "next read … flags – long distance",
        "format": "sam",
        "output_subdir": "36",
        "generator_func": "sam_bam_generators.generate_sam_24",
        "params": {"special_reference": "large_ref.fa"},
        "description": "Reads with mates mapped far apart",
        "output_files": ["alignment.sam", "large_ref.fa"]
    },
    {
        "id": "37",
        "name": "Short intron / splice (N)",
        "format": "sam",
        "output_subdir": "37",
        "generator_func": "sam_bam_generators.generate_sam_25",
        "params": {},
        "description": "Read with N cigar operator for splicing",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "38",
        "name": "Long intron / splice (N)",
        "format": "sam",
        "output_subdir": "38",
        "generator_func": "sam_bam_generators.generate_sam_26",
        "params": {},
        "description": "Read with long N cigar operator",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "39",
        "name": "Empty read (all bases deleted)",
        "format": "sam",
        "output_subdir": "39",
        "generator_func": "sam_bam_generators.generate_sam_27",
        "params": {},
        "description": "Read with all bases deleted (D in CIGAR)",
        "output_files": ["alignment.sam"]
    },
    {
        "id": "40",
        "name": "Empty read (all bases softclipped)",
        "format": "sam",
        "output_subdir": "40",
        "generator_func": "sam_bam_generators.generate_sam_28",
        "params": {},
        "description": "Read with all bases softclipped (S in CIGAR)",
        "output_files": ["alignment.sam"]
    },
    {
        "id": "41",
        "name": "Empty read (all bases hardclipped)",
        "format": "sam",
        "output_subdir": "41",
        "generator_func": "sam_bam_generators.generate_sam_29",
        "params": {},
        "description": "Read with all bases hardclipped (H in CIGAR)",
        "output_files": ["alignment.sam"]
    },
    {
        "id": "42",
        "name": "Empty read (no nucleotides in read / * in sam)",
        "format": "sam",
        "output_subdir": "42",
        "generator_func": "sam_bam_generators.generate_sam_30",
        "params": {},
        "description": "Read with SEQ and QUAL set to '*'",
        "output_files": ["alignment.sam"]
    },
    {
        "id": "43",
        "name": "Quality scores absent",
        "format": "sam",
        "output_subdir": "43",
        "generator_func": "sam_bam_generators.generate_sam_31",
        "params": {},
        "description": "Read with QUAL set to '*'",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "44",
        "name": "Optional tags",
        "format": "sam",
        "output_subdir": "44",
        "generator_func": "sam_bam_generators.generate_sam_32",
        "params": {},
        "description": "Read with various optional tags",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "45",
        "name": "Read groups",
        "format": "sam",
        "output_subdir": "45",
        "generator_func": "sam_bam_generators.generate_sam_33",
        "params": {},
        "description": "Reads with different read group tags",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "46",
        "name": "Reverse Complement (different + same) – short distance",
        "format": "sam",
        "output_subdir": "46",
        "generator_func": "sam_bam_generators.generate_sam_34",
        "params": {},
        "description": "Pairs with different reverse complement flag combinations",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "47",
        "name": "Reverse Complement (different + same) – long distance",
        "format": "sam",
        "output_subdir": "47",
        "generator_func": "sam_bam_generators.generate_sam_35",
        "params": {"special_reference": "large_ref.fa"},
        "description": "Pairs with reverse flags and large separation",
        "output_files": ["alignment.sam", "large_ref.fa"]
    },
    {
        "id": "48",
        "name": "Reverse Complement (different + same) – unmapped",
        "format": "sam",
        "output_subdir": "48",
        "generator_func": "sam_bam_generators.generate_sam_36",
        "params": {},
        "description": "Unmapped pairs with different reverse flags",
        "output_files": ["alignment.sam"]
    },
    {
        "id": "49",
        "name": "Reverse Complement (different + same) – half mapped",
        "format": "sam",
        "output_subdir": "49",
        "generator_func": "sam_bam_generators.generate_sam_37",
        "params": {},
        "description": "Half-mapped pairs with different reverse flags",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "50",
        "name": "Circular reference",
        "format": "sam",
        "output_subdir": "50",
        "generator_func": "sam_bam_generators.generate_sam_38",
        "params": {},
        "description": "Read overlapping circular reference boundary",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "51",
        "name": "(BAM input) - Generate BAM",
        "format": "bam",
        "output_subdir": "51",
        "generator_func": "sam_bam_generators.generate_sam_39",
        "params": {},
        "description": "Basic BAM file generation",
        "output_files": ["alignment.bam", "simple_ref.fa"]
    },
    {
        "id": "52",
        "name": "(BAM output) - Generate SAM",
        "format": "sam",
        "output_subdir": "52",
        "generator_func": "sam_bam_generators.generate_sam_40",
        "params": {},
        "description": "Test BAM output capability",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "53",
        "name": "(CRAM input) - Generate CRAM",
        "format": "cram",
        "output_subdir": "53",
        "generator_func": "sam_bam_generators.generate_sam_41",
        "params": {},
        "description": "Basic CRAM file generation",
        "output_files": ["alignment.cram", "simple_ref.fa"]
    },
    {
        "id": "54",
        "name": "(CRAM output) - Generate SAM",
        "format": "sam",
        "output_subdir": "54",
        "generator_func": "sam_bam_generators.generate_sam_42",
        "params": {},
        "description": "Test CRAM output capability",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },

]
