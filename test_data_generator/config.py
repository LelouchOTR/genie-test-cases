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
        "params": {"special_reference": "large_ref.fa"},
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
        "id": "SAM_06",
        "name": "Mapped read pair – enclosed + TLEN",
        "format": "sam",
        "output_subdir": "SAM_06",
        "generator_func": "sam_bam_generators.generate_sam_06",
        "params": {},
        "description": "Paired reads where one read encloses the other",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_07",
        "name": "Mapped read pair – overlapping + TLEN",
        "format": "sam",
        "output_subdir": "SAM_07",
        "generator_func": "sam_bam_generators.generate_sam_07",
        "params": {},
        "description": "Paired reads with overlapping alignment",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_08",
        "name": "Mapped read pair – no overlapping + TLEN",
        "format": "sam",
        "output_subdir": "SAM_08",
        "generator_func": "sam_bam_generators.generate_sam_08",
        "params": {},
        "description": "Paired reads with no overlap between mates",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_09",
        "name": "Mapped read pair – long distance + TLEN",
        "format": "sam",
        "output_subdir": "SAM_09",
        "generator_func": "sam_bam_generators.generate_sam_09",
        "params": {},
        "description": "Paired reads mapped far apart (>1M bases)",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_10",
        "name": "Mapped read pair – different reference + TLEN",
        "format": "sam",
        "output_subdir": "SAM_10",
        "generator_func": "sam_bam_generators.generate_sam_10",
        "params": {},
        "description": "Paired reads mapped to different reference sequences",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_11",
        "name": "Secondary alignment",
        "format": "sam",
        "output_subdir": "SAM_11",
        "generator_func": "sam_bam_generators.generate_sam_11",
        "params": {},
        "description": "Read with primary and secondary alignments",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_12",
        "name": "Supplementary / chimeric alignment",
        "format": "sam",
        "output_subdir": "SAM_12",
        "generator_func": "sam_bam_generators.generate_sam_12",
        "params": {},
        "description": "Chimeric read with supplementary alignments",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_13",
        "name": "Base substitution (M, =, X)",
        "format": "sam",
        "output_subdir": "SAM_13",
        "generator_func": "sam_bam_generators.generate_sam_13",
        "params": {},
        "description": "Read with M, =, X in CIGAR string",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_14",
        "name": "Base insertion",
        "format": "sam",
        "output_subdir": "SAM_14",
        "generator_func": "sam_bam_generators.generate_sam_14",
        "params": {},
        "description": "Read with I in CIGAR string",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_15",
        "name": "Base deletion",
        "format": "sam",
        "output_subdir": "SAM_15",
        "generator_func": "sam_bam_generators.generate_sam_15",
        "params": {},
        "description": "Read with D in CIGAR string",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_16",
        "name": "Softclips",
        "format": "sam",
        "output_subdir": "SAM_16",
        "generator_func": "sam_bam_generators.generate_sam_16",
        "params": {},
        "description": "Read with S in CIGAR string",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_17",
        "name": "Padding (P)",
        "format": "sam",
        "output_subdir": "SAM_17",
        "generator_func": "sam_bam_generators.generate_sam_17",
        "params": {},
        "description": "Read with P in CIGAR string",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_18",
        "name": "Hardclips",
        "format": "sam",
        "output_subdir": "SAM_18",
        "generator_func": "sam_bam_generators.generate_sam_18",
        "params": {},
        "description": "Read with H in CIGAR string",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_19",
        "name": "PCR duplicate flag",
        "format": "sam",
        "output_subdir": "SAM_19",
        "generator_func": "sam_bam_generators.generate_sam_19",
        "params": {},
        "description": "Reads with different PCR duplicate flag values",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_20",
        "name": "Paired end – different flags per mate",
        "format": "sam",
        "output_subdir": "SAM_20",
        "generator_func": "sam_bam_generators.generate_sam_20",
        "params": {},
        "description": "Pairs with different flag combinations per mate",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_21",
        "name": "next read … flags – unmapped",
        "format": "sam",
        "output_subdir": "SAM_21",
        "generator_func": "sam_bam_generators.generate_sam_21",
        "params": {},
        "description": "Unmapped read with mate flags set",
        "output_files": ["alignment.sam"]
    },
    {
        "id": "SAM_22",
        "name": "next read … flags – half mapped",
        "format": "sam",
        "output_subdir": "SAM_22",
        "generator_func": "sam_bam_generators.generate_sam_22",
        "params": {},
        "description": "Mapped read with unmapped mate flags",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_23",
        "name": "next read … flags – short distance",
        "format": "sam",
        "output_subdir": "SAM_23",
        "generator_func": "sam_bam_generators.generate_sam_23",
        "params": {},
        "description": "Properly paired reads with nearby mates",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_24",
        "name": "next read … flags – long distance",
        "format": "sam",
        "output_subdir": "SAM_24",
        "generator_func": "sam_bam_generators.generate_sam_24",
        "params": {},
        "description": "Reads with mates mapped far apart",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_25",
        "name": "Short intron / splice (N)",
        "format": "sam",
        "output_subdir": "SAM_25",
        "generator_func": "sam_bam_generators.generate_sam_25",
        "params": {},
        "description": "Read with N cigar operator for splicing",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_26",
        "name": "Long intron / splice (N)",
        "format": "sam",
        "output_subdir": "SAM_26",
        "generator_func": "sam_bam_generators.generate_sam_26",
        "params": {},
        "description": "Read with long N cigar operator",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_27",
        "name": "Empty read (all bases deleted)",
        "format": "sam",
        "output_subdir": "SAM_27",
        "generator_func": "sam_bam_generators.generate_sam_27",
        "params": {},
        "description": "Read with all bases deleted (D in CIGAR)",
        "output_files": ["alignment.sam"]
    },
    {
        "id": "SAM_28",
        "name": "Empty read (all bases softclipped)",
        "format": "sam",
        "output_subdir": "SAM_28",
        "generator_func": "sam_bam_generators.generate_sam_28",
        "params": {},
        "description": "Read with all bases softclipped (S in CIGAR)",
        "output_files": ["alignment.sam"]
    },
    {
        "id": "SAM_29",
        "name": "Empty read (all bases hardclipped)",
        "format": "sam",
        "output_subdir": "SAM_29",
        "generator_func": "sam_bam_generators.generate_sam_29",
        "params": {},
        "description": "Read with all bases hardclipped (H in CIGAR)",
        "output_files": ["alignment.sam"]
    },
    {
        "id": "SAM_30",
        "name": "Empty read (no nucleotides in read / * in sam)",
        "format": "sam",
        "output_subdir": "SAM_30",
        "generator_func": "sam_bam_generators.generate_sam_30",
        "params": {},
        "description": "Read with SEQ and QUAL set to '*'",
        "output_files": ["alignment.sam"]
    },
    {
        "id": "SAM_31",
        "name": "Quality scores absent",
        "format": "sam",
        "output_subdir": "SAM_31",
        "generator_func": "sam_bam_generators.generate_sam_31",
        "params": {},
        "description": "Read with QUAL set to '*'",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_32",
        "name": "Optional tags",
        "format": "sam",
        "output_subdir": "SAM_32",
        "generator_func": "sam_bam_generators.generate_sam_32",
        "params": {},
        "description": "Read with various optional tags",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_33",
        "name": "Read groups",
        "format": "sam",
        "output_subdir": "SAM_33",
        "generator_func": "sam_bam_generators.generate_sam_33",
        "params": {},
        "description": "Reads with different read group tags",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_34",
        "name": "Reverse Complement (different + same) – short distance",
        "format": "sam",
        "output_subdir": "SAM_34",
        "generator_func": "sam_bam_generators.generate_sam_34",
        "params": {},
        "description": "Pairs with different reverse complement flag combinations",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_35",
        "name": "Reverse Complement (different + same) – long distance",
        "format": "sam",
        "output_subdir": "SAM_35",
        "generator_func": "sam_bam_generators.generate_sam_35",
        "params": {},
        "description": "Pairs with reverse flags and large separation",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_36",
        "name": "Reverse Complement (different + same) – unmapped",
        "format": "sam",
        "output_subdir": "SAM_36",
        "generator_func": "sam_bam_generators.generate_sam_36",
        "params": {},
        "description": "Unmapped pairs with different reverse flags",
        "output_files": ["alignment.sam"]
    },
    {
        "id": "SAM_37",
        "name": "Reverse Complement (different + same) – half mapped",
        "format": "sam",
        "output_subdir": "SAM_37",
        "generator_func": "sam_bam_generators.generate_sam_37",
        "params": {},
        "description": "Half-mapped pairs with different reverse flags",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_38",
        "name": "Circular reference",
        "format": "sam",
        "output_subdir": "SAM_38",
        "generator_func": "sam_bam_generators.generate_sam_38",
        "params": {},
        "description": "Read overlapping circular reference boundary",
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
    {
        "id": "SAM_40",
        "name": "(BAM output) - Generate SAM",
        "format": "sam",
        "output_subdir": "SAM_40",
        "generator_func": "sam_bam_generators.generate_sam_40",
        "params": {},
        "description": "Test BAM output capability",
        "output_files": ["alignment.sam", "simple_ref.fa"]
    },
    {
        "id": "SAM_41",
        "name": "(CRAM input) - Generate CRAM",
        "format": "cram",
        "output_subdir": "SAM_41",
        "generator_func": "sam_bam_generators.generate_sam_41",
        "params": {},
        "description": "Basic CRAM file generation",
        "output_files": ["alignment.cram", "simple_ref.fa"]
    },
    {
        "id": "SAM_42",
        "name": "(CRAM output) - Generate SAM",
        "format": "sam",
        "output_subdir": "SAM_42",
        "generator_func": "sam_bam_generators.generate_sam_42",
        "params": {},
        "description": "Test CRAM output capability",
        "output_files": ["alignment.sam", "simple_ref.fa"]
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
        "id": "FASTQ_04",
        "name": "Paired End - mates different length",
        "format": "fastq", 
        "output_subdir": "FASTQ_04",
        "generator_func": "fastq_generators.generate_fastq_04",
        "params": {},
        "description": "Paired-end reads with different lengths per mate",
        "output_files": ["reads_1.fq", "reads_2.fq"]
    },
    {
        "id": "FASTQ_05",
        "name": "GZIP compressed input",
        "format": "fastq",
        "output_subdir": "FASTQ_05",
        "generator_func": "fastq_generators.generate_fastq_05",
        "params": {},
        "description": "FASTQ input compressed with gzip",
        "output_files": ["reads.fq.gz"]
    },
    {
        "id": "FASTQ_06",
        "name": "GZIP compressed output",
        "format": "fastq",
        "output_subdir": "FASTQ_06",
        "generator_func": "fastq_generators.generate_fastq_06",
        "params": {},
        "description": "Test tool's ability to output gzipped FASTQ",
        "output_files": ["reads.fq"]
    },
    {
        "id": "FASTQ_07",
        "name": "Paired End - different read names", 
        "format": "fastq",
        "output_subdir": "FASTQ_07",
        "generator_func": "fastq_generators.generate_fastq_07",
        "params": {},
        "description": "Paired-end files with non-matching read names",
        "output_files": ["reads_1.fq", "reads_2.fq"]
    },
    {
        "id": "FASTQ_08",
        "name": "Full quality score range",
        "format": "fastq",
        "output_subdir": "FASTQ_08",
        "generator_func": "fastq_generators.generate_fastq_08",
        "params": {},
        "description": "All valid Phred+33 quality scores",
        "output_files": ["reads.fq"]
    },
    {
        "id": "FASTQ_09",
        "name": "ACGT nucleotide only",
        "format": "fastq",
        "output_subdir": "FASTQ_09",
        "generator_func": "fastq_generators.generate_fastq_09",
        "params": {},
        "description": "Reads containing only ACGT bases",
        "output_files": ["reads.fq"]
    },
    {
        "id": "FASTQ_10",
        "name": "IUPAC ambiguity codes",
        "format": "fastq",
        "output_subdir": "FASTQ_10",
        "generator_func": "fastq_generators.generate_fastq_10",
        "params": {},
        "description": "Reads with IUPAC ambiguity characters",
        "output_files": ["reads.fq"]
    },
    {
        "id": "FASTQ_11",
        "name": "Special read name characters",
        "format": "fastq",
        "output_subdir": "FASTQ_11",
        "generator_func": "fastq_generators.generate_fastq_11",
        "params": {},
        "description": "Read names with special characters",
        "output_files": ["reads.fq"]
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
