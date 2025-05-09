from pathlib import Path
from . import utils  # Use relative import
import gzip
import shutil
import string
from tqdm import tqdm

# === FASTQ Generators ===

def generate_fastq_01(output_dir: Path, **kwargs):
    """FASTQ_01: Single End – constant read length"""
    file_path = output_dir / "reads.fastq"
    content = ""
    content += utils.create_fastq_entry("read1_const", "ACGTACGTACGT", "!!!!!!!!!!!!") # len 12
    content += utils.create_fastq_entry("read2_const", "CCCCCCCCCCCC", "############") # len 12
    content += utils.create_fastq_entry("read3_const", "GGGGGGGGGGGG", "$$$$$$$$$$$$") # len 12
    with open(file_path, "w") as f:
        f.write(content)

def generate_fastq_02(output_dir: Path, **kwargs):
    """FASTQ_02: Single End – variable read length"""
    file_path = output_dir / "reads.fastq"
    content = ""
    content += utils.create_fastq_entry("read1_var", "ACGTACGTACGT", "!!!!!!!!!!!!") # len 12
    content += utils.create_fastq_entry("read2_var", "CCCCCCCC", "########")       # len 8
    content += utils.create_fastq_entry("read3_var", "GGGGGGGGGGGGGG", "$$$$$$$$$$$$$$") # len 14
    with open(file_path, "w") as f:
        f.write(content)

def generate_fastq_03(output_dir: Path, **kwargs):
    """FASTQ_03: Paired End - both mates same length"""
    file1_path = output_dir / "reads_1.fastq"
    file2_path = output_dir / "reads_2.fastq"
    content1 = ""
    content2 = ""
    read_len = 15

    content1 += utils.create_fastq_entry("pair1_same/1", "A"*read_len, "!"*read_len)
    content2 += utils.create_fastq_entry("pair1_same/2", "T"*read_len, "#"*read_len)
    content1 += utils.create_fastq_entry("pair2_same/1", "C"*read_len, "$"*read_len)
    content2 += utils.create_fastq_entry("pair2_same/2", "G"*read_len, "%"*read_len)
    content1 += utils.create_fastq_entry("pair3_same/1", "N"*read_len, "&"*read_len)
    content2 += utils.create_fastq_entry("pair3_same/2", "A"*read_len, "'"*read_len)

    with open(file1_path, "w") as f1, open(file2_path, "w") as f2:
        f1.write(content1)
        f2.write(content2)

def generate_fastq_04(output_dir: Path, **kwargs):
    """FASTQ_04: Paired End - mates different length"""
    file1_path = output_dir / "reads_1.fastq"
    file2_path = output_dir / "reads_2.fastq"
    content1 = ""
    content2 = ""
    len1 = 12
    len2 = 10

    content1 += utils.create_fastq_entry("pair1_diff/1", "A"*len1, "!"*len1)
    content2 += utils.create_fastq_entry("pair1_diff/2", "T"*len2, "#"*len2)
    content1 += utils.create_fastq_entry("pair2_diff/1", "C"*len1, "$"*len1)
    content2 += utils.create_fastq_entry("pair2_diff/2", "G"*len2, "%"*len2)
    content1 += utils.create_fastq_entry("pair3_diff/1", "N"*len1, "&"*len1)
    content2 += utils.create_fastq_entry("pair3_diff/2", "A"*len2, "'"*len2)

    with open(file1_path, "w") as f1, open(file2_path, "w") as f2:
        f1.write(content1)
        f2.write(content2)

def generate_fastq_05(output_dir: Path, **kwargs):
    """FASTQ_05: (fastq.gz input)"""
    # Generate the content first
    file_path_fq = output_dir / "reads.fastq.temp"
    content = ""
    content += utils.create_fastq_entry("read1_gz", "ACGTACGTACGT", "!!!!!!!!!!!!")
    content += utils.create_fastq_entry("read2_gz", "CCCCCCCCCCCC", "############")
    content += utils.create_fastq_entry("read3_gz", "GGGGGGGGGGGG", "$$$$$$$$$$$$")
    with open(file_path_fq, "w") as f:
        f.write(content)

    # Gzip the file
    file_path_gz = output_dir / "reads.fastq.gz"
    with open(file_path_fq, 'rb') as f_in, gzip.open(file_path_gz, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    file_path_fq.unlink() # Remove temporary file

def generate_fastq_06(output_dir: Path, **kwargs):
    """FASTQ_06: (fastq.gz output) - Generate plain fastq"""
    # This case expects the *tool under test* to produce gz output.
    # So we just provide a standard fastq file as input.
    generate_fastq_01(output_dir, **kwargs) # Reuse generator for constant length

def generate_fastq_07(output_dir: Path, **kwargs):
    """FASTQ_07: Paired End - different read names"""
    file1_path = output_dir / "reads_1.fastq"
    file2_path = output_dir / "reads_2.fastq"
    content1 = ""
    content2 = ""
    read_len = 10

    # Mates have completely different base names, no /1 or /2 convention
    content1 += utils.create_fastq_entry("pair1_read_A", "A"*read_len, "!"*read_len)
    content2 += utils.create_fastq_entry("pair1_read_B", "T"*read_len, "#"*read_len)
    content1 += utils.create_fastq_entry("pair2_forward", "C"*read_len, "$"*read_len)
    content2 += utils.create_fastq_entry("pair2_reverse", "G"*read_len, "%"*read_len)
    content1 += utils.create_fastq_entry("pair3_left", "G"*read_len, "&"*read_len)
    content2 += utils.create_fastq_entry("pair3_rightSide", "T"*read_len, "'"*read_len)

    with open(file1_path, "w") as f1, open(file2_path, "w") as f2:
        f1.write(content1)
        f2.write(content2)

def generate_fastq_08(output_dir: Path, **kwargs):
    """FASTQ_08: Alphabet test Quality score range"""
    file_path = output_dir / "reads.fastq"
    # Phred+33 ASCII range: ! (33) to ~ (126)
    qual_chars = "".join([chr(i) for i in range(33, 127)])
    # Sequence needs to match length
    sequence = "A" * len(qual_chars)
    content = utils.create_fastq_entry("qual_range_read", sequence, qual_chars)
    with open(file_path, "w") as f:
        f.write(content)

def generate_fastq_09(output_dir: Path, **kwargs):
    """FASTQ_09: Alphabet test nucleotide alphabet ACTG"""
    file_path = output_dir / "reads.fastq"
    # Strictly ACGT, repeating
    sequence = "ACGT" * 10
    quality = "!" * len(sequence)
    content = utils.create_fastq_entry("acgt_read", sequence, quality)
    with open(file_path, "w") as f:
        f.write(content)

def generate_fastq_10(output_dir: Path, **kwargs):
    """FASTQ_10: Alphabet test nucleotide alphabet IUPAC"""
    file_path = output_dir / "reads.fastq"
    # Standard IUPAC ambiguity codes + N
    sequence = "ACGTRYSWKMBDHVN" * 2
    quality = "!" * len(sequence)
    content = utils.create_fastq_entry("iupac_read", sequence, quality)
    with open(file_path, "w") as f:
        f.write(content)

def generate_fastq_11(output_dir: Path, **kwargs):
    """FASTQ_11: Read name alphabet test"""
    file_path = output_dir / "reads.fastq"
    # Allowed: Any # tqdm.writeable ASCII except whitespace chars like space, tab
    # Let's try a good range, excluding space (chr 32)
    read_name = "".join([chr(i) for i in range(33, 127)])
    sequence = "ACGT"
    quality = "!!!!"
    content = utils.create_fastq_entry(read_name, sequence, quality)
    with open(file_path, "w") as f:
        f.write(content)

def generate_fastq_12(output_dir: Path, **kwargs):
    """FASTQ_12: Paired End - unequal file lengths"""
    file1_path = output_dir / "reads_1.fastq"
    file2_path = output_dir / "reads_2.fastq"
    content1 = ""
    content2 = ""
    read_len = 10

    # Two normal pairs
    content1 += utils.create_fastq_entry("unequal_pair1/1", "A"*read_len, "!"*read_len)
    content2 += utils.create_fastq_entry("unequal_pair1/2", "T"*read_len, "#"*read_len)
    content1 += utils.create_fastq_entry("unequal_pair2/1", "C"*read_len, "$"*read_len)
    content2 += utils.create_fastq_entry("unequal_pair2/2", "G"*read_len, "%"*read_len)
    # Extra read in file 1 (only) to make unequal lengths
    content1 += utils.create_fastq_entry("extra_read/1", "N"*read_len, "&"*read_len)

    with open(file1_path, "w") as f1, open(file2_path, "w") as f2:
        f1.write(content1)
        f2.write(content2)
