from pathlib import Path
import pysam
from . import utils # Relative import
from tqdm import tqdm
from colorama import Fore

# Helper function for setting common paired flags
def set_paired_flags(segment: pysam.AlignedSegment, is_read1: bool):
    segment.is_paired = True
    segment.is_read1 = is_read1
    segment.is_read2 = not is_read1

# Helper function for setting mate info
def set_mate_info(segment: pysam.AlignedSegment, header: pysam.AlignmentHeader,
                  mate_ref_name: str | None, mate_start: int,
                  mate_is_unmapped: bool, mate_is_reverse: bool):
    segment.mate_is_unmapped = mate_is_unmapped
    segment.mate_is_reverse = mate_is_reverse
    
    if mate_ref_name is not None and not mate_is_unmapped:
        try:
            ref_id = header.references.index(mate_ref_name)
        except ValueError:
            raise
        
        segment.next_reference_id = ref_id
        segment.next_reference_start = mate_start
    else:
        segment.next_reference_id = segment.reference_id if segment.reference_id is not None else -1
        segment.next_reference_start = segment.reference_start if segment.reference_start is not None else -1


# === SAM/BAM/CRAM Generators ===

def generate_sam_01(output_dir: Path, **kwargs):
    """SAM_01: Unmapped read – single end"""
    file_path = output_dir / "alignment.sam"
    header = utils.get_default_sam_header()

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        a = pysam.AlignedSegment()
        a.query_name = "unmapped_se_1"
        a.query_sequence = "AGCTAGCTAGCT"
        a.query_qualities = pysam.qualitystring_to_array("!!!!!!!!!!!!")
        a.flag = 4 # 0x4 = unmapped
        # Ensure not marked as paired
        a.is_paired = False
        samfile.write(a)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")

def generate_sam_02(output_dir: Path, **kwargs):
    """SAM_02: Unmapped pair"""
    file_path = output_dir / "alignment.sam"
    header = utils.get_default_sam_header()

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Common flags for unmapped pair: Paired, Read Unmapped, Mate Unmapped
        # 0x1 (paired) + 0x4 (unmapped) + 0x8 (mate unmapped) = 13
        # Plus 0x40 (read 1) or 0x80 (read 2)
        flag_r1 = 13 + 0x40 # = 77
        flag_r2 = 13 + 0x80 # = 141

        # Pair 1 (R1 first)
        r1 = pysam.AlignedSegment()
        r1.query_name = "unmapped_pair_1"
        r1.query_sequence = "AAAA" * 3
        r1.query_qualities = pysam.qualitystring_to_array("!!!!"*3)
        r1.flag = flag_r1
        samfile.write(r1)
        r2 = pysam.AlignedSegment()
        r2.query_name = "unmapped_pair_1"
        r2.query_sequence = "TTTT" * 3
        r2.query_qualities = pysam.qualitystring_to_array("####"*3)
        r2.flag = flag_r2
        samfile.write(r2)

        # Pair 2 (R2 first)
        r2 = pysam.AlignedSegment()
        r2.query_name = "unmapped_pair_2"
        r2.query_sequence = "GGGG" * 3
        r2.query_qualities = pysam.qualitystring_to_array("$$$$"*3)
        r2.flag = flag_r2 # R2 flags
        samfile.write(r2)
        r1 = pysam.AlignedSegment()
        r1.query_name = "unmapped_pair_2"
        r1.query_sequence = "CCCC" * 3
        r1.query_qualities = pysam.qualitystring_to_array("%%%%"*3)
        r1.flag = flag_r1 # R1 flags
        samfile.write(r1)

        # Pair 3 (R1 first)
        r1 = pysam.AlignedSegment()
        r1.query_name = "unmapped_pair_3"
        r1.query_sequence = "ACGT" * 3
        r1.query_qualities = pysam.qualitystring_to_array("&&&&"*3)
        r1.flag = flag_r1
        samfile.write(r1)
        r2 = pysam.AlignedSegment()
        r2.query_name = "unmapped_pair_3"
        r2.query_sequence = "TGCA" * 3
        r2.query_qualities = pysam.qualitystring_to_array("****"*3)
        r2.flag = flag_r2
        samfile.write(r2)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")

def generate_sam_03(output_dir: Path, **kwargs):
    """SAM_03: Half-mapped read pair"""
    file_path = output_dir / "alignment.sam"
    # Use special large reference specified in test case params
    ref_path = utils.copy_reference_to_output(output_dir, ref_name="large_ref.fa")
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    map_pos = 10
    read_len = 12

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Pair 1: R1 mapped, R2 unmapped (R1 first)
        r1 = pysam.AlignedSegment()
        r1.query_name = "half_mapped_1"
        r1.query_sequence = "A" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("!" * read_len)
        set_paired_flags(r1, is_read1=True)
        r1.reference_id = ref_id
        r1.reference_start = map_pos
        r1.mapping_quality = 60
        r1.cigarstring = f"{read_len}M"
        # Mate is unmapped, R2, same ref, pos=0 (convention)
        set_mate_info(r1, header, mate_ref_name="*", mate_start=0, mate_is_unmapped=True, mate_is_reverse=False)
        # Flags: 0x1 (paired) + 0x8 (mate unmapped) + 0x40 (R1) = 73
        r1.flag = 73
        samfile.write(r1)

        r2 = pysam.AlignedSegment()
        r2.query_name = "half_mapped_1"
        r2.query_sequence = utils.reverse_complement(r1.query_sequence)  # Custom reverse complement
        r2.query_qualities = pysam.qualitystring_to_array("#" * read_len)
        set_paired_flags(r2, is_read1=False)
        r2.is_unmapped = True
        # Mate is mapped, R1, same ref, pos=map_pos
        set_mate_info(r2, header, mate_ref_name="*", mate_start=0, mate_is_unmapped=True, mate_is_reverse=False)
        # Flags: 0x1 (paired) + 0x4 (read unmapped) + 0x80 (R2) = 133
        r2.flag = 133
        samfile.write(r2)

        # Pair 2: R1 unmapped, R2 mapped (R2 first)
        r2 = pysam.AlignedSegment()
        r2.query_name = "half_mapped_2"
        r2.query_sequence = "G" * read_len
        r2.query_qualities = pysam.qualitystring_to_array("$" * read_len)
        set_paired_flags(r2, is_read1=False)
        r2.reference_id = ref_id
        r2.reference_start = map_pos + 50
        r2.mapping_quality = 60
        r2.cigarstring = f"{read_len}M"
        # Mate is unmapped, R1, same ref, pos=0 (convention)
        set_mate_info(r2, header, mate_ref_name=ref_name, mate_start=-1, mate_is_unmapped=True, mate_is_reverse=False)
        # Flags: 0x1 (paired) + 0x8 (mate unmapped) + 0x80 (R2) = 137
        r2.flag = 137
        samfile.write(r2)

        r1 = pysam.AlignedSegment()
        r1.query_name = "half_mapped_2"
        r1.query_sequence = "C" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("%" * read_len)
        set_paired_flags(r1, is_read1=True)
        r1.is_unmapped = True
        # Mate is mapped, R2, same ref, pos=map_pos+50
        set_mate_info(r1, header, mate_ref_name=ref_name, mate_start=map_pos + 50, mate_is_unmapped=False,
mate_is_reverse=False)
        # Flags: 0x1 (paired) + 0x4 (read unmapped) + 0x40 (R1) = 69
        r1.flag = 69
        samfile.write(r1)

        # Pair 3: R1 mapped, R2 unmapped (R1 first) - Mate on different chrom
        r1 = pysam.AlignedSegment()
        r1.query_name = "half_mapped_3"
        r1.query_sequence = "N" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("&" * read_len)
        set_paired_flags(r1, is_read1=True)
        r1.reference_id = ref_id
        r1.reference_start = map_pos + 100
        r1.mapping_quality = 60
        r1.cigarstring = f"{read_len}M"
        # Mate is unmapped, R2, different ref ('ref2'), pos=0 (convention)
        set_mate_info(r1, header, mate_ref_name="ref2", mate_start=-1, mate_is_unmapped=True, mate_is_reverse=False)
        r1.flag = 73 # Same flags as pair 1 R1
        samfile.write(r1)

        r2 = pysam.AlignedSegment()
        r2.query_name = "half_mapped_3"
        r2.query_sequence = "A" * read_len
        r2.query_qualities = pysam.qualitystring_to_array("'" * read_len)
        set_paired_flags(r2, is_read1=False)
        r2.is_unmapped = True
        # Mate is mapped, R1, different ref ('ref1'), pos=map_pos+100
        set_mate_info(r2, header, mate_ref_name=ref_name, mate_start=map_pos + 100, mate_is_unmapped=False,
mate_is_reverse=False)
        r2.flag = 133 # Same flags as pair 1 R2
        samfile.write(r2)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_04(output_dir: Path, **kwargs):
    """SAM_04: Mapped read single end"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    read_len = 12

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Read 1 (matches ref1 start)
        with pysam.FastaFile(str(ref_path)) as fasta:
            seq = fasta.fetch("ref1", 0, read_len).upper()
        r1 = pysam.AlignedSegment()
        r1.query_name = "mapped_se_1"
        r1.query_sequence = seq
        r1.query_qualities = pysam.qualitystring_to_array("!" * read_len)
        r1.reference_id = header.references.index("ref1")
        r1.reference_start = 0
        r1.mapping_quality = 60
        r1.cigarstring = f"{read_len}M"
        r1.flag = 0 # Basic mapped SE read
        samfile.write(r1)

        # Read 2 (matches ref1 elsewhere, reverse strand)
        with pysam.FastaFile(str(ref_path)) as fasta:
            seq = fasta.fetch("ref1", 50, 50 + read_len).upper()
        r2 = pysam.AlignedSegment()
        r2.query_name = "mapped_se_2"
        r2.query_sequence = seq.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]
        r2.query_qualities = pysam.qualitystring_to_array("#" * read_len)
        r2.reference_id = header.references.index("ref1")
        r2.reference_start = 50
        r2.mapping_quality = 60
        r2.cigarstring = f"{read_len}M"
        r2.flag = 16 # 0x10 = reverse strand
        samfile.write(r2)

        # Read 3 (matches ref2 start)
        with pysam.FastaFile(str(ref_path)) as fasta:
            seq = fasta.fetch("ref2", 0, read_len).upper()
        r3 = pysam.AlignedSegment()
        r3.query_name = "mapped_se_3" 
        r3.query_sequence = seq
        r3.query_qualities = pysam.qualitystring_to_array("$" * read_len)
        r3.reference_id = header.references.index("ref2")
        r3.reference_start = 0
        r3.mapping_quality = 60
        r3.cigarstring = f"{read_len}M"
        r3.flag = 0
        samfile.write(r3)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_05(output_dir: Path, **kwargs):
    """SAM_05: Mapped read pair – same position + TLEN"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    map_pos = 20
    read_len = 15

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        r1 = pysam.AlignedSegment()
        r1.query_name = "samepos_pair"
        r1.query_sequence = "A" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("!" * read_len)
        set_paired_flags(r1, is_read1=True)
        r1.reference_id = ref_id
        r1.reference_start = map_pos
        r1.mapping_quality = 60
        r1.cigarstring = f"{read_len}M"
        # Mate info: R2, mapped, same ref, same pos, not reversed
        set_mate_info(r1, header, mate_ref_name=ref_name, mate_start=map_pos, mate_is_unmapped=False,
mate_is_reverse=False)
        # TLEN: For same position, TLEN is +/- read_len depending on strand. Assume FR pair -> TLEN = read_len
        r1.template_length = read_len
        # Flags: 0x1 (paired) + 0x2 (proper pair) + 0x40 (R1) = 67
        r1.flag = 67
        samfile.write(r1)

        r2 = pysam.AlignedSegment()
        r2.query_name = "samepos_pair"
        r2.query_sequence = "T" * read_len
        r2.query_qualities = pysam.qualitystring_to_array("#" * read_len)
        set_paired_flags(r2, is_read1=False)
        r2.reference_id = ref_id
        r2.reference_start = map_pos
        r2.mapping_quality = 60
        r2.cigarstring = f"{read_len}M"
        # Mate info: R1, mapped, same ref, same pos, not reversed
        set_mate_info(r2, header, mate_ref_name=ref_name, mate_start=map_pos, mate_is_unmapped=False,
mate_is_reverse=False)
        # TLEN: For same position, TLEN is +/- read_len. Assume FR pair -> TLEN = -read_len
        r2.template_length = -read_len
        # Flags: 0x1 (paired) + 0x2 (proper pair) + 0x80 (R2) = 131
        r2.flag = 131
        samfile.write(r2)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_06(output_dir: Path, **kwargs):
    """SAM_06: Mapped read pair – enclosed + TLEN"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    r1_start = 30
    r1_len = 50 # R1: 30 - 80
    r2_start = 40
    r2_len = 20 # R2: 40 - 60 (enclosed by R1)

    # Assume FR orientation: R1 forward, R2 reverse
    # TLEN = pos(R2_rightmost) - pos(R1_leftmost) + 1 = (r2_start + r2_len - 1) - r1_start + 1 = 59 - 30 + 1 = 30
    tlen = (r2_start + r2_len - 1) - r1_start + 1

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        r1 = pysam.AlignedSegment()
        r1.query_name = "enclosed_pair"
        r1.query_sequence = "A" * r1_len
        r1.query_qualities = pysam.qualitystring_to_array("!" * r1_len)
        set_paired_flags(r1, is_read1=True)
        r1.reference_id = ref_id
        r1.reference_start = r1_start
        r1.mapping_quality = 60
        r1.cigarstring = f"{r1_len}M"
        # Mate info: R2, mapped, same ref, r2_start, IS reversed
        set_mate_info(r1, header, ref_name, r2_start, False, True)
        r1.template_length = tlen
        # Flags: 0x1 (paired) + 0x2 (proper pair) + 0x20 (mate reverse) + 0x40 (R1) = 99
        r1.flag = 99
        samfile.write(r1)

        r2 = pysam.AlignedSegment()
        r2.query_name = "enclosed_pair"
        r2.query_sequence = "T" * r2_len
        r2.query_qualities = pysam.qualitystring_to_array("#" * r2_len)
        set_paired_flags(r2, is_read1=False)
        r2.reference_id = ref_id
        r2.reference_start = r2_start
        r2.mapping_quality = 60
        r2.cigarstring = f"{r2_len}M"
        r2.is_reverse = True # R2 is reverse
        # Mate info: R1, mapped, same ref, r1_start, NOT reversed
        set_mate_info(r2, header, ref_name, r1_start, False, False)
        r2.template_length = -tlen
        # Flags: 0x1 (paired) + 0x2 (proper pair) + 0x10 (read reverse) + 0x80 (R2) = 147
        r2.flag = 147
        samfile.write(r2)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_07(output_dir: Path, **kwargs):
    """SAM_07: Mapped read pair – overlapping + TLEN"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    r1_start = 50
    r1_len = 40 # R1: 50 - 90
    r2_start = 70
    r2_len = 40 # R2: 70 - 110 (overlaps R1)

    # Calculate TLEN based on actual positions
    tlen = (r2_start + r2_len - 1) - r1_start + 1

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        r1 = pysam.AlignedSegment()
        r1.query_name = "overlap_pair"
        r1.query_sequence = "C" * r1_len
        r1.query_qualities = pysam.qualitystring_to_array("!" * r1_len)
        set_paired_flags(r1, is_read1=True)
        r1.reference_id = ref_id
        r1.reference_start = r1_start
        r1.mapping_quality = 60
        r1.cigarstring = f"{r1_len}M"
        # Mate info: R2, mapped, same ref, r2_start, IS reversed
        set_mate_info(r1, header, mate_ref_name=ref_name, mate_start=r2_start, mate_is_unmapped=False,
mate_is_reverse=True)
        r1.template_length = tlen
        # Flags: 0x1 (paired) + 0x2 (proper pair) + 0x20 (mate reverse) + 0x40 (R1) = 99
        r1.flag = 99
        samfile.write(r1)

        r2 = pysam.AlignedSegment()
        r2.query_name = "overlap_pair"
        r2.query_sequence = "G" * r2_len
        r2.query_qualities = pysam.qualitystring_to_array("#" * r2_len)
        set_paired_flags(r2, is_read1=False)
        r2.reference_id = ref_id
        r2.reference_start = r2_start
        r2.mapping_quality = 60
        r2.cigarstring = f"{r2_len}M"
        r2.is_reverse = True # R2 is reverse
        # Mate info: R1, mapped, same ref, r1_start, NOT reversed
        set_mate_info(r2, header, mate_ref_name=ref_name, mate_start=r1_start, mate_is_unmapped=False,
mate_is_reverse=False)
        r2.template_length = -tlen
        # Flags: 0x1 (paired) + 0x2 (proper pair) + 0x10 (read reverse) + 0x80 (R2) = 147
        r2.flag = 147
        samfile.write(r2)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_08(output_dir: Path, **kwargs):
    """SAM_08: Mapped read pair – no overlapping + TLEN"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    r1_start = 10
    r1_len = 30 # R1: 10 - 40
    r2_start = 60
    r2_len = 30 # R2: 60 - 90 (no overlap)

    # Assume FR orientation: R1 forward, R2 reverse
    # TLEN = pos(R2_rightmost) - pos(R1_leftmost) + 1 = (r2_start + r2_len - 1) - r1_start + 1 = 89 - 10 + 1 = 80
    tlen = (r2_start + r2_len - 1) - r1_start + 1

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        r1 = pysam.AlignedSegment()
        r1.query_name = "no_overlap_pair"
        r1.query_sequence = "A" * r1_len
        r1.query_qualities = pysam.qualitystring_to_array("!" * r1_len)
        set_paired_flags(r1, is_read1=True)
        r1.reference_id = ref_id
        r1.reference_start = r1_start
        r1.mapping_quality = 60
        r1.cigarstring = f"{r1_len}M"
        set_mate_info(r1, header, mate_ref_name=ref_name, mate_start=r2_start, mate_is_unmapped=False,
mate_is_reverse=True)
        r1.template_length = tlen
        r1.flag = 99 # Paired, Proper, Mate Reverse, R1
        samfile.write(r1)

        r2 = pysam.AlignedSegment()
        r2.query_name = "no_overlap_pair"
        r2.query_sequence = "T" * r2_len
        r2.query_qualities = pysam.qualitystring_to_array("#" * r2_len)
        set_paired_flags(r2, is_read1=False)
        r2.reference_id = ref_id
        r2.reference_start = r2_start
        r2.mapping_quality = 60
        r2.cigarstring = f"{r2_len}M"
        r2.is_reverse = True
        set_mate_info(r2, header, mate_ref_name=ref_name, mate_start=r1_start, mate_is_unmapped=False,
mate_is_reverse=False)
        r2.template_length = -tlen
        r2.flag = 147 # Paired, Proper, Read Reverse, R2
        samfile.write(r2)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_09(output_dir: Path, **kwargs):
    """SAM_09: Mapped read pair - long distance + TLEN (>1M bases)"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir, ref_name="large_ref.fa")
    
    with pysam.FastaFile(str(ref_path)) as fasta:
        references = fasta.references
        lengths = fasta.lengths
        ref_name = references[0]
    
    header = pysam.AlignmentHeader.from_references(references, lengths)
    ref_id = header.references.index(ref_name)
    ref_length = header.lengths[ref_id]
    
    r1_start = 5
    r2_start = ref_length - 100
    read_len = 100
    tlen = (r2_start + read_len) - r1_start

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Read 1 (leftmost)
        r1 = pysam.AlignedSegment(header)
        r1.query_name = "long_dist_pair"
        r1.query_sequence = "A" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("!" * read_len)
        r1.reference_id = ref_id
        r1.reference_start = r1_start
        r1.mapping_quality = 60
        r1.cigarstring = f"{read_len}M"
        set_mate_info(r1, header, ref_name, r2_start, False, True)
        r1.template_length = tlen
        r1.flag = 0x1 | 0x40 | 0x20

        # Read 2 (rightmost)
        r2 = pysam.AlignedSegment(header)
        r2.query_name = "long_dist_pair"
        r2.query_sequence = "T" * read_len
        r2.query_qualities = pysam.qualitystring_to_array("#" * read_len)
        r2.reference_id = ref_id
        r2.reference_start = r2_start
        r2.mapping_quality = 60
        r2.cigarstring = f"{read_len}M"
        r2.is_reverse = True
        set_mate_info(r2, header, ref_name, r1_start, False, False)
        r2.template_length = -tlen
        r2.flag = 0x1 | 0x80 | 0x10

        samfile.write(r1)
        samfile.write(r2)

def generate_sam_10(output_dir: Path, **kwargs):
    """SAM_10: Mapped read pair – different reference + TLEN"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref1_name = "ref1"
    ref1_id = header.references.index(ref1_name)
    ref2_name = "ref2"
    ref2_id = header.references.index(ref2_name)
    r1_start = 10
    r1_len = 25
    r2_start = 20
    r2_len = 25

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        r1 = pysam.AlignedSegment()
        r1.query_name = "diff_ref_pair"
        r1.query_sequence = "A" * r1_len
        r1.query_qualities = pysam.qualitystring_to_array("!" * r1_len)
        set_paired_flags(r1, is_read1=True)
        r1.reference_id = ref1_id # Maps to ref1
        r1.reference_start = r1_start
        r1.mapping_quality = 60
        r1.cigarstring = f"{r1_len}M"
        # Mate info: R2, mapped, ref2, r2_start, assume not reversed
        set_mate_info(r1, header, mate_ref_name=ref2_name, mate_start=r2_start, mate_is_unmapped=False,
mate_is_reverse=False)
        # TLEN is usually 0 when mates are on different chromosomes
        r1.template_length = 0
        # Flags: 0x1 (paired) + 0x40 (R1). Not proper pair (0x2).
        r1.flag = 65
        samfile.write(r1)

        r2 = pysam.AlignedSegment()
        r2.query_name = "diff_ref_pair"
        r2.query_sequence = "T" * r2_len
        r2.query_qualities = pysam.qualitystring_to_array("#" * r2_len)
        set_paired_flags(r2, is_read1=False)
        r2.reference_id = ref2_id # Maps to ref2
        r2.reference_start = r2_start
        r2.mapping_quality = 60
        r2.cigarstring = f"{r2_len}M"
        # Mate info: R1, mapped, ref1, r1_start, assume not reversed
        set_mate_info(r2, header, mate_ref_name=ref1_name, mate_start=r1_start, mate_is_unmapped=False,
mate_is_reverse=False)
        r2.template_length = 0
        # Flags: 0x1 (paired) + 0x80 (R2). Not proper pair (0x2).
        r2.flag = 129
        samfile.write(r2)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_11(output_dir: Path, **kwargs):
    """SAM_11: Secondary alignment"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    read_len = 30
    seq = "C" * read_len
    qual = "!" * read_len

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Primary alignment
        p = pysam.AlignedSegment()
        p.query_name = "secondary_test_read"
        p.query_sequence = seq
        p.query_qualities = pysam.qualitystring_to_array(qual)
        p.reference_id = ref_id
        p.reference_start = 10
        p.mapping_quality = 60
        p.cigarstring = f"{read_len}M"
        # Flags: 0 (basic mapped SE read, primary)
        p.flag = 0
        samfile.write(p)

        # Secondary alignment 1
        s1 = pysam.AlignedSegment()
        s1.query_name = "secondary_test_read" # Same name
        s1.query_sequence = seq
        s1.query_qualities = pysam.qualitystring_to_array(qual)
        s1.reference_id = ref_id
        s1.reference_start = 50 # Different position
        s1.mapping_quality = 20 # Lower MAPQ typical
        s1.cigarstring = f"{read_len}M"
        # Flags: 0x100 (secondary)
        s1.flag = 256
        samfile.write(s1)

        # Secondary alignment 2 (different CIGAR)
        s2 = pysam.AlignedSegment()
        s2.query_name = "secondary_test_read"
        s2.query_sequence = seq
        s2.query_qualities = pysam.qualitystring_to_array(qual)
        s2.reference_id = ref_id
        s2.reference_start = 90
        s2.mapping_quality = 10
        s2.cigarstring = f"5S{read_len-5}M" # Add softclip
        # Flags: 0x100 (secondary)
        s2.flag = 256
        samfile.write(s2)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_12(output_dir: Path, **kwargs):
    """SAM_12: Supplementary / chimeric alignment"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref1_name = "ref1"
    ref1_id = header.references.index(ref1_name)
    ref2_name = "ref2"
    ref2_id = header.references.index(ref2_name)

    # Simulate a read split across two references
    read_name = "chimeric_read"
    seq_part1 = "A" * 20
    seq_part2 = "T" * 25
    full_seq = seq_part1 + seq_part2
    qual = "!" * len(full_seq)

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Alignment part 1 (representative, primary)
        p1 = pysam.AlignedSegment()
        p1.query_name = read_name
        p1.query_sequence = full_seq # Full sequence usually on primary non-suppl.
        p1.query_qualities = pysam.qualitystring_to_array(qual)
        p1.reference_id = ref1_id
        p1.reference_start = 10
        p1.mapping_quality = 60
        # CIGAR shows part1 matches, part2 softclipped
        p1.cigarstring = f"{len(seq_part1)}M{len(seq_part2)}S"
        # Flags: 0 (basic mapped SE, primary, non-supplementary)
        p1.flag = 0
        # SA tag points to supplementary alignment
        p1.set_tag('SA', f"{ref2_name},{50},+,{len(seq_part1)}S{len(seq_part2)}M,{60},0;", 'Z')
        samfile.write(p1)

        # Alignment part 2 (supplementary)
        s1 = pysam.AlignedSegment()
        s1.query_name = read_name
        s1.query_sequence = full_seq # Often full seq here too, or sometimes '*'
        s1.query_qualities = pysam.qualitystring_to_array(qual) # Or '*'
        s1.reference_id = ref2_id # Different reference
        s1.reference_start = 50
        s1.mapping_quality = 60
        # CIGAR shows part1 softclipped, part2 matches
        s1.cigarstring = f"{len(seq_part1)}S{len(seq_part2)}M"
        # Flags: 0x800 (supplementary)
        s1.flag = 2048
        samfile.write(s1)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_13(output_dir: Path, **kwargs):
    """SAM_13: Base substitution (M, =, X)"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    ref_path = utils.copy_reference_to_output(output_dir)
    with pysam.FastaFile(str(ref_path)) as fasta:
        ref_seq = fasta.fetch(ref_name, 10, 40).upper() # ACGTACGTACGTACGTACGTACGTACGTAC

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Read with M, =, X
        r1 = pysam.AlignedSegment()
        r1.query_name = "subst_read_M_eq_X"
        # Create sequence: match first 5 (=), mismatch next 5 (X), match next 5 (=), mismatch last 5 (X)
        # Ref: ACGTA CGTAC GTACG TACGT ACGTA CGTAC
        # Seq: ACGTA TTTTT GTACG AAAAA ACGTA TTTTT
        seq = ref_seq[0:5] + "TTTTT" + ref_seq[10:15] + "AAAAA" + ref_seq[20:25] + "TTTTT"
        r1.query_sequence = seq
        r1.query_qualities = pysam.qualitystring_to_array("!" * len(seq))
        r1.reference_id = ref_id
        r1.reference_start = 10
        r1.mapping_quality = 60
        # CIGAR: 5=5X5=5X5=5X (or use M: 5=5M5=5M5=5M)
        # Using =/X requires MD/NM tags usually. Let's use M for simplicity here.
        # If strict =/X is needed, CIGAR would be e.g. 5=5X5=5X5=5X
        r1.cigarstring = "5=5X5=5X5=5X" # Explicit matches and mismatches
        r1.flag = 0
        # Add NM tag (Number of Mismatches): 5+5+5 = 15
        r1.set_tag('NM', 10, 'i')
        # Add MD tag (String for mismatching positions): 5^CGTAC5^ACGTA5^CGTAC5
        r1.set_tag('MD', '5^TTTTT5^AAAAA5^TTTTT5', 'Z')
        samfile.write(r1)

        # Read with just M (can represent match or mismatch)
        r2 = pysam.AlignedSegment()
        r2.query_name = "subst_read_M"
        # Seq: ACGTA TTTTT GTACG AAAAA ACGTA TTTTT (same as above)
        r2.query_sequence = seq
        r2.query_qualities = pysam.qualitystring_to_array("#" * len(seq))
        r2.reference_id = ref_id
        r2.reference_start = 50 # Different location
        r2.mapping_quality = 60
        r2.cigarstring = f"{len(seq)}M" # Use M for the whole length
        r2.flag = 0
        r2.set_tag('NM', 10, 'i') # Still need NM/MD if using M for mismatches
        r2.set_tag('MD', '5^TTTTT5^AAAAA5^TTTTT5', 'Z')
        samfile.write(r2)


    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_14(output_dir: Path, **kwargs):
    """SAM_14: Base insertion"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    with pysam.FastaFile(str(ref_path)) as fasta:
        ref_seq_part1 = fasta.fetch(ref_name, 10, 20).upper() # 10 bases
        ref_seq_part2 = fasta.fetch(ref_name, 20, 30).upper() # 10 bases
    insertion = "NNN" # 3 bases inserted

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        r1 = pysam.AlignedSegment()
        r1.query_name = "insertion_read"
        # Sequence contains the insertion
        r1.query_sequence = ref_seq_part1 + insertion + ref_seq_part2
        r1.query_qualities = pysam.qualitystring_to_array("!" * len(r1.query_sequence))
        r1.reference_id = ref_id
        r1.reference_start = 10
        r1.mapping_quality = 60
        # CIGAR: 10 matches, 3 insertion, 10 matches
        r1.cigarstring = f"{len(ref_seq_part1)}M{len(insertion)}I{len(ref_seq_part2)}M"
        r1.flag = 0
        samfile.write(r1)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_15(output_dir: Path, **kwargs):
    """SAM_15: Base deletion"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    with pysam.FastaFile(str(ref_path)) as fasta:
        ref_seq_part1 = fasta.fetch(ref_name, 10, 20).upper() # 10 bases
        # Deletion of 3 bases from ref: ref pos 20, 21, 22
        ref_seq_part2 = fasta.fetch(ref_name, 23, 33).upper() # 10 bases

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        r1 = pysam.AlignedSegment()
        r1.query_name = "deletion_read"
        # Sequence does NOT contain the deleted bases
        r1.query_sequence = ref_seq_part1 + ref_seq_part2
        r1.query_qualities = pysam.qualitystring_to_array("!" * len(r1.query_sequence))
        r1.reference_id = ref_id
        r1.reference_start = 10
        r1.mapping_quality = 60
        # CIGAR: 10 matches, 3 deletion, 10 matches
        r1.cigarstring = f"{len(ref_seq_part1)}M3D{len(ref_seq_part2)}M"
        r1.flag = 0
        samfile.write(r1)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_16(output_dir: Path, **kwargs):
    """SAM_16: Softclips"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    with pysam.FastaFile(str(ref_path)) as fasta:
        ref_match_seq = fasta.fetch(ref_name, 20, 40).upper() # 20 bases match
    clip5_seq = "NNNNN" # 5 bases clipped at start
    clip3_seq = "YYYYY" # 5 bases clipped at end

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Read 1: 5' softclip
        r1 = pysam.AlignedSegment()
        r1.query_name = "softclip_5prime"
        # Sequence includes clipped bases
        r1.query_sequence = clip5_seq + ref_match_seq
        r1.query_qualities = pysam.qualitystring_to_array("!" * len(r1.query_sequence))
        r1.reference_id = ref_id
        r1.reference_start = 20 # Maps starting at pos 20
        r1.mapping_quality = 60
        r1.cigarstring = f"{len(clip5_seq)}S{len(ref_match_seq)}M"
        r1.flag = 0
        samfile.write(r1)

        # Read 2: 3' softclip
        r2 = pysam.AlignedSegment()
        r2.query_name = "softclip_3prime"
        r2.query_sequence = ref_match_seq + clip3_seq
        r2.query_qualities = pysam.qualitystring_to_array("#" * len(r2.query_sequence))
        r2.reference_id = ref_id
        r2.reference_start = 20
        r2.mapping_quality = 60
        r2.cigarstring = f"{len(ref_match_seq)}M{len(clip3_seq)}S"
        r2.flag = 0
        samfile.write(r2)

        # Read 3: Both 5' and 3' softclip
        r3 = pysam.AlignedSegment()
        r3.query_name = "softclip_both"
        r3.query_sequence = clip5_seq + ref_match_seq + clip3_seq
        r3.query_qualities = pysam.qualitystring_to_array("$" * len(r3.query_sequence))
        r3.reference_id = ref_id
        r3.reference_start = 20
        r3.mapping_quality = 60
        r3.cigarstring = f"{len(clip5_seq)}S{len(ref_match_seq)}M{len(clip3_seq)}S"
        r3.flag = 0
        samfile.write(r3)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_17(output_dir: Path, **kwargs):
    """SAM_17: Padding (P)"""
    # Padding is rare in mapping, more common in multiple alignments.
    # It consumes neither read nor reference.
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    with pysam.FastaFile(str(ref_path)) as fasta:
        ref_seq1 = fasta.fetch(ref_name, 10, 20).upper() # 10 bases
        ref_seq2 = fasta.fetch(ref_name, 20, 30).upper() # 10 bases

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        r1 = pysam.AlignedSegment()
        r1.query_name = "padding_read"
        # Sequence length only includes M/I/S/=/X ops
        r1.query_sequence = ref_seq1 + ref_seq2
        r1.query_qualities = pysam.qualitystring_to_array("!" * len(r1.query_sequence))
        r1.reference_id = ref_id
        r1.reference_start = 10
        r1.mapping_quality = 60
        # CIGAR: 10 matches, 2 padding, 10 matches
        r1.cigarstring = f"{len(ref_seq1)}M2P{len(ref_seq2)}M"
        r1.flag = 0
        samfile.write(r1)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_18(output_dir: Path, **kwargs):
    """SAM_18: Hardclips"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    with pysam.FastaFile(str(ref_path)) as fasta:
        ref_match_seq = fasta.fetch(ref_name, 20, 40).upper() # 20 bases match
    clip5_len = 5 # 5 bases hardclipped at start
    clip3_len = 5 # 5 bases hardclipped at end

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Read 1: 5' hardclip
        r1 = pysam.AlignedSegment()
        r1.query_name = "hardclip_5prime"
        # Sequence does NOT include clipped bases
        r1.query_sequence = ref_match_seq
        r1.query_qualities = pysam.qualitystring_to_array("!" * len(r1.query_sequence))
        r1.reference_id = ref_id
        r1.reference_start = 20 # Maps starting at pos 20
        r1.mapping_quality = 60
        r1.cigarstring = f"{clip5_len}H{len(ref_match_seq)}M"
        r1.flag = 0
        samfile.write(r1)

        # Read 2: 3' hardclip
        r2 = pysam.AlignedSegment()
        r2.query_name = "hardclip_3prime"
        r2.query_sequence = ref_match_seq
        r2.query_qualities = pysam.qualitystring_to_array("#" * len(r2.query_sequence))
        r2.reference_id = ref_id
        r2.reference_start = 20
        r2.mapping_quality = 60
        r2.cigarstring = f"{len(ref_match_seq)}M{clip3_len}H"
        r2.flag = 0
        samfile.write(r2)

        # Read 3: Both 5' and 3' hardclip
        r3 = pysam.AlignedSegment()
        r3.query_name = "hardclip_both"
        r3.query_sequence = ref_match_seq
        r3.query_qualities = pysam.qualitystring_to_array("$" * len(r3.query_sequence))
        r3.reference_id = ref_id
        r3.reference_start = 20
        r3.mapping_quality = 60
        r3.cigarstring = f"{clip5_len}H{len(ref_match_seq)}M{clip3_len}H"
        r3.flag = 0
        samfile.write(r3)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_19(output_dir: Path, **kwargs):
    """SAM_19: PCR duplicate flag"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    read_len = 10

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Read 1: Not duplicate
        r1 = pysam.AlignedSegment()
        r1.query_name = "pcr_dup_1_not_dup"
        r1.query_sequence = "A" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("!" * read_len)
        r1.reference_id = ref_id
        r1.reference_start = 10
        r1.mapping_quality = 60
        r1.cigarstring = f"{read_len}M"
        r1.flag = 0 # Not duplicate (0x400 is clear)
        samfile.write(r1)

        # Read 2: Is duplicate
        r2 = pysam.AlignedSegment()
        r2.query_name = "pcr_dup_2_is_dup"
        r2.query_sequence = "C" * read_len
        r2.query_qualities = pysam.qualitystring_to_array("#" * read_len)
        r2.reference_id = ref_id
        r2.reference_start = 30
        r2.mapping_quality = 60
        r2.cigarstring = f"{read_len}M"
        r2.flag = 0x400 # Set duplicate flag
        # r2.is_duplicate = True # Alternative
        samfile.write(r2)

        # Read 3: Paired, R1 is duplicate, R2 is not (example)
        r3_1 = pysam.AlignedSegment()
        r3_1.query_name = "pcr_dup_3_mixed_pair"
        r3_1.query_sequence = "G" * read_len
        r3_1.query_qualities = pysam.qualitystring_to_array("$" * read_len)
        set_paired_flags(r3_1, is_read1=True)
        r3_1.reference_id = ref_id
        r3_1.reference_start = 50
        r3_1.mapping_quality = 60
        r3_1.cigarstring = f"{read_len}M"
        set_mate_info(r3_1, header, mate_ref_name=ref_name, mate_start=80, mate_is_unmapped=False, mate_is_reverse=True)
        r3_1.template_length = 40 # Approx
        # Flags: Paired, Proper, Mate Reverse, R1, Duplicate
        r3_1.flag = 0x1 + 0x2 + 0x20 + 0x40 + 0x400 # = 1123
        samfile.write(r3_1)

        r3_2 = pysam.AlignedSegment()
        r3_2.query_name = "pcr_dup_3_mixed_pair"
        r3_2.query_sequence = "T" * read_len
        r3_2.query_qualities = pysam.qualitystring_to_array("%" * read_len)
        set_paired_flags(r3_2, is_read1=False)
        r3_2.reference_id = ref_id
        r3_2.reference_start = 80
        r3_2.mapping_quality = 60
        r3_2.cigarstring = f"{read_len}M"
        r3_2.is_reverse = True
        set_mate_info(r3_2, header, mate_ref_name=ref_name, mate_start=50, mate_is_unmapped=False,
mate_is_reverse=False)
        r3_2.template_length = -40
        # Flags: Paired, Proper, Read Reverse, R2 (NO duplicate)
        r3_2.flag = 0x1 + 0x2 + 0x10 + 0x80 # = 147
        samfile.write(r3_2)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_20(output_dir: Path, **kwargs):
    """SAM_20: Paired end – different flags per mate"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    read_len = 15

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Pair 1: R1 normal, R2 QC Fail
        r1_start, r2_start = 10, 40
        tlen = (r2_start + read_len -1) - r1_start + 1
        r1 = pysam.AlignedSegment()
        r1.query_name = "diff_flags_pair1"
        r1.query_sequence = "A" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("!" * read_len)
        set_paired_flags(r1, is_read1=True)
        r1.reference_id = ref_id
        r1.reference_start = r1_start
        r1.mapping_quality = 60
        r1.cigarstring = f"{read_len}M"
        set_mate_info(r1, header, mate_ref_name=ref_name, mate_start=r2_start, mate_is_unmapped=False,
mate_is_reverse=True)
        r1.template_length = tlen
        # Flags: Paired, Proper, Mate Reverse, R1 (Normal: 99)
        r1.flag = 99
        samfile.write(r1)

        r2 = pysam.AlignedSegment()
        r2.query_name = "diff_flags_pair1"
        r2.query_sequence = "T" * read_len
        r2.query_qualities = pysam.qualitystring_to_array("#" * read_len)
        set_paired_flags(r2, is_read1=False)
        r2.reference_id = ref_id
        r2.reference_start = r2_start
        r2.mapping_quality = 60
        r2.cigarstring = f"{read_len}M"
        r2.is_reverse = True
        set_mate_info(r2, header, mate_ref_name=ref_name, mate_start=r1_start, mate_is_unmapped=False,
mate_is_reverse=False)
        r2.template_length = -tlen
        # Flags: Paired, Proper, Read Reverse, R2, QC Fail (0x200)
        r2.flag = 147 + 0x200 # = 659
        # r2.is_qcfail = True # Alternative
        samfile.write(r2)

        # Pair 2: R1 Duplicate, R2 normal
        r1_start, r2_start = 70, 100
        tlen = (r2_start + read_len -1) - r1_start + 1
        r1 = pysam.AlignedSegment()
        r1.query_name = "diff_flags_pair2"
        r1.query_sequence = "C" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("$" * read_len)
        set_paired_flags(r1, is_read1=True)
        r1.reference_id = ref_id
        r1.reference_start = r1_start
        r1.mapping_quality = 60
        r1.cigarstring = f"{read_len}M"
        set_mate_info(r1, header, mate_ref_name=ref_name, mate_start=r2_start, mate_is_unmapped=False,
mate_is_reverse=True)
        r1.template_length = tlen
        # Flags: Paired, Proper, Mate Reverse, R1, Duplicate (0x400)
        r1.flag = 99 + 0x400 # = 1123
        # r1.is_duplicate = True # Alternative
        samfile.write(r1)

        r2 = pysam.AlignedSegment()
        r2.query_name = "diff_flags_pair2"
        r2.query_sequence = "G" * read_len
        r2.query_qualities = pysam.qualitystring_to_array("%" * read_len)
        set_paired_flags(r2, is_read1=False)
        r2.reference_id = ref_id
        r2.reference_start = r2_start
        r2.mapping_quality = 60
        r2.cigarstring = f"{read_len}M"
        r2.is_reverse = True
        set_mate_info(r2, header, mate_ref_name=ref_name, mate_start=r1_start, mate_is_unmapped=False,
mate_is_reverse=False)
        r2.template_length = -tlen
        # Flags: Paired, Proper, Read Reverse, R2 (Normal: 147)
        r2.flag = 147
        samfile.write(r2)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_21(output_dir: Path, **kwargs):
    """SAM_21: next read … flags – unmapped"""
    # Assumption: Show an unmapped read where mate flags are set (e.g., mate unmapped)
    file_path = output_dir / "alignment.sam"
    header = utils.get_default_sam_header()

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Unmapped SE read, but marked as paired, mate also unmapped
        r1 = pysam.AlignedSegment()
        r1.query_name = "unmapped_with_mate_flags"
        r1.query_sequence = "A" * 10
        r1.query_qualities = pysam.qualitystring_to_array("!" * 10)
        r1.flag = 0x1 + 0x4 + 0x8 # Paired, Unmapped, Mate Unmapped
        # Mate position often set to self if unmapped
        r1.next_reference_id = -1  # *
        r1.next_reference_start = 0
        r1.mate_is_unmapped = True
        r1.mate_is_reverse = False
        samfile.write(r1)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Note: Assumed test requires unmapped read with mate flags set.")

def generate_sam_22(output_dir: Path, **kwargs):
    """SAM_22: next read … flags – half mapped"""
    # Assumption: Show a mapped read where mate flags indicate mate is unmapped.
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Mapped read, but mate is unmapped
        r1 = pysam.AlignedSegment()
        r1.query_name = "mapped_mate_unmapped"
        r1.query_sequence = "C" * 12
        r1.query_qualities = pysam.qualitystring_to_array("!" * 12)
        r1.reference_id = ref_id
        r1.reference_start = 20
        r1.mapping_quality = 60
        r1.cigarstring = "12M"
        r1.is_paired = True
        r1.mate_is_unmapped = True # Key flag
        # Mate position often set to self if mate unmapped
        set_mate_info(r1, header, mate_ref_name=ref_name, mate_start=20, mate_is_unmapped=True, mate_is_reverse=False)
        # Flags: 0x1 (paired) + 0x8 (mate unmapped) = 9
        r1.flag = 9
        samfile.write(r1)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")
    # tqdm.write(f"{Fore.GREEN}  Note: Assumed test requires mapped read with mate_unmapped flag set.")

def generate_sam_23(output_dir: Path, **kwargs):
    """SAM_23: next read … flags – short distance"""
    # Assumption: Show a mapped read where mate flags indicate mate is mapped nearby.
    # This is essentially a standard proper pair, focusing on the mate flags.
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    r1_start, r2_start = 30, 70
    read_len = 20
    tlen = (r2_start + read_len - 1) - r1_start + 1

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # R1 of a proper pair
        r1 = pysam.AlignedSegment()
        r1.query_name = "mapped_mate_nearby"
        r1.query_sequence = "G" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("!" * read_len)
        r1.reference_id = ref_id
        r1.reference_start = r1_start
        r1.mapping_quality = 60
        r1.cigarstring = f"{read_len}M"
        r1.is_paired = True
        r1.is_proper_pair = True
        r1.is_read1 = True
        r1.mate_is_reverse = True # Mate flags set correctly
        r1.next_reference_id = ref_id
        r1.next_reference_start = r2_start
        r1.template_length = tlen
        # Flags: 0x1 + 0x2 + 0x20 + 0x40 = 99
        r1.flag = 99
        samfile.write(r1)
        # (Mate R2 not strictly needed for this test, but good practice)
        r2 = pysam.AlignedSegment()
        r2.query_name = "mapped_mate_nearby"
        r2.query_sequence = "T" * read_len
        r2.query_qualities = pysam.qualitystring_to_array("#" * read_len)
        r2.reference_id = ref_id
        r2.reference_start = r2_start
        r2.mapping_quality = 60
        r2.cigarstring = f"{read_len}M"
        r2.is_paired = True
        r2.is_proper_pair = True
        r2.is_read2 = True
        r2.is_reverse = True
        r2.mate_is_reverse = False # Mate flags set correctly
        r2.next_reference_id = ref_id
        r2.next_reference_start = r1_start
        r2.template_length = -tlen
        # Flags: 0x1 + 0x2 + 0x10 + 0x80 = 147
        r2.flag = 147
        samfile.write(r2)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")
    # tqdm.write(f"{Fore.GREEN}  Note: Assumed test requires mapped read with mate flags indicating nearby mate.")


def generate_sam_24(output_dir: Path, **kwargs):
    """SAM_24: next read … flags – long distance"""
    # Assumption: Show a mapped read where mate flags indicate mate is mapped far away.
    # Similar to SAM_09, but focus is on the flags of one read.
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    r1_start = 5
    r1_len = 20
    r2_start = 130 # Far away on short ref
    r2_len = 20
    tlen = (r2_start + r2_len - 1) - r1_start + 1 # Large TLEN

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # R1 of a long-distance pair
        r1 = pysam.AlignedSegment()
        r1.query_name = "mapped_mate_far"
        r1.query_sequence = "N" * r1_len
        r1.query_qualities = pysam.qualitystring_to_array("!" * r1_len)
        r1.reference_id = ref_id
        r1.reference_start = r1_start
        r1.mapping_quality = 60
        r1.cigarstring = f"{r1_len}M"
        r1.is_paired = True
        # May not be 'proper pair' depending on definition
        # r1.is_proper_pair = False
        r1.is_read1 = True
        r1.mate_is_reverse = True # Mate flags set correctly
        r1.next_reference_id = ref_id
        r1.next_reference_start = r2_start # Far position
        r1.template_length = tlen # Large TLEN
        # Flags: 0x1 + 0x20 + 0x40 = 97 (Paired, Mate Reverse, R1 - NOT proper)
        r1.flag = 97
        samfile.write(r1)
        # (Mate R2 not strictly needed for this test)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")
    # tqdm.write(f"{Fore.GREEN}  Note: Assumed test requires mapped read with mate flags indicating far away mate.")

def generate_sam_25(output_dir: Path, **kwargs):
    """SAM_25: Short intron / splice (N)"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    exon1_len = 25
    intron_len = 50 # Short intron
    exon2_len = 25
    start_pos = 10
    ref_path = utils.copy_reference_to_output(output_dir)
    with pysam.FastaFile(str(ref_path)) as fasta:
        ref_exon1 = fasta.fetch(ref_name, start_pos, start_pos + exon1_len).upper()
        ref_exon2 = fasta.fetch(ref_name, start_pos + exon1_len + intron_len, 
                              start_pos + exon1_len + intron_len + exon2_len).upper()

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        r1 = pysam.AlignedSegment()
        r1.query_name = "short_splice_read"
        # Sequence contains only exons
        r1.query_sequence = ref_exon1 + ref_exon2
        r1.query_qualities = pysam.qualitystring_to_array("!" * len(r1.query_sequence))
        r1.reference_id = ref_id
        r1.reference_start = start_pos
        r1.mapping_quality = 60
        # CIGAR: Exon1 match, Intron skip, Exon2 match
        r1.cigarstring = f"{exon1_len}M{intron_len}N{exon2_len}M"
        r1.set_tag('SA', f"ref1:{start_pos + exon1_len + intron_len}+{exon2_len}M,60,0;", 'Z')
        r1.flag = 0
        samfile.write(r1)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_26(output_dir: Path, **kwargs):
    """SAM_26: Long intron / splice (N)"""
    # Note: Our reference is short. We simulate a long intron skip relative
    # to the reference size, but it won't be biologically 'long'.
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    exon1_len = 15
    intron_len = 100 # Long relative to ref1 (len 160)
    exon2_len = 15
    start_pos = 5
    with pysam.FastaFile(str(ref_path)) as fasta:
        ref_exon1 = fasta.fetch(ref_name, 5, 5 + 15).upper()
        ref_exon2_start = 5 + 15 + 100
        ref_exon2 = fasta.fetch(ref_name, ref_exon2_start, ref_exon2_start +
                                15).upper()

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        r1 = pysam.AlignedSegment()
        r1.query_name = "long_splice_read"
        r1.query_sequence = ref_exon1 + ref_exon2
        r1.query_qualities = pysam.qualitystring_to_array("!" * len(r1.query_sequence))
        r1.reference_id = ref_id
        r1.reference_start = start_pos
        r1.mapping_quality = 60
        r1.cigarstring = f"{exon1_len}M{intron_len}N{exon2_len}M"
        r1.set_tag('SA', f"ref1,85,+,25M,60,0;", 'Z')
        r1.flag = 0
        samfile.write(r1)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_27(output_dir: Path, **kwargs):
    """SAM_27: Empty read (all bases deleted)"""
    # This is unusual. Typically implies an alignment placeholder.
    # Usually unmapped.
    file_path = output_dir / "alignment.sam"
    header = utils.get_default_sam_header()
    deletion_len = 50

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        r1 = pysam.AlignedSegment()
        r1.query_name = "empty_deleted_read"
        r1.query_sequence = "*" # No sequence
        r1.query_qualities = None # No qualities represented by '*' in SAM
        r1.flag = 4 # Unmapped
        # CIGAR might be tricky. If it represents a deletion relative to a
        # reference alignment attempt, it might be e.g. 50D.
        # If completely unmapped placeholder, CIGAR might be '*'
        # Let's assume it represents a failed alignment attempt with deletion
        r1.reference_id = header.references.index("ref1") # Nominal ref
        r1.reference_start = 10 # Nominal start
        r1.mapping_quality = 0
        r1.cigarstring = f"{deletion_len}D" # All deleted CIGAR
        samfile.write(r1)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")

def generate_sam_28(output_dir: Path, **kwargs):
    """SAM_28: Empty read (all bases softclipped)"""
    # Read sequence exists but doesn't map anywhere. Usually unmapped.
    file_path = output_dir / "alignment.sam"
    header = utils.get_default_sam_header()
    read_len = 50

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        r1 = pysam.AlignedSegment()
        r1.query_name = "empty_softclipped_read"
        r1.query_sequence = "N" * read_len # Sequence must be present
        r1.query_qualities = pysam.qualitystring_to_array("!" * read_len) # Qualities must be present
        r1.reference_id = 0  # Must have valid reference for mapped read
        r1.reference_start = 0
        r1.mapping_quality = 60
        r1.cigarstring = f"{read_len}S"
        r1.flag = 0  # Clear unmapped flag (0x4)
        samfile.write(r1)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")

def generate_sam_29(output_dir: Path, **kwargs):
    """SAM_29: Empty read (all bases hardclipped)"""
    # Similar to all deleted, sequence is gone. Usually unmapped.
    file_path = output_dir / "alignment.sam"
    header = utils.get_default_sam_header()
    hardclip_len = 50

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        r1 = pysam.AlignedSegment()
        r1.query_name = "empty_hardclipped_read"
        r1.query_sequence = "*" # No sequence
        r1.query_qualities = None # No qualities
        r1.flag = 4 # Unmapped
        r1.reference_id = -1
        r1.reference_start = -1
        r1.mapping_quality = 0
        r1.cigarstring = f"{hardclip_len}H" # All hardclipped
        samfile.write(r1)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")

def generate_sam_30(output_dir: Path, **kwargs):
    """SAM_30: Empty read (no nucleotides in read / * in sam)"""
    # Explicitly setting SEQ and QUAL to '*'. Usually unmapped.
    file_path = output_dir / "alignment.sam"
    header = utils.get_default_sam_header()

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        r1 = pysam.AlignedSegment()
        r1.query_name = "empty_star_seq_qual_read"
        r1.query_sequence = "*" # Explicitly '*'
        r1.query_qualities = None # Explicitly '*'
        r1.flag = 4 # Unmapped
        r1.reference_id = -1
        r1.reference_start = -1
        r1.mapping_quality = 0
        r1.cigarstring = "*" # CIGAR is also often '*' here
        samfile.write(r1)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")

def generate_sam_31(output_dir: Path, **kwargs):
    """SAM_31: Quality scores absent"""
    # Sequence present, but QUAL is '*'.
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    read_len = 20

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        r1 = pysam.AlignedSegment()
        r1.query_name = "no_qual_read"
        r1.query_sequence = "A" * read_len # Sequence must be present
        r1.query_qualities = None # Set QUAL to '*'
        r1.reference_id = ref_id
        r1.reference_start = 10
        r1.mapping_quality = 60
        r1.cigarstring = f"{read_len}M"
        r1.flag = 0
        samfile.write(r1)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_32(output_dir: Path, **kwargs):
    """SAM_32: Optional tags"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    read_len = 15

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Read 1: Various tags
        r1 = pysam.AlignedSegment()
        r1.query_name = "optional_tags_read"
        # Create a sequence with a mismatch
        with pysam.FastaFile(str(ref_path)) as fasta:
            ref_seq = fasta.fetch(ref_name, 10, 10 + 15).upper()
        seq = ref_seq[:5] + "N" + ref_seq[6:] # Mismatch at pos 5 (index)
        r1.query_sequence = seq
        r1.query_qualities = pysam.qualitystring_to_array("!" * read_len)
        r1.reference_id = ref_id
        r1.reference_start = 10
        r1.mapping_quality = 60
        r1.cigarstring = f"{read_len}M" # Use M
        r1.flag = 0
        # Add tags
        r1.set_tag('NM', 1, 'i') # Edit distance
        r1.set_tag('MD', '5A9', 'Z') # Mismatch string (pos 5 was A)
        r1.set_tag('AS', 50, 'i') # Alignment score
        r1.set_tag('XS', 40, 'i') # Secondary alignment score (example)
        r1.set_tag('RG', 'readgroup1', 'Z') # Read group
        r1.set_tag('BC', 'ACGT', 'Z') # Barcode sequence
        r1.set_tag('MI', 'molecule1', 'Z') # Molecular identifier
        r1.set_tag('ZZ', 'custom_tag_value', 'Z') # Custom tag
        samfile.write(r1)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_33(output_dir: Path, **kwargs):
    """SAM_33: Read groups"""
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    # Get base header and add RG lines
    header_dict = utils.get_default_sam_header().to_dict()
    header_dict['RG'] = [
        {'ID': 'rg1', 'SM': 'sample1', 'LB': 'lib1', 'PL': 'ILLUMINA', 'PU': 'unit1'},
        {'ID': 'rg2', 'SM': 'sample2', 'LB': 'lib2', 'PL': 'ONT', 'PU': 'unit2'},
    ]
    header = pysam.AlignmentHeader.from_dict(header_dict)

    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    read_len = 10

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Read 1 in rg1
        r1 = pysam.AlignedSegment()
        r1.query_name = "read_group_1"
        r1.query_sequence = "A" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("!" * read_len)
        r1.reference_id = ref_id
        r1.reference_start = 10
        r1.mapping_quality = 60
        r1.cigarstring = f"{read_len}M"
        r1.flag = 0
        r1.set_tag('RG', 'rg1', 'Z')
        samfile.write(r1)

        # Read 2 in rg2
        r2 = pysam.AlignedSegment()
        r2.query_name = "read_group_2"
        r2.query_sequence = "C" * read_len
        r2.query_qualities = pysam.qualitystring_to_array("#" * read_len)
        r2.reference_id = ref_id
        r2.reference_start = 30
        r2.mapping_quality = 60
        r2.cigarstring = f"{read_len}M"
        r2.flag = 0
        r2.set_tag('RG', 'rg2', 'Z')
        samfile.write(r2)

        # Read 3 in rg1 again
        r3 = pysam.AlignedSegment()
        r3.query_name = "read_group_3"
        r3.query_sequence = "G" * read_len
        r3.query_qualities = pysam.qualitystring_to_array("$" * read_len)
        r3.reference_id = ref_id
        r3.reference_start = 50
        r3.mapping_quality = 60
        r3.cigarstring = f"{read_len}M"
        r3.flag = 0
        r3.set_tag('RG', 'rg1', 'Z')
        samfile.write(r3)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_34(output_dir: Path, **kwargs):
    """SAM_34: Reverse Complement (different + same) – short distance"""
    # Assumption: Show pairs where one or both reads are reverse complemented (0x10/0x20 flags)
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    read_len = 20

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Pair 1: Standard FR (R1 Fwd, R2 Rev) - Already covered (e.g. SAM_08)
        r1_start, r2_start = 10, 50
        tlen = (r2_start + read_len - 1) - r1_start + 1
        r1 = pysam.AlignedSegment()
        r1.query_name = "revcomp_pair_FR"
        r1.query_sequence = "A" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("!" * read_len)
        set_paired_flags(r1, is_read1=True)
        r1.reference_id = ref_id; r1.reference_start = r1_start
        r1.mapping_quality = 60; r1.cigarstring = f"{read_len}M"
        set_mate_info(r1, header, ref_name, r2_start, False, True) # Mate is Rev
        r1.template_length = tlen; r1.flag = 99 # Proper, Mate Rev, R1
        samfile.write(r1)
        r2 = pysam.AlignedSegment()
        r2.query_name = "revcomp_pair_FR"
        r2.query_sequence = "T" * read_len
        r2.query_qualities = pysam.qualitystring_to_array("#" * read_len)
        set_paired_flags(r2, is_read1=False); r2.is_reverse = True # Read is Rev
        r2.reference_id = ref_id; r2.reference_start = r2_start
        r2.mapping_quality = 60; r2.cigarstring = f"{read_len}M"
        set_mate_info(r2, header, ref_name, r1_start, False, False) # Mate is Fwd
        r2.template_length = -tlen; r2.flag = 147 # Proper, Read Rev, R2
        samfile.write(r2)

        # Pair 2: RF (R1 Rev, R2 Fwd)
        r1_start, r2_start = 80, 120
        tlen = (r2_start + read_len - 1) - r1_start + 1 # R2 is now leftmost
        r1 = pysam.AlignedSegment()
        r1.query_name = "revcomp_pair_RF"
        r1.query_sequence = "C" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("$" * read_len)
        set_paired_flags(r1, is_read1=True); r1.is_reverse = True # Read is Rev
        r1.reference_id = ref_id; r1.reference_start = r1_start
        r1.mapping_quality = 60; r1.cigarstring = f"{read_len}M"
        set_mate_info(r1, header, ref_name, r2_start, False, False) # Mate is Fwd
        r1.template_length = tlen; r1.flag = 83 # Proper, Read Rev, R1
        samfile.write(r1)
        r2 = pysam.AlignedSegment()
        r2.query_name = "revcomp_pair_RF"
        r2.query_sequence = "G" * read_len
        r2.query_qualities = pysam.qualitystring_to_array("%" * read_len)
        set_paired_flags(r2, is_read1=False) # Read is Fwd
        r2.reference_id = ref_id; r2.reference_start = r2_start
        r2.mapping_quality = 60; r2.cigarstring = f"{read_len}M"
        set_mate_info(r2, header, ref_name, r1_start, False, True) # Mate is Rev
        r2.template_length = -tlen; r2.flag = 163 # Proper, Mate Rev, R2
        samfile.write(r2)

        # Pair 3: FF (Both Fwd - not proper pair)
        r1_start, r2_start = 10, 50 # Reuse positions
        tlen = 0 # TLEN undefined / 0 for non-FR/RF
        r1 = pysam.AlignedSegment()
        r1.query_name = "revcomp_pair_FF"
        r1.query_sequence = "N" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("&" * read_len)
        set_paired_flags(r1, is_read1=True) # Read Fwd
        r1.reference_id = ref_id; r1.reference_start = r1_start
        r1.mapping_quality = 60; r1.cigarstring = f"{read_len}M"
        set_mate_info(r1, header, ref_name, r2_start, False, False) # Mate Fwd
        r1.template_length = tlen; r1.flag = 65 # Paired, R1 (Not proper, Mate not Rev)
        samfile.write(r1)
        r2 = pysam.AlignedSegment()
        r2.query_name = "revcomp_pair_FF"
        r2.query_sequence = "Y" * read_len
        r2.query_qualities = pysam.qualitystring_to_array("'" * read_len)
        set_paired_flags(r2, is_read1=False) # Read Fwd
        r2.reference_id = ref_id; r2.reference_start = r2_start
        r2.mapping_quality = 60; r2.cigarstring = f"{read_len}M"
        set_mate_info(r2, header, ref_name, r1_start, False, False) # Mate Fwd
        r2.template_length = tlen; r2.flag = 129 # Paired, R2 (Not proper, Mate not Rev)
        samfile.write(r2)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")
    # tqdm.write(f"{Fore.GREEN}  Note: Assumed test requires pairs with different is_reverse/mate_is_reverse flag combinations.")

def generate_sam_35(output_dir: Path, **kwargs):
    """SAM_35: Reverse Complement (different + same) – long distance"""
    # Assumption: Similar to SAM_34 but with large separation / TLEN
    # Combine SAM_09 logic with SAM_34 logic
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    read_len = 15
    r1_start = 5
    r2_start = 140 # Far apart

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Pair 1: FR, long distance
        tlen = (r2_start + read_len - 1) - r1_start + 1
        r1 = pysam.AlignedSegment()
        r1.query_name = "revcomp_long_FR"
        r1.query_sequence = "A" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("!" * read_len)
        set_paired_flags(r1, is_read1=True)
        r1.reference_id = ref_id; r1.reference_start = r1_start
        r1.mapping_quality = 60; r1.cigarstring = f"{read_len}M"
        set_mate_info(r1, header, ref_name, r2_start, False, True) # Mate is Rev
        r1.template_length = tlen
        # Flag: Paired, Mate Rev, R1. Maybe not 'proper' (0x2) due to distance.
        r1.flag = 97 # 0x1 + 0x20 + 0x40
        samfile.write(r1)
        r2 = pysam.AlignedSegment()
        r2.query_name = "revcomp_long_FR"
        r2.query_sequence = "T" * read_len
        r2.query_qualities = pysam.qualitystring_to_array("#" * read_len)
        set_paired_flags(r2, is_read1=False); r2.is_reverse = True # Read is Rev
        r2.reference_id = ref_id; r2.reference_start = r2_start
        r2.mapping_quality = 60; r2.cigarstring = f"{read_len}M"
        set_mate_info(r2, header, ref_name, r1_start, False, False) # Mate is Fwd
        r2.template_length = -tlen
        # Flag: Paired, Read Rev, R2. Maybe not 'proper'.
        r2.flag = 145 # 0x1 + 0x10 + 0x80
        samfile.write(r2)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")
    # tqdm.write(f"{Fore.GREEN}  Note: Assumed test requires pairs with different reverse flags and large separation.")


def generate_sam_36(output_dir: Path, **kwargs):
    """SAM_36: Reverse Complement (different + same) – unmapped"""
    # Assumption: Show unmapped pairs with different reverse flag combinations.
    file_path = output_dir / "alignment.sam"
    header = utils.get_default_sam_header()
    read_len = 10

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Pair 1: R1 Fwd, R2 Rev (but both unmapped)
        r1 = pysam.AlignedSegment()
        r1.query_name = "revcomp_unmapped_FRlike"
        r1.query_sequence = "A" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("!" * read_len)
        set_paired_flags(r1, is_read1=True) # R1 Fwd
        r1.is_unmapped = True; r1.mate_is_unmapped = True
        r1.mate_is_reverse = True # Mate Rev flag set
        r1.flag = 0x1 + 0x4 + 0x8 + 0x20 + 0x40 # = 109
        samfile.write(r1)
        r2 = pysam.AlignedSegment()
        r2.query_name = "revcomp_unmapped_FRlike"
        r2.query_sequence = "T" * read_len
        r2.query_qualities = pysam.qualitystring_to_array("#" * read_len)
        set_paired_flags(r2, is_read1=False); r2.is_reverse = True # R2 Rev
        r2.is_unmapped = True; r2.mate_is_unmapped = True
        r2.mate_is_reverse = False # Mate Fwd flag set
        r2.flag = 0x1 + 0x4 + 0x8 + 0x10 + 0x80 # = 141 (Note: This is same as basic unmapped R2)
        # Let's recalculate R2 flag: Paired(1) + Unmapped(4) + MateUnmapped(8) + ReadRev(16) + R2(128) = 157
        r2.flag = 157
        samfile.write(r2)

        # Pair 2: Both Fwd (unmapped)
        r1 = pysam.AlignedSegment()
        r1.query_name = "revcomp_unmapped_FFlike"
        r1.query_sequence = "C" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("$" * read_len)
        set_paired_flags(r1, is_read1=True) # R1 Fwd
        r1.is_unmapped = True; r1.mate_is_unmapped = True
        r1.mate_is_reverse = False # Mate Fwd flag set
        r1.flag = 0x1 + 0x4 + 0x8 + 0x40 # = 77
        samfile.write(r1)
        r2 = pysam.AlignedSegment()
        r2.query_name = "revcomp_unmapped_FFlike"
        r2.query_sequence = "G" * read_len
        r2.query_qualities = pysam.qualitystring_to_array("%" * read_len)
        set_paired_flags(r2, is_read1=False) # R2 Fwd
        r2.is_unmapped = True; r2.mate_is_unmapped = True
        r2.mate_is_reverse = False # Mate Fwd flag set
        r2.flag = 0x1 + 0x4 + 0x8 + 0x80 # = 133
        samfile.write(r2)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Note: Assumed test requires unmapped pairs with different reverse flags.")

def generate_sam_37(output_dir: Path, **kwargs):
    """SAM_37: Reverse Complement (different + same) – half mapped"""
    # Assumption: Show half-mapped pairs with different reverse flag combinations.
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    map_pos = 20
    read_len = 15

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        # Pair 1: R1 mapped Fwd, R2 unmapped Rev
        r1 = pysam.AlignedSegment()
        r1.query_name = "revcomp_half_1"
        r1.query_sequence = "A" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("!" * read_len)
        set_paired_flags(r1, is_read1=True) # R1 Fwd
        r1.reference_id = ref_id; r1.reference_start = map_pos
        r1.mapping_quality = 60; r1.cigarstring = f"{read_len}M"
        r1.mate_is_unmapped = True
        r1.mate_is_reverse = True # Mate Rev flag set
        # Flags: Paired(1) + MateUnmapped(8) + MateRev(32) + R1(64) = 105
        r1.flag = 105
        samfile.write(r1)
        r2 = pysam.AlignedSegment()
        r2.query_name = "revcomp_half_1"
        r2.query_sequence = "T" * read_len
        r2.query_qualities = pysam.qualitystring_to_array("#" * read_len)
        set_paired_flags(r2, is_read1=False); r2.is_reverse = True # R2 Rev
        r2.is_unmapped = True
        r2.mate_is_unmapped = False # Mate is mapped
        r2.mate_is_reverse = False # Mate is Fwd
        set_mate_info(r2, header, ref_name, map_pos, False, False)
        # Flags: Paired(1) + Unmapped(4) + ReadRev(16) + R2(128) = 149
        r2.flag = 149
        samfile.write(r2)

        # Pair 2: R1 mapped Rev, R2 unmapped Fwd
        map_pos = 60
        r1 = pysam.AlignedSegment()
        r1.query_name = "revcomp_half_2"
        r1.query_sequence = "C" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("$" * read_len)
        set_paired_flags(r1, is_read1=True); r1.is_reverse = True # R1 Rev
        r1.reference_id = ref_id; r1.reference_start = map_pos
        r1.mapping_quality = 60; r1.cigarstring = f"{read_len}M"
        r1.mate_is_unmapped = True
        r1.mate_is_reverse = False # Mate Fwd flag set
        # Flags: Paired(1) + MateUnmapped(8) + ReadRev(16) + R1(64) = 89
        r1.flag = 89
        samfile.write(r1)
        r2 = pysam.AlignedSegment()
        r2.query_name = "revcomp_half_2"
        r2.query_sequence = "G" * read_len
        r2.query_qualities = pysam.qualitystring_to_array("%" * read_len)
        set_paired_flags(r2, is_read1=False) # R2 Fwd
        r2.is_unmapped = True
        r2.mate_is_unmapped = False # Mate is mapped
        r2.mate_is_reverse = True # Mate is Rev
        set_mate_info(r2, header, ref_name, map_pos, False, True)
        # Flags: Paired(1) + Unmapped(4) + MateRev(32) + R2(128) = 165
        r2.flag = 165
        samfile.write(r2)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")
    # tqdm.write(f"{Fore.GREEN}  Note: Assumed test requires half-mapped pairs with different reverse flags.")

def generate_sam_38(output_dir: Path, **kwargs):
    """SAM_38: Circular reference"""
    # This requires the aligner to handle circularity. We simulate a read
    # that *could* wrap around if the reference were circular.
    # Pysam itself doesn't enforce circularity; the test tool would need to.
    file_path = output_dir / "alignment.sam"
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    ref_len = header.get_reference_length(ref_name) # 160
    read_len = 30
    # Place read near the end, e.g., starting at ref_len - 10
    start_pos = ref_len - 10 # = 150

    with pysam.AlignmentFile(str(file_path), "w", header=header) as samfile:
        r1 = pysam.AlignedSegment()
        r1.query_name = "circular_ref_read"
        # Sequence would need to match ref[150:160] + ref[0:20] if circular
        # We just put a placeholder sequence here.
        r1.query_sequence = "A" * read_len
        r1.query_qualities = pysam.qualitystring_to_array("!" * read_len)
        r1.reference_id = ref_id
        r1.reference_start = start_pos # Starts near end
        r1.mapping_quality = 60
        # CIGAR would likely be split if aligner handled circularity,
        # e.g., 10M followed by alignment at start.
        # For this test data, just show a simple match starting near end.
        r1.cigarstring = f"{read_len}M" # This will run off the end of linear ref
        r1.flag = 0
        # Add a note that the tool needs to interpret this against a circular ref
        r1.set_tag("ZC", "circular_test", "Z")
        samfile.write(r1)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")
    # tqdm.write(f"{Fore.GREEN}  Note: Read placed near end of linear reference. Tool under test must handle circularity.")

def generate_sam_39(output_dir: Path, **kwargs):
    """SAM_39: (bam input) - Generate BAM"""
    file_path = output_dir / "alignment.bam" # Output BAM
    ref_path = utils.copy_reference_to_output(output_dir)
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    read_len = 12

    # Use 'wb' mode for BAM output
    with pysam.AlignmentFile(str(file_path), "wb", header=header) as bamfile:
        # Add a simple mapped read (similar to generate_sam_04)
        r1 = pysam.AlignedSegment()
        r1.query_name = "bam_input_read_1"
        with pysam.FastaFile(str(ref_path)) as fasta:
            seq = fasta.fetch(ref_name, 0, read_len).upper()
            r1.query_sequence = seq
        r1.query_qualities = pysam.qualitystring_to_array("!" * read_len)
        r1.reference_id = ref_id
        r1.reference_start = 0
        r1.mapping_quality = 60
        r1.cigarstring = f"{read_len}M"
        r1.flag = 0
        bamfile.write(r1)

        r2 = pysam.AlignedSegment()
        r2.query_name = "bam_input_read_2"
        with pysam.FastaFile(str(ref_path)) as fasta:
            seq = fasta.fetch(ref_name, 50, 50 + read_len).upper()
            r2.query_sequence = seq
        r2.query_qualities = pysam.qualitystring_to_array("#" * read_len)
        r2.reference_id = ref_id
        r2.reference_start = 50
        r2.mapping_quality = 60
        r2.cigarstring = f"{read_len}M"
        r2.flag = 0
        bamfile.write(r2)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_40(output_dir: Path, **kwargs):
    """SAM_40: (bam output) - Generate SAM input for testing BAM output"""
    # This case expects the *tool under test* to produce bam output.
    # So we provide a standard sam file as input. Reference also needed.
    # Reuse SAM_04 generator.
    generate_sam_04(output_dir, **kwargs)
    # tqdm.write(f"{Fore.GREEN}  Note: Generated standard alignment.sam and copied reference for testing bam output capability.")
    # tqdm.write(f"{Fore.GREEN}  Note: Generated standard alignment.sam for testing bam output capability.")

def generate_sam_41(output_dir: Path, **kwargs):
    """SAM_41: (cram input) - Generate CRAM"""
    file_path = output_dir / "alignment.cram" # Output CRAM
    ref_path = utils.copy_reference_to_output(output_dir) # Reference is MANDATORY for CRAM
    header = utils.get_default_sam_header()
    ref_name = "ref1"
    ref_id = header.references.index(ref_name)
    read_len = 12

    # Use 'wc' mode for CRAM output, provide reference path
    with pysam.AlignmentFile(str(file_path), "wc", header=header, reference_filename=str(ref_path)) as cramfile:
        # Add a simple mapped read (similar to generate_sam_04)
        r1 = pysam.AlignedSegment()
        r1.query_name = "cram_input_read_1"
        with pysam.FastaFile(str(ref_path)) as fasta:
            seq = fasta.fetch(ref_name, 0, read_len).upper()
            r1.query_sequence = seq
        r1.query_qualities = pysam.qualitystring_to_array("!" * read_len)
        r1.reference_id = ref_id
        r1.reference_start = 0
        r1.mapping_quality = 60
        r1.cigarstring = f"{read_len}M"
        r1.flag = 0
        cramfile.write(r1)

        r2 = pysam.AlignedSegment()
        r2.query_name = "cram_input_read_2"
        with pysam.FastaFile(str(ref_path)) as fasta:
            seq = fasta.fetch(ref_name, 50, 50 + read_len).upper()
        r2.query_sequence = seq
        r2.query_qualities = pysam.qualitystring_to_array("#" * read_len)
        r2.reference_id = ref_id
        r2.reference_start = 50
        r2.mapping_quality = 60
        r2.cigarstring = f"{read_len}M"
        r2.flag = 0
        cramfile.write(r2)

    # tqdm.write(f"{Fore.GREEN}   {file_path}")
    # tqdm.write(f"{Fore.GREEN}  Copied reference: {ref_path}")

def generate_sam_42(output_dir: Path, **kwargs):
    """SAM_42: (cram output) - Generate SAM"""
    # This case expects the *tool under test* to produce cram output.
    # So we just provide a standard sam file as input. Reference also needed.
    # Reuse SAM_04 generator.
    generate_sam_04(output_dir, **kwargs)
    # tqdm.write(f"{Fore.GREEN}  Note: Generated standard alignment.sam and copied reference for testing cram output capability.")
