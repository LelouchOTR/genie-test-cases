from pathlib import Path
import pysam
from . import utils # Relative import

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
        segment.next_reference_id = header.references.index(mate_ref_name)
        segment.next_reference_start = mate_start
    else:
        # If mate is unmapped, standard practice is to set mate position info
        # to the same as the read itself (or 0 if this read is also unmapped)
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

    print(f"  Generated: {file_path}")

[Rest of the SAM/BAM generator functions would go here...]
