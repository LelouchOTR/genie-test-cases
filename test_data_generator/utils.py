import os
import pysam
from pathlib import Path
import shutil

# Define the path to the reference file relative to this utils file

def _generate_large_reference_file(ref_path: Path):
    """Generates the large_ref.fa file and its index."""
    # Create repetitive sequence to make a large file
    repeat_unit = "ACGT" * 250  # 1000 bp repeat unit
    num_repeats = 1000          # 1000 * 1000 = 1,000,000 bp
    description = ">large_ref\n" # Simplified description

    try:
        ref_path.parent.mkdir(parents=True, exist_ok=True) # Ensure parent dir exists
        with open(ref_path, "w") as f:
            f.write(description)
            current_len = 0
            line_len = 80 # Standard FASTA line length
            for _ in range(num_repeats):
                f.write(repeat_unit)
                current_len += len(repeat_unit)
                # Add newlines for standard FASTA formatting
                if current_len >= line_len:
                    f.write("\n")
                    current_len = 0
            # Ensure the file ends with a newline if needed
            if current_len > 0:
                 f.write("\n")

        # Create FASTA index silently
        pysam.faidx(str(ref_path))
    except Exception as e:
        # If generation fails, raise an error to stop the process clearly.
        raise RuntimeError(f"Failed to generate large reference {ref_path}: {e}") from e

def ensure_reference_exists(ref_name: str):
    """
    Checks if a reference file exists in the reference source directory.
    If not, generates it (currently only supports 'large_ref.fa').
    """
    ref_dir = _THIS_DIR / "reference"
    ref_path = ref_dir / ref_name
    ref_index_path = ref_path.with_suffix(ref_path.suffix + '.fai')

    # Check if both the FASTA and its index exist
    if not ref_path.exists() or not ref_index_path.exists():
        if ref_name == "large_ref.fa":
            # Generation happens here if needed. No print statements to keep tqdm clean.
            _generate_large_reference_file(ref_path)
        else:
            # Raise an error if a non-large reference is missing, as we don't auto-generate others.
             raise FileNotFoundError(f"Required reference file '{ref_name}' not found in {ref_dir} and cannot be auto-generated.")
_THIS_DIR = Path(__file__).parent
REFERENCE_FASTA_PATH = _THIS_DIR / "reference" / "simple_ref.fa"

def ensure_dir(dir_path: Path):
    """Creates a directory if it doesn't exist."""
    dir_path.mkdir(parents=True, exist_ok=True)

def write_readme(output_dir: Path, case_config: dict):
    """Creates a README.md file in the output directory describing the test case."""
    readme_path = output_dir / "README.md"
    content = f"""## Test Case: {case_config['name']}

**Description:** {case_config['description']}

**Format:** {case_config['format']}

**Generated Files:**
"""
    for filename in case_config.get('output_files', ['N/A']):
        content += f"- `{filename}`\n"

    if 'notes' in case_config:
        content += f"\n**Notes:**\n{case_config['notes']}\n"

    with open(readme_path, 'w') as f:
        f.write(content)

def create_fastq_entry(read_id: str, sequence: str, quality: str, comment: str = "") -> str:
    """Formats a single FASTQ entry."""
    if len(sequence) != len(quality):
        raise ValueError(f"Sequence length ({len(sequence)}) and quality length ({len(quality)}) must match for read {read_id}")
    comment_str = f" {comment}" if comment else ""
    return f"@{read_id}{comment_str}\n{sequence}\n+\n{quality}\n"

def get_default_sam_header() -> pysam.AlignmentHeader:
    """Creates a basic pysam AlignmentHeader using the default reference."""
    ensure_reference_exists("simple_ref.fa")
    
    # Create header directly from the reference FASTA
    with pysam.FastaFile(str(REFERENCE_FASTA_PATH)) as fasta:
        references = [(name, fasta.get_reference_length(name)) 
                     for name in fasta.references]
    
    header_dict = {
        'HD': {'VN': '1.6', 'SO': 'unsorted'},
        'SQ': [{'SN': name, 'LN': length} for name, length in references]
    }
    return pysam.AlignmentHeader.from_dict(header_dict)

def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.
    
    Args:
        seq: Input DNA sequence string
        
    Returns:
        Reverse complemented sequence
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                 'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
                 'N': 'N', 'n': 'n'}
    return ''.join([complement.get(base, 'N') for base in reversed(seq)])

def copy_reference_to_output(output_dir: Path, ref_name: str = "simple_ref.fa") -> Path:
    """Copies the specified reference FASTA to the output directory."""
    src_path = _THIS_DIR / "reference" / ref_name
    dest_path = output_dir / src_path.name
    
    # Copy main FASTA file
    shutil.copy(str(src_path), str(dest_path))
    
    # Also copy the index if it exists
    for ext in ['.fai', '.amb', '.ann', '.bwt', '.pac', '.sa']:
        index_src = src_path.with_suffix(src_path.suffix + ext)
        if index_src.exists():
            shutil.copy(str(index_src), str(dest_path.with_suffix(dest_path.suffix + ext)))
    
    return dest_path
