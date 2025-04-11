import os
import pysam
from pathlib import Path
import shutil

# Define the path to the reference file relative to this utils file
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
    if not REFERENCE_FASTA_PATH.exists():
        raise FileNotFoundError(f"Reference FASTA not found at: {REFERENCE_FASTA_PATH}")
    
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
