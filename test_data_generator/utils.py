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

def copy_reference_to_output(output_dir: Path) -> Path:
    """Copies the default reference FASTA to the output directory."""
    dest_path = output_dir / REFERENCE_FASTA_PATH.name
    shutil.copy(str(REFERENCE_FASTA_PATH), str(dest_path))
    # Also copy the index if it exists
    ref_fai = REFERENCE_FASTA_PATH.with_suffix(REFERENCE_FASTA_PATH.suffix + '.fai')
    if ref_fai.exists():
         shutil.copy(str(ref_fai), str(dest_path.with_suffix(dest_path.suffix + '.fai')))
    return dest_path
