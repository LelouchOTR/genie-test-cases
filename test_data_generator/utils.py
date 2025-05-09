import os
import pysam
from pathlib import Path
import shutil
import textwrap
from tqdm import tqdm

# Define the path to the reference file relative to this utils file

def _generate_large_reference_file(ref_path: Path) -> None:
    """Generates a large reference FASTA file (large_ref.fa) and its index.

    Creates a repetitive sequence (~1Mbp ACGT repeats) to efficiently generate
    a large reference file and then indexes it using pysam.faidx.

    Args:
        ref_path: The full path where the reference FASTA file should be created.
                  The index file (.fai) will be created alongside it.

    Raises:
        RuntimeError: If the file generation or indexing fails.
    """
    # Create repetitive sequence to make a large file
    repeat_unit = "ACGT" * 250  # 1000 bp repeat unit
    num_repeats = 2000          # 2000 * 1000 = 2,000,000 bp (~2 Mbp)
    description = ">large_ref\n" # Simplified description

    try:
        ref_path.parent.mkdir(parents=True, exist_ok=True) # Ensure parent dir exists
        with open(ref_path, "w") as f:
            f.write(description)
            current_len = 0
            FASTA_LINE_LENGTH = 80  # Standard FASTA line length
            line_len = FASTA_LINE_LENGTH
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

def ensure_reference_exists(ref_name: str) -> None:
    """Checks if a reference FASTA and its index exist, generating if needed.

    Verifies the presence of the specified reference FASTA file and its '.fai'
    index in the 'reference/' directory. If `ref_name` is 'large_ref.fa' and
    it or its index is missing, it calls `_generate_large_reference_file`.
    For other missing reference files, it raises an error.

    Args:
        ref_name: The base name of the reference file (e.g., "simple_ref.fa").

    Raises:
        FileNotFoundError: If a required reference (other than 'large_ref.fa')
                           or its index is not found.
        RuntimeError: If generation of 'large_ref.fa' fails.
    """
    ref_dir = _THIS_DIR / "reference"
    ref_path = ref_dir / ref_name
    ref_index_path = ref_path.with_suffix(ref_path.suffix + '.fai')

    # Check if both the FASTA and its index exist
    if not ref_path.exists() or not ref_index_path.exists():
        if ref_name == "large_ref.fa":
            if not ref_path.exists():
                tqdm.write(f"Generating large reference file: {ref_path}...")
                _generate_large_reference_file(ref_path)  # This also creates the index
            elif not ref_index_path.exists():
                # Should be created by _generate_large_reference_file, but index might have been deleted
                tqdm.write(f"Generating index for existing large reference file: {ref_path}...")
                try:
                    pysam.faidx(str(ref_path))
                except Exception as e:
                    raise RuntimeError(f"Failed to generate index for {ref_path}: {e}") from e
        else:
            # For non-large references, check if the FASTA exists first
            if not ref_path.exists():
                raise FileNotFoundError(f"Required reference file '{ref_name}' not found in {ref_dir}.")
            # If FASTA exists but index doesn't, try to create the index
            elif not ref_index_path.exists():
                tqdm.write(f"Index file '{ref_index_path.name}' not found. Attempting to create index for {ref_name}...")
                try:
                    pysam.faidx(str(ref_path))
                    tqdm.write(f"Successfully created index: {ref_index_path.name}")
                except Exception as e:
                    raise RuntimeError(f"Failed to create index for {ref_name}: {e}") from e


_THIS_DIR = Path(__file__).parent
REFERENCE_FASTA_PATH = _THIS_DIR / "reference" / "simple_ref.fa"

def ensure_dir(dir_path: Path) -> None:
    """Creates a directory if it doesn't exist, including parent directories.

    Args:
        dir_path: The pathlib.Path object representing the directory to create.
    """
    dir_path.mkdir(parents=True, exist_ok=True)

def write_readme(output_dir: Path, case_config: dict) -> None:
    """Creates a README.md file in the output directory describing the test case.

    Uses information from the test case configuration dictionary to generate
    a descriptive README file.

    Args:
        output_dir: The directory where the README.md file will be written.
        case_config: The dictionary containing configuration details for the
                     specific test case (from config.TEST_CASES).
    """
    readme_path = output_dir / "README.md"
    # Use textwrap.dedent to allow cleaner indentation in the source code
    content = textwrap.dedent(f"""\
        ## Test Case: {case_config['name']}

        **Description:** {case_config['description']}

        **Format:** {case_config['format']}
    """)

    # Manually add the generated files list, ensuring correct formatting
    content += "\n**Generated Files:**\n"
    for filename in case_config.get('output_files', ['N/A']):
        content += f"- `{filename}`\n"

    if 'notes' in case_config:
        # Add notes with proper spacing
        content += f"\n**Notes:**\n{case_config['notes']}\n"

    with open(readme_path, 'w') as f:
        f.write(content)

def create_fastq_entry(read_id: str, sequence: str, quality: str, comment: str = "") -> str:
    """Formats a single FASTQ entry as a string.

    Args:
        read_id: The read identifier (without the starting '@').
        sequence: The nucleotide sequence.
        quality: The Phred quality score string. Must be same length as sequence.
        comment: Optional comment to append after the read ID.

    Returns:
        A string formatted as a complete 4-line FASTQ entry.

    Raises:
        ValueError: If the sequence and quality strings have different lengths.
    """
    if len(sequence) != len(quality):
        raise ValueError(f"Sequence length ({len(sequence)}) and quality length ({len(quality)}) must match for read {read_id}")
    comment_str = f" {comment}" if comment else ""
    return f"@{read_id}{comment_str}\n{sequence}\n+\n{quality}\n"

def get_default_sam_header() -> pysam.AlignmentHeader:
    """Creates a basic pysam AlignmentHeader using the default 'simple_ref.fa'.

    Ensures 'simple_ref.fa' exists, then reads its sequence names and lengths
    to construct a minimal SAM header dictionary suitable for pysam.

    Returns:
        A pysam.AlignmentHeader object based on 'simple_ref.fa'.
    """
    ensure_reference_exists("simple_ref.fa")
    
    # Create header directly from the reference FASTA
    with pysam.FastaFile(str(REFERENCE_FASTA_PATH)) as fasta:
        references = [(name, fasta.get_reference_length(name)) 
                     for name in fasta.references]
    
    # Get the absolute path URI for the reference file
    ref_uri = REFERENCE_FASTA_PATH.absolute().as_uri()
    header_dict = {
        'HD': {'VN': '1.6', 'SO': 'unsorted'},
        'SQ': [{'SN': name, 'LN': length, 'UR': ref_uri} for name, length in references]
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
    """Copies reference FASTA and generates .fai index in output directory.
    
    Copies the specified reference FASTA to the output directory and ensures
    a .fai index exists. Also copies any BWA-style index files if present.

    Args:
        output_dir: Destination directory for the reference files
        ref_name: Base name of reference FASTA file (default: simple_ref.fa)

    Returns:
        Path to copied reference FASTA in output directory

    Raises:
        FileNotFoundError: If source reference FASTA doesn't exist
        RuntimeError: If .fai index generation fails
    """
    src_path = _THIS_DIR / "reference" / ref_name
    dest_path = output_dir / src_path.name
    
    # Copy FASTA file
    shutil.copy(str(src_path), str(dest_path))
    
    # Copy BWA-style indexes if present
    bwa_exts = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    for ext in bwa_exts:
        if (src_idx := src_path.with_suffix(src_path.suffix + ext)).exists():
            shutil.copy(str(src_idx), str(dest_path.with_suffix(dest_path.suffix + ext)))
    
    # Generate .fai index if missing
    dest_fai = dest_path.with_suffix(dest_path.suffix + '.fai')
    if not dest_fai.exists():
        try:
            pysam.faidx(str(dest_path))
        except Exception as e:
            raise RuntimeError(f"Failed to generate .fai index for {dest_path}: {e}") from e
    
    return dest_path
