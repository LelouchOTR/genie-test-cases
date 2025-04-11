import sys
from pathlib import Path
import pysam

def main():
    """Generate a large reference genome file (>1M bases) for testing."""
    ref_path = Path(__file__).parent / "large_ref.fa"
    
    # Create repetitive sequence to make a large file
    repeat_unit = "ACGT" * 250  # 1000 bp repeat unit
    num_repeats = 1000          # 1000 * 1000 = 1,000,000 bp
    description = ">large_ref circular=true\n"
    
    with open(ref_path, "w") as f:
        f.write(description)
        for _ in range(num_repeats):
            f.write(repeat_unit)
            
    print(f"Generated large reference: {ref_path}")
    print("Creating FASTA index...")
    pysam.faidx(str(ref_path))

if __name__ == "__main__":
    main()
