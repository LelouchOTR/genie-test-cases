import os
import sys
from pathlib import Path
from .config import TEST_CASES, BASE_OUTPUT_DIR
from . import utils # Import utils to access helper functions if needed directly

def main():
    """Generates all configured test data files."""
    print("Starting test data generation...")
    base_output_path = Path(BASE_OUTPUT_DIR)
    print(f"Output directory: {base_output_path.resolve()}")

    generated_count = 0
    for i, case_config in enumerate(TEST_CASES):
        print(f"\nProcessing Case {i+1}/{len(TEST_CASES)}: {case_config['name']} ({case_config['id']})")

        # Construct output path: base_dir / format / subdir
        output_dir = base_output_path / case_config['format'] / case_config['output_subdir']

        # Create the directory
        utils.ensure_dir(output_dir)
        print(f"  Ensured directory: {output_dir}")

        # Call the specific generator function
        generator_func = case_config.get('generator_func')
        if generator_func:
            try:
                # Pass output directory and any specific params
                params = case_config.get('params', {})
                generator_func(output_dir, **params)
                generated_count += 1
            except Exception as e:
                print(f"  ERROR generating data for {case_config['id']}: {e}")
                # Optionally: raise e # Stop execution on error
                continue # Continue to next case
        else:
            print(f"  WARNING: No generator function defined for {case_config['id']}.")
            continue # Skip README generation if no data was attempted

        # Write the README file
        try:
            utils.write_readme(output_dir, case_config)
            print(f"  Generated: {output_dir / 'README.md'}")
        except Exception as e:
            print(f"  ERROR generating README for {case_config['id']}: {e}")


    print(f"\nTest data generation complete. Generated data for {generated_count} cases.")

if __name__ == "__main__":
    # Ensure the script is run from the directory containing generate_tests.py
    # or adjust paths accordingly if run from elsewhere.
    os.chdir(_THIS_DIR) # Change working dir to script dir for consistency
    main()
