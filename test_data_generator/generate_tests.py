import os
import sys
from pathlib import Path
from .config import TEST_CASES, BASE_OUTPUT_DIR
from . import utils # Import utils to access helper functions if needed directly

def main():
    """Generates all configured test data files."""
    print(Fore.CYAN + "\nStarting test data generation...")
    base_output_path = Path(BASE_OUTPUT_DIR)
    print(Fore.YELLOW + f"Output directory: {base_output_path.resolve()}")

    success_count = 0
    error_count = 0
    total_cases = len(TEST_CASES)
    
    # Create progress bar
    with tqdm(TEST_CASES, unit="case", 
             bar_format="{l_bar}%s{bar}%s{r_bar}" % (Fore.BLUE, Fore.RESET)) as pbar:
        for case_config in pbar:
            case_id = case_config['id']
            pbar.set_description(f"{Fore.WHITE}Processing {case_id}")

        # Construct output path: base_dir / format / subdir
        output_dir = base_output_path / case_config['format'] / case_config['output_subdir']

        # Create the directory
        utils.ensure_dir(output_dir)
        print(f"  Ensured directory: {output_dir}")

        # Load and call the specific generator function
        func_path = case_config.get('generator_func')
        if func_path:
            try:
                # Split "module.function" into components
                module_name, func_name = func_path.rsplit('.', 1)
                # Dynamically import the module
                module = __import__(f'test_data_generator.{module_name}',
                                  fromlist=[func_name])
                # Get the actual function object
                generator_func = getattr(module, func_name)
                
                # Call it with params
                params = case_config.get('params', {})
                generator_func(output_dir, **params)
                success_count += 1
            except Exception as e:
                error_count += 1
                print(f"{Fore.RED}  ERROR in {case_id}: {str(e)}")
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


    # Final summary
    print(Fore.GREEN + f"\n✅ Successfully generated {success_count}/{total_cases} cases")
    if error_count > 0:
        print(Fore.RED + f"❌ Failed to generate {error_count}/{total_cases} cases")
    print(Style.RESET_ALL)

if __name__ == "__main__":
    # Ensure the script is run from the directory containing generate_tests.py
    # or adjust paths accordingly if run from elsewhere.
    os.chdir(Path(__file__).parent) # Change working dir to script dir for consistency
    main()
