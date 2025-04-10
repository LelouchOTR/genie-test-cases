import os
import sys
from pathlib import Path
from tqdm import tqdm
from colorama import Fore, Back, Style, init
init(autoreset=True)  # Auto-reset colors after each print
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
    
    # Create a single progress bar
    with tqdm(total=total_cases, unit="case", 
             bar_format="{l_bar}%s{bar}%s{r_bar}" % (Fore.BLUE, Fore.RESET),
             desc=Fore.WHITE + "Generating") as pbar:
        for case_config in TEST_CASES:
            case_id = case_config['id']
            
            # Update progress bar description with current case
            pbar.set_description(f"{Fore.WHITE}Processing {case_id}")
            
            output_dir = base_output_path / case_config['format'] / case_config['output_subdir']
            utils.ensure_dir(output_dir)

            try:
                func_path = case_config.get('generator_func')
                if func_path:
                    module_name, func_name = func_path.rsplit('.', 1)
                    module = __import__(f'test_data_generator.{module_name}',
                                      fromlist=[func_name])
                    generator_func = getattr(module, func_name)
                    
                    # Run the generator
                    params = case_config.get('params', {})
                    generator_func(output_dir, **params)
                    
                    # Write README after successful generation
                    utils.write_readme(output_dir, case_config)
                    success_count += 1
                else:
                    error_count += 1
                    tqdm.write(f"{Fore.RED}  WARNING: No generator for {case_id}")
                
            except Exception as e:
                error_count += 1
                tqdm.write(f"{Fore.RED}  ERROR in {case_id}: {str(e)}")
                
            finally:
                # Update progress bar by 1 step regardless of success/failure
                pbar.update(1)

    # Final summary
    tqdm.write(Fore.GREEN + f"\n✅ Success: {success_count}/{total_cases}")
    if error_count > 0:
        tqdm.write(Fore.RED + f"❌ Errors: {error_count}/{total_cases}")
    print(Style.RESET_ALL)

if __name__ == "__main__":
    # Ensure the script is run from the directory containing generate_tests.py
    # or adjust paths accordingly if run from elsewhere.
    os.chdir(Path(__file__).parent) # Change working dir to script dir for consistency
    main()
