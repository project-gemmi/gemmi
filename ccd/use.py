#!/usr/bin/env python
#
# This script reads monomer definitions from a gzipped CIF file (components.cif.gz)
# and extracts specified monomer blocks into separate CIF files.
#
# This script reads monomer definitions from a gzipped CIF file (components.cif.gz)
# and extracts specified monomer blocks into separate CIF files.
#
# For example, to extract the 'ARG' and 'ATP' monomers:
# python use.py ARG ATP
#
# This will create two files: ./orig/ARG.cif and ./orig/ATP.cif.

import sys
import os
import gzip
import subprocess

# Programmatically add the gemmi build path to sys.path
# This allows 'import gemmi' to succeed without an explicit PYTHONPATH env var.
gemmi_build_path = os.path.join(os.path.expanduser("~"), "gemmi", "gemmi", "build", "py")
if gemmi_build_path not in sys.path:
    sys.path.insert(0, gemmi_build_path)

try:
    import gemmi
except ImportError:
    print("Error: The 'gemmi' module could not be imported.")
    print("Please ensure that the 'gemmi' library is installed at:")
    print(f"  {gemmi_build_path}")
    print("Or that it's otherwise available in your Python environment.")
    sys.exit(1)


def run_acedrg(code: str, input_file: str, output_dir: str):
    """
    Runs the Acedrg command for a given monomer.

    Args:
        code: The monomer code (e.g., 'ATP').
        input_file: Path to the input CIF file (e.g., 'orig/ATP.cif').
        output_dir: The directory for Acedrg output (e.g., 'acedrg').
    """
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            print(f"Created directory: {output_dir}")
        except OSError as e:
            print(f"Error creating directory {output_dir}: {e}")
            return  # Stop if we can't create the output directory

    output_prefix = os.path.join(output_dir, code)
    #command = ["acedrg", "-c", input_file, "-o", output_prefix]
    command = ["acedrg", "--typeOut", "-c", input_file, "-o", output_prefix]
    print(f"\nRunning: {' '.join(command)}")

    try:
        # Execute the command, capturing output and checking for errors
        _ = subprocess.run(
            command, check=True, capture_output=True, text=True, encoding='utf-8'
        )
        print(f"Acedrg finished successfully for {code}.")
        # Optionally, uncomment the following lines to see Acedrg's output
        # if _.stdout:
        #     print("Acedrg Log:\n" + _.stdout)
        # if _.stderr:
        #     print("Acedrg Errors:\n" + _.stderr)

    except FileNotFoundError:
        print("\nError: 'acedrg' command not found.")
        print("Please ensure Acedrg is installed and its location is in your system's PATH.")
        # Exit the script because Acedrg is a required part of the process
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"\nError running Acedrg for {code}:")
        print(f"Command failed: {' '.join(command)}")
        print(f"Return code: {e.returncode}")
        print(f"--- Acedrg STDOUT ---\n{e.stdout}")
        print(f"--- Acedrg STDERR ---\n{e.stderr}")
    except Exception as e:
        print(f"\nAn unexpected error occurred while running Acedrg for {code}: {e}")


def extract_and_process_monomers(
    monomer_codes: list[str], input_cif: str, orig_dir: str, acedrg_dir: str
):
    """
    Reads a CIF file, finds monomer blocks, writes them to files, and runs Acedrg.

    Args:
        monomer_codes: A list of monomer codes to extract.
        input_cif: Path to the input gzipped CIF file.
        orig_dir: The directory to save the initial monomer .cif files.
        acedrg_dir: The directory for the Acedrg output.
    """
    if not os.path.exists(orig_dir):
        try:
            os.makedirs(orig_dir)
            print(f"Created directory: {orig_dir}")
        except OSError as e:
            print(f"Error creating directory {orig_dir}: {e}")
            sys.exit(1)

    try:
        with gzip.open(input_cif, "rt", encoding="utf-8") as gz_file:
            cif_content = gz_file.read()
        doc = gemmi.cif.read_string(cif_content)
    except FileNotFoundError:
        print(f"Error: Input file not found at '{input_cif}'")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading or parsing '{input_cif}': {e}")
        sys.exit(1)

    for code in monomer_codes:
        block = doc.find_block(code)
        if block:
            new_doc = gemmi.cif.Document()
            new_doc.add_copied_block(block)

            output_filename = f"{code}.cif"
            output_path = os.path.join(orig_dir, output_filename)

            try:
                new_doc.write_file(output_path)
                print(f"Successfully wrote {output_path}")
                # After writing the file, run Acedrg
                run_acedrg(code, output_path, acedrg_dir)
            except Exception as e:
                print(f"Error writing file or running Acedrg for monomer '{code}': {e}")
        else:
            print(f"\nWarning: Monomer '{code}' not found in '{input_cif}'")


def main():
    """
    Main function to parse command-line arguments and initiate extraction and processing.
    """
    if len(sys.argv) < 2:
        print(f"Usage: python {sys.argv[0]} MONOMER1 MONOMER2 ...")
        print("Example: python {sys.argv[0]} ARG ATP")
        sys.exit(1)

    monomer_codes = [arg.upper() for arg in sys.argv[1:]]

    input_cif_file = "components.cif.gz"
    orig_output_dir = "orig"
    acedrg_output_dir = "acedrg"

    extract_and_process_monomers(monomer_codes, input_cif_file, orig_output_dir, acedrg_output_dir)


if __name__ == "__main__":
    main()
