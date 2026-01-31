#!/bin/bash -eu

# Usage:
#   ./try.sh ARG ATP        # process orig/ARG.cif, orig/ATP.cif
#   ./try.sh orig/a/*.cif   # process all files in orig/a/

#:rm -f gemmi/*.cif
#../run-tests.sh G

# Group files by subdirectory for proper output-dir handling
declare -A subdir_files

for x in "$@"; do
  if [[ "$x" == *.cif ]]; then
    # It's a file path - extract subdir
    subdir=$(dirname "$x" | sed 's|^orig||; s|^/||')
    [[ -z "$subdir" ]] && subdir="."
  else
    # It's a monomer code
    x="orig/${x^^}.cif"
    subdir="."
  fi
  subdir_files["$subdir"]+="$x "
done

# Run gemmi drg for each subdirectory group
for subdir in "${!subdir_files[@]}"; do
  read -ra files <<< "${subdir_files[$subdir]}"
  if [[ "$subdir" == "." ]]; then
    outdir="gemmi"
  else
    outdir="gemmi/$subdir"
  fi
  mkdir -p "$outdir"
  echo "Processing ${#files[@]} files -> $outdir"
  ../build/gemmi drg --typeOut -v --output-dir="$outdir" "${files[@]}" 2>>drg.log ||:
done

# Compare each file
for x in "$@"; do
  if [[ "$x" == *.cif ]]; then
    code=$(basename "$x" .cif)
    subdir=$(dirname "$x" | sed 's|^orig||; s|^/||')
  else
    code=${x^^}
    subdir=""
  fi

  if [[ -n "$subdir" ]]; then
    acedrg_file="acedrg/$subdir/$code.cif"
    gemmi_file="gemmi/$subdir/$code.cif"
  else
    acedrg_file="acedrg/$code.cif"
    gemmi_file="gemmi/$code.cif"
  fi

  echo "========================${code}========================"
  ../build/gemmi mondiff "$acedrg_file" "$gemmi_file" ||:
done
echo "Script ended successfully!"
