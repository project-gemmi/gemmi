#!/bin/bash -eu

# Usage:
#   ./try.sh ARG ATP        # process orig/ARG.cif, orig/ATP.cif
#   ./try.sh orig/a/*.cif   # process all files in orig/a/

#:rm -f gemmi/*.cif
#../run-tests.sh G

files=()
for x in "$@"; do
  if [[ "$x" == *.cif ]]; then
    # It's a file path
    files+=("$x")
  else
    # It's a monomer code
    code=${x^^}
    files+=("orig/$code.cif")
  fi
done

# Run gemmi drg on all files
../build/gemmi drg --typeOut -v --output-dir=gemmi "${files[@]}" 2>drg.log ||:

# Compare each file
for f in "${files[@]}"; do
  if [[ "$f" == *.cif ]]; then
    # Extract code and subdir from path like orig/a/ARG.cif
    code=$(basename "$f" .cif)
    subdir=$(dirname "$f" | sed 's|^orig||; s|^/||')
    if [[ -n "$subdir" ]]; then
      acedrg_file="acedrg/$subdir/$code.cif"
      gemmi_file="gemmi/$subdir/$code.cif"
    else
      acedrg_file="acedrg/$code.cif"
      gemmi_file="gemmi/$code.cif"
    fi
  else
    code=${f^^}
    acedrg_file="acedrg/$code.cif"
    gemmi_file="gemmi/$code.cif"
  fi
  echo "========================${code}========================"
  ../build/gemmi mondiff "$acedrg_file" "$gemmi_file" ||:
done
echo "Script ended successfully!"
