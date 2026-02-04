#!/bin/bash -eu

# Usage:
#   ./try.sh ARG ATP        # process orig/ARG.cif, orig/ATP.cif
#   ./try.sh orig/a/*.cif   # process all files in orig/a/

#:rm -f gemmi/*.cif
#../run-tests.sh G

# Convert args to absolute input files and find unique subdirs
declare -a all_files
declare -a subdirs
for x in "$@"; do
  if [[ "$x" == *.cif ]]; then
    all_files+=("$x")
    subdir=$(dirname "$x" | sed 's|^orig||; s|^/||')
    [[ -z "$subdir" ]] && subdir="."
  else
    all_files+=("orig/${x^^}.cif")
    subdir="."
  fi
  # Add to subdirs if not already present
  if [[ ! " ${subdirs[*]:-} " =~ " $subdir " ]]; then
    subdirs+=("$subdir")
  fi
done

# Run gemmi drg for each subdirectory group
for subdir in "${subdirs[@]}"; do
  if [[ "$subdir" == "." ]]; then
    outdir="gemmi"
    pattern="orig/*.cif"
  else
    outdir="gemmi/$subdir"
    pattern="orig/$subdir/*.cif"
  fi
  # Filter files matching this subdir
  files=()
  for f in "${all_files[@]}"; do
    f_subdir=$(dirname "$f" | sed 's|^orig||; s|^/||')
    [[ -z "$f_subdir" ]] && f_subdir="."
    if [[ "$f_subdir" == "$subdir" ]]; then
      files+=("$f")
    fi
  done
  mkdir -p "$outdir"
  echo "Processing ${#files[@]} files -> $outdir"
  if ! ../build/gemmi drg  --output-dir="$outdir" "${files[@]}" 2> >(tee -a drg.log >&2); then
    echo "ERROR: gemmi drg failed. Check drg.log for details." >&2
    exit 1
  fi
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
  ../build/gemmi mondiff --only=at "$acedrg_file" "$gemmi_file" ||:
done
echo "Script ended successfully!"
