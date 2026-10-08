#!/usr/bin/env bash
###############################################################################
# 00_tool_versions.sh                                                 (server)
#
# The version of every tool the analysis used on the server, for the methods.
# Each tool is looked for on PATH and then in every conda environment, so
# nothing needs activating first. One line per place a tool was found:
#   tool <tab> version <tab> conda env (or PATH) <tab> path <tab> first line it printed
#
# USAGE (from the repository root, on the server):
#   bash 01_methylation_server/00_tool_versions.sh | tee tool_versions_server.tsv
###############################################################################
set -uo pipefail
[[ -f config/paths.sh ]] && source config/paths.sh
MSUITE2=${MSUITE2:-msuite2}
TOOLS=(msuite2 bowtie2 hisat2 samtools bedtools metilene defiant findMotifsGenome.pl annotatePeaks.pl
       featureCounts rmats.py shiba.py STAR python3 R xz)

first_line() { tr -d '\r' | grep -av '^[[:space:]]*$' | head -1 | cut -c1-120; }
to() { if command -v timeout >/dev/null 2>&1; then timeout 20 "$@"; else "$@"; fi; }   # no tool may hang the script

version_of() {                      # $1 tool, $2 path
  local t=$1 p=$2 v=""
  case "$t" in
    msuite2)       v=$(grep -a -m1 -oiE 'version[^0-9]{0,12}[0-9]+\.[0-9]+(\.[0-9]+)?' "$p" 2>/dev/null | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?')
                   [[ -z "$v" ]] && v=$(to "$p" -v 2>&1 | grep -aoE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1) ;;
    bowtie2|hisat2|STAR) v=$(to "$p" --version 2>&1 | grep -aoE '[0-9]+\.[0-9]+(\.[0-9]+[a-z]?)?' | head -1) ;;
    samtools)      v=$(to "$p" --version 2>&1 | head -1 | awk '{print $2}') ;;
    bedtools)      v=$(to "$p" --version 2>&1 | awk '{print $2}') ;;
    metilene)      v=$(to "$p" --version 2>&1; to "$p" -h 2>&1; to "$p" 2>&1) ; v=$(grep -aoiE '(version|v)[: ]*[0-9]+\.[0-9]+[-.0-9]*' <<< "$v" | head -1)
                   # the binary prints no version: look for one compiled into it, or in its folder name
                   [[ -z "$v" ]] && v=$(strings "$p" 2>/dev/null | grep -m1 -oiE 'metilene[ _-]*v?[0-9]+\.[0-9]+[-.0-9]*')
                   [[ -z "$v" ]] && v=$(readlink -f "$p" | grep -oiE 'v?[0-9]+\.[0-9]+[-.0-9]*' | tail -1) ;;
    defiant)       v=$(to "$p" -h 2>&1; to "$p" --version 2>&1; to "$p" -v 2>&1) ; v=$(grep -aoiE '(version|v)[: ]*[0-9]+\.[0-9]+[.0-9]*' <<< "$v" | head -1)
                   [[ -z "$v" ]] && v=$(strings "$p" 2>/dev/null | grep -m1 -oiE '(defiant|version)[ _:-]*v?[0-9]+\.[0-9]+[.0-9]*')
                   [[ -z "$v" ]] && v=$(readlink -f "$p" | grep -oiE 'v?[0-9]+\.[0-9]+[.0-9]*' | tail -1) ;;
    findMotifsGenome.pl|annotatePeaks.pl)
                   local d; d=$(dirname "$(readlink -f "$p")")
                   v=$(awk -F'\t' '$1 == "homer" {print $2; exit}' "$d/../config.txt" "$d/config.txt" 2>/dev/null) ;;
    featureCounts) v=$(to "$p" -v 2>&1 | grep -aoE 'v[0-9.]+' | head -1) ;;
    rmats.py)      v=$(to "$p" --version 2>&1 | first_line) ;;
    shiba.py)      v=$(cat "$(dirname "$(readlink -f "$p")")/VERSION" 2>/dev/null) ;;
    python3)       v=$("$p" --version 2>&1 | awk '{print $2}') ;;
    R)             v=$("$p" --version 2>&1 | head -1 | awk '{print $3}') ;;
    xz)            v=$("$p" --version 2>&1 | head -1 | awk '{print $NF}') ;;
  esac
  echo "${v:-not parsed}"
}

raw_of() {                          # the first line a tool prints, to check the parsed version by eye
  case "$1" in
    findMotifsGenome.pl|annotatePeaks.pl|shiba.py|msuite2) echo "-" ;;
    defiant)  to "$2" -h 2>&1 | first_line ;;
    metilene) to "$2" 2>&1 | first_line ;;
    *)        to "$2" --version 2>&1 | first_line ;;
  esac
}

report() {                          # $1 tool, $2 path, $3 where
  printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$(version_of "$1" "$2")" "$3" "$2" "$(raw_of "$1" "$2")"
}

printf 'tool\tversion\tfound_in\tpath\tfirst_line\n'
printf 'host\t%s\t-\t-\t%s\n' "$(hostname)" "$(date '+%F %T %Z')"
ENVS=()
if command -v conda >/dev/null 2>&1; then
  while IFS= read -r e; do ENVS+=("$e"); done < <(conda env list 2>/dev/null | awk '$1 !~ /^#/ && NF {print $NF}')
fi
for t in "${TOOLS[@]}"; do
  found=0
  if [[ $t == msuite2 && -x "$MSUITE2" ]]; then report "$t" "$MSUITE2" "fixed path"; found=1; fi
  if p=$(command -v "$t" 2>/dev/null); then report "$t" "$p" "PATH"; found=1; fi
  for e in ${ENVS[@]+"${ENVS[@]}"}; do
    if [[ -x "$e/bin/$t" && "$e/bin/$t" != "${p:-}" ]]; then report "$t" "$e/bin/$t" "conda:$(basename "$e")"; found=1; fi
  done
  (( found )) || printf '%s\tnot found\t-\t-\t-\n' "$t"
done
