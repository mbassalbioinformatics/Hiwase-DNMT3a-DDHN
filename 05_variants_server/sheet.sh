#!/usr/bin/env bash
###############################################################################
# sheet.sh                         sourced by the scripts in 05_variants_server/
#
# Reads config/paths.sh and the sample sheet (config/samples.tsv) and defines:
#   ALL                    every analysed library (methylation = analysed), sheet order
#   CASES                  the cases with Mutect2 roles, sheet order
#   M2_NORMALS[case]       Mutect2 normals of a case (variant_role = normal; the first is the baseline)
#   M2_TUMOURS[case]       Mutect2 tumours of a case (variant_role = tumour)
#   CNV_PANEL              the copy-number panel of normals (cnv_panel = yes)
#   MATCHED[library]       the library whose heterozygous sites it is segmented with:
#                          the first normal of the same person (none for that normal itself)
###############################################################################
_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
[[ -s "$_ROOT/config/paths.sh" ]] || { echo "ERROR: $_ROOT/config/paths.sh missing (copy config/paths.example.sh)"; exit 1; }
source "$_ROOT/config/paths.sh"
SHEET=${SHEET:-$_ROOT/config/samples.tsv}
[[ -s "$SHEET" ]] || { echo "ERROR: $SHEET missing (copy config/samples.example.tsv)"; exit 1; }
mapfile -t _ROWS < <(awk -F'\t' '
  function v(x) { return x == "" ? "-" : x }
  NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i
            split("sample case person methylation variant_role cnv_panel", need, " ")
            for (k in need) if (!(need[k] in h)) { print "MISSING " need[k]; exit 1 }
            next }
  $h["methylation"] == "analysed" { print v($h["sample"]), v($h["case"]), v($h["person"]), v($h["variant_role"]), v($h["cnv_panel"]) }' "$SHEET")
[[ "${_ROWS[0]:-}" == MISSING* ]] && { echo "ERROR: $SHEET has no column ${_ROWS[0]#MISSING }"; exit 1; }
ALL=(); CASES=(); CNV_PANEL=()
declare -A M2_NORMALS=() M2_TUMOURS=() MATCHED=() _FIRST=()
for _r in "${_ROWS[@]}"; do
  read -r _s _c _p _role _panel <<< "$_r"
  ALL+=("$_s")
  if [[ "$_role" == normal || "$_role" == tumour ]] && [[ " ${CASES[*]} " != *" $_c "* ]]; then CASES+=("$_c"); fi
  case "$_role" in
    normal) M2_NORMALS[$_c]+="${M2_NORMALS[$_c]:+ }$_s"; [[ -n "${_FIRST[$_p]:-}" ]] || _FIRST[$_p]=$_s ;;
    tumour) M2_TUMOURS[$_c]+="${M2_TUMOURS[$_c]:+ }$_s" ;;
  esac
  [[ "$_panel" == yes ]] && CNV_PANEL+=("$_s")
done
for _r in "${_ROWS[@]}"; do
  read -r _s _c _p _role _panel <<< "$_r"
  [[ -n "${_FIRST[$_p]:-}" && "${_FIRST[$_p]}" != "$_s" ]] && MATCHED[$_s]=${_FIRST[$_p]}
done
(( ${#ALL[@]} )) || { echo "ERROR: no analysed library in $SHEET"; exit 1; }
unset _r _s _c _p _role _panel _ROWS _FIRST _ROOT
