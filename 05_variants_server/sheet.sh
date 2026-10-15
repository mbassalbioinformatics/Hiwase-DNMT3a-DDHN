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
# Needs bash >= 4.4.
###############################################################################
_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ ! -s "$_ROOT/config/paths.sh" ]]; then echo "ERROR: $_ROOT/config/paths.sh missing (copy config/paths.example.sh)"; exit 1; fi
source "$_ROOT/config/paths.sh"
SHEET=${SHEET:-$_ROOT/config/samples.tsv}
if [[ ! -s "$SHEET" ]]; then echo "ERROR: $SHEET missing (copy config/samples.example.tsv)"; exit 1; fi
mapfile -t _ROWS < <(awk -F'\t' '
  function v(x) { return x == "" ? "-" : x }
  NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i
            n = split("sample case person methylation variant_role cnv_panel", need, " ")
            for (k = 1; k <= n; k++) if (!(need[k] in h)) { print "MISSING " need[k]; exit 1 }
            next }
  $h["methylation"] == "analysed" { print v($h["sample"]), v($h["case"]), v($h["person"]), v($h["variant_role"]), v($h["cnv_panel"]) }' "$SHEET")
if [[ "${_ROWS[0]:-}" == MISSING* ]]; then echo "ERROR: $SHEET has no column ${_ROWS[0]#MISSING }"; exit 1; fi
if (( ${#_ROWS[@]} == 0 )); then echo "ERROR: no analysed library in $SHEET"; exit 1; fi
ALL=(); CASES=(); CNV_PANEL=()
declare -gA M2_NORMALS=() M2_TUMOURS=() MATCHED=() _FIRST=()
for _r in "${_ROWS[@]}"; do
  read -r _s _c _p _role _panel <<< "$_r"
  ALL+=("$_s")
  if [[ "$_role" == normal || "$_role" == tumour ]] && [[ " ${CASES[*]:-} " != *" $_c "* ]]; then CASES+=("$_c"); fi
  case "$_role" in
    normal) M2_NORMALS[$_c]+="${M2_NORMALS[$_c]:+ }$_s"
            if [[ -z "${_FIRST[$_p]:-}" ]]; then _FIRST[$_p]=$_s; fi ;;
    tumour) M2_TUMOURS[$_c]+="${M2_TUMOURS[$_c]:+ }$_s" ;;
  esac
  if [[ "$_panel" == yes ]]; then CNV_PANEL+=("$_s"); fi
done
for _r in "${_ROWS[@]}"; do
  read -r _s _c _p _role _panel <<< "$_r"
  if [[ -n "${_FIRST[$_p]:-}" && "${_FIRST[$_p]}" != "$_s" ]]; then MATCHED[$_s]=${_FIRST[$_p]}; fi
done
unset _r _s _c _p _role _panel _ROWS _FIRST _ROOT
