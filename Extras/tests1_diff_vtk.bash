#! /bin/bash
# Compare with tolerances selected VTK outputs from SMARDDA/PFC Test Deck 1
# against their regression reference files.
#
# This scripts needs to know the location of two folder:
#   HSAUX  - the smiteraux1 directory where tests.bash has been runtwo directories
#   HSAOLD - the directory containing the test-deck reference VTK files.
# 
# The first argument to this script sets HSAUX.
# The second argument to this script optionally sets HSAOLD.
#
# Defaults:
#  HSAUX  = ../smiteraux1' (relative to thisscript) 
#  HSAOLD  = $HSAUX/References

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

# Get the smiteraux test directory
if [ -n "$1" ]; then
	HSAUX="$1"
elif [ -z "${HSAUX}" ]; then
	HSAUX="$(dirname "${SCRIPT_DIR}")/smiteraux1"
fi

# Get the smiteraux reference directory
if [ -n "$2" ]; then
	HSAOLD="$2"
elif [ -z "${HSAOLD}" ]; then
	HSAOLD="${HSAUX}/References"
fi

if [ ! -d "${HSAUX}" ]; then
	echo "Test directory \$HSAUX=${HSAUX} does not exist, try different \$HSAUX "
	echo "Script quitting"
	exit 1
fi

if [ ! -d "${HSAOLD}" ]; then
	echo "Reference directory \$HSAOLD=${HSAOLD} does not exist, try different \$HSAOLD "
	echo "Script quitting"
	exit 1
fi

do_diff() {
	file_extension="${2##*.}"
	if [ "${file_extension}" = "log" ]; then
		"${SCRIPT_DIR}/compare_logs" "$1" "$2"
	else
		diff -bu "$1" "$2" >/dev/null
	fi
}

# Recursively check all VTK files in these directories
DIRS=(
	"Test-EQ3-inrshad1-inres1"
	"Test-EQX40-shado14r-resf14r"
	"Test-EQX40-shado14r-resf14r-3g"
	"Test-MNEW8-sha12-res12"
	"Val-F11-fw264t-fw4t"
	"Val-F11-fw264t-fw4t-1l"
	"Val-vde-bot-baf"
	"Val-vde-botsy-baf"
	"Val-HR-bott-baff"
	"Val-smanal-mjg"
)

VTK_FILES=()

for dir in "${DIRS[@]}"; do

	found_files=$(cd "${HSAOLD}" && find "${dir}" -type f -name "*.vtk")

	for f in ${found_files}; do
		VTK_FILES+=("${f}")
	done
	
done


echo
echo "Detailed VTK comparison report:"

python3 "$(dirname "${SCRIPT_DIR}")/tools/vtk_folder_compare.py" \
    "${HSAOLD}" "${HSAUX}" \
	--files "${VTK_FILES[@]}"

exit $?
