#! /bin/bash
# Compare with tolerances selected VTK outputs from SMARDDA/PFC Test Deck 2
# against their regression reference files.
#
# This scripts needs to know the location of two folders:
#   HSAUX  - the smiteraux2 directory where tests.bash has been run.
#   HSAOLD - the directory containing the test-deck reference VTK files.
# 
# The first argument to this script optionally sets HSAUX.
# The second argument to this script optionally sets HSAOLD.
#
# Defaults:
#  HSAUX  = ../smiteraux2' (relative to this script) 
#  HSAOLD  = $HSAUX/References

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

# Get the smiteraux2 test directory
if [ -n "$1" ]; then
	HSAUX="$1"
elif [ -z "${HSAUX}" ]; then
	HSAUX="$(dirname "${SCRIPT_DIR}")/smiteraux2"
fi

# Get the smiteraux2 reference directory
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

# Recursively check all VTK files in these directories
DIRS=(
	"DEMO-double-null"
	"ITER-valid"
	"JET-ant"
	"JET-radn"
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
