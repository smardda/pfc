#! /bin/bash
# diff all testdeck outputs from SMARDDA/PFC test-deck
#
# This scripts needs to know the location of two directories: the path to
# 'smiteraux1' (HSAUX; where the 'tests.bash' script has been run), and the path
# to the test-deck reference outputs (HSAOLD).
# The first argument to this script sets HSAUX, the second sets HSAOLD.
# If no arguments are passed, HSAUX is set to '../smiteraux1' (relative to this
# script's directory), and HSAOLD is set to '../smiteraux-results/smiteraux1'.

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

# Get the smiteraux test directory
if [ -n "$1" ]; then
	HSAUX="$1"
elif [ -z "${HSAUX}" ]; then
	HSAUX="$(dirname "${SCRIPT_DIR}")/smiteraux"
fi
# Get the smiteraux reference directory
if [ -n "$2" ]; then
	HSAOLD="$2"
elif [ -z "${HSAOLD}" ]; then
	HSAOLD="$(dirname "${HSAUX}")/smiteraux-results/smiteraux1"
fi

if [ ! -d "${HSAUX}" ]; then
	echo "smiteraux directory \$HSAUX=${HSAUX} does not exist, try different \$HSAUX "
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
	elif [ "${file_extension}" = "out" ]; then
		"${SCRIPT_DIR}/compare_out" "$1" "$2"
	elif [ "${file_extension}" = "hds" ]; then
		"${SCRIPT_DIR}/compare_hds" "$1" "$2"
	elif [ "${file_extension}" = "ps" ]; then
		# For .ps files, just check the corresponding file is generated.
		# The file content depends on the version of gnuplot being used, so
		# diff-ing the files is not practical.
		[ -f "$1" ] #
	else
		diff -bu "$1" "$2" >/dev/null
	fi
}

FILES=(
	"Data/Geometry/beancan.log"
	"Data/Geometry/beancansy.log"
	"Data/Geometry/dome_outrefl_skeleton_july2012.log"
	"Data/Geometry/dome_outrefl_skeleton_july2012sy.log"
	"Data/Geometry/fw18.log"
	"Data/Geometry/fw18sy.log"
	"Data/Geometry/FW2.log"
	"Data/Geometry/FW3.log"
	"Data/Geometry/FW4.log"
	"Data/Geometry/FW5.log"
	"Data/Geometry/FW6.log"
	"Data/Geometry/.STOP-F11-fw264t-fw4t"
	"Data/Geometry/.STOP-vde-bot-baf"
	"Data/Geometry/.STOP-vde-botsy-baf"
	"Data/Geometry/wdfd03_shaped_mesh10mm.log"
	"Data/Geometry/wdpdd06_mesh10mm_onlyfront.log"
	"Data/Geometry/wdpdd06_mesh30mm.log"
	"Data/VTK/baf.log"
	"Data/VTK/baf.vtk"
	"Data/VTK/beancansy.vtk"
	"Data/VTK/beancan.vtk"
	"Data/VTK/bot.log"
	"Data/VTK/botsy.log"
	"Data/VTK/botsy.vtk"
	"Data/VTK/bot.vtk"
	"Data/VTK/dome_outrefl_skeleton_july2012sy.vtk"
	"Data/VTK/dome_outrefl_skeleton_july2012.vtk"
	"Data/VTK/fw18sy.vtk"
	"Data/VTK/fw18.vtk"
	"Data/VTK/fw264t.log"
	"Data/VTK/fw264t.vtk"
	"Data/VTK/FW2.vtk"
	"Data/VTK/FW3.vtk"
	"Data/VTK/fw4t.log"
	"Data/VTK/fw4t.vtk"
	"Data/VTK/FW4.vtk"
	"Data/VTK/FW5.vtk"
	"Data/VTK/FW6.vtk"
	"Data/VTK/.STOP-F11-fw264t-fw4t"
	"Data/VTK/.STOP-vde-bot-baf"
	"Data/VTK/.STOP-vde-botsy-baf"
	"Data/VTK/wdfd03_shaped_mesh10mm.vtk"
	"Data/VTK/wdpdd06_mesh10mm_onlyfront.vtk"
	"Data/VTK/wdpdd06_mesh30mm_tfm1.log"
	"Data/VTK/wdpdd06_mesh30mm_tfm1.vtk"
	"Data/VTK/wdpdd06_mesh30mm.vtk"
)
# Recursively check all files in these directories
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
for dir in "${DIRS[@]}"; do
	found_files=$(cd "${HSAOLD}" && find "${dir}" -type f)
	for f in ${found_files}; do
		FILES+=("${f}")
	done
done

diff_count=0
VTK_FILES=()

for file in "${FILES[@]}"; do
	echo "${file}"

	if [[ "${file}" == *.vtk ]]; then
		VTK_FILES+=("${file}")
	fi

	if ! do_diff "${HSAUX}/${file}" "${HSAOLD}/${file}"; then
		diff_count=$((diff_count + 1))
		diff -q "${HSAUX}/${file}" "${HSAOLD}/${file}"
	fi
done
echo "Compared ${#FILES[@]} files. Found ${diff_count} with differences."

if [ "${diff_count}" -ne 0 ]; then
	echo
	echo "Detailed VTK comparison report:"
	python3 "$(dirname "${SCRIPT_DIR}")/tools/vtk_folder_compare.py" \
	"${HSAOLD}" "${HSAUX}" \
	--files "${VTK_FILES[@]}"
fi

exit ${diff_count}
