#!/usr/bin/env bash

# Generate a '.deb' installer for SMARDDA.
# This script requires the 'fakeroot' package to be installed:
#  'apt install fakeroot'

set -o errexit  # exit early on any error

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

OUT_DIR="."
VERSION=""
while getopts "o:v:" opt; do
  case "$opt" in
    o) OUT_DIR="${OPTARG}"
      ;;
    v) VERSION="${OPTARG}"
      ;;
    *) echo "make-deb: vnvalid flag '${opt}'."; exit 1
  esac
done

if [ -z "${VERSION}" ]; then
    echo "make-deb: you must specify a package version using the '-v' option"
    exit 1
fi
if [ ! -d "${OUT_DIR}" ]; then
    echo "make-deb: no such directory '${OUT_DIR}' to write deb to"
    exit 1
fi

NAME="smardda"
ARCH="$(dpkg --print-architecture)"
PACKAGING_DIR="/tmp"
AUTHOR="UKAEA"
ISSUES_URL="https://git.ccfe.ac.uk/smardda/smardda-pfc-new"

PACKAGING_ROOT="${PACKAGING_DIR}/${NAME}_${VERSION}_${ARCH}"
mkdir -p "${PACKAGING_ROOT}"
INSTALL_SCRIPT="$(realpath "${SCRIPT_DIR}"/../exec/smardda-install)"
${INSTALL_SCRIPT} --prefix "${PACKAGING_ROOT}"

# Build the control file
DEBIAN_DIR="${PACKAGING_ROOT}/DEBIAN"
mkdir -p "${DEBIAN_DIR}"
echo "Package: ${NAME}
Version: ${VERSION}
Maintainer: ${AUTHOR}
Section: science, research
Priority: optional
Depends: libgfortran5
Architecture: ${ARCH}
Bugs: ${ISSUES_URL}
Description: ${NAME}
 Fieldline tracing kernel of SMITER = SMARDDA for ITER" > "${DEBIAN_DIR}/control"
chmod 644 "${DEBIAN_DIR}/control"

# Make deb
mkdir -p "${OUT_DIR}"
fakeroot dpkg-deb --build "${PACKAGING_ROOT}" "${OUT_DIR}"
