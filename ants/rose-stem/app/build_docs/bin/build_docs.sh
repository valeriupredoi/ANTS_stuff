#!/bin/bash -l
set -eu
echo "Building the documentation"
cd ${DOCSDIR}
make clean
make clean-apidoc
make apidoc
make html
echo "Documentation built. "
echo "Run 'gio open ${DOCSDIR}/build/html/index.html' to view"
