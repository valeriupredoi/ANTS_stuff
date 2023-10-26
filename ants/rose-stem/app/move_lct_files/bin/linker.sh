#!/bin/bash -l
# Convenience wrapper for using in rose bunch to symmlink various lct files
set -eux
LINKME=$1

LISTINGS=`ls ${ROSE_DATA}/${LINKME}`

for filename in ${LISTINGS}; do 
    ln -sf ${ROSE_DATA}/${LINKME}/${filename} ${ROSE_DATA}/${LINKME}_${filename}
done
