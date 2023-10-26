#!/usr/bin/env bash
set -eu

# CYLC7 compatible suite log dumping script for trac template completion.

export DB_LOCATION=${CYLC_SUITE_RUN_DIR}/log/db
export OUTPUT_LOG=${CYLC_SUITE_RUN_DIR}/trac.log

# User info
date > $OUTPUT_LOG
echo "-----" >> $OUTPUT_LOG
echo $USER ran suite located at: $CYLC_SUITE_RUN_DIR >> $OUTPUT_LOG
echo "-----" >> $OUTPUT_LOG

# Verson info
echo "=== Revision information ===" >> $OUTPUT_LOG

echo "{{{" >> $OUTPUT_LOG
cat ${CYLC_SUITE_RUN_DIR}/log/rose-suite-run.version >> $OUTPUT_LOG
echo "}}}" >> $OUTPUT_LOG

# Generate test status table
echo "=== Test Results ===" >> $OUTPUT_LOG
echo " || '''task''' || '''status''' || " >> $OUTPUT_LOG
sqlite3 -separator " || " $DB_LOCATION "select '', name, status, '' from task_states" >> $OUTPUT_LOG
