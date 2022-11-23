#!/bin/bash
# Execute command and compare wall time to a reference
set -e

SCRIPT_PATH=$(dirname "${BASH_SOURCE}")
SCRIPT_PATH=$(cd "${SCRIPT_PATH}" && pwd)

# YAML file with reference wall times
REFERENCE_FILE="${SCRIPT_PATH}"/reference.yml

# Reference key
REF="$1"
# Command to time
CMD="${@:2}"
# Logfile
LOG=$(mktemp)

# Run and time
/usr/bin/time -f %e ${CMD} 2>&1 | tee "$LOG"
RET=${PIPESTATUS[0]}
test "$RET" != 0 && exit "$RET"
TIME=$(tail -n1 "$LOG")

# Read reference
REF_TIME=$(shyaml get-value $REF < $REFERENCE_FILE)

if (( $(echo "$TIME > $REF_TIME" | bc) )); then
  echo "TOO LONG: $TIME > $REF_TIME"
  exit 1
else
  echo "OKAY: $TIME < $REF_TIME"
fi
