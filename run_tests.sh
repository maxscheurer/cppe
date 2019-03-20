#!/bin/bash

PART="$1"
shift

if [ "$PART" == "pytest" ]; then
	python3 -m pytest "$@"
else
	echo "Unknown part: $PART, known are: testdata pytest gdb_pytest valgrind_pytest" >&2
	exit 1
fi
