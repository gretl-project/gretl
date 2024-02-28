#!/bin/bash

# Collection of helper functions needed for shell scripts executing test files.
# This file will be sourced by ./run_<TYPE>.sh files

function log_msg_start {
	echo ""
	echo "*******************************************************************"
    echo "INFO: Start time executing $1 for functions: $(date '+ %F %T')"
 }

function log_msg_end {
	echo ""
	echo "*******************************************************************"
    echo "INFO: Finished executing $1 for functions: $(date '+ %F %T')"
}

function delete_file {
	find . -iname $1 -exec rm -rf {} \;
	echo "INFO: Removed file: $1."
}

function count_number_of_inp_files {
	N=$(ls -f *.inp | wc -l)

	if [ "$N" = 0 ]; then
    	printf "INFO: No files found. Finish.\n"
    	exit 0
	fi

	printf "INFO: Found %d files in current directory for computation.\n" "$N"
}

function execute_inp_files {
	# Loop over inp-files in current directory and execute.

	failed_files=()
	error_code=0

	for f in *.inp; do
		echo "$(date '+ %F %T')"
		echo `basename $f`

    	gretlcli -b -e -q $f #> /dev/null 2>&1

		if [ "$?" -ne 0 ]; then
			# Print 'Error', update status variable, and log the failed script
			echo -e " [\e[0;31mERROR\e[0m]"
			error_code=1
			failed_files+=("$f")
		else
			# Print 'OK' if the script succeeded
			echo -e " [\e[0;32mOK\e[0m]"
		fi
    	echo "--------------------------------------"
	done
}





