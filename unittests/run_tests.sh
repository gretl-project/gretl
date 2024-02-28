#!/bin/bash

# Store filenames of scripts which error
fails=()

# Store the current directory path
HERE=`pwd`

# Initialize status variable
exitcode=0

showHelp() {
# `cat << EOF` This means that cat should stop reading when EOF is detected
cat << EOF
Usage: ./run_tests.sh --OPTION

--help           Display help

--practice       Run 'practice' set of tests

--commands       Run set of tests for gretl commands

--functions      Run set of tests for gretl functions

--fundamentals   Run set of tests testing fundamentals

--plots          Run set of tests testing plots

--all            Run all scripts

EOF
# EOF is found above and hence cat command stops reading.
# This is equivalent to echo but much neater when printing out.
}

# Function to execute general tasks
execGeneral() {
	# Store current directory
    cwd=$(pwd)

    # Change to the specified directory or exit if it fails
    cd $1 || exit

	# Execute functions; $1 refers to <DIRECTORY> and $2 to <TYPE>
	log_msg_start $2
	delete_file session.inp
	delete_file string_table.txt
	count_number_of_inp_files
	execute_inp_files $2
	log_msg_end $2

	# Return to the original directory or exit if it fails
	cd $cwd || exit
}

# Function to execute practice tests
execPractice() {
    DIRECTORY="./practice_scripts"
    TYPE="practice_scripts"

    # Source the helper script
    source "./helper.sh"

    # Execute general tasks
    execGeneral $DIRECTORY $TYPE
}

# Similar functions for commands, functions, fundamentals, and plots tests
execCommands() {
    DIRECTORY="./test_scripts/commands"
    TYPE="commands"

    source "./helper.sh"
    execGeneral $DIRECTORY $TYPE
}

execFunctions() {
    DIRECTORY="./test_scripts/functions"
    TYPE="unit-tests"

    source "./helper.sh"
    execGeneral $DIRECTORY $TYPE
}

execFundamentals() {
    DIRECTORY="./test_scripts/fundamentals"
    TYPE="fundamentals"

    source "./helper.sh"
    execGeneral $DIRECTORY $TYPE
}

execPlots() {
    DIRECTORY="./test_scripts/plots"
    TYPE="plots"

    source "./helper.sh"
    execGeneral $DIRECTORY $TYPE
}

# Main script execution starts here
echo "$1"
if [[ "$1" == "--help" ]] ; then
    showHelp
elif [[ "$1" == "--practice" ]] ; then
    execPractice
    exitcode=$error_code
    fails+=("$failed_files")
elif [[ "$1" == "--commands" ]] ; then
    execCommands
    exitcode=$error_code
    fails+=("$failed_files")
elif [[ "$1" == "--functions" ]] ; then
    execFunctions
    exitcode=$error_code
    fails+=("$failed_files")
elif [[ "$1" == "--fundamentals" ]] ; then
    execFundamentals
    exitcode=$error_code
    fails+=("$failed_files")
elif [[ "$1" == "--plots" ]] ; then
    execPlots
    exitcode=$error_code
    fails+=("$failed_files")
elif [[ "$1" == "--all" ]] ; then
    execFundamentals
    exitcode=$((exitcode + $error_code))
    fails+=("$failed_files")

    execFunctions
    exitcode=$((exitcode + $error_code))
    fails+=("$failed_files")

	execCommands
	exitcode=$((exitcode + $error_code))
	fails+=("$failed_files")

	execPractice
	exitcode=$((exitcode + $error_code))
	fails+=("$failed_files")

	execPlots
	exitcode=$((exitcode + $error_code))
	fails+=("$failed_files")
else
    showHelp
    exit 1
fi

# Print summary of failed files
echo -n "Number of failed unit-tests = $exitcode"

printf "\n=================================================\n"
if [ "$exitcode" -gt "0" ]; then
	printf "Summary: Scripts which failed for $1:\n"
	echo "------------------------------"

	printf 'File: \t ./%s\n' "${fails[@]}"
else
	printf "*** INFO: No errors found ***\n"
	# Exit with success status
fi
printf "=================================================\n"

# Exit with error code
exit $exitcode

