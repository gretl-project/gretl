#!/bin/bash

# Remove the 'fails' file if it exists
rm -f fails

# Store the current directory path
HERE=`pwd`

# Initialize status variable
my_status=0

# Start testing scripts in the 'practice_scripts' directory
echo "*** practice scripts ***"
for f in `find ./practice_scripts -name "*.inp"` ; do
   # Print the name of the script being tested
   echo -n `basename $f`

   # Run the script with gretlcli in batch and quiet mode
   gretlcli -b -q -e $f > /dev/null 2>&1

   # Check if the script failed
   if [ $? != 0 ] ; then
      # Print 'Failed', update status variable, and log the failed script
      echo " [\e[0;31mFailed\e[0m]"
      my_status=1
      echo $f >> $HERE/fails
   else
      # Print 'OK' if the script succeeded
      echo -e " [\e[0;32mOK\e[0m]"
   fi
done

# Start testing scripts in the 'commands', 'functions', and 'fundamentals' directories
for d in commands functions fundamentals ; do
   echo "*** $d ***"
   cd ./test_scripts/$d
   for f in `find . -name "*.inp"` ; do
      echo -n `basename $f`
      gretlcli -b -q -e $f > /dev/null 2>&1
      if [ $? != 0 ] ; then
           echo -e " [\e[0;31mFailed\e[0m]"
           my_status=1
           echo $f >> $HERE/fails
      else
           echo -e " [\e[0;32mOK\e[0m]"
      fi
   done
   # Return to the original directory
   cd $HERE
done

# If there were any failures, print the names of the failed scripts
if test -f fails ; then
   echo "Failed script(s):"
   cat fails
fi

# Exit with the status code. 0 if all scripts passed, 1 if any script failed.
exit $my_status
