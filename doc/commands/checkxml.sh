#!/bin/sh

XMLLINT=$(which xmllint 2> /dev/null)

if [ "x$XMLLINT" = x ] ; then
  # OK, we wont insist that this check can be done
  echo "xmllint not found, check against DTD skipped"
  exit 0
fi

for lang in en it pt ; do
  if xmllint --noout --dtdvalid gretl_commands.dtd gretl_commands_${lang}.xml ; then
    echo "gretl_commands_${lang}.xml validated OK"
  else
    echo "*** gretl_commands_${lang}.xml does not validate ***"
    exit 1
  fi
done

for lang in en it pt ; do
  if xmllint --noout --dtdvalid gretl_functions.dtd gretl_functions_${lang}.xml ; then
    echo "gretl_functions_${lang}.xml validated OK"
  else
    echo "*** gretl_functions_${lang}.xml does not validate ***"
    exit 1
  fi
done

exit 0
