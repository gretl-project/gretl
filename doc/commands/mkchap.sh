#!/bin/sh
cat commands.xml | fgrep -v '<?x' >> ../chapters/cmdlist.xml
