#!/bin/sh

recode utf8..windows-1252 < toc.hhc > toc.ansi && \
  mv toc.ansi toc.hhc
