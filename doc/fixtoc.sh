#!/bin/sh

recode -x: u8..cp1252 < toc.hhc > toc.ansi && \
  mv toc.ansi toc.hhc
