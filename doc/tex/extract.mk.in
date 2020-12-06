CC = gcc -Wall -W -O2

scriptbits: extract_scripts
	./extract_scripts > $@

extract_scripts: extract_scripts.c
	$(CC) -o $@ $<
