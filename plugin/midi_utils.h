#ifndef UTIL_H
#define UTIL_H

int write_var_len (long val, FILE *fp);
int delta_time (int beat, int nticks, FILE *fp);
int write_be_long (unsigned long val, FILE *fp);
int write_be_24 (int val, FILE *fp);
int write_be_short (unsigned short val, FILE *fp);
int write_midi_byte (unsigned char c, FILE *fp);
int write_midi_meta (int val, FILE *fp);
int write_midi_eot (FILE *fp);
void write_midi_header (int format, int tracks, int ticks,
			FILE *fp);

const char *get_patch_name (int pnum);

#endif
