/* MIDI utility functions for gretl audio plugin */

#include <stdio.h>

#include "miditypes.h"
#include "midi_utils.h"

int write_var_len (long val, FILE *fp)
{
    long buf = val & 0x7f;
    int wrote = 0;

    while ((val >>= 7) > 0) {
	buf <<= 8;
	buf |= 0x80;
	buf += (val & 0x7f);
    }

    while (1) {
	putc(buf, fp);
	wrote++;
	if (buf & 0x80) buf >>= 8;
	else break;
    }

    return wrote;
}

int delta_time (double beat, int nticks, FILE *fp)
{
    double dval = beat * nticks;
    long val = (long) dval;

    return write_var_len(val, fp); 
}

int write_be_long (unsigned long val, FILE *fp)
{
    unsigned char c[4];

    c[0] = (val >> 24) & 0xFFL; 
    c[1] = (val >> 16) & 0xFFL;
    c[2] = (val >> 8) & 0xFFL;
    c[3] = val & 0xFFL;
    
    return fwrite(c, 1, 4, fp);
}

int write_be_24 (int val, FILE *fp)
{
    unsigned char c[3];

    c[0] = (val >> 16) & 0xFFL;
    c[1] = (val >> 8) & 0xFFL;
    c[2] = val & 0xFFL;
    
    return fwrite(c, 1, 3, fp);
}

int write_be_short (unsigned short val, FILE *fp)
{
    unsigned char c[2];

    c[0] = (val >> 8) & 0xFF;
    c[1] = val & 0xFF;

    return fwrite(c, 1, 2, fp);
}

int write_midi_byte (unsigned char c, FILE *fp) 
{
    putc(c, fp);

    return 1;
}

int write_midi_meta (int val, FILE *fp)
{
    unsigned char c = val;

    putc(MIDI_META, fp);
    putc(c, fp);

    return 2;
}

int write_midi_eot (FILE *fp)
{
    unsigned char eot[] = { 0xff, 0x2f, 0x00 };

    return fwrite(eot, 1, 3, fp);
}

void write_midi_header (int format, int tracks, int ticks,
			FILE *fp)
{
    const char *hdr = "MThd";

    fwrite(hdr, 1, 4, fp);
    write_be_long(6, fp);
    write_be_short(format, fp);
    write_be_short(tracks, fp);
    write_be_short(ticks, fp);
}

const char *get_patch_name (int pnum)
{
    switch (pnum) {
    case PC_GRAND:
	return "Acoustic Grand Piano";
    case PC_ELEC_PIANO_2:
	return "Electric Piano 2";
    case PC_HARPSICHORD:
	return "Harpsichord";
    case PC_MUSIC_BOX:
	return "Music Box";
    case PC_HARP:
	return "Harp";
    case PC_TRUMPET:
	return "Trumpet";
    case PC_ALTO_SAX:
	return "Alto Sax";
    case PC_BASSOON:
	return "Bassoon";
    case PC_MARIMBA:
	return "Marimba";
    case PC_CELLO:
	return "Cello";
    default:
	return "Unknown";
    }
}


