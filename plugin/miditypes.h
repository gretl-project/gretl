#ifndef MIDITYPES_H
#define MIDITYPES_H

#define MIDI_META            0xff

#define MIDI_SEQUENCE_NUMBER 0x00
#define MIDI_TEXT_EVENT      0x01
#define MIDI_COPYRIGHT       0x02
#define MIDI_TRACK_NAME      0x03
#define MIDI_INSTR_NAME      0x04
#define MIDI_PROGRAM_CHANGE  0xc0

#define MIDI_PAN_MSB         0x0a
#define MIDI_PAN_LSB         0x2a

#define MIDI_END_TRACK       0x2f

#define MIDI_TEMPO           0x51
#define MIDI_TIME_SIG        0x58

#define MIDI_NOTE_ON         0x90

#define PC_GRAND             0
#define PC_ELEC_PIANO_2      5
#define PC_HARPSICHORD       6
#define PC_MUSIC_BOX        10
#define PC_MARIMBA          12
#define PC_ACCORDION        21
#define PC_GUITAR           22
#define PC_GTR_OVERDRIV     29
#define PC_CELLO            42
#define PC_HARP             46
#define PC_TRUMPET          56
#define PC_ALTO_SAX         65
#define PC_BASSOON          70
#define PC_CLARINET         71

#endif /* MIDITYPES_H */
