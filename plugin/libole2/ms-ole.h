/**
 * ms-ole.h: MS Office OLE support for Gnumeric
 *
 * Authors:
 *    Michael Meeks (michael@imaginator.com)
 *    Arturo Tena (arturo@directmail.org)
 *
 * Copyright 1998-2001 Ximian, Inc., Arturo Tena
 *
 * Adapted for gretl by Allin Cottrell
 *
 **/

#ifndef MS_OLE_H
#define MS_OLE_H

/* Allin Cottrell modifications here */

#ifndef _WIN32
# include <fcntl.h>	/* for mode_t */
# include <sys/types.h>
#else
# include <sys/types.h>
typedef /* unsigned */ long caddr_t;
#endif

#include <glib.h>

typedef enum {
    MS_OLE_ERR_OK,
    MS_OLE_ERR_EXIST,
    MS_OLE_ERR_INVALID,
    MS_OLE_ERR_FORMAT,
    MS_OLE_ERR_PERM,
    MS_OLE_ERR_MEM,
    MS_OLE_ERR_SPACE,
    MS_OLE_ERR_NOTEMPTY,
    MS_OLE_ERR_BADARG
} MsOleErr;

typedef enum {
    MsOleSeekSet,
    MsOleSeekCur,
    MsOleSeekEnd
} MsOleSeek;

typedef enum  {
    MsOleStorageT = 1,
    MsOleStreamT  = 2,
    MsOleRootT    = 5
} MsOleType;

typedef guint32 MsOlePos;
typedef gint32  MsOleSPos;

typedef struct _MsOle       MsOle;
typedef struct _MsOleStat   MsOleStat;
typedef struct _MsOleStream MsOleStream;

struct _MsOleStat {
    MsOleType type;
    MsOlePos  size;
};

struct _MsOleStream {
    MsOlePos size;
    gint       (*read_copy)	(MsOleStream *stream,
				 guint8 *ptr,
				 MsOlePos length);
    guint8 *	(*read_ptr)	(MsOleStream *stream,
				 MsOlePos length);
    MsOleSPos	(*lseek)	(MsOleStream *stream,
				 MsOleSPos bytes,
				 MsOleSeek type);
    MsOlePos	(*tell)		(MsOleStream *stream);
    MsOlePos	(*write)	(MsOleStream *stream,
				 guint8 *ptr,
				 MsOlePos length);
    /**
     * Private.
     **/
    enum {
	MsOleSmallBlock,
	MsOleLargeBlock
    } type;
    MsOle    *file;
    void     *pps;      /* Straight PPS */
    GArray   *blocks;	/* A list of the blocks in the file
                           if NULL: no file */
    MsOlePos position;	/* Current offset into file.
                           Points to the next byte to read */
};

MsOleErr ms_ole_open    (MsOle **fs, const char *path);
void     ms_ole_destroy (MsOle **fs);

#define MS_OLE_GET_GUINT8(p)  (*((const guint8 *)(p) + 0))
#define MS_OLE_GET_GUINT16(p) (guint16)(*((const guint8 *)(p)+0) |        \
					(*((const guint8 *)(p)+1)<<8))
#define MS_OLE_GET_GUINT32(p) (guint32)(*((const guint8 *)(p)+0) |        \
					(*((const guint8 *)(p)+1)<<8) |   \
					(*((const guint8 *)(p)+2)<<16) |  \
					(*((const guint8 *)(p)+3)<<24))
#define MS_OLE_GET_GUINT64(p) (MS_OLE_GET_GUINT32(p) | \
			       (((guint32)MS_OLE_GET_GUINT32((const guint8 *)(p)+4))<<32))

MsOleErr ms_ole_stream_open_workbook (MsOleStream ** const stream,
				      MsOle *fs);

MsOleErr ms_ole_stream_close (MsOleStream ** const stream);

MsOlePos ms_ole_stream_position (const MsOleStream *s);

#endif	/* MS_OLE_H */
