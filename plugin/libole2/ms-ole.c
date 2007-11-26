/**
 * ms-ole.c: MS Office OLE support for Gnumeric
 *
 * Authors:
 *    Michael Meeks (michael@ximian.com)
 *    Arturo Tena   (arturo@directmail.org)
 *
 * Copyright 1998-2000 Helix Code, Inc., Arturo Tena
 **/

#include <stdio.h>

#if defined(WIN32)
# include "winconfig.h"
#elif defined (HAVE_CONFIG_H)
# include "config.h"
#endif

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

#include <sys/stat.h>	/* for struct stat */
#include <fcntl.h>

#include <assert.h>
#include <ctype.h>
#include <glib.h>
#include <string.h>

#ifdef HAVE_MMAP
# include <sys/mman.h>
#endif

#include "ms-ole.h"

#ifndef HAVE_UNISTD_H
# define S_IRUSR 0000400
# define S_IWUSR 0000200
# define S_IRGRP 0000040
# define S_IWGRP 0000020
# define _S_ISREG(m) (((m)&0170000) == 0100000)
# define S_ISREG(m) _S_ISREG(m)
# define O_NONBLOCK 0x4000
#endif

#ifndef PROT_READ
# define PROT_READ 0x1
#endif

#ifndef PROT_WRITE
# define PROT_WRITE 0x2
#endif

#ifndef MAP_FAILED
/* Someone needs their head examining - BSD ? */
#	define MAP_FAILED ((void *)-1)
#endif

#if !defined(MAP_SHARED) || !defined(HAVE_MMAP)
/* Only define this where mmap() is not supported */
#	define MAP_SHARED 0
#endif

/* Implementational detail - not for global header */
#define OLE_DEBUG 0
#define OLE_CHAIN_DEBUG 0

/* FIXME tenix add ADD_BBD_LIST_BLOCK where it should be used) */
#define ADD_BBD_LIST_BLOCK   0xfffffffc       /* -4 */
#define SPECIAL_BLOCK        0xfffffffd       /* -3 (BBD_LIST BLOCK) */
#define END_OF_CHAIN         0xfffffffe       /* -2 */
#define UNUSED_BLOCK         0xffffffff       /* -1 */

/* FIXME tenix laola reads this from the header */
#define BB_BLOCK_SIZE     512
#define SB_BLOCK_SIZE      64

/**
 * Structure describing an OLE file
 **/
struct _MsOle {
    int               ref_count;
    gboolean          ole_mmap;
    guint8           *mem;
    guint32           length;
    MsOleSysWrappers *syswrap;

    char              mode;
    int               file_des;
    int               dirty;
    GArray           *bb;      /* Big  blocks status  */
    GArray           *sb;      /* Small block status  */
    GArray           *sbf;     /* The small block file */
    guint32           num_pps; /* Count of number of property sets */
    GList            *pps;     /* Property Storage -> struct _PPS, always 1 valid entry or NULL */
    /* if memory mapped */
    GPtrArray        *bbattr;  /* Pointers to block structures */
    /* end if memory mapped */
};


#define BLOCK_COUNT(f) (((f)->length + BB_BLOCK_SIZE - 1) / BB_BLOCK_SIZE)


/**
 * Default system calls wrappers
 **/

static int
open2_wrap (const char *pathname, int flags)
{
#ifdef O_BINARY
    return open (pathname, flags | O_BINARY);
#else
    return open (pathname, flags);
#endif
}

static int
open3_wrap (const char *pathname, int flags, mode_t mode)
{
#ifdef O_BINARY
    return open (pathname, flags | O_BINARY, mode);
#else
    return open (pathname, flags, mode);
#endif
}

static ssize_t
read_wrap (int fd, void *buf, size_t count)
{
    return read (fd, buf, count);
}

static int
close_wrap (int fd)
{
    return close (fd);
}

static ssize_t
write_wrap (int fd, const void *buf, size_t count)
{
    return write (fd, (void *)buf, count);
}

static off_t
lseek_wrap (int fd, off_t offset, int whence)
{
    return lseek (fd, offset, whence);
}

static int
isregfile_wrap (int fd)
{
    struct stat st;

    if (fstat (fd, &st))
	return 0;

    return S_ISREG (st.st_mode);
}

static int
getfilesize_wrap (int fd, guint32 *size)
{
    struct stat st;

    if (fstat (fd, &st))
	return -1;

    *size = st.st_size;
    return 0;
}

#ifdef HAVE_MMAP
static void *
mmap_wrap (void *start, size_t length, int prot,
	   int flags, int fd, off_t offset)
{
    return mmap (start, length, prot, flags, fd, offset);
}

static int
munmap_wrap (void *start, size_t length)
{
    return munmap (start, length);
}
#endif

static MsOleSysWrappers ms_ole_default_wrappers = {
    open2_wrap,
    open3_wrap,
    read_wrap,
    close_wrap,
    write_wrap,
    lseek_wrap,
    isregfile_wrap,
    getfilesize_wrap,

#ifdef HAVE_MMAP
    mmap_wrap,
    munmap_wrap
#else
    NULL,
    NULL
#endif
};

static void
take_wrapper_functions (MsOle *f, MsOleSysWrappers *wrappers)
{
    if (wrappers == NULL)
	f->syswrap = &ms_ole_default_wrappers;
    else
	f->syswrap = wrappers;
}

typedef guint32 PPS_IDX ;

#if OLE_DEBUG > 0
/* Very grim, but quite necessary */
#       define ms_array_index(a,b,c) (b)my_array_hack ((a), sizeof(b), (c))

static guint32
my_array_hack (GArray *a, guint s, guint32 idx)
{
    g_assert (a != NULL);
    g_assert (idx >= 0);
    g_assert (idx < a->len);
    g_assert (s == 4);
    return ((guint32 *)a->data)[idx];
}
#else
/* Far far faster... */
#       define ms_array_index(a,b,c) g_array_index (a, b, c)
#endif


typedef guint32 BLP;	/* Block pointer */


#define BB_THRESHOLD   0x1000

#define PPS_ROOT_INDEX    0
#define PPS_BLOCK_SIZE 0x80
#define PPS_END_OF_CHAIN 0xffffffff

typedef struct _PPS PPS;

#define PPS_SIG 0x13579753
#define IS_PPS(p) (((PPS *)(p))->sig == PPS_SIG)

struct _PPS {
    int      sig;
    char    *name;
    GList   *children;
    PPS     *parent;
    guint32  size;
    BLP      start;
    MsOleType type;
    PPS_IDX  idx; /* Only used on write */
};

#define BB_R_PTR(f,b) ((f)->ole_mmap ? ((f)->mem + ((b) + 1) * BB_BLOCK_SIZE) : \
				       (get_block_ptr (f, b, FALSE)))

#define GET_SB_R_PTR(f,b) (BB_R_PTR(f, g_array_index ((f)->sbf, BLP, (b)/(BB_BLOCK_SIZE/SB_BLOCK_SIZE))) \
			   + (((b)%(BB_BLOCK_SIZE/SB_BLOCK_SIZE))*SB_BLOCK_SIZE))

#define MAX_CACHED_BLOCKS  32

typedef struct {
    guint32  blk;
    gboolean dirty;
    int      usage;
    guint8   *data;
} BBBlkAttr;

static BBBlkAttr *
bb_blk_attr_new (guint32 blk)
{
    BBBlkAttr *attr = g_new (BBBlkAttr, 1);
    attr->blk   = blk;
    attr->dirty = FALSE;
    attr->usage = 0;
    attr->data  = 0;
    return attr;
}

static guint8 *
get_block_ptr (MsOle *f, BLP b, gboolean forwrite)
{
    BBBlkAttr *attr, *tmp, *min;
    size_t offset;
    guint32 i, blks;

    g_assert (f);
    g_assert (b < f->bbattr->len);

    /* Have we cached it ? */
    attr = g_ptr_array_index (f->bbattr, b);
    g_assert (attr);
    g_assert (attr->blk == b);

    if (attr->data) {
	attr->usage++;
	if (forwrite)
	    attr->dirty = TRUE;
	return attr->data;
    }

    /* LRU strategy */
    min  = NULL;
    blks = 0;
    for (i = 0;i<f->bbattr->len;i++) {
	tmp = g_ptr_array_index (f->bbattr, i);
	if (tmp->data) {
	    blks++;
	    if (!min)
		min = tmp;
	    else if (tmp->usage < min->usage)
		min = tmp;
	}
	tmp->usage = (guint32)tmp->usage*0.707;
    }
    if (blks < MAX_CACHED_BLOCKS)
	min = 0;

    g_assert (!attr->data);
    if (min) {
	g_assert (min->data);
#if OLE_DEBUG > 2
	g_print ("Replacing cache block %d with %d\n", min->blk, b);
#endif
	attr->data  = min->data;
	min->data   = 0;
	min->usage  = 0;
    } else
	attr->data = g_new (guint8, BB_BLOCK_SIZE);

    offset = (b+1)*BB_BLOCK_SIZE;
    f->syswrap->lseek (f->file_des, offset, SEEK_SET);
    f->syswrap->read (f->file_des, attr->data, BB_BLOCK_SIZE);
    attr->usage = 1;
    attr->dirty = forwrite;

    return attr->data;
}


/* This is a list of big blocks which contain a flat description of all blocks
   in the file. Effectively inside these blocks is a FAT of chains of other BBs,
   so the theoretical max size = 128 BB Fat blocks, thus = 128*512*512/4 blocks
   ~= 8.4MBytes */
/* FIXME tenix the max size would actually be 109*512*512/4 + 512 blocks ~=
   7MBytes if we don't take in count the additional Big Block Depot lists.
   Number of additional lists is in header:0x48, the location of the first
   additional list is in header:0x44, the location of the second additional
   list is at the very end of the first additional list and so on, the last
   additional list have at the end a END_OF_CHAIN.
   Each additional list can address 128*512/4*512 blocks ~= 8MBytes */
/* The number of Big Block Descriptor (fat) Blocks */
#define GET_NUM_BBD_BLOCKS(f)   (MS_OLE_GET_GUINT32 ((f)->mem + 0x2c))
#define SET_NUM_BBD_BLOCKS(f,n) (MS_OLE_SET_GUINT32 ((f)->mem + 0x2c, (n)))
/* The block locations of the Big Block Descriptor Blocks */
#define MAX_SIZE_BBD_LIST           109
/* FIXME tenix next is broken with big files */
#define GET_BBD_LIST(f,i)           (MS_OLE_GET_GUINT32 ((f)->mem + 0x4c + (i)*4))
/* FIXME tenix next is broken with big files */
#define SET_BBD_LIST(f,i,n)         (MS_OLE_SET_GUINT32 ((f)->mem + 0x4c + (i)*4, (n)))
#define NEXT_BB(f,n)                (g_array_index ((f)->bb, BLP, n))
#define NEXT_SB(f,n)                (g_array_index ((f)->sb, BLP, n))
/* Additional Big Block Descriptor (fat) Blocks */
#define MAX_SIZE_ADD_BBD_LIST       127
#define GET_NUM_ADD_BBD_LISTS(f)   (MS_OLE_GET_GUINT32 ((f)->mem + 0x48))
#define GET_FIRST_ADD_BBD_LIST(f)  (MS_OLE_GET_GUINT32 ((f)->mem + 0x44))

/* Get the start block of the root directory ( PPS ) chain */
#define GET_ROOT_STARTBLOCK(f)   (MS_OLE_GET_GUINT32 ((f)->mem + 0x30))
#define SET_ROOT_STARTBLOCK(f,i) (MS_OLE_SET_GUINT32 ((f)->mem + 0x30, i))
/* Get the start block of the SBD chain */
#define GET_SBD_STARTBLOCK(f)    (MS_OLE_GET_GUINT32 ((f)->mem + 0x3c))
#define SET_SBD_STARTBLOCK(f,i)  (MS_OLE_SET_GUINT32 ((f)->mem + 0x3c, i))


/* NB it is misleading to assume that Microsofts linked lists link correctly.
   It is not the case that pps_next(f, pps_prev(f, n)) = n ! For the final list
   item there are no valid links. Cretins. */
#define PPS_GET_NAME_LEN(p)   (MS_OLE_GET_GUINT16 (p + 0x40))
#define PPS_SET_NAME_LEN(p,i) (MS_OLE_SET_GUINT16 (p + 0x40, (i)))
#define PPS_GET_PREV(p)   ((PPS_IDX) MS_OLE_GET_GUINT32 (p + 0x44))
#define PPS_GET_NEXT(p)   ((PPS_IDX) MS_OLE_GET_GUINT32 (p + 0x48))
#define PPS_GET_DIR(p)    ((PPS_IDX) MS_OLE_GET_GUINT32 (p + 0x4c))
#define PPS_SET_PREV(p,i) ((PPS_IDX) MS_OLE_SET_GUINT32 (p + 0x44, i))
#define PPS_SET_NEXT(p,i) ((PPS_IDX) MS_OLE_SET_GUINT32 (p + 0x48, i))
#define PPS_SET_DIR(p,i)  ((PPS_IDX) MS_OLE_SET_GUINT32 (p + 0x4c, i))
/* These get other interesting stuff from the PPS record */
#define PPS_GET_STARTBLOCK(p)      ( MS_OLE_GET_GUINT32 (p + 0x74))
#define PPS_GET_SIZE(p)            ( MS_OLE_GET_GUINT32 (p + 0x78))
#define PPS_GET_TYPE(p) ((MsOleType)( MS_OLE_GET_GUINT8 (p + 0x42)))
#define PPS_SET_STARTBLOCK(p,i)    ( MS_OLE_SET_GUINT32 (p + 0x74, i))
#define PPS_SET_SIZE(p,i)          ( MS_OLE_SET_GUINT32 (p + 0x78, i))
#define PPS_SET_TYPE(p,i)          ( MS_OLE_SET_GUINT8  (p + 0x42, i))

/* Try to mark the Big Block "b" as as unused if it is marked as "c", in the
   FAT "f". */
#define TRY_MARK_UNUSED_BLOCK(f,block,mark) {                                    \
        if (g_array_index ((f), BLP, (block)) != (mark)) {			 \
	        g_warning ("Tried to mark as unused the block %d which has %d\n",\
                   (block), g_array_index ((f), BLP, (block)));                  \
        } else { g_array_index ((f), BLP, (block)) = UNUSED_BLOCK; } }

/* FIXME: This needs proper unicode support ! current support is a guess */
/* Length is in bytes == 1/2 the final text length */
/* NB. Different from biff_get_text, looks like a bug ! */
static char *
pps_get_text (guint8 *ptr, int length)
{
    int lp;
    char *ans;
    guint16 c;
    guint8 *inb;

    length = (length+1)/2;

    if (length <= 0 ||
	length > (PPS_BLOCK_SIZE/4)) {
#if OLE_DEBUG > 0
	g_print ("Nulled name of length %d\n", length);
#endif
	return 0;
    }

    ans = g_malloc(length + 1);

    inb = ptr;
    for (lp = 0; lp < length; lp++) {
	c = MS_OLE_GET_GUINT16 (inb);
	ans [lp] = (char) c;
	inb += 2;
    }
    ans [lp] = 0;

    return ans;
}

/*
 * get_next_block:
 * @f:   the file handle
 * @blk: an index into the big block fat
 *
 * Return value: the block index of the BBD block.
 */
static BLP
get_next_block (MsOle *f, BLP blk, gboolean *err)
{
    BLP bbd = GET_BBD_LIST (f, blk / (BB_BLOCK_SIZE / 4));

    if (bbd > BLOCK_COUNT (f)) {
	*err = TRUE;
	return 0;
    } else
	*err = FALSE;

    return MS_OLE_GET_GUINT32 (BB_R_PTR (f, bbd) +
			       4 * (blk % (BB_BLOCK_SIZE / 4)));
}

/* Builds the FAT */
static int
read_bb (MsOle *f)
{
    /* FIXME tenix may be later we wish to split this function */
    guint32  numbbd;
    BLP      lp;
    guint32  num_add_bbd_lists;
    BLP      missing_lps;
    BLP      missing_bbds;
    guint32  visited_add_bbd_list;
    BLP tmp;
    BLP bbd;

    g_return_val_if_fail (f, 0);
    g_return_val_if_fail (f->mem, 0);

    f->bb   = g_array_new (FALSE, FALSE, sizeof(BLP));
    numbbd  = GET_NUM_BBD_BLOCKS  (f);

    /* Add BBD's that live in the BBD list */
    for (lp = 0; (lp < BLOCK_COUNT (f) - 1) &&
	     (lp < MAX_SIZE_BBD_LIST * BB_BLOCK_SIZE / 4); lp++) {
	gboolean err;

	tmp = get_next_block (f, lp, &err);
	if (err)
	    return 0;

	g_array_append_val (f->bb, tmp);
    }

    /* Add BBD's that live in the additional BBD lists */
    num_add_bbd_lists = GET_NUM_ADD_BBD_LISTS (f);
    if (num_add_bbd_lists > 0) {
	if (lp != MAX_SIZE_BBD_LIST * BB_BLOCK_SIZE / 4)
	    return 0;

	visited_add_bbd_list = GET_FIRST_ADD_BBD_LIST (f);
	missing_lps = BLOCK_COUNT (f) - 1
	    - MAX_SIZE_BBD_LIST*BB_BLOCK_SIZE/4;
	for (lp = 0; lp < missing_lps; lp++) {
	    if ((lp != 0) && !(lp%(MAX_SIZE_ADD_BBD_LIST*
				   (BB_BLOCK_SIZE/4)))) {
		/* This lp lives in the next add bbd list */
		visited_add_bbd_list = MS_OLE_GET_GUINT32 (
							   BB_R_PTR(f,visited_add_bbd_list)
							   +4*MAX_SIZE_ADD_BBD_LIST);
		if (visited_add_bbd_list == END_OF_CHAIN) {
		    if (lp + 1 != missing_lps) {
			/* FIXME tenix error */
		    }
		}
	    }

	    /* tmp here means the number of one block that
	       belongs to the fat */
	    bbd = MS_OLE_GET_GUINT32 (BB_R_PTR (f, visited_add_bbd_list) + 4*((lp/(BB_BLOCK_SIZE/4))%MAX_SIZE_ADD_BBD_LIST));
	    tmp = MS_OLE_GET_GUINT32 (BB_R_PTR(f,bbd) +
				      4 * (lp % (BB_BLOCK_SIZE / 4)));
	    g_array_append_val (f->bb, tmp);
	}
	/* FIXME tenix do we check if we have visited all lp's but
	   there are more additional lists? */
    }

    /* Mark the bbd list blocks as unused */
    for (lp = 0; lp < MIN (numbbd, MAX_SIZE_BBD_LIST); lp++) {
	TRY_MARK_UNUSED_BLOCK (f->bb, GET_BBD_LIST(f,lp),
			       SPECIAL_BLOCK);
    }
    if (num_add_bbd_lists > 0) {
	visited_add_bbd_list = GET_FIRST_ADD_BBD_LIST (f);
	TRY_MARK_UNUSED_BLOCK (f->bb, visited_add_bbd_list,
			       ADD_BBD_LIST_BLOCK);
	missing_bbds = numbbd - MAX_SIZE_BBD_LIST;
	for (lp = 0; lp < missing_bbds; lp++) {
	    if ((lp != 0) && !(lp % (MAX_SIZE_ADD_BBD_LIST))) {
		/* This lp lives in the next add bbd list */
		visited_add_bbd_list = MS_OLE_GET_GUINT32 (
							   BB_R_PTR(f,visited_add_bbd_list)
							   + 4*MAX_SIZE_ADD_BBD_LIST);
		if (visited_add_bbd_list == END_OF_CHAIN) {
		    if (lp + 1 != missing_lps) {
			/* FIXME tenix error */
		    }
		}
		TRY_MARK_UNUSED_BLOCK (f->bb,
				       visited_add_bbd_list,
				       ADD_BBD_LIST_BLOCK);
	    }

	    bbd = MS_OLE_GET_GUINT32 (BB_R_PTR(f, visited_add_bbd_list) + 4*(lp%MAX_SIZE_ADD_BBD_LIST));
	    TRY_MARK_UNUSED_BLOCK (f->bb, bbd, SPECIAL_BLOCK);
	}
    }

    g_assert (f->bb->len < BLOCK_COUNT (f));

    return 1;
}

static guint8 *
get_pps_ptr (MsOle *f, PPS_IDX i, gboolean forwrite)
{
    int lp;
    BLP blk = GET_ROOT_STARTBLOCK (f);

    lp = i/(BB_BLOCK_SIZE/PPS_BLOCK_SIZE);
    while (lp && blk != END_OF_CHAIN) {
	if (blk == SPECIAL_BLOCK ||
	    blk == UNUSED_BLOCK) {
	    g_warning ("Duff block in root chain\n");
	    return 0;
	}
	lp--;
	blk = NEXT_BB (f, blk);
    }
    if (blk == END_OF_CHAIN) {
	g_warning ("Serious error finding pps %d\n", i);
	return 0;
    }

#if OLE_DEBUG > 0
    g_print ("get_pps_ptr: blk = %d\n", blk);
#endif

    return BB_R_PTR(f, blk) + (i%(BB_BLOCK_SIZE/PPS_BLOCK_SIZE))*PPS_BLOCK_SIZE;
}

static gint
pps_compare_func (PPS *a, PPS *b)
{
    g_return_val_if_fail (a, 0);
    g_return_val_if_fail (b, 0);
    g_return_val_if_fail (a->name, 0);
    g_return_val_if_fail (b->name, 0);

    return g_strcasecmp (b->name, a->name);
}

static void
pps_decode_tree (MsOle *f, PPS_IDX p, PPS *parent)
{
    PPS    *pps;
    guint8 *mem;

    if (p == PPS_END_OF_CHAIN)
	return;

    pps           = g_new (PPS, 1);
    pps->sig      = PPS_SIG;
    mem           = get_pps_ptr (f, p, FALSE);
    if (!mem) {
	g_warning ("Serious directory error %d\n", p);
	f->pps = NULL;
	return;
    }
#if OLE_DEBUG > 0
    g_print ("pps_decode_tree: mem (offset)= %#8.8x\n", mem - f->mem);
#endif
    pps->name     = pps_get_text  (mem, PPS_GET_NAME_LEN(mem));
    pps->type     = PPS_GET_TYPE  (mem);
    pps->size     = PPS_GET_SIZE  (mem);
    pps->children = NULL;
    pps->parent   = parent;
    pps->idx      = 0;
    if (!pps->name) { /* Make safe */
	g_print ("how odd: blank named file in directory\n");
	g_free (pps);
	return;
    }

    f->num_pps++;

    if (parent) {
#if OLE_DEBUG > 0
	g_print ("Inserting '%s' into '%s'\n", pps->name, parent->name);
#endif
	parent->children = g_list_insert_sorted (parent->children, pps,
						 (GCompareFunc)pps_compare_func);
    }
    else {
#if OLE_DEBUG > 0
	g_print ("Setting root to '%s'\n", pps->name);
#endif
	f->pps = g_list_append (0, pps);
    }

    if (PPS_GET_NEXT(mem) != PPS_END_OF_CHAIN)
	pps_decode_tree (f, PPS_GET_NEXT(mem), parent);

    if (PPS_GET_PREV(mem) != PPS_END_OF_CHAIN)
	pps_decode_tree (f, PPS_GET_PREV(mem), parent);

    if (PPS_GET_DIR (mem) != PPS_END_OF_CHAIN)
	pps_decode_tree (f, PPS_GET_DIR(mem), pps);

    pps->start   = PPS_GET_STARTBLOCK (mem);

#if OLE_DEBUG > 1
    g_print ("PPS decode : '%s'\n", pps->name?pps->name:"Null");
#endif
    return;
}

static int
read_pps (MsOle *f)
{
    PPS *pps;
    g_return_val_if_fail (f, 0);

    f->num_pps = 0;
    pps_decode_tree (f, PPS_ROOT_INDEX, NULL);

    if (!f->pps || g_list_length (f->pps) < 1 ||
	g_list_length (f->pps) > 1) {
	g_warning ("Invalid root chain\n");
	return 0;
    } else if (!f->pps->data) {
	g_warning ("No root entry\n");
	return 0;
    }

    /* Fiddle root, perhaps our get_text is broken */
    /* perhaps it is just an MS oddity in coding */
    pps = f->pps->data;
    if (pps->name)
	g_free (pps->name);
    pps->name = g_strdup ("Root Entry");

    { /* Free up the root chain */
	BLP blk, last;
	last = blk = GET_ROOT_STARTBLOCK (f);
	while (blk != END_OF_CHAIN) {
	    last = blk;
	    blk = NEXT_BB (f, blk);
	    g_array_index (f->bb, BLP, last) = UNUSED_BLOCK;
	}
    }

    if (!f->pps) {
	g_warning ("Root directory too small\n");
	return 0;
    }
    return 1;
}

static void
destroy_pps (GList *l)
{
    GList *tmp;

    if (!l)
	return;

    for (tmp = l; tmp; tmp = g_list_next (tmp)) {
	PPS *pps = tmp->data;
	if (pps->name) {
	    g_free (pps->name);
	    pps->name = NULL;
	}
	destroy_pps (pps->children);
	pps->children = NULL;
	g_free (pps);
	pps = NULL;
    }
    g_list_free (l);
}


static int
read_sb (MsOle *f)
{
    BLP ptr;
    int lastidx, idx;
    PPS *root;

    g_return_val_if_fail (f, 0);
    g_return_val_if_fail (f->pps, 0);

    root = f->pps->data;
    g_return_val_if_fail (root, 0);

    f->sbf = g_array_new (FALSE, FALSE, sizeof(BLP));
    f->sb  = g_array_new (FALSE, FALSE, sizeof(BLP));

    /* List of big blocks in SB file */
    ptr = root->start;
#if OLE_DEBUG > 0
    g_print ("Starting Small block file at %d\n", root->start);
#endif
    while (ptr != END_OF_CHAIN) {
	if (ptr == UNUSED_BLOCK ||
	    ptr == SPECIAL_BLOCK) {
	    g_warning ("Corrupt small block file: serious error, "
		       "invalid block in chain\n");
	    g_array_free (f->sbf, TRUE);
	    f->sbf = 0;
	    return 0;
	}
	g_array_append_val (f->sbf, ptr);
	ptr = NEXT_BB (f, ptr);
    }

    /* Description of small blocks */
    lastidx = -1;
    idx     = 0;
    ptr = GET_SBD_STARTBLOCK (f);

    if (f->sbf->len == 0 && ptr != END_OF_CHAIN) {
#if 0
	g_warning ("No small block file, but small block depot start block exists!: "
		   "ignore depot, since there's no small block files after all.\n");
#endif
	ptr = END_OF_CHAIN;
    }

    while (ptr != END_OF_CHAIN) {
	guint32 lp;
	if (ptr == UNUSED_BLOCK ||
	    ptr == SPECIAL_BLOCK) {
	    g_warning ("Corrupt file descriptor: serious error, "
		       "invalid block in chain\n");
	    g_array_free (f->sb, TRUE);
	    f->sb = 0;
	    return 0;
	}
	for (lp = 0;lp<BB_BLOCK_SIZE/4;lp++) {
	    BLP p = MS_OLE_GET_GUINT32 (BB_R_PTR(f, ptr) + lp*4);
	    g_array_append_val (f->sb, p);

	    if (p != UNUSED_BLOCK)
		lastidx = idx;
	    idx++;
	}
	ptr = NEXT_BB (f, ptr);
    }
    if (lastidx>0)
	g_array_set_size (f->sb, lastidx+1);

    if (f->sbf->len * BB_BLOCK_SIZE < f->sb->len*SB_BLOCK_SIZE) {
	g_warning ("Not enough small block file for descriptors\n"
		   "sbf->len == %d, sb->len == %d\n", f->sbf->len,
		   f->sb->len);
	return 0;
    }

    return 1;
}

static int
ms_ole_setup (MsOle *f)
{
    if (!f->ole_mmap) {
	guint32 i;
	f->bbattr = g_ptr_array_new ();
	for (i = 0; i < BLOCK_COUNT (f); i++)
	    g_ptr_array_add (f->bbattr, bb_blk_attr_new (i));
    }

    if (read_bb  (f) &&
	read_pps (f) &&
	read_sb  (f)) {
	return 1;
    }
    return 0;
}

static MsOle *
ms_ole_new ()
{
    MsOle *f = g_new0 (MsOle, 1);

    f->mem    = (guint8 *)0xdeadbeef;
    f->length = 0;
    f->mode   = 'r';
    f->bb     = 0;
    f->bbattr = 0;
    f->sb     = 0;
    f->sbf    = 0;
    f->pps    = 0;
    f->dirty  = 0;

    return f;
}


/**
 * ms_ole_ref:
 * @fs: filesystem object.
 *
 * Increment by one the count of references to the filesystem.
 **/
void
ms_ole_ref (MsOle *fs)
{
    g_return_if_fail (fs != NULL);
    fs->ref_count++;
}


/**
 * ms_ole_unref:
 * @fs: filesystem object.
 *
 * Decrement by one the count of references to the filesystem.
 **/
void
ms_ole_unref (MsOle *fs)
{
    g_return_if_fail (fs != NULL);
    fs->ref_count--;
}


/**
 * ms_ole_open_vfs:
 * @fs: filesystem object.
 * @path: path to the filesystem-in-the file on the actual filesystem.
 * @try_mmap: TRUE if try to mmap(2) the filesystem-in-a-file,
 *            instead of opening.
 * @wrappers: system functions wrappers, %NULL if standard functions are used.
 *
 * Opens the filesystem-in-the-file @path and creates the filesystem object @fs.
 *
 * Return value: a #MsOleErr code.
 **/
MsOleErr
ms_ole_open_vfs (MsOle **fs, const char *name,
		 gboolean try_mmap,
		 MsOleSysWrappers *wrappers)
{
    int prot = PROT_READ | PROT_WRITE;
    MsOle *f;
    int file;

    if (!fs)
	return MS_OLE_ERR_BADARG;

#if OLE_DEBUG > 0
    g_print ("New OLE file '%s'\n", name);
#endif

    f = *fs = ms_ole_new ();
    take_wrapper_functions (f, wrappers);

    /* Allin Cottrell modification: we only need to read */

    f->file_des = file = f->syswrap->open2 (name, O_RDONLY);
    f->ref_count = 0;
    f->mode = 'r';
	
    if ((file == -1) || !(f->syswrap->isregfile (file))) {
	g_warning ("No such file '%s'\n", name);
	g_free (f) ;
	*fs = NULL;
	return MS_OLE_ERR_EXIST;
    }

    if (f->syswrap->getfilesize (file, &(f->length) )) {
	g_warning ("Couldn't get the size of file '%s'\n", name);
	f->syswrap->close (file) ;
	g_free (f);
	*fs = NULL;
	return MS_OLE_ERR_EXIST;
    }

    if (f->length <= 0x4c) { /* Bad show */
#if OLE_DEBUG > 0
	g_warning ("File '%s' too short\n", name);
#endif
	f->syswrap->close (file) ;
	g_free (f) ;
	*fs = NULL;
	return MS_OLE_ERR_FORMAT;
    }

    if (try_mmap && f->syswrap->mmap) {
	f->mem = f->syswrap->mmap(0, f->length, prot, MAP_SHARED, file, 0);
	if (!f->mem || (caddr_t) f->mem == (caddr_t) MAP_FAILED) {
	    g_warning ("I can't mmap that file, falling back to slower method");
	    f->mem = NULL;
	} else
	    f->ole_mmap = TRUE;
    } else {
	f->mem = NULL;
#if OLE_DEBUG > 0
	g_warning ("I won't mmap that file, using a slower method\n");
#endif
    }


    if (f->mem == NULL) {
	f->ole_mmap = FALSE;
	f->mem = g_new (guint8, BB_BLOCK_SIZE);

	if (!f->mem ||
	    f->syswrap->read (file, f->mem, BB_BLOCK_SIZE) == -1) {
	    g_warning ("Error reading header\n");
	    f->syswrap->close (file) ;
	    g_free (f);
	    *fs = NULL;
	    return MS_OLE_ERR_EXIST;
	}
    }

    if (MS_OLE_GET_GUINT32 (f->mem    ) != 0xe011cfd0 ||
	MS_OLE_GET_GUINT32 (f->mem + 4) != 0xe11ab1a1) {
	g_warning("Failed OLE2 magic number %x %x\n",
		  MS_OLE_GET_GUINT32(f->mem), MS_OLE_GET_GUINT32(f->mem+4));
	ms_ole_destroy (fs);
	return MS_OLE_ERR_FORMAT;
    }

    {
	unsigned short bbs = MS_OLE_GET_GUINT16 (f->mem + 0x1e);
	unsigned short sbs = MS_OLE_GET_GUINT16 (f->mem + 0x20);

	if ((1 << bbs) != BB_BLOCK_SIZE)
	    g_warning ("Big-block-size mismatch [%d] -- expect trouble.", bbs);

	if ((1 << sbs) != SB_BLOCK_SIZE)
	    g_warning ("Small-block-size mismatch [%d] -- expect trouble.", sbs);
    }

#if 0
    if (f->length % BB_BLOCK_SIZE)
	g_warning ("Warning file '%s': %d bytes, non-integer number of blocks\n",
		   name, f->length);
#endif

    if (!ms_ole_setup (f)) {
	g_warning ("'%s' : duff file !\n", name);
	ms_ole_destroy (fs);
	return MS_OLE_ERR_FORMAT;
    }

    g_assert (f->bb->len < BLOCK_COUNT (f));

#if OLE_DEBUG > 0
    g_print ("New OLE file '%s'\n", name);
#endif
    /* If writing then when destroy commit it */
    return MS_OLE_ERR_OK;
}


/**
 * ms_ole_destroy:
 * @fs: filesystem object.
 *
 * Closes the filesystem @fs and truncates any free blocks.
 **/
void
ms_ole_destroy (MsOle **ptr)
{
    MsOle *f = *ptr;

#if OLE_DEBUG > 0
    g_print ("FIXME: should truncate to remove unused blocks\n");
#endif
    if (f) {
	if (f->ref_count != 0)
	    g_warning ("Unclosed files exist on this OLE stream\n");

	if (f->mem == (void *)0xdeadbeef)
	    f->mem = NULL;
	else if (f->ole_mmap) {
#ifdef HAVE_MMAP
	    munmap (f->mem, f->length);
#else
	    g_warning ("Unmapping while we dont have mmap call");
#endif
	} else {
	    if (f->bbattr) {
		guint32 i;
		for (i = 0; i < f->bbattr->len; i++) {
		    BBBlkAttr *attr = g_ptr_array_index (f->bbattr, i);
		    g_free (attr->data);
		    attr->data = NULL;
		    g_free (attr);
		}
		f->bbattr = NULL;
	    }

	    if (f->dirty) {
		f->syswrap->lseek (f->file_des, 0, SEEK_SET);
		f->syswrap->write (f->file_des, f->mem,
				   BB_BLOCK_SIZE);
	    }
	    g_free (f->mem);
	    f->mem = NULL;
	}

	destroy_pps (f->pps);
	f->pps = NULL;

	f->syswrap->close (f->file_des);
	g_free (f);

#if OLE_DEBUG > 0
	g_print ("Closing OLE file\n");
#endif
    }
    *ptr = NULL;
}

static MsOlePos
tell_pos (MsOleStream *s)
{
    return s->position;
}

/**
 * ms_ole_lseek:
 * @s: stream object.
 * @bytes: number of bytes to set the stream pointer.
 * @type: relative from where the stream pointer will be set.
 *
 * Set the stream pointer for @s as many as @bytes bytes according to @type.
 *
 * Return value: the new position of the stream pointer.
 **/
static MsOleSPos
ms_ole_lseek (MsOleStream *s, MsOleSPos bytes, MsOleSeek type)
{
    /* FIXME tenix improve limits detection: avoid gint vs guint limits */
    MsOleSPos newpos;

    g_return_val_if_fail (s, -1);

    if (type == MsOleSeekSet)
	newpos = bytes;
    else if (type == MsOleSeekCur)
	newpos = s->position + bytes;
    else
	newpos = s->size + bytes;

    if (newpos > s->size || newpos < 0) {
	g_warning ("Invalid seek");
	return -1;
    }
    s->position = newpos;
    return newpos;
}


/*
 *  Returns:
 *  NULL    - on error
 */
static guint8*
ms_ole_read_ptr_bb (MsOleStream *s, MsOlePos length)
{
    int blklen;
    guint8 *ans;
    guint32 len = length;
    int blockidx = s->position / BB_BLOCK_SIZE;

    g_return_val_if_fail (s, NULL);

    if (!s->blocks || blockidx >= s->blocks->len) {
	g_warning ("Reading from NULL file\n");
	return NULL;
    }

    blklen = BB_BLOCK_SIZE - s->position % BB_BLOCK_SIZE;

    if (len > blklen && !s->file->ole_mmap)
	return NULL;

    while (len > blklen) {
	len -= blklen;
	blklen = BB_BLOCK_SIZE;
	if (blockidx >= (s->blocks->len - 1)
	    || (ms_array_index (s->blocks, BLP, blockidx)
		!= blockidx + 1))
	    return NULL;
	blockidx++;
    }
    /* Straight map, simply return a pointer */
    ans = BB_R_PTR (s->file, ms_array_index (s->blocks, BLP,
					     s->position / BB_BLOCK_SIZE))
	+ s->position % BB_BLOCK_SIZE;
    ms_ole_lseek (s, length, MsOleSeekCur);

    return ans;
}


/*
 *  Returns:
 *  NULL    - on error
 */
static guint8*
ms_ole_read_ptr_sb (MsOleStream *s, MsOlePos length)
{
    int blklen;
    guint8 *ans;
    guint32 len = length;
    int blockidx = s->position / SB_BLOCK_SIZE;

    g_return_val_if_fail (s, NULL);

    if (!s->blocks || blockidx >= s->blocks->len) {
	g_warning ("Reading from NULL file\n");
	return NULL;
    }

    blklen = SB_BLOCK_SIZE - s->position % SB_BLOCK_SIZE;

    if (len > blklen && !s->file->ole_mmap)
	return NULL;

    while (len > blklen) {
	len -= blklen;
	blklen = SB_BLOCK_SIZE;
	if (blockidx >= (s->blocks->len - 1)
	    || (ms_array_index (s->blocks, BLP, blockidx)
		!= blockidx + 1))
	    return NULL;
	blockidx++;
    }
    /* Straight map, simply return a pointer */
    ans = GET_SB_R_PTR (s->file, ms_array_index (s->blocks, BLP,
						 s->position / SB_BLOCK_SIZE))
	+ s->position%SB_BLOCK_SIZE;
    ms_ole_lseek (s, length, MsOleSeekCur);

    return ans;
}


/*
 *  Returns:
 *  zero    - on error
 *  no zero - on success
 */
static gint
ms_ole_read_copy_bb (MsOleStream *s, guint8 *ptr, MsOlePos length)
{
    guint8 *src;
    int offset = s->position % BB_BLOCK_SIZE;
    int blkidx = s->position / BB_BLOCK_SIZE;

    g_return_val_if_fail (s, 0);
    g_return_val_if_fail (ptr, 0);

    if (!s->blocks) {
	g_warning ("Reading from NULL file\n");
	return 0;
    }

    while (length > 0) {
	BLP block;
	int cpylen = BB_BLOCK_SIZE - offset;
	if (cpylen > length)
	    cpylen = length;

	if (s->position + cpylen > s->size
	    || blkidx == s->blocks->len) {
#if OLE_DEBUG > 0
	    g_print ("Trying 2 to read beyond end of stream %d+%d %d\n",
		     s->position, cpylen, s->size);
#endif
	    return 0;
	}
	g_assert (blkidx < s->blocks->len);
	block = ms_array_index (s->blocks, BLP, blkidx);
	src = BB_R_PTR (s->file, block) + offset;

	memcpy (ptr, src, cpylen);
	ptr    += cpylen;
	length -= cpylen;

	offset = 0;

	blkidx++;
	s->position += cpylen;
    }

    return 1;
}


/*
 *  Returns:
 *  zero    - on error
 *  no zero - on success
 */
static gint
ms_ole_read_copy_sb (MsOleStream *s, guint8 *ptr, MsOlePos length)
{
    int offset = s->position % SB_BLOCK_SIZE;
    int blkidx = s->position / SB_BLOCK_SIZE;
    guint8 *src;

    g_return_val_if_fail (s, 0);
    g_return_val_if_fail (ptr, 0);

    if (!s->blocks) {
	g_warning ("Reading from NULL file\n");
	return 0;
    }

    while (length > 0) {
	int cpylen = SB_BLOCK_SIZE - offset;
	BLP block;
	if (cpylen>length)
	    cpylen = length;
	if (s->position + cpylen > s->size || blkidx == s->blocks->len) {
#if OLE_DEBUG > 0
	    g_print ("Trying 3 to read beyond end of stream %d+%d %d\n",
		     s->position, cpylen, s->size);
#endif
	    return 0;
	}
	g_assert (blkidx < s->blocks->len);
	block = ms_array_index (s->blocks, BLP, blkidx);
	src = GET_SB_R_PTR (s->file, block) + offset;

	memcpy (ptr, src, cpylen);
	ptr += cpylen;
	length -= cpylen;

	offset = 0;

	blkidx++;
	s->position += cpylen;
    }

    return 1;
}

/**
 * pps_create:
 * @f: ole file handle.
 * @p: returned pps.
 * @parent: parent pps.
 * @name: its name.
 * @type: the type.
 *
 * Creates a storage or stream.
 *
 * Return value: error status.
 **/
static MsOleErr
pps_create (MsOle *f, GList **p, GList *parent, const char *name,
	    MsOleType type)
{
    PPS *pps, *par;

    if (!p || !parent || !parent->data || !name) {
	g_warning ("duff arguments to pps_create");
	return MS_OLE_ERR_BADARG;
    }

    pps  = g_new (PPS, 1);
    if (!pps) {
	return MS_OLE_ERR_MEM;
    }

    pps->sig      = PPS_SIG;
    pps->name     = g_strdup (name);
    pps->type     = type;
    pps->size     = 0;
    pps->start    = END_OF_CHAIN;
    pps->children = NULL;
    pps->parent   = parent->data;

    par = (PPS *)parent->data;
    par->children = g_list_insert_sorted (par->children, pps,
					  (GCompareFunc)pps_compare_func);
    *p = g_list_find (par->children, pps);
    f->num_pps++;

    return MS_OLE_ERR_OK;
}


/**
 * find_in_pps:
 * @l: the parent storage chain element.
 *
 * Find the right Stream ... John 4:13-14 ...
 * in a storage
 *
 * Return value: %NULL if not found or pointer to the child list
 **/
static GList *
find_in_pps (GList *l, const char *name)
{
    PPS   *pps;
    GList *cur;

    g_return_val_if_fail (l != NULL, NULL);
    g_return_val_if_fail (l->data != NULL, NULL);
    pps = l->data;
    g_return_val_if_fail (IS_PPS (pps), NULL);

    if (pps->type == MsOleStorageT ||
	pps->type == MsOleRootT) {
	cur = pps->children;
    } else {
	g_warning ("trying to enter a stream '%s'",
		   pps->name? pps->name : "no name");
	return NULL;
    }

    for ( ;cur ; cur = g_list_next (cur)) {
	PPS *pps = cur->data;
	g_return_val_if_fail (IS_PPS (pps), NULL);

	if (!pps->name)
	    continue;

	if (!g_strcasecmp (pps->name, name))
	    return cur;
    }
    return NULL;
}


/**
 * path_to_pps:
 * @pps:  pointer to pps to return value in.
 * @f:    ole file hande.
 * @path: path to find.
 * @file: file to find in path.
 * @create_if_not_found: create the pps with the given path if not found.
 *
 * Locates a stream or storage with the given path.
 *
 * Return value: a #MsOleErr code.
 **/
static MsOleErr
path_to_pps (PPS **pps, MsOle *f, const char *path,
	     const char *file,
	     gboolean create_if_not_found)
{
    guint     lp;
    gchar   **dirs;
    GList    *cur, *parent;

    g_return_val_if_fail (f != NULL, MS_OLE_ERR_BADARG);
    g_return_val_if_fail (path != NULL, MS_OLE_ERR_BADARG);

    dirs = g_strsplit (path, "/", -1);
    g_return_val_if_fail (dirs != NULL, MS_OLE_ERR_BADARG);

    parent = cur = f->pps;

    for (lp = 0; dirs[lp]; lp++) {
	if (dirs[lp][0] == '\0' || !cur) {
	    g_free (dirs[lp]);
	    continue;
	}

	parent = cur;

	cur = find_in_pps (parent, dirs[lp]);

	if (!cur && create_if_not_found &&
	    pps_create (f, &cur, parent, dirs[lp], MsOleStorageT) !=
	    MS_OLE_ERR_OK)
	    cur = NULL;
	/* else carry on not finding them before dropping out */

	g_free (dirs[lp]);
    }
    g_free (dirs);

    if (!cur || !cur->data)
	return MS_OLE_ERR_EXIST;

    if (file[0] == '\0') { /* We just want a directory */
	*pps = cur->data;
	g_return_val_if_fail (IS_PPS (cur->data), MS_OLE_ERR_INVALID);
	return MS_OLE_ERR_OK;
    }

    parent = cur;
    cur = find_in_pps (parent, file);

    /* now the file */
    if (!cur) {
	if (create_if_not_found) {
	    MsOleErr result;
	    result = pps_create (f, &cur, parent, file,
				 MsOleStreamT);
	    if (result == MS_OLE_ERR_OK) {
		*pps = cur->data;
		g_return_val_if_fail (IS_PPS (cur->data),
				      MS_OLE_ERR_INVALID);
		return MS_OLE_ERR_OK;
	    } else
		return result;
	}
	return MS_OLE_ERR_EXIST;
    }

    if (cur && cur->data) {
	*pps = cur->data;
	g_return_val_if_fail (IS_PPS (cur->data), MS_OLE_ERR_INVALID);
	return MS_OLE_ERR_OK;
    }

    return MS_OLE_ERR_EXIST;
}


/**
 * ms_ole_directory:
 * @names: array where the names are storesd, it's %NULL ended.
 * @fs: filesystem object.
 * @dirpath: directory path.
 *
 * Gets the names of the streams and directories in the directory @dirpath.
 *
 * Returns: a #MsOleErr code.
 **/
MsOleErr
ms_ole_directory (char ***names, MsOle *f, const char *path)
{
    char    **ans;
    PPS      *pps;
    MsOleErr  result;
    GList    *l;
    int       lp;

    g_return_val_if_fail (f != NULL, MS_OLE_ERR_BADARG);
    g_return_val_if_fail (path != NULL, MS_OLE_ERR_BADARG);

    if ((result = path_to_pps (&pps, f, path, "", FALSE)) !=
	MS_OLE_ERR_OK)
	return result;

    if (!pps)
	return MS_OLE_ERR_INVALID;

    l   = pps->children;
    ans = g_new (char *, g_list_length (l) + 1);

    lp = 0;
    for (; l; l = g_list_next (l)) {
	pps = (PPS *)l->data;

	if (!pps->name)
	    continue;

	ans[lp] = g_strdup (pps->name);
	lp++;
    }
    ans[lp] = NULL;

    *names = ans;
    return MS_OLE_ERR_OK;
}

/**
 * ms_ole_stream_open:
 * @stream: stream object.
 * @fs: filesystem object.
 * @dirpath: directory of the stream.
 * @name: stream name.
 * @mode: mode of opening stream.
 *
 * Opens the stream in @dirpath with the name @name and creates the stream
 * object @stream. If @mode is '%r' it opens read only, and if it is '%w'
 * it opens for write only.
 *
 * Return value: a #MsOleErr code.
 **/
MsOleErr
ms_ole_stream_open (MsOleStream ** const stream, MsOle *f,
		    const char *path, const char *fname, char mode)
{
    PPS         *p;
    MsOleStream *s;
    int lp, panic = 0;
    MsOleErr     result;

    if (!stream)
	return MS_OLE_ERR_BADARG;
    *stream = NULL;

    if (!path || !f)
	return MS_OLE_ERR_BADARG;

    if (mode == 'w' && f->mode != 'w') {
	g_print ("Opening stream '%c' when file is '%c' only\n",
		 mode, f->mode);
	return MS_OLE_ERR_PERM;
    }

    if ((result = path_to_pps (&p, f, path, fname, (mode == 'w'))) !=
	MS_OLE_ERR_OK)
	return result;

    if (p->type != MsOleStreamT)
	return MS_OLE_ERR_INVALID;

    s           = g_new0 (MsOleStream, 1);
    s->file     = f;
    s->pps      = p;
    s->position = 0;
    s->size     = p->size;
    s->blocks   = NULL;

#if OLE_DEBUG > 0
    g_print ("Parsing blocks\n");
#endif
    if (s->size >= BB_THRESHOLD) {
	BLP b = p->start;

	s->read_copy = ms_ole_read_copy_bb;
	s->read_ptr  = ms_ole_read_ptr_bb;
	s->lseek     = ms_ole_lseek;
	s->tell      = tell_pos;
	s->write     = NULL;

	s->blocks    = g_array_new (FALSE, FALSE, sizeof(BLP));
	s->type      = MsOleLargeBlock;
	for (lp = 0; !panic && (lp < (s->size + BB_BLOCK_SIZE - 1) / BB_BLOCK_SIZE); lp++) {
	    g_array_append_val (s->blocks, b);
#if OLE_DEBUG > 2
	    g_print ("Block %d\n", b);
#endif
	    if (b == END_OF_CHAIN ||
		b == SPECIAL_BLOCK ||
		b == UNUSED_BLOCK) {

		g_warning ("Panic: broken stream, truncating to block %d\n", lp);
		s->size = (lp-1)*BB_BLOCK_SIZE;
		panic   = 1;

#if OLE_DEBUG > 0
		if (b == END_OF_CHAIN)
		    g_warning ("Warning: bad file length in '%s'\n", p->name);
		else if (b == SPECIAL_BLOCK)
		    g_warning ("Warning: special block in '%s'\n", p->name);
		else if (b == UNUSED_BLOCK)
		    g_warning ("Warning: unused block in '%s'\n", p->name);
#endif
	    } else
		b = NEXT_BB (f, b);
	}
	if (b != END_OF_CHAIN) {
	    BLP next;
	    g_warning ("Panic: extra unused blocks on end of '%s', %x wiping it\n",
		       p->name, b);
	    while (b != END_OF_CHAIN &&
		   b != UNUSED_BLOCK &&
		   b != SPECIAL_BLOCK &&
		   b < f->bb->len) {
		next = NEXT_BB (f, b);
		g_array_index (f->bb, BLP, b) = END_OF_CHAIN;
		b = next;
	    }
	}
    } else {
	BLP b = p->start;

	s->read_copy = ms_ole_read_copy_sb;
	s->read_ptr  = ms_ole_read_ptr_sb;
	s->lseek     = ms_ole_lseek;
	s->tell      = tell_pos;
	s->write     = NULL;

	if (s->size > 0)
	    s->blocks = g_array_new (FALSE, FALSE, sizeof(BLP));
	else
	    s->blocks = NULL;

	s->type = MsOleSmallBlock;

	for (lp = 0; !panic && (lp < (s->size + SB_BLOCK_SIZE - 1) / SB_BLOCK_SIZE); lp++) {
	    g_array_append_val (s->blocks, b);
#if OLE_DEBUG > 0
	    g_print ("Block %d\n", b);
#endif
	    if (b == END_OF_CHAIN ||
		b == SPECIAL_BLOCK ||
		b == UNUSED_BLOCK) {

		g_warning ("Panic: broken stream, truncating to block %d\n", lp);
		s->size = (lp-1)*SB_BLOCK_SIZE;
		panic   = 1;
#if OLE_DEBUG > 0
		if (b == END_OF_CHAIN)
		    g_warning ("Warning: bad file length in '%s'\n", p->name);
		else if (b == SPECIAL_BLOCK)
		    g_warning ("Warning: special block in '%s'\n", p->name);
		else if (b == UNUSED_BLOCK)
		    g_warning ("Warning: unused block in '%s'\n", p->name);
#endif
	    } else
		b = NEXT_SB (f, b);
	}
	if (b != END_OF_CHAIN) {
	    BLP next;
	    g_warning ("Panic: extra unused blocks on end of '%s', wiping it\n",
		       p->name);
	    while (b != END_OF_CHAIN &&
		   b != UNUSED_BLOCK &&
		   b != SPECIAL_BLOCK &&
		   b < f->sb->len) {
		next = NEXT_SB (f, b);
		g_array_index (f->sb, BLP, b) = END_OF_CHAIN;
		b = next;
	    }
	    if (b != END_OF_CHAIN)
		g_warning ("Panic: even more serious block error\n");
	}
    }
    *stream = s;
    ms_ole_ref (s->file);

    return MS_OLE_ERR_OK;
}

/**
 * ms_ole_stream_close:
 * @stream: stream object to be closed.
 *
 * Closes the @stream.
 *
 * Return value: a #MsOleErr code.
 **/
MsOleErr
ms_ole_stream_close (MsOleStream ** const s)
{
    if (*s) {
	if ((*s)->file && (*s)->file->mode == 'w')
	    ((PPS *)(*s)->pps)->size = (*s)->size;

	if ((*s)->blocks)
	    g_array_free ((*s)->blocks, TRUE);

	ms_ole_unref ((*s)->file);

	g_free (*s);
	*s = NULL;

	return MS_OLE_ERR_OK;
    }
    return MS_OLE_ERR_BADARG;
}
