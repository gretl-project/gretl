/* macros for file selector */

#define IS_DAT_ACTION(i) (i == SAVE_DATA || \
                          i == SAVE_GZDATA || \
                          i == SAVE_BIN1 || \
                          i == SAVE_BIN2 || \
                          i == OPEN_DATA)

#define OPEN_DATA_ACTION(i)  (i == OPEN_DATA || \
                              i == OPEN_CSV || \
	                      i == OPEN_BOX || \
                              i == OPEN_GNUMERIC || \
	                      i == OPEN_EXCEL || \
                              i == OPEN_DES) 

#define APPEND_DATA_ACTION(i) (i == APPEND_CSV || \
                               i == APPEND_GNUMERIC || \
                               i == APPEND_EXCEL)

#define SAVE_TEX_ACTION(i) (i == SAVE_TEX_TAB || \
                            i == SAVE_TEX_EQ || \
	                    i == SAVE_TEX_TAB_FRAG || \
                            i == SAVE_TEX_EQ_FRAG)

#define SAVE_DATA_ACTION(i) (i >= SAVE_DATA && i < END_SAVE_DATA)

#define SAVE_GRAPH_ACTION(i) (i == SAVE_GNUPLOT || \
                              i == SAVE_THIS_GRAPH || \
                              i == SAVE_LAST_GRAPH || \
                              i == SAVE_BOXPLOT_EPS || \
                              i == SAVE_BOXPLOT_PS || \
                              i == SAVE_BOXPLOT_XPM)
