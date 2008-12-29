#ifndef GRETL_ENUMS_H
#define GRETL_ENUMS_H

enum extra_cmds {
    RELABEL = NC,
    VSETMISS,
    GSETMISS,
    SMPLDUM,
    SMPLBOOL,
    SMPLRAND,
    MARKERS,
    STORE_MODEL,
    VAR_SUMMARY,
    ALL_SUMMARY,
    ALL_CORR,
    GENR_RANDOM,
    SEED_RANDOM,
    ONLINE,
    EXPORT,
    MEANTEST2,
    MODEL_GENR,
    GR_PLOT,
    GR_XY,
    GR_IMP,
    GR_DUMMY,
    GR_XYZ,
    GR_BOX,
    GR_NBOX,
    GR_3D,
    COMPACT,
    EXPAND, 
    COEFFINT,
    COVAR,
    STAT_TABLE,
    H_TEST,
    TRAMO,
    X12A,
    VIEW_SERIES,
    VIEW_MODEL,
    VIEW_LOG,
    VIEW_DATA,
    VIEW_SCRIPT,
    VIEW_CODEBOOK,
    VIEW_MODELTABLE,
    VIEW_FUNC_INFO,
    VIEW_FUNC_CODE,
    VIEW_FILE,
    DATA_REPORT,
    SCRIPT_OUT,
    CONSOLE,
    EDIT_HEADER,
    EDIT_SCRIPT,
    EDIT_NOTES,
    EDIT_FUNC_CODE,
    EDIT_GP,
    EDIT_R,
    CLI_HELP,
    GUI_HELP,
    CLI_HELP_EN,
    GUI_HELP_EN,
    FUNCS_HELP,
    MODELTABLE,
    GRAPHPAGE,
    LMTEST_LOGS,
    LMTEST_SQUARES,
    LMTEST_WHITE,
    LMTEST_WHITE_NOX,
    LMTEST_BP,
    LMTEST_BPK,
    LMTEST_BG,
    LMTEST_GROUPWISE,
    KERNEL_DENSITY,
    CREATE_DATASET,
    HCCME,
    VAR_IRF,
    VAR_DECOMP,
    IRF_BOOT,
    HTEST, 
    HTESTNP,
    MODEL_RESTR,
    SYS_RESTR,
    VECM_RESTR,
    ELLIPSE,
    LAGS_DIALOG,
    COPY_FORMATS,
    MINIBUF,
    VLAGSEL,
    VAROMIT,
    DEFINE_LIST,
    DEFINE_MATRIX,
    PANEL_WLS,
    PANEL_MODE,
    TSPLOTS,
    ITERATIONS,
    CUSUMSQ,
    PANEL_B,
    IMPORT,
    BOOTSTRAP,
    TRANSPOS,
    DATASORT,
    WORKDIR,
    DFGLS,
    CORC,
    HILU,
    PWE,
    AFR, /* actual, fitted, residual */
    GPT_ADDLINE,
    GUI_CMD_MAX
};

enum file_ops {
    OPEN_DATA = GUI_CMD_MAX + 1, /* don't collide with extra_cmds */
    OPEN_RATS_DB,
    OPEN_PCGIVE_DB,
    OPEN_SCRIPT,
    OPEN_CSV,
    APPEND_DATA,
    APPEND_CSV,
    OPEN_ASCII,
    APPEND_ASCII,
    OPEN_OCTAVE,
    APPEND_OCTAVE,
    OPEN_GNUMERIC,
    APPEND_GNUMERIC,
    OPEN_XLS,
    APPEND_XLS,
    OPEN_WF1,
    APPEND_WF1,
    OPEN_DTA,
    APPEND_DTA,
    OPEN_SAV,
    APPEND_SAV,
    OPEN_JMULTI,
    APPEND_JMULTI,
    OPEN_ODS,
    APPEND_ODS,
    OPEN_SESSION,
    OPEN_MARKERS,
    END_OPEN,      /* marker for end of file open section */
    SAVE_DATA,
    SAVE_DATA_AS,
    SAVE_DBDATA,
    EXPORT_OCTAVE,
    EXPORT_R,
    EXPORT_CSV,
    EXPORT_DAT,
    EXPORT_JM,
    COPY_CSV,
    END_SAVE_DATA,  /* marker for end of data-saving section */
    SAVE_TEX,
    SAVE_RTF,
    SAVE_TEXT,
    SAVE_SCRIPT,
    SAVE_OUTPUT,
    SAVE_SESSION,
    SAVE_GNUPLOT,
    SAVE_GP_CMDS,
    SAVE_R_CMDS,
    SAVE_CONSOLE,
    SAVE_FUNCTIONS,
    SAVE_FUNCTIONS_AS,
    SAVE_BOOT_DATA,
    SET_PROG,
    SET_DIR,
    SET_WDIR,
    FILE_OP_MAX
};

#define SAVE_DATA_ACTION(i) (i >= SAVE_DATA && i < END_SAVE_DATA)

enum browser_codes {
    TEXTBOOK_DATA = FILE_OP_MAX + 1, /* don't collide with file_ops enum */
    PS_FILES,
    FUNC_FILES,
    REMOTE_FUNC_FILES,
    NATIVE_DB,
    REMOTE_DB,
    NATIVE_SERIES,
    RATS_SERIES,
    PCGIVE_SERIES,
    REMOTE_SERIES,
    MAINWIN
};

enum pref_tabs {
    TAB_NONE = 0,
    TAB_MAIN,
    TAB_DBS,
    TAB_PROGS,
    TAB_VCV,
    TAB_MAN,
    TAB_MAX
};

enum clipstuff {
    TARGET_UTF8_STRING,
    TARGET_STRING,
    TARGET_TEXT,
    TARGET_COMPOUND_TEXT,
    TARGET_RTF
};

enum data_status_flags {
    HAVE_DATA     = 1 << 0,
    BOOK_DATA     = 1 << 1,
    USER_DATA     = 1 << 2,
    IMPORT_DATA   = 1 << 3,
    GUI_DATA      = 1 << 4,
    MODIFIED_DATA = 1 << 5,
    GZIPPED_DATA  = 1 << 6
};

enum register_data_flags {
    DATAFILE_OPENED = 1,
    OPENED_VIA_CLI,
    OPENED_VIA_SESSION,
    DATA_APPENDED,
    NULLDATA_STARTED
};

enum drag_types {
    GRETL_FILENAME,
    GRETL_DBSERIES_PTR,
    GRETL_MODEL_PTR,
    GRETL_REMOTE_DB_PTR,
    GRETL_REMOTE_FNPKG_PTR,
    GRETL_GRAPH_FILE
};

enum file_lists {
    FILE_LIST_DATA,
    FILE_LIST_SESSION,
    FILE_LIST_SCRIPT,
    FILE_LIST_WDIR,
};

enum varclick_actions {
    VARCLICK_NONE,
    VARCLICK_INSERT_ID,
    VARCLICK_INSERT_NAME,
    VARCLICK_INSERT_TEXT
};

enum font_selections {
    FIXED_FONT_SELECTION,
    APP_FONT_SELECTION
};

enum calc_functions {
    CALC_PVAL,
    CALC_DIST,
    CALC_TEST,
    CALC_NPTEST,
    CALC_GRAPH,
    CALC_GRAPH_ADD,
    CALC_RAND,
    CALC_MAX
};

enum auto_dummies {
    TS_DUMMIES,
    PANEL_UNIT_DUMMIES,
    PANEL_TIME_DUMMIES
};

#define MULTI_FORMAT_ENABLED(c) (c == SUMMARY || c == VAR_SUMMARY || \
                                 c == ALL_SUMMARY || c == AFR || \
	                         c == CORR || c == ALL_CORR || \
	                         c == FCAST || c == COEFFINT || \
	                         c == COVAR || c == VIEW_MODEL || \
                                 c == VIEW_MODELTABLE || c == VAR || c == VECM || \
                                 c == VAR_IRF || c == VAR_DECOMP)

#endif /* GRETL_ENUMS_H */
