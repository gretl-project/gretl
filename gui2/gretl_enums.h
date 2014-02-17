#ifndef GRETL_ENUMS_H
#define GRETL_ENUMS_H

enum extra_cmds {
    RELABEL = NC,
    VSETMISS,
    GSETMISS,
    SMPLDUM,
    SMPLBOOL,
    SMPLRAND,
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
    GR_FBOX,
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
    VIEW_PKG_INFO,
    VIEW_PKG_CODE,
    VIEW_BUNDLE,
    VIEW_FILE,
    DATA_REPORT,
    SCRIPT_OUT,
    CONSOLE,
    EDIT_HEADER,
    EDIT_SCRIPT,
    EDIT_NOTES,
    EDIT_PKG_CODE,
    EDIT_PKG_SAMPLE,
    EDIT_GP,
    EDIT_R,
    EDIT_OX,
    EDIT_OCTAVE,
    EDIT_PYTHON,
    EDIT_X12A,
    EDIT_MAX,
    CLI_HELP,
    GUI_HELP,
    CLI_HELP_EN,
    GUI_HELP_EN,
    FUNCS_HELP,
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
    AFR, /* actual, fitted, residual */
    GPT_ADDLINE,
    GPT_CURVE,
    OLOGIT,
    MLOGIT,
    OPROBIT,
    REPROBIT,
    IV_LIML,
    IV_GMM,
    BFGS_CONFIG,
    COUNTMOD,
    BWFILTER,
    POLYWEIGHTS,
    X12AHELP,
    LOESS,
    NADARWAT,
    SSHEET,
    CLUSTER,
    GUI_FUNCS,
    MENU_ATTACH,
    GUI_CMD_MAX
};

#define help_role(r) (r >= CLI_HELP && r <= FUNCS_HELP)

enum file_ops {
    OPEN_DATA = GUI_CMD_MAX + 1, /* don't collide with extra_cmds */
    OPEN_RATS_DB,
    OPEN_PCGIVE_DB,
    OPEN_SCRIPT,
    APPEND_DATA,
    OPEN_SESSION,
    OPEN_MARKERS,
    OPEN_LABELS,
    OPEN_BARS,
    OPEN_GFN,
    END_OPEN,        /* marker for end of file open section */
    AUTO_SAVE_DATA,
    SAVE_DATA,
    SAVE_DATA_AS,
    EXPORT_OCTAVE,
    EXPORT_R,
    EXPORT_CSV,
    EXPORT_DAT,
    EXPORT_JM,
    EXPORT_DB,
    EXPORT_GDT,
    EXPORT_GDTB,
    COPY_CSV,
    END_SAVE_DATA,  /* marker for end of data-saving section */
    SAVE_TEX,
    SAVE_RTF,
    SAVE_TEXT,
    SAVE_SCRIPT,
    SAVE_OUTPUT,
    SAVE_SESSION,
    SAVE_GNUPLOT,
    SAVE_GRAPHIC,
    SAVE_GP_CMDS,
    SAVE_R_CMDS,
    SAVE_OX_CMDS,
    SAVE_OCTAVE_CMDS,
    SAVE_PYTHON_CMDS,
    SAVE_CONSOLE,
    SAVE_CMD_LOG,
    SAVE_FUNCTIONS,
    SAVE_FUNCTIONS_AS,
    SAVE_BOOT_DATA,
    SAVE_MARKERS,
    SAVE_LABELS,
    SAVE_GFN_SPEC,
    END_SAVE_OTHER, /* marker for end of other user-file saving */
    EDIT_FUNCTIONS,
    SET_PROG,
    SET_DIR,
    SET_WDIR,
    SET_FDIR,
    SET_DBDIR,
    SET_OTHER,
    SAVE_DATA_PKG,
    SAVE_REMOTE_DB,
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
    REMOTE_DATA_PKGS,
    NATIVE_SERIES,
    RATS_SERIES,
    PCGIVE_SERIES,
    REMOTE_SERIES,
    REMOTE_ADDONS,
    MAINWIN
};

#define BROWSER_ROLE(r) (r >= TEXTBOOK_DATA && r < MAINWIN)

enum pref_tabs {
    TAB_NONE = 0,
    TAB_MAIN,
    TAB_NET,
    TAB_PROGS,
    TAB_MPI,
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
    GZIPPED_DATA  = 1 << 6,
    SESSION_DATA  = 1 << 7
};

enum register_data_flags {
    DATAFILE_OPENED = 1,
    OPENED_VIA_CLI,
    OPENED_VIA_SESSION,
    DATA_APPENDED,
    NULLDATA_STARTED,
    DATA_PASTED
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
    CALC_PLOT,
    CALC_MAX
};

enum auto_dummies {
    TS_DUMMIES,
    PANEL_UNIT_DUMMIES,
    PANEL_TIME_DUMMIES
};

#endif /* GRETL_ENUMS_H */
