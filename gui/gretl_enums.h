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
    TDISAGG,
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
    VIEW_PKG_SAMPLE,
    VIEW_BUNDLE,
    VIEW_DBNOMICS,
    VIEW_FILE,
    VIEW_DOC,
    VIEW_BIBITEM,
    VIEW_DBSEARCH,
    DATA_REPORT,
    SCRIPT_OUT,
    FNCALL_OUT,
    CONSOLE,
    EDIT_HEADER,
    EDIT_HANSL,
    EDIT_NOTES,
    EDIT_PKG_CODE,
    EDIT_PKG_SAMPLE,
    EDIT_PKG_HELP,
    EDIT_PKG_GHLP,
    EDIT_GP,
    EDIT_R,
    EDIT_OX,
    EDIT_OCTAVE,
    EDIT_PYTHON,
    EDIT_STATA,
    EDIT_JULIA,
    EDIT_DYNARE,
    EDIT_LPSOLVE,
    EDIT_X12A,
    EDIT_SPEC,
    EDIT_MAX,
    CMD_HELP,
    GUI_HELP,
    FUNC_HELP,
    CMD_HELP_EN,
    GUI_HELP_EN,
    FUNC_HELP_EN,
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
    ALAGSEL,
    VLAGSEL,
    VAROMIT,
    DEFINE_LIST,
    DEFINE_MATRIX,
    PANEL_WLS,
    FE_LOGISTIC,
    PANEL_MODE,
    TS_TO_PANEL,
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
    REGLS,
    REGLS_ADV,
    REGLS_PLOTSEL,
    BWFILTER,
    POLYWEIGHTS,
    EMAFILTER,
    X12AHELP,
    MAILHELP,
    LOESS,
    NADARWAT,
    SSHEET,
    CLUSTER,
    GUI_FUNCS,
    MENU_ATTACH,
    DAILY_PURGE,
    BUILD_PKG,
    PKG_FILES,
    PKG_DEPS,
    EDITOR,
    MIDAS_LIST,
    MIDAS_PARM,
    PKGHELP,
    DBNHELP,
    MAPHELP,
    MDHELP,
    KALMAN,
    GUI_CMD_MAX
};

#define help_role(r) (r >= CMD_HELP && r <= FUNC_HELP_EN)
#define english_help_role(r) (r >= CMD_HELP_EN && r <= FUNC_HELP_EN)

#define editing_hansl(r) (r == EDIT_HANSL || \
			  r == EDIT_PKG_CODE ||	\
			  r == EDIT_PKG_SAMPLE)

#define editing_alt_script(r) (r >= EDIT_R && r <= EDIT_LPSOLVE)

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
    OPEN_SPEC,
    OPEN_ANY,
    OPEN_MAP,
    UPLOAD_PKG,
    INSTALL_PKG,
    END_OPEN,        /* marker for end of file open section */
    AUTO_SAVE_DATA,
    SAVE_DATA,
    SAVE_DATA_AS,
    EXPORT_OCTAVE,
    EXPORT_R,
    EXPORT_CSV,
    EXPORT_DAT,
    EXPORT_DTA,
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
    SAVE_STATA_CMDS,
    SAVE_JULIA_CODE,
    SAVE_DYNARE_CODE,
    SAVE_LPSOLVE_CODE,
    SAVE_SPEC_FILE,
    SAVE_X13_SPC,
    SAVE_HELP_TEXT,
    SAVE_CONSOLE,
    SAVE_CMD_LOG,
    SAVE_FUNCTIONS,
    SAVE_FUNCTIONS_AS,
    SAVE_BOOT_DATA,
    SAVE_MARKERS,
    SAVE_LABELS,
    SAVE_GFN_SPEC,
    SAVE_GFN_ZIP,
    SAVE_MAP,
    WRITE_MAP,
    END_SAVE_OTHER, /* marker for end of other user-file saving */
    EDIT_FUNCTIONS,
    SET_PROG,
    SET_DIR,
    SET_WDIR,
    SET_FDIR,
    SET_DBDIR,
    SET_OTHER,
    SELECT_PDF,
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
    DBNOMICS_TOP,
    DBNOMICS_DB,
    REMOTE_DATA_PKGS,
    NATIVE_SERIES,
    RATS_SERIES,
    PCGIVE_SERIES,
    REMOTE_SERIES,
    DBNOMICS_SERIES,
    ADDONS_FILES,
    PKG_REGISTRY,
    MAINWIN
};

#define BROWSER_ROLE(r) (r >= TEXTBOOK_DATA && r < MAINWIN)

enum pref_tabs {
    TAB_NONE = 0,
    TAB_MAIN,
    TAB_PROGS,
    TAB_EDITOR,
    TAB_NET,
    TAB_VCV,
#ifdef HAVE_MPI
    TAB_MPI,
#endif
    TAB_MAX
};

enum clipstuff {
    TARGET_UTF8_STRING,
    TARGET_STRING,
    TARGET_TEXT,
    TARGET_COMPOUND_TEXT,
    TARGET_RTF,
    TARGET_SVG,
    TARGET_EMF,
    TARGET_EPS,
    TARGET_PDF,
    TARGET_PNG,
    TARGET_HTM
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
    DATAFILE_OPENED    = 1 << 0,
    OPENED_VIA_CLI     = 1 << 1,
    OPENED_VIA_SESSION = 1 << 2,
    DATA_APPENDED      = 1 << 3,
    NULLDATA_STARTED   = 1 << 4,
    DATA_PASTED        = 1 << 5,
    FOCUS_CONSOLE      = 1 << 6
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
    FILE_LIST_GFN,
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
    PANEL_TIME_DUMMIES,
    DISCRETE_DUMMIES
};

enum script_output_policies {
    OUTPUT_POLICY_UNSET,
    OUTPUT_POLICY_REPLACE,
    OUTPUT_POLICY_APPEND,
    OUTPUT_POLICY_NEW_WINDOW
};

enum icon_size {
    ICON_SIZE_AUTO,
    ICON_SIZE_SMALL,
    ICON_SIZE_MEDIUM
};

#endif /* GRETL_ENUMS_H */
