#ifndef FORMATTER_H
#define FORMATTER_H

typedef struct _ARGLIST ARGLIST;
typedef struct _OPTION OPTION;
typedef struct _GUI_ACCESS GUI_ACCESS;
typedef struct _COMMAND COMMAND;
typedef struct _COMMANDLIST COMMANDLIST;

struct _ARGLIST {
    int n_args;
    char **args;
};

struct _OPTION {
    char *flag;
    char *effect;
};

struct _GUI_ACCESS {
    char *menu_path;
    char *other_access;
};



struct _COMMAND {
    char *name;
    char *xref;
    char *section;
    ARGLIST arglist;
    ARGLIST optargs;
    int n_options;
    OPTION *options;
    int n_examples;
    char **examples;
    int n_descrip_paras;
    char **descrip;
    char *optnotes;
    char *addendum;
    GUI_ACCESS gui_access;
};

void xml_format_command (COMMAND *cmd, FILE *fp);

#endif /* FORMATTER_H */
