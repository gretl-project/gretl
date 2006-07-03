/* For stuff shared by library.c with cmdstack.c and session.c */

#ifndef LIB_PRIVATE_H
#define LIB_PRIVATE_H

int check_specific_command (char *s);

char *get_lib_cmdline (void);

CMD *get_lib_cmd (void);

void lib_cmd_destroy_context (void);

void lib_modelspec_free (void);

#endif /* LIB_PRIVATE_H */
