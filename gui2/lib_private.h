/* For stuff shared by library.c with cmdstack.c and session.c */

#ifndef LIB_PRIVATE_H
#define LIB_PRIVATE_H

int parse_lib_command (void);

char *get_lib_cmdline (void);

void lib_cmd_destroy_context (void);

#endif /* LIB_PRIVATE_H */
