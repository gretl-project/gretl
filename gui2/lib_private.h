/* For stuff shared by library,c with cmdstack.c */

#ifndef LIB_PRIVATE_H
#define LIB_PRIVATE_H

int check_specific_command (char *s);

char *get_lib_cmdline (void);

CMD *get_lib_cmd (void);

#endif /* LIB_PRIVATE_H */
