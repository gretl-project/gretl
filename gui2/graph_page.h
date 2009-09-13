#ifndef GRAPH_PAGE_H
#define GRAPH_PAGE_H

int display_graph_page (void);

void clear_graph_page (int on_exit);

int graph_page_add_file (const char *fname);

int graph_page_get_n_graphs (void);

int in_graph_page (const char *fname);

int save_graph_page (const char *fname);

#endif /* GRAPH_PAGE_H */
