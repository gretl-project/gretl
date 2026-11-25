#ifndef PANEL_MARKERS_H
#define PANEL_MARKERS_H

/* Offer a dialog for constructing "nice" observation markers
   (e.g. "FR:2026") for a panel dataset from a unit-identifier
   series and a period-identifier series.
*/

void panel_markers_dialog (const int *ulist, const int *tlist);

int get_candidate_lists (int **ulist, int **tlist);

int panel_markers_choice (void);

#endif /* PANEL_MARKERS_H */
