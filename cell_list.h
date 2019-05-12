#ifndef CELL_LIST_H
#define CELL_LIST_H

int init_clist(Atoms* atm, Sim *sim, Box *bx, double rCut);

int cell_index_sim(double x, double y, double z, Sim *sim);
// return number of cell by real{x;y;z} coordintates and sim record

void free_clist(Sim *sim);

#endif // CELL_LIST_H
