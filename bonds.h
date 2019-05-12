#ifndef BONDS_H
#define BONDS_H

//int init_bonds(Sim *sim, int **nbonds, int ***bonds, int ***neigs/*, Sim *sim, Box *box*/, int def_nbonds);

int read_bondlist(Atoms *atm, Field *field);

void save_bondlist(char *fname, Field *field);

int read_bond(int id, FILE *f, Field *field);


int autobonding(Atoms *atm, Field *field, Sim *sim, Box *bx, int Naut, double rc, FILE *f, int *nbonds, int **bonds, int **neigs, int mxbond);

int read_linkage(FILE *f, int Nlnk, Field *field, int nBonds);
double bond_iter(double r2, Bond *bnd, double &eng);
double bond_iter_r(double r2, double &r, Bond *bnd, double &eng);

int create_bond(int at1, int at2, int type, Atoms *atm, Field *field);

void exec_bondlist(Atoms *atm, Field *field, Sim *sim, Box *bx);

//void bonding(Atoms *atm, Spec *sp, VdW ***vdws, Sim *sim, Box *bx, Bond *btps, Field *field);


int free_bonds(Sim *sim, int **nbonds, int ***bonds, int ***neigs/*, Sim *sim, Box *box*/, int nAtm);

#endif // BONDS_H
