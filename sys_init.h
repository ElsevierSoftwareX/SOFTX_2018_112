#ifndef SYS_INIT_H
#define SYS_INIT_H

//create atoms struct

//Atoms* init_atoms(int n);
Field* init_field(Sim *sim);
Atoms* init_atoms_and_box(Field *field, Sim *sim, Box *&box);

int init_md(Atoms *atm, Field *field, Sim *sim, Box *bx);
void prepare_md(Atoms *atm, Sim *sim, Field *field);
void info_md(Sim *sim); // return information about simulation

void free_md(Atoms *atm, Field *field, Sim *sim);
int init_neighbors(Atoms *atm, Sim *sim);
void free_neighbors(Atoms *atm, Sim *sim);
void free_field(Field *field);
void free_atoms(Atoms **atm);


#endif  /* SYS_INIT_H */
