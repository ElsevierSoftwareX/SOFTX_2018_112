// ejump.h
#ifndef EJUMP_H
#define EJUMP_H

void init_ejump(Atoms *atm, Field *field, Sim *sim);
int ejump(Atoms *atm, Spec *sp, VdW ***vdws, Sim *sim, Box *bx);
// electron jumps routine, return the number of electron hops
int ejump_min(Atoms *atm, Spec *sp, VdW ***vdws, Sim *sim, Box *bx);
// electron jumps routine (chose minimal energy ), return the number of electron hops
int ejump_metr(Atoms *atm, Spec *sp, VdW ***vdws, Sim *sim, Box *bx);
int free_ejump(Sim *sim, Field *field);

#endif  /* EJUMP_H */
