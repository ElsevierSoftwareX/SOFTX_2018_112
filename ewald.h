// ewald.h
//  заголовочный файл модул€, работающего с алгоритмом Ёвальда
#ifndef EWALD_H
#define EWALD_H

//simultaneous calculation of sinus and cosinus of the same angle
void sincos(double arg, double &s, double &c);

void init_ewald(Ewald *ew, Box *bx, Sim *sim, Atoms *atm);
double ewald_const(Atoms *atm, Spec *sp, Sim *sim, Box *box);
// return constant part of Columbic potential energy via Ewald method
//   (!) need to be recalculated only then volume or summ(q) are changed

double ewald_rec(Atoms *atm, Spec *sp, Box *bx, Sim *sim, Ewald *ew);
// calculate reciprocal part of Ewald summ and corresponding forces

double coul_iter(double r2, double &r, double chprd, double alpha, double &eng);
// real part of Ewald for each pair

void free_ewald(Ewald *ew, Atoms *atm);

#endif  /* EWALD_H */
