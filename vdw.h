#ifndef VDW_H
#define VDW_H

//constants:
//  names for user (in species.txt)
const char lj_name[5] = "lnjs";  // lennard-jones
const char bh_name[5] = "bkhm";  // buckingham          U = A exp(-r/ro) - C/r^6
const char CuCl_name[5] = "p746";  // potential 7-4-6   U = A/r^7 - B/r^4 - C/r^6
const char BHM_name[5] = "bmhs"; // Born-Mayer-Huggins U = Aexp[B(s-r)] - C/r^6 - D/r^8

// identifers
const int  lj_type = 0;
const int  bh_type = 1;
const int  CuCl_type = 2;
const int  BHM_type = 3;

// names[type]
const char vdw_names[4][20] = {"Lenard-Jones", "Buckingham", "CuCl(7-4-6)", "Born-Mayer-Huggins"};


int prepare_vdw(char *name, VdW &pp);
// prepare VdW parameters after reading from INPUT file

double vdw_iter(double r2, VdW *vdw, double &eng);
// calculate energy and return force of vdw iteraction
//  r2 - square of distance

// pair potential functions:
double fer_lj(double r2, double &r, VdW *vdw, double &eng);
double fe_lj(double r2, VdW *vdw, double &eng);
double e_lj(double r2, VdW *vdw);
double er_lj(double r2, double r, VdW *vdw);
double fer_buckingham(double r2, double &r, VdW *vdw, double &eng);
double fe_buckingham(double r2, VdW *vdw, double &eng);
double e_buckingham(double r2, VdW *vdw);
double er_buckingham(double r2, double r, VdW *vdw);
double fer_bhm(double r2, double &r, VdW *vdw, double &eng);
double fe_bhm(double r2, VdW *vdw, double &eng);
double e_bhm(double r2, VdW *vdw);
double er_bhm(double r2, double r, VdW *vdw);
double fer_746(double r2, double &r, VdW *vdw, double &eng);
double fe_746(double r2, VdW *vdw, double &eng);
double e_746(double r2, VdW *vdw);
double er_746(double r2, double r, VdW *vdw);

#endif /* VDW_H */
