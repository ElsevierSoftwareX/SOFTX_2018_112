//MODULE vdw.cpp
//with HEADER vdw.h CONTAINS CONST FOR VDW READING AND APPLYING
#include <math.h>
#include <string.h>
#include <stdio.h>   //! temp, for loggin and debuging


#include "const.h"
#include "dataStruct.h"  // Sim, Box, Atoms ....
#include "vdw.h"

// Here and next fer_, fe_ and e_ functions must contain the same code
//   prefexis fer_, fe_ and e_ means force and energy calculation by r also, force and energy and energy
double fer_lj(double r2, double &r, VdW *vdw, double &eng)
// calculate force and energy by Lennard-Jones pair potential: U = 4e[(s/r)^12 - (s/r)^6]
//  r2 - square of distance, vdw - parameters (the same for next PP functions)
{
    double r2i = 1.0 / r2;
    double sr2 = vdw->p1 * r2i;
    double sr6 = sr2 * sr2 * sr2;

    eng += vdw->p0 * sr6 * (sr6 - 1.0);
    return vdw->p2 * r2i * sr6 * (2.0 * sr6 - 1.0);
}

double fe_lj(double r2, VdW *vdw, double &eng)
// calculate force and energy by Lennard-Jones pair potential: U = 4e[(s/r)^12 - (s/r)^6]
//  r2 - square of distance, vdw - parameters (the same for next PP functions)
{
    double r2i = 1.0 / r2;
    double sr2 = vdw->p1 * r2i;
    double sr6 = sr2 * sr2 * sr2;

    eng += vdw->p0 * sr6 * (sr6 - 1.0);
    return vdw->p2 * r2i * sr6 * (2.0 * sr6 - 1.0);
}

double e_lj(double r2, VdW *vdw)
// calculate energy by Lennard-Jones pair potential
{
    double r2i = 1.0 / r2;
    double sr2 = vdw->p1 * r2i;
    double sr6 = sr2 * sr2 * sr2;

    return vdw->p0 * sr6 * (sr6 - 1.0);
}

double er_lj(double r2, double r, VdW *vdw)
// calculate energy by Lennard-Jones pair potential
{
    double r2i = 1.0 / r2;
    double sr2 = vdw->p1 * r2i;
    double sr6 = sr2 * sr2 * sr2;

    return vdw->p0 * sr6 * (sr6 - 1.0);
}

double fer_buckingham(double r2, double &r, VdW *vdw, double &eng)
// calculate force and energy (with r) by Buckingham pair potential: U = A exp(-r/ro) - C/r^6
{
   double r2i = 1.0 / r2;
   double r4i = r2i * r2i;
   if (r == 0.0)    // calculate if unkonwn (zero) and use otherwise
     r = sqrt(r2);

   eng += vdw->p0 * exp(-r/vdw->p1) - vdw->p2 * r4i * r2i;
   return vdw->p0 * exp(-r/vdw->p1) / r / vdw->p1 - 6.0 * vdw->p2 * r4i * r4i;
}

double fe_buckingham(double r2, VdW *vdw, double &eng)
// calculate force and energy by Buckingham pair potential: U = A exp(-r/ro) - C/r^6
{
   double r2i = 1.0 / r2;
   double r  = sqrt(r2);
   double r4i = r2i * r2i;

   eng += vdw->p0 * exp(-r/vdw->p1) - vdw->p2 * r4i * r2i;
   return vdw->p0 * exp(-r/vdw->p1) / r / vdw->p1 - 6.0 * vdw->p2 * r4i * r4i;
}

double e_buckingham(double r2, VdW *vdw)
// calculate energy by Buckingham pair potential: U = A exp(-r/ro) - C/r^6
{
   double r2i = 1.0 / r2;
   double r  = sqrt(r2);
   double r4i = r2i * r2i;

   return vdw->p0 * exp(-r/vdw->p1) - vdw->p2 * r4i * r2i;
}

double er_buckingham(double r2, double r, VdW *vdw)
// calculate energy by Buckingham pair potential: U = A exp(-r/ro) - C/r^6
{
   double r2i = 1.0 / r2;
   double r4i = r2i * r2i;

   return vdw->p0 * exp(-r/vdw->p1) - vdw->p2 * r4i * r2i;
}

double fer_bhm(double r2, double &r, VdW *vdw, double &eng)
// calculate force and energy (with r) by Born-Huggins-Maier pair potential: U = Aexp[B(s-r)] - C/r^6 - D/r^8
{
   double r2i = 1.0 / r2;
   double r4i = r2i * r2i;
   if (r == 0.0)    // calculate if unkonwn (zero) and use otherwise
     r = sqrt(r2);

   eng += vdw->p0 * exp(vdw->p1*(vdw->p2 - r)) - vdw->p3 * r4i * r2i - vdw->p4 * r4i * r4i;
   return vdw->p0 * vdw->p1 * exp(vdw->p1*(vdw->p2 - r)) / r  - 6.0 * vdw->p3 * r4i * r4i - 8.0 * vdw->p4 * r4i * r4i * r2i;
}

double fe_bhm(double r2, VdW *vdw, double &eng)
// calculate force and energy by Born-Huggins-Maier pair potential: U = Aexp[B(s-r)] - C/r^6 - D/r^8
{
   double r2i = 1.0 / r2;
   double r  = sqrt(r2);
   double r4i = r2i * r2i;

   eng += vdw->p0 * exp(vdw->p1*(vdw->p2 - r)) - vdw->p3 * r4i * r2i - vdw->p4 * r4i * r4i;
   return vdw->p0 * vdw->p1 * exp(vdw->p1*(vdw->p2 - r)) / r  - 6.0 * vdw->p3 * r4i * r4i - 8.0 * vdw->p4 * r4i * r4i * r2i;
}

double e_bhm(double r2, VdW *vdw)
// calculate energy by Born-Huggins-Maier pair potential: U = Aexp[B(s-r)] - C/r^6 - D/r^8
{
   double r2i = 1.0 / r2;
   double r  = sqrt(r2);
   double r4i = r2i * r2i;

   return vdw->p0 * exp(vdw->p1*(vdw->p2 - r)) - vdw->p3 * r4i * r2i - vdw->p4 * r4i * r4i;
}

double er_bhm(double r2, double r, VdW *vdw)
// calculate energy by Born-Huggins-Maier pair potential: U = Aexp[B(s-r)] - C/r^6 - D/r^8
{
   double r2i = 1.0 / r2;
   double r4i = r2i * r2i;

   return vdw->p0 * exp(vdw->p1*(vdw->p2 - r)) - vdw->p3 * r4i * r2i - vdw->p4 * r4i * r4i;
}

double fer_746(double r2, double &r, VdW *vdw, double &eng)
// calculate force and energy (with r) by "746" pair potential (Staddford et.al.[]): U = A/r^7 - B/r^4 - C/r^6
{
   double r2i = 1.0 / r2;
   double r4i = r2i * r2i;
   double ri;
   if (r == 0.0)    // calculate from (r2i) if r unkonwn (zero) and calculate from r otherwise
     ri = sqrt(r2i);
   else
     ri = 1.0 / r;

   eng += r4i * (vdw->p0 * r2i * ri - vdw->p1 - vdw->p2 * r2i);
   return r4i * r2i * (7.0 * vdw->p0 * r2i * ri - 4.0 * vdw->p1 - 6.0 * vdw->p2 * r2i);
}


double fe_746(double r2, VdW *vdw, double &eng)
// calculate force and energy by "746" pair potential (Staddford et.al.[]): U = A/r^7 - B/r^4 - C/r^6
{
   double r2i = 1.0 / r2;
   double ri  = sqrt(r2i);
   double r4i = r2i * r2i;

   eng += r4i * (vdw->p0 * r2i * ri - vdw->p1 - vdw->p2 * r2i);
   return r4i * r2i * (7.0 * vdw->p0 * r2i * ri - 4.0 * vdw->p1 - 6.0 * vdw->p2 * r2i);
}

double e_746(double r2, VdW *vdw)
// calculate energy by "746" pair potential (Staddford et.al.[]): U = A/r^7 - B/r^4 - C/r^6
{
   double r2i = 1.0 / r2;
   double ri  = sqrt(r2i);
   double r4i = r2i * r2i;

   return r4i * (vdw->p0 * r2i * ri - vdw->p1 - vdw->p2 * r2i);
}

double er_746(double r2, double r, VdW *vdw)
// calculate energy by "746" pair potential (Staddford et.al.[]): U = A/r^7 - B/r^4 - C/r^6
{
   double r2i = 1.0 / r2;
   double ri  = 1.0 / r;
   double r4i = r2i * r2i;

   return r4i * (vdw->p0 * r2i * ri - vdw->p1 - vdw->p2 * r2i);
}
// END PP functions definition

int prepare_vdw(char *name, VdW &pp)
// prepare VdW parameters after reading from INPUT file
{
    double x;
    int res = 1;

    //fscanf(f, " %8s %8s %8s %lf %lf %lf", aname, bname, vdwnm, &pp.r2cut, &pp.p0, &pp.p1);
    //printf("res[%d]. a=%s; p0=%e; p1=%e p2=%e;\n", res, name, pp.p0, pp.p1, pp.p2);


    pp.r2cut *= r_scale; // convert external length units to MD units
    // user enter radii, but prog need sqaure of radii
    pp.r2cut = pp.r2cut * pp.r2cut;

    //DEFINE THE TYPE OF THE PAIR POTENTIAL AND DO PREPARATION:
    if (strcmp(name, lj_name) == 0)  // U = 4e[(s/r)^12 - (s/r)^6]
      {
         //printf("res[%d]lenard jones detected\n", res);
         pp.type = lj_type;

         //all parameters need to be converted into MD units
         pp.p0 *= 4 * E_scale;  // 4*epsilon in L-J pair potential

         pp.p1 *= r_scale;       // sigma
         pp.p1 = pp.p1 * pp.p1;  // sigma^2  in L-J pair potential

         pp.p2 = 6 * pp.p0;    // 24*epsilon for force calculation

         //p0 = 4e;  p1 = s^2;  p2 = 24e;
         pp.feng_r = fer_lj;
         pp.feng = fe_lj;
         pp.eng = e_lj;
         pp.eng_r = er_lj;

         pp.p3 = 0.0;
         pp.p4 = 0.0;
      }
    else if (strcmp(name, bh_name) == 0)  // U = A exp(-r/ro) - C/r^6
      {
         pp.type = bh_type;

         pp.p0 *= E_scale; // A  [eV]
         pp.p1 *= r_scale; // ro  [A]

         x = r_scale * r_scale;
         x = x * x * x; // r^6
         pp.p2 *= E_scale * x; // C [eV*A^6]

         pp.feng_r = fer_buckingham;
         pp.feng = fe_buckingham;
         pp.eng = e_buckingham;
         pp.eng_r = er_buckingham;

         pp.p3 = 0.0;
         pp.p4 = 0.0;
      }
    else if (strcmp(name, BHM_name) == 0)  // U = Aexp[B(s-r)] - C/r^6 - D/r^8
      {
         pp.type = BHM_type;

         pp.p0 *= E_scale; // A  [eV]
         pp.p1 /= r_scale; // B  [1/A]
         pp.p2 *= r_scale; // sigma  [A]

         x = r_scale * r_scale;
         x = x * x * x; // r^6
         pp.p3 *= E_scale * x; // C [eV*A^6]
         pp.p4 *= E_scale * x * r_scale * r_scale; // D [eV*A^8]

         pp.feng_r = fer_bhm;
         pp.feng = fe_bhm;
         pp.eng = e_bhm;
         pp.eng_r = er_bhm;
      }
    else if (strcmp(name, CuCl_name) == 0)  // U = A/r^7 - B/r^4 - C/r^6
      {
         pp.type = CuCl_type;

         x = r_scale * r_scale;
         x = x * x; // r_scale^4

         pp.p1 *= E_scale * x;  // B

         x = x * r_scale * r_scale; // r_scale^6;
         pp.p2 *= E_scale * x;  // C

         x = x * r_scale;
         pp.p0 *= E_scale * x;  // A

         //p0 = A; p1 = B;  p2 = C;
         pp.feng_r = fer_746;
         pp.feng = fe_746;
         pp.eng = e_746;
         pp.eng_r = er_746;

         pp.p3 = 0.0;
         pp.p4 = 0.0;

      }
    else
      res = 0;

    //printf("res[%d]\n", res);
    return res;
}
//end 'read_vdw' function


double vdw_iter(double r2, VdW *vdw, double &eng)
// calculate energy and return force of vdw iteraction  (   Fx/dx = -(1/r)*dU(r)/dr   )
//  r2 - square of distance
{
   double r2i, sr2, sr6;
   double ri, r4i;
   double r;

   switch (vdw->type)
    {
       case lj_type:   // U = 4e[(s/r)^12 - (s/r)^6]
         r2i = 1.0 / r2;
         sr2 = vdw->p1 * r2i;
         sr6 = sr2 * sr2 * sr2;

         eng += vdw->p0 * sr6 * (sr6 - 1.0);
         return vdw->p2 * r2i * sr6 * (2.0 * sr6 - 1.0);
         //break; /// break после return не имеет смысла

       case bh_type:    // U = A exp(-r/ro) - C/r^6
         r2i = 1.0 / r2;
         r  = sqrt(r2);
         r4i = r2i * r2i;

         eng += vdw->p0 * exp(-r/vdw->p1) - vdw->p2 * r4i * r2i;
         return vdw->p0 * exp(-r/vdw->p1) / r / vdw->p1 - 6.0 * vdw->p2 * r4i * r4i;

       case BHM_type:    // U = Aexp[B(s-r)] - C/r^6 - D/r^8
         r2i = 1.0 / r2;
         r  = sqrt(r2);
         r4i = r2i * r2i;

         eng += vdw->p0 * exp(vdw->p1*(vdw->p2 - r)) - vdw->p3 * r4i * r2i - vdw->p4 * r4i * r4i;
         return vdw->p0 * vdw->p1 * exp(vdw->p1*(vdw->p2 - r)) / r  - 6.0 * vdw->p3 * r4i * r4i - 8.0 * vdw->p4 * r4i * r4i * r2i;


       case CuCl_type:   // U = A/r^7 - B/r^4 - C/r^6
         r2i = 1.0 / r2;
         ri  = sqrt(r2i);
         r4i = r2i * r2i;
         //double sr2 = vdw->p1 * r2i;
         //double sr6 = sr2 * sr2 * sr2;

         eng += r4i * (vdw->p0 * r2i * ri - vdw->p1 - vdw->p2 * r2i);
         return r4i * r2i * (7.0 * vdw->p0 * r2i * ri - 4.0 * vdw->p1 - 6.0 * vdw->p2 * r2i);
         //break; /// break после return не имеет смысла
    }
}
// end 'vdw_iter' function


