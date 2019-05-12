// UNIT FOR EWALD ALGORITMS
//   ...also SINCOS function is also described here
#include <stdlib.h>  // malloc, alloc, rand, NULL
#include <math.h>   // log, sqrt
//#include <stdio.h>   //!временно, для вывода тестов

#include "const.h"
#include "dataStruct.h"  // Sim, Box, Atoms ....
#include "utils.h"  // int_size, pointer_size, etc...
#include "ewald.h"

void sincos(double arg, double &s, double &c)
// function of simultaneous calculation sin and cos (faster then sin and cos separately)
{
    asm("fsincos" : [c]"=t"(c), [s]"=u"(s) : "0"(arg));
}

void init_ewald(Ewald *ew, Box *bx, Sim *sim, Atoms *atm)
// create arrays for ewald recipoal space calculation:
//   el - exp(i 2pi x/a * l); em - exp (i 2 pi y/b * m); en - exp (i 2pi z/c * n)
//   -c - cosinus (or real) part; -s - sinus (or imaginary) part
//   lm = el * em; ck = q * el * em * en

{
   int i;
   //double **arr, **arr1;
   int Nat = atm->nAt;
   int kx = sim->kx;
   int ky = sim->ky;
   int kz = sim->kz;

   //constants definition:
   ew->alpha = sim->alpha;
   ew->scale = 2 * twopi * bx->rvol * Fcoul_scale / sim->eps;
   ew->scale2 = 2 * ew->scale;
   ew->mr4a2 = -0.25 / ew->alpha / ew->alpha; // -1/(4a^2)

   ew->rkcut = kx * bx->ip1;
   //! найти её один раз и вперёд - наименьшая среди ka * ipa
   if (ew->rkcut > ky * bx->ip2)
     ew->rkcut = ky * bx->ip2;
   //printf("ip2=%f  2pi=%f rkcut=%f\n", bx->ip2, twopi, ew->rkcut);
   if (ew->rkcut > kz * bx->ip3)
     ew->rkcut = kz * bx->ip3;
   //printf("ip3=%f  2pi=%f rkcut=%f\n", bx->ip3, twopi, ew->rkcut);
   ew->rkcut *= twopi * 1.05; // according to DL_POLY source
   ew->rkcut2 = ew->rkcut * ew->rkcut;
   //printf("ip1=%f  2pi=%f rkcut=%f rkcut2=%f\n", bx->ip1, twopi, ew->rkcut, ew->rkcut2);


   // create elc and els arrays: only 2 elements as we will fill [0..N][0] elements by [0..N][k] elements
   ew->elc = (double**)malloc(Nat * pointer_size);
   ew->els = (double**)malloc(Nat * pointer_size);
   //arr = *elc;
   //arr1 = *els;
   for (i = 0; i < Nat; i++)
    {
      ew->elc[i] = (double*)malloc(2 * double_size);
      ew->els[i] = (double*)malloc(2 * double_size);
    }

   // create emc and ems arrays for every ky
   ew->emc = (double**)malloc(Nat * pointer_size);
   ew->ems = (double**)malloc(Nat * pointer_size);
   //arr = *emc;
   //arr1 = *ems;
   for (i = 0; i < Nat; i++)
    {
      ew->emc[i] = (double*)malloc(ky * double_size);
      ew->ems[i] = (double*)malloc(ky * double_size);

      // constant part:  exp(i * 2pi * 0 * m)
      ew->emc[i][0] = 1.0;
      ew->ems[i][0] = 0.0;
    }

   // create enc and ens arrays for every kz
   ew->enc = (double**)malloc(Nat * pointer_size);
   ew->ens = (double**)malloc(Nat * pointer_size);
   //arr = *enc;
   //arr1 = *ens;
   for (i = 0; i < Nat; i++)
    {
      ew->enc[i] = (double*)malloc(kz * double_size);
      ew->ens[i] = (double*)malloc(kz * double_size);

      // constant part:  exp(i * 2pi * 0 * n)
      ew->enc[i][0] = 1.0;
      ew->ens[i][0] = 0.0;
    }

   // create lmc, lms, ckc, cks buffer-arrays:
   ew->lmc = (double*)malloc(Nat * double_size);
   ew->lms = (double*)malloc(Nat * double_size);
   ew->ckc = (double*)malloc(Nat * double_size);
   ew->cks = (double*)malloc(Nat * double_size);
}
// end 'init_ewald' function

double ewald_const(Atoms *atm, Spec *sp, Sim *sim, Box *box)
// return constant part of Columbic potential energy via Ewald method
//   (!) need to be recalculated only then volume or summ(q) are changed
{
   int i;
   double q;
   double sq = 0.0;
   double eng = 0.0;

   for (i = 0; i < atm->nAt; i++)
     {
        q = sp[atm->types[i]].charge;
        sq += q;
        eng += q*q;
     }

   eng *= (-1.0) * sim->alpha / sqrtpi;  // constant part of Ewald summation
   q = -0.5 * pi * (sq * sq / sim->alpha / sim->alpha) * box->rvol;

   return Fcoul_scale * (eng + q) / sim->eps;
}
// end 'ewald_const' function

double ewald_rec(Atoms *atm, Spec *sp, Box *bx, Sim *sim, Ewald *ew)
// calculate reciprocal part of Ewald summ and corresponding forces
{
   int i, l, m, n;
   int mmin = 0;
   int nmin = 1;
   double x, ch; // temporary variables (for complex *=,  for charge)
   double rkx, rky, rkz, rk2, akk, eng, sumC, sumS;
   int Nat = atm->nAt;
   int kx = sim->kx;
   int ky = sim->ky;
   int kz = sim->kz;
   // double rkcut2 = sim->rkcut2;  //! В DL_POLY это вычисляемая величина
   //printf("ewald_rec Nat=%d kx=%d ky=%d kz=%d  rkut2=%f\n", Nat, kx, ky, kz, ew->rkcut2);

   eng = 0.0;
   //! тут всё верно для прямоугольной геометрии. Если ячейка будет кривая, код нужно править
   for (i = 0; i < Nat; i++)
     {
        // exp (i 2pi * 0 * l) for em- and en- arrays this step omitted as they set in 'init_ewald'
        //   el- arrays need to refresh (according cycle by l)
        ew->elc[i][0] = 1.0;
        ew->els[i][0] = 0.0;

        // exp (ikr)
        sincos(twopi * atm->xs[i] * bx->ra, ew->els[i][1], ew->elc[i][1]);
        sincos(twopi * atm->ys[i] * bx->rb, ew->ems[i][1], ew->emc[i][1]);
        sincos(twopi * atm->zs[i] * bx->rc, ew->ens[i][1], ew->enc[i][1]);
        //printf("a=%f  sin:%f  cos: %f\n", twopi * atm[i].x * bx->ra, els[i][1], elc[i][1]);
     }

    // fil exp(iky) array by complex multiplication
    for (l = 2; l < ky; l++)
      for (i = 0; i < Nat; i++)
        {
           ew->emc[i][l] = ew->emc[i][l-1] * ew->emc[i][1] - ew->ems[i][l-1] * ew->ems[i][1];
           ew->ems[i][l] = ew->ems[i][l-1] * ew->emc[i][1] + ew->emc[i][l-1] * ew->ems[i][1];
        }

    // fil exp(ikz) array by complex multiplication
    for (l = 2; l < kz; l++)
      for (i = 0; i < Nat; i++)
        {
           ew->enc[i][l] = ew->enc[i][l-1] * ew->enc[i][1] - ew->ens[i][l-1] * ew->ens[i][1];
           ew->ens[i][l] = ew->ens[i][l-1] * ew->enc[i][1] + ew->enc[i][l-1] * ew->ens[i][1];
        }

    // MAIN LOOP OVER K-VECTORS:
    for (l = 0; l < kx; l++)
      {
         rkx = l * twopi * bx->ra; // only for rect geometry!
         // move exp(ikx[l]) to ikx[0] for memory saving (ikx[i>1] are not used)
         if (l == 1)
           for (i = 0; i < Nat; i++)
             {
                ew->elc[i][0] = ew->elc[i][1];
                ew->els[i][0] = ew->els[i][1];
             }
         else if (l > 1)
           for (i = 0; i < Nat; i++)
             {
                // exp(ikx[0]) = exp(ikx[0]) * exp(ikx[1])
                x = ew->elc[i][0];
                ew->elc[i][0] = x * ew->elc[i][1] - ew->els[i][0] * ew->els[i][1];
                ew->els[i][0] = ew->els[i][0] * ew->elc[i][1] + x * ew->els[i][1];
             }

         //ky - loop:
         for (m = mmin; m < ky; m++)
           {
              rky = m * twopi * bx->rb;
              //fil temp arrays for keeping e^ikx * e^iky
              if (m >= 0)
                for (i = 0; i < Nat; i++)
                  {
                     ew->lmc[i] = ew->elc[i][0] * ew->emc[i][m] - ew->els[i][0] * ew->ems[i][m];
                     ew->lms[i] = ew->els[i][0] * ew->emc[i][m] + ew->ems[i][m] * ew->elc[i][0];
                  }
              else // for negative ky give complex adjustment to positive ky:
                for (i = 0; i < Nat; i++)
                  {
                     ew->lmc[i] = ew->elc[i][0] * ew->emc[i][-m] + ew->els[i][0] * ew->ems[i][-m];
                     ew->lms[i] = ew->els[i][0] * ew->emc[i][-m] - ew->ems[i][-m] * ew->elc[i][0];
                  }

              //kz - loop:
              for (n = nmin; n < kz; n++)
                {
                   rkz = n * twopi * bx->rc;
                   // rk2 = (2pi * l / a)^2 + (2pi * m / b)^2 + (2pi * n / c)^2   !only for rectangular geometry!
                   //rk2 = twopi2 * (l * l * bx->ra2 + m * m * bx->rb2 + n * n * bx->rc2);
                   rk2 = rkx * rkx + rky * rky + rkz * rkz;
                   //! у нас cuttof и rk2 возможно в разных единицах измерения, надо это провентилировать
                   //printf("rk2 * rkcut2 :  %f  *  %f\n", rk2, ew->rkcut2);
                   if (rk2 < ew->rkcut2) // cutoff
                     {
                        // calculate summ(ikr*q[iAt])
                        sumC = 0; sumS = 0;
                        if (n >= 0)
                          for (i = 0; i < Nat; i++)
                            {
                               ch = sp[atm->types[i]].charge;

                               ew->ckc[i] = ch * (ew->lmc[i] * ew->enc[i][n] - ew->lms[i] * ew->ens[i][n]);
                               ew->cks[i] = ch * (ew->lms[i] * ew->enc[i][n] + ew->lmc[i] * ew->ens[i][n]);

                               sumC += ew->ckc[i];
                               sumS += ew->cks[i];
                            }
                        else // for negative kz give complex adjustment to positive kz:
                          for (i = 0; i < Nat; i++)
                            {
                               ch = sp[atm->types[i]].charge;

                               ew->ckc[i] = ch * (ew->lmc[i] * ew->enc[i][-n] + ew->lms[i] * ew->ens[i][-n]);
                               ew->cks[i] = ch * (ew->lms[i] * ew->enc[i][-n] - ew->lmc[i] * ew->ens[i][-n]);

                               sumC += ew->ckc[i];
                               sumS += ew->cks[i];
                            }

                        //energy and force calculation!
                        akk = exp(rk2*ew->mr4a2) / rk2;
                        eng += akk * (sumC * sumC + sumS * sumS);
                        //printf("akk=%f, sumC=%f,  sumS=%f\n", akk, sumC, sumS);

                        for (i = 0; i < Nat; i++)
                          {
                             //! rkx = 2pi * l / a - посмотреть что быстрее, вводить rkx, rky и rkz или сделать как щас
                             x = akk * (ew->cks[i] * sumC - ew->ckc[i] * sumS);
                             // in DLPOLY there is a factor  4 * pi / V / (4piee0):
                             x *= ew->scale2;
                             atm->fxs[i] += rkx * x;
                             atm->fys[i] += rky * x;
                             atm->fzs[i] += rkz * x;
                          }

                        //printf("eng=%f, \n", eng);
                     }
                } // end n-loop (over kz-vectors)

              nmin = 1 - kz;

           } // end m-loop (over ky-vectors)

         mmin = 1 - ky;

      }  // end l-loop (over kx-vectors)


    //printf("rvol=%f, coul-scale=%f,  eng=%f\n", bx->rvol, Fcoul_scale, eng);

    //! надо ещё добавить эти постоянные члены
    //printf("ewald_rec=%f\n", scale * eng);
    return ew->scale * eng;
}
// end 'ewald_rec' function

double coul_iter(double r2, double &r, double chprd, double alpha, double &eng)
// calculate energy and return force of real part Coulombic iteraction via Ewald procedure
//  r2 - square of distance, chi, chj - charges of i and j particles
//  chprd - production of charges
//  eng - for saving energy
//! тут надо ещё ввести эпсилон в закон кулона
{
   //double r;
   double ar; //alpha * r
   double erfcar; // erfc(alpha * r);
   double kqq = chprd * Fcoul_scale; // q[i]*q[j]*1/4pie0;

   //brute force calc:
   if (r == 0)
     r = sqrt(r2);  // if r is unknown, calculate it
   //! надо предусмотреть вариант, когда r2 неизвестно
   ar = alpha * r;
   erfcar = erfc(ar);

   //printf("r2=%f  ar=%f  effcar=%f  E=%4.2E,    F=%4.2E\n", r2, ar, erfcar, kqq * erfcar /  r, kqq / r / r2 * (erfcar + 2 * ar / sqrtpi * exp(-ar*ar)));

   eng += kqq * erfcar /  r;
   //return 0;
   return kqq / r / r2 * (erfcar + 2 * ar / sqrtpi * exp(-ar*ar));
}
// end 'coul_iter' function

void free_ewald(Ewald *ew, Atoms *atm)
// free memory from ewald algorithm
//   el - exp(i 2pi x/a * l); em - exp (i 2 pi y/b * m); en - exp (i 2pi z/c * n)
//   -c - cosinus (or real) part; -s - sinus (or imaginary) part
//   lm = el * em; ck = q * el * em * en
{
   int i;
   //double **arr, **arr1;

   // free elc and els arrays: only 2 elements as we will fill [0..N][0] elements by [0..N][k] elements
   for (i = 0; i < atm->nAt; i++)
     {
        delete[] ew->elc[i];
        delete[] ew->els[i];
        delete[] ew->emc[i];
        delete[] ew->ems[i];
        delete[] ew->enc[i];
        delete[] ew->ens[i];
     }

   delete[] ew->elc;
   delete[] ew->els;
   delete[] ew->emc;
   delete[] ew->ems;
   delete[] ew->enc;
   delete[] ew->ens;

   // free lmc, lms, ckc, cks buffer-arrays:
   delete[] ew->lmc;
   delete[] ew->lms;
   delete[] ew->ckc;
   delete[] ew->cks;
}
// end 'free_ewald' function
