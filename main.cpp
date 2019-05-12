//  //! стоит обратить внимание
//#define STDOUT 1  // вывод на экран промежуточных результатов

#include <stdio.h>   // FILE
#include <stdlib.h>  // malloc, alloc, rand, NULL
#include <time.h>   // time
#include <math.h>   // log, sqrt, fabs
//#include <strings.h> // strcpy, strcat
//#include <conio.h>
//#include <cstring> // strcpy strcat

#include "const.h"      // pi, 2pi, physical const, etc...
#include "dataStruct.h"  // Sim, Box, Atoms ....
#include "sys_init.h"   // init_md, free_md
#include "temperature.h"  // termostates etc
#include "utils.h"  // int_size, pointer_size, etc...
#include "out_md.h" // md output
#include "vdw.h"     // vdw_iter and vdw constants
#include "rdf.h"    // function for work with rdf
#include "ewald.h"   // ewald procedure
#include "box.h"     // rect_periodic, rect_box, box_prop
#include "cell_list.h"  // cell list method
#include "ejump.h"  // electron transfer functions
#include "bonds.h"
#include "angles.h"

void clear_force(Atoms *atm, Spec *sp, Sim *sim, Box *bx)
// set forces to zero (or according to external electric fieled) and
//     number of pairs(nPs), total probability (tPbs) for e-jump procedure
{
  int i, j;
  int N = atm->nAt;

  for (i = 0; i < N; i++)
    {
      //clear neighbors list
      //! может быть это не всегда и нужно
      for (j = 0; j < sim->nNbors[i]; j++)
        sim->tnbors[i][j] = 0;
      sim->nNbors[i] = 0;

      //external constant electric field
      //! ввести в бокс, например напряжённость dU/dr, чтобы не вычислять каждй раз
      atm->fxs[i] = -sp[atm->types[i]].charge * sim->Ux * bx->ra;   // F = -q dU/dr
      atm->fys[i] = 0.0;
      atm->fzs[i] = 0.0;

      // for e-jump procedure
      //sim->nPairs[i] = 0;
      //tPbs[i] = 0;

    }
}
// end 'clear_force' function

void save_neigh(int iat, int jat, double r, double r2, Sim *sim)
// save a neighbour and its type in neighbors arrays
{
   int Ni = 0;
   int Nj = 0;

   if (r == 0.0)
     r = sqrt(r2);

   if ((sim->nNbors[iat] < sim->maxNbors) || (sim->nNbors[jat] < sim->maxNbors))
    {
         Ni = sim->nNbors[iat];
         Nj = sim->nNbors[jat];

         sim->nbors[iat][Ni] = jat;
         sim->nbors[jat][Nj] = iat;
         sim->distances[iat][Ni] = r;
         sim->distances[jat][Nj] = r;

         // a neighbor in VdW range
         if (r2 <= sim->mxRvdw2)
           {
              sim->tnbors[iat][Ni] |= 1 << 1;  // the first bit is reserved for VdW neighbour
              sim->tnbors[jat][Nj] |= 1 << 1;
           }

         // a neighbor in eJump range
         if (r2 <= sim->r2Elec)
           {
              sim->tnbors[iat][Ni] |= 1 << 2;  // the second bit is reserved for ejump neighbour
              sim->tnbors[jat][Nj] |= 1 << 2;
           }

         // a neighbor in eJump range
         if (r2 <= sim->r2Bond)
           {
              sim->tnbors[iat][Ni] |= 1 << 3;  // the third bit is reserved for can form bond
              sim->tnbors[jat][Nj] |= 1 << 3;
           }

         /*if ((iat<0)||(iat > sim->nAt))
           printf("iat=%d\n", iat);
         if ((jat<0)||(jat > sim->nAt))
           printf("jat=%d\n", jat);
         */
         sim->nNbors[iat]++;
         sim->nNbors[jat]++;
    }
   else
     printf("WARNING[111]: the maximal number of neighbors is reached\n");

}

void pair_coul_vdw(int i, int j, Atoms *atm, Field *field, Box *bx, Sim *sim)
// calculate energy and determine forces for ion pairs in case of VdW and Coul forces only
{
   double r, r2, f, r2b;//, eng_vdw, eng_real;
   double dx = atm->xs[i] - atm->xs[j];
   double dy = atm->ys[i] - atm->ys[j];
   double dz = atm->zs[i] - atm->zs[j];
   int it, jt, bt;
   int link, sec;    // id of atom type, which can link and the second atom
   VdW *vdw;
   Spec* sp = field->species;

   // periodic boundary for rectangular geometry
   rect_periodic(dx, dy, dz, bx);  //! need be replaced for NOT rectangular geomery
   r2 = dx*dx + dy*dy + dz*dz;
   if (r2 <= sim->r2Max)
   {
      it = atm->types[i];
      jt = atm->types[j];
      vdw = field->vdws[it][jt];

      r = 0.0;
      // eng_vdw = 0.0; eng_real = 0.0;
      f = 0.0;
      if (sp[it].charged)  //! ещё одна возможная оптимизация, составить массив произведений зарядов [i][j]
        if (sp[jt].charged) //! хотя не факт, что перемещение по двумерному массиву будет быстрее, чем два IF
          f += coul_iter(r2, r, sp[it].charge * sp[jt].charge, sim->alpha, sim->engEwaldReal);

      if (vdw != NULL)
        if (r2 <= vdw->r2cut)
          //f += vdw -> feng(r2, vdw, eng);
          f += vdw -> feng_r(r2, r, vdw, sim->engVdW);

      //! verify bond creating!
      //! для оптимизации мб поставить сюда проверку на самый дальний радиус связывания?
      if (sp[it].canBond || sp[jt].canBond)
      {
        if (sp[it].canBond)
        {
           link = it;
           sec = jt;
        }
        else
        {
           link = jt;
           sec = it;
        }

        bt = sp[link].bondKeys[sec];
        //! может опять же для оптимизации сделать связь с нулевым индексом пустой и символизирующей отсутствие, тогда не прийдется вычитать потом
        if (bt) // zero bt means no bond
         if (atm->parents[i] != j)    //! этого раньше не было, но по логике ведь должно быть?
          if (atm->parents[j] != i)  // exclude double bonding
          {
             //if (sec != 3)
             //  printf("try to link wrong types(sec): %d %d\n", link, sec);
             //if ((link != 5)&&(link != 8))
             //  printf("try to link wrong types(link): %d %d\n", link, sec);

             //! вообще в bondKeys должно хранится
             bt--;  // as bt=0 reserved as no bond
             r2b = field->bdata[bt].rMax2;
             if (r2 <= r2b)
             {
                  create_bond(i, j, bt, atm, field);
                //printf("aft create bond\n");
             }
          }
      }
      //end bond creating

      if (isnan(f))
        {printf("force between pair %d and %d  is NAN\n", i, j);exit(1);}
      else
        if (f*f > 1e10)
          {
              printf("force between pair %s[%d] and %s[%d] (r2=%f) is colossal(%f)\n", sp[it].name, i, sp[jt].name, j, r2, f);
              //exit(1);
              return;
          }

      // forces update
      atm->fxs[i] += f * dx;
      atm->fxs[j] -= f * dx;
      atm->fys[i] += f * dy;
      atm->fys[j] -= f * dy;
      atm->fzs[i] += f * dz;
      atm->fzs[j] -= f * dz;
      //return eng;
   }
   //else return 0.0; // too far to interact
}

void pair_coul_vdw_lst(int i, int j, Atoms *atm, Field *field, Box *bx, Sim *sim)
// calculate energy and determine forces for ion pairs in case of VdW and Coul forces only (with neighbors list filling)
{
   double r, r2, f;//, eng;
   double dx = atm->xs[i] - atm->xs[j];
   double dy = atm->ys[i] - atm->ys[j];
   double dz = atm->zs[i] - atm->zs[j];
   int it, jt;
   int Ni, Nj; //! я не понимаю почему, но если это не объявить тут, программа вылетает
   VdW *vdw;
   Spec* sp = field->species;

   // periodic boundary for rectangular geometry
   rect_periodic(dx, dy, dz, bx);  //! need be replaced for NOT rectangular geomery
   r2 = dx*dx + dy*dy + dz*dz;
   if (r2 <= sim->r2Max)
   {
      it = atm->types[i];
      jt = atm->types[j];

      r = 0.0;
      //eng = 0.0;
      f = 0.0;
      if (sp[it].charged)  //! ещё одна возможная оптимизация, составить массив произведений зарядов [i][j]
        if (sp[jt].charged) //! хотя не факт, что перемещение по двумерному массиву будет быстрее, чем два IF
          f += coul_iter(r2, r, sp[it].charge * sp[jt].charge, sim->alpha, sim->engEwaldReal);

      vdw = field->vdws[it][jt];
      if (vdw != NULL)
        if (r2 <= vdw->r2cut)
          //f += vdw -> feng(r2, vdw, eng);
          f += vdw -> feng_r(r2, r, vdw, sim->engVdW);

      //list saving
      save_neigh(i, j, r, r2, sim);
      /*
      if ((sim->nNbors[i] < sim->maxNbors) && (sim->nNbors[j] < sim->maxNbors))
      {
         Ni = sim->nNbors[i];
         Nj = sim->nNbors[j];

         sim->nbors[i][Ni] = j;
         sim->nbors[j][Nj] = i;
         sim->distances[i][Ni] = r;
         sim->distances[j][Nj] = r;
         sim->nNbors[i]++;
         sim->nNbors[j]++;
      }
      else
        printf("WARNING[111]: the maximal number of neighbors is reached\n");
      */


      // forces update
      atm->fxs[i] += f * dx;
      atm->fxs[j] -= f * dx;
      atm->fys[i] += f * dy;
      atm->fys[j] -= f * dy;
      atm->fzs[i] += f * dz;
      atm->fzs[j] -= f * dz;
      //return eng;
   }
   //else return 0.0; // too far to interact
}

void pair_onlytwoatomic(int i, int j, Atoms *atm, Field *field, Box *bx, Sim *sim)
// calculate energy and determine forces for ion pairs in the case of only_twoatomic optimization
//   return potential energy
{
   double r, r2, f;//, eng;
   double dx = atm->xs[i] - atm->xs[j];
   double dy = atm->ys[i] - atm->ys[j];
   double dz = atm->zs[i] - atm->zs[j];
   int it, jt;
   int Ni, Nj; //! я не понимаю почему, но если это не объявить тут, программа вылетает
   VdW *vdw;
   Bond *bnd;
   Spec* sp = field->species;

   // periodic boundary for rectangular geometry
   rect_periodic(dx, dy, dz, bx);  //! need be replaced for NOT rectangular geomery
   r2 = dx*dx + dy*dy + dz*dz;
   if (r2 <= sim->r2Max)
   {
      it = atm->types[i];
      jt = atm->types[j];

      r = 0.0;
      //eng = 0.0;
      f = 0.0;
      if (sp[it].charged)  //! ещё одна возможная оптимизация, составить массив произведений зарядов [i][j]
        if (sp[jt].charged) //! хотя не факт, что перемещение по двумерному массиву будет быстрее, чем два IF
          f += coul_iter(r2, r, sp[it].charge * sp[jt].charge, sim->alpha, sim->engEwaldReal);

      if (sim->bondAtoms[i] != j) // atoms are not bonded
      {
        vdw = field->vdws[it][jt];
        if (vdw != NULL)
            if (r2 <= vdw->r2cut)
                //f += vdw -> feng(r2, vdw, eng);
                f += vdw -> feng_r(r2, r, vdw, sim->engVdW);
      }
      else // atoms are bonded, specific intramolecular potential
      {
         bnd = &field->bdata[sim->bondTypes[i]];
         if (r2 > bnd->rMax2)
           {
              // break the bond
              sim->bondAtoms[i] = -1;
              sim->bondAtoms[j] = -1;
              if (atm->types[i] == bnd->spec1) // change atom types, according to bond parameters
              {
                atm->types[i] = bnd->spec1br;
                atm->types[j] = bnd->spec2br;
              }
              else
              {
                atm->types[i] = bnd->spec2br;
                atm->types[j] = bnd->spec1br;
              }

              if (sp[atm->types[i]].canBond)
                {
                   if (sim->nCanBondAtm == sim->maxCanBondAtm)
                     printf("WARNING[115] the maximal number of atoms which can form bond is reached!\n");
                   else
                     {
                        sim->canBondAtms[sim->nCanBondAtm] = i;
                        sim->nCanBondAtm++;
                     }
                }
              if (sp[atm->types[j]].canBond)
                {
                   if (sim->nCanBondAtm == sim->maxCanBondAtm)
                     printf("WARNING[115] the maximal number of atoms which can form bond is reached!\n");
                   else
                     {
                        sim->canBondAtms[sim->nCanBondAtm] = j;
                        sim->nCanBondAtm++;
                     }
                }

              sim->nBndBr++;    // inc the number of breaking bonds
              sim->nBonds--;
              sim->degFree++;
              sim->revDegFree = 1.0 / sim->degFree; //! надо бы это вынести за цикл, ну какой смысл это каждую разрушенную связь перевычислять

              // vdw interaction:
              vdw = field->vdws[it][jt];
              if (vdw != NULL)
                  if (r2 <= vdw->r2cut)
                      //f += vdw -> feng(r2, vdw, eng);
                      f += vdw -> feng_r(r2, r, vdw, sim->engVdW);

              // system energy change
              //! надо куда-то сохранять изменение энергии?
              //eng -= bnd->energy;

           }
         else
           f += bond_iter_r(r2, r, bnd, sim->engBond);
      }

      //list saving
      save_neigh(i, j, r, r2, sim);
      /*
      if ((sim->nNbors[i] < sim->maxNbors) && (sim->nNbors[j] < sim->maxNbors))
      {
         Ni = sim->nNbors[i];
         Nj = sim->nNbors[j];

         sim->nbors[i][Ni] = j;
         sim->nbors[j][Nj] = i;
         sim->distances[i][Ni] = r;
         sim->distances[j][Nj] = r;
         sim->nNbors[i]++;
         sim->nNbors[j]++;
      }
      else
        printf("WARNING[111]: the maximal number of neighbors is reached\n");
      */

      // forces update
      atm->fxs[i] += f * dx;
      atm->fxs[j] -= f * dx;
      atm->fys[i] += f * dy;
      atm->fys[j] -= f * dy;
      atm->fzs[i] += f * dz;
      atm->fzs[j] -= f * dz;
      //return eng;
   }
   //else return 0.0; // too far to interact
}

void cell_list_coul_vdw(Atoms *atm, Field *field, Box *bx, Sim *sim)
// function for calculation of vdw and real part of ewald iteraction via cell_list algorithm
//   return potential energy
{
   int iC, jC; // cell indexes
   int i, j; // atom indexes
   //double eng = 0.0; // potential energy
//   int itNumb = 0;
//   int itNumb1;

   for (iC = 0; iC < sim->nHead; iC++)
   {
      i = sim->chead[iC];
      //itNumb1 = 0;
      //printf("iC=%d  i=%d  \n", iC, i, j);

      // loop by all atoms in the cell
      while (i >= 0)
      {
         // atoms in the same cell
         j = sim->clist[i];
         while (j >= 0)
         {
            //eng += pair_coul_vdw_lst(i, j, atm, vdws, sp, bx, sim, btps, barr);
            /*eng += */sim->pair(i, j, atm, field, bx, sim);
            j = sim->clist[j];
         }

         // atoms in neighbors cells
         for (jC = 0; jC < sim->nHeadNeig[iC]; jC++)
         {
            j = sim->chead[sim->lstHNeig[iC][jC]];
            while (j >= 0)
            {
               //if ((i < 0)||(j < 0)||(i >= atm->nAt)||(j >= atm->nAt))
               //  printf("cell list error!\n");
               //printf("iC=%d  i=%d  j=%d\n", iC, i, j);
               //! щас здесь всегда с сохранением соседеней, но иногда можно и без них
               //eng += pair_coul_vdw_lst(i, j, atm, vdws, sp, bx, sim, btps, barr);
               /*eng += */sim->pair(i, j, atm, field, bx, sim);
               j = sim->clist[j];
               //itNumb++;
            }
         }
         //! во всех книгах пишут, что нужна проверка (j < i), чтобы избежать повторного начисления
         //!   но, кажется, это лишено смысла, в соседних ячейках частицы явно разные
         //!   а в своей если считать первый j как clist[i] то вроде бы обход осущестлвяется только по одной паре
         i = sim->clist[i];
      } // end loop by i
      //printf("cell[%d] nPair = %d\n", iC, itNumb1);
   }

   //printf("number of pair examinated=%d  energy=%f\n", itNumb, eng);
   //return eng;
}
// end 'cell_list_coul_vdw' function

void all_pair(Atoms *atm, Field *field, Box *bx, Sim *sim)
// function for calculation of vdw and real part of ewald iteraction
{
   int i, j;
   int N = atm->nAt;

   for (i = 0; i < N-1; i++)
     for (j = i + 1; j < N; j++)
     {
       //eng += pair_coul_vdw(i, j, atm, vdws, sp, bx, sim, btps, barr);
       sim->pair(i, j, atm, field, bx, sim);
     }
}
// end 'all_pair' function

void integrate1(Atoms *atm, Spec *spec, Sim *sim, Box *box, double &chit, double &conint)
// the first part of verlet integrator
//  first variant: X += V*dt + F/m * dt^2/2;  V += F/m * dt/2
//  second variant - наоборот, сначала V, потом X: V += F/m * dt/2 X += V*dt + F/m * dt^2/2;
{
  int i, t;
  // if no variable timestep:
  //static double tSt2 = sim->tSt * sim->tSt / 2;   //! проверить все static правильно ли я использую
  static double tSt = sim->tSt;
  static double tStHalf = tSt * 0.5;
  double rmass;
  //double rMASSxTstHalf;
  double charge;

  // Hoover thermostat applying
  if (sim->nvt)
    nvtscale(atm, sim, chit, conint);

#pragma vector always
  for (i = 0; i < atm->nAt; i++)
  {
    t = atm->types[i];
    rmass = spec[t].rmass;
    charge = spec[t].charge;

    //the first stage of velocity update
    atm->vxs[i] += rmass * atm->fxs[i] * tStHalf;
    atm->vys[i] += rmass * atm->fys[i] * tStHalf;
    atm->vzs[i] += rmass * atm->fzs[i] * tStHalf;

    //update coordinates:
    /*
    atm[i].x += atm[i].vx * tSt + rmass * atm[i].fx * tSt2;
    atm[i].y += atm[i].vy * tSt + rmass * atm[i].fy * tSt2;
    atm[i].z += atm[i].vz * tSt + rmass * atm[i].fz * tSt2;
    */
    atm->xs[i] += atm->vxs[i] * tSt;
    atm->ys[i] += atm->vys[i] * tSt;
    atm->zs[i] += atm->vzs[i] * tSt;

    //box applying
    //! only for rectangular!
    if (!rect_box(atm, i, spec, box))
      {
         break;
      }

    //external field
    //! use Ux/ra as const
    sim->engElecField += charge * atm->xs[i] * sim->Ux * box->ra;
  }
}
// end integrate1() function

void integrate1_lst(Atoms *atm, Spec *spec, Sim *sim, Box *box, double &chit, double &conint)
// the first part of verlet integrator (with cell-list implementation)
//  first variant: X += V*dt + F/m * dt^2/2;  V += F/m * dt/2
//  second variant - наоборот, сначала V, потом X: V += F/m * dt/2 X += V*dt + F/m * dt^2/2;
{
  int i, t, c;
  // if no variable timestep:
  //static double tSt2 = sim->tSt * sim->tSt / 2;   //! проверить все static правильно ли я использую
  static double tSt = sim->tSt;
  static double tStHalf = tSt * 0.5;
  double rmass;
  double charge;

  // Hoover thermostat applying
  if (sim->nvt)
    nvtscale(atm, sim, chit, conint);

  // Clear C_LIST
  for (i = 0; i < sim->nHead; i++)
    sim->chead[i] = -1;

#pragma vector always
  for (i = 0; i < atm->nAt; i++)
  {
    t = atm->types[i];
    rmass = spec[t].rmass;
    charge = spec[t].charge;

    //the first stage of velocity update
    atm->vxs[i] += rmass * atm->fxs[i] * tStHalf;
    atm->vys[i] += rmass * atm->fys[i] * tStHalf;
    atm->vzs[i] += rmass * atm->fzs[i] * tStHalf;

    //update coordinates:
    atm->xs[i] += atm->vxs[i] * tSt;
    atm->ys[i] += atm->vys[i] * tSt;
    atm->zs[i] += atm->vzs[i] * tSt;

    //box applying
    //! only for rectangular!
    rect_box(atm, i, spec, box);

    //save the atom in cell list
    c = cell_index_sim(atm->xs[i], atm->ys[i], atm->zs[i], sim);
    sim->clist[i] = sim->chead[c];
    sim->chead[c] = i;

    //external field
    //! use Ux/ra as const
    sim->engElecField += charge * atm->xs[i] * sim->Ux * box->ra;
  }
}
// end integrate1_lst() function

void integrate2(Atoms *atm, Spec *sp, Sim* sim, int tScale, double &chit, double &conint)
// the second part of verlet integrator
//   with termostat applying
//   return kinetic energy (mv2/2)
{
  double k;
  int i, ind;
  int N = atm->nAt;
  double tSt = sim->tSt;
  // if no variable timestep:
  static double tStHalf = tSt * 0.5;
  double tempA = 0.0;  // temperature
  //double sigma = sqrt(sim->tTemp);
  //double kinE;
  //static double tStNu = tSt * nu;
  double rmass;
  //double rMASSxTstHalf;
/*
  double *vxs = atm->vxs;
  double *vys = atm->vys;
  double *vzs = atm->vzs;
  double *fxs = atm->fxs;
  double *fys = atm->fys;
  double *fzs = atm->fzs;
  int *types = atm->types;
*/

#pragma vector always
  for (i = 0; i < N; i++)
  {
    //the second stage of velocity update
    ind = atm->types[i];

    //if (isnan(atm->vxs[i]))
    //    printf("vxs[%d] is nan\n", i);

    rmass = sp[ind].rmass;
    atm->vxs[i] += rmass * atm->fxs[i] * tStHalf;
    atm->vys[i] += rmass * atm->fys[i] * tStHalf;
    atm->vzs[i] += rmass * atm->fzs[i] * tStHalf;

    //if (isnan(atm->fxs[i]))
    //    {printf("int2: fxs[%d] is nan\n", i);break;}

    /*
    rMASSxTstHalf = sp[ind].rmass * tStHalf;
    atm->vxs[i] += atm->fxs[i] * rMASSxTstHalf;
    atm->vys[i] += atm->fys[i] * rMASSxTstHalf;
    atm->vzs[i] += atm->fzs[i] * rMASSxTstHalf;
    */
    /*
    vxs[i] += fxs[i] * rMASSxTstHalf;
    vys[i] += fys[i] * rMASSxTstHalf;
    vzs[i] += fzs[i] * rMASSxTstHalf;
    */

    //printf("atm[%d].vx=%f\n", i, atm[i].vx);
    //во френкелевском примере не умножалось на массу, а потом делилось на 3/m
    //видимо для случая, когда все частицы - одинаковы
    // я сделал это более логично- сначала домножал, а потом разделил на 3

    tempA += (atm->vxs[i] * atm->vxs[i] + atm->vys[i] * atm->vys[i] + atm->vzs[i] * atm->vzs[i]) * sp[ind].mass;
    //tempA += (vxs[i] * vxs[i] + vys[i] * vys[i] + vzs[i] * vzs[i]) * sp[ind].mass;
  }

  sim->engKin = 0.5 * tempA; // summ(mv^2/2)
  //tempA = tempA * r3kB / N;         //  T = <mv^2>/3/kB

  //T-scale:
  if (tScale)
    {
       k = sqrt(sim->tKin / sim->engKin);
       for (i = 0; i < N; i++)
         {
            atm->vxs[i] *= k;
            atm->vys[i] *= k;
            atm->vzs[i] *= k;
         }
       sim->engKin = sim->tKin;
    }
  //printf("aft T-scale: %d tempA=%f engKin=%f vx[0]=%f\n", tScale, tempA, sim->engKin, atm->vxs[0]);

  //applying Nose-Hoover thermostat:
  //! я не уверен нужно ли его после T-Scale или до. И вообще, как они взаимодействуют
  if (sim->nvt)
    {
       //printf("call nvtscale(atm, sim, %f, %f)\n", chit, conint);
       sim->engKin = nvtscale(atm, sim, chit, conint);
       //! в ДЛПОЛИ есть ещё строчка consv = conint + 0.5 qmass * chit^2
    }

  //applying Andersen termostat
/*
  tempA /= (3 * N); // instaneous tempertaure
  for (i = 0; i < N; i++)
    if (rand01() < tStNu)
    {
      kinE -= (atm[i].vx * atm[i].vx + atm[i].vy * atm[i].vy + atm[i].vz * atm[i].vz) * sp[atm[i].type].mass / 0.5;
      atm[i].vx = gauss(sigma, l0);
      atm[i].vy = gauss(sigma, l0);
      atm[i].vz = gauss(sigma, l0);
      kinE += (atm[i].vx * atm[i].vx + atm[i].vy * atm[i].vy + atm[i].vz * atm[i].vz) * sp[atm[i].type].mass / 0.5;
    }
*/
}
// end 'integrate2()' function

int main(int argc, char *argv[])
{

   int i, j;
   int iSt = 0;  // nuber of steps, index of step
   double tSim = 0.0; // time of simulation
   int mxJump = 0;      // (c) maximal number of eJump procedure per timestep
   int tJump = 0;       // (c) the total number of eJump procedure
   char revfname[32];
   char c;       // for press anykey
   int nRDFout;
   double **rdf;

   // for Nose-Hoover:
   double nu = 0.005;
   double l0 = 0.05;
   double chit = 0.0;
   double conint = 0.0;

   printf("azTotMD by Raskovalov Anton\n");
   int start_time = time(NULL);

   //SYSTEM INITIALISATION
   Sim *sim = (Sim*)malloc(sizeof(Sim));
   Field *field = init_field(sim);
   Box *box;
   Atoms *atoms = init_atoms_and_box(field, sim, box);
   if (!init_md(atoms, field, sim, box))
     {
        printf("FATAL ERROR: SYSTEM CAN'T BE INITIALIZED\n");
        scanf("%c", &c);
        return 0;
     }
   prepare_md(atoms, sim, field);

   init_neighbors(atoms, sim);  //! он может быть и не всегда нужен
   Ewald *ewald = (Ewald*)malloc(sizeof(Ewald));
   init_ewald(ewald, box, sim, atoms);
   init_rdf(&rdf, sim, field, box, nRDFout);
   init_ejump(atoms, field, sim);

   // variable functions:
   void (*integrator1)(Atoms *atm, Spec *spec, Sim *sim, Box *box, double &chit, double &conint);
   void (*forcefield)(Atoms *atm, Field *field, Box *bx, Sim *sim);
   int (*ejumper)(Atoms *atm, Spec *sp, VdW ***vdws, Sim *sim, Box *bx);

   srand(time(NULL));   // ranomize for Metropolis scheme

   if (sim->gaussian)
     gauss_temp(atoms, field->species, sim);

   //OUTPUT FILES
   FILE *jf;        // electron jump file

   FILE *info_file = fopen("info.dat", "w");
   info_header(info_file, atoms, field->species);

   FILE *stat_file = fopen("stat.dat", "w");
   stat_header(stat_file, sim, field->species);

   FILE *hf = fopen("hist.dat", "w");
   history_header(hf);

   FILE *msd_file = fopen("msd.dat", "w");
   msd_header(msd_file, sim, field);

   if (sim->eJump)
   {
     jf = fopen("jumps.dat", "w");
     fprintf(jf, "time\tstep\ttot\tpX\tnX\tp\tn");
     for (i = 0; i < field->nSpec; i++)
       if ((field->species[i].donacc >> 0) & 1)
         for (j = 0; j < field->nSpec; j++)
           if ((field->species[j].donacc >> 1) & 1)
             fprintf(jf, "\t%s->%s", field->species[i].name, field->species[j].name);

     fprintf(jf, "\n");

     info_md(sim);
     if (sim->ejtype == ejt_min)
       ejumper = ejump_min;
     else if (sim->ejtype == ejt_eq)
       ejumper = ejump;
     else
       ejumper = ejump_metr;
   }

   if (sim->eJump || sim->only_twoatomic)
     {
        sim->pair = pair_coul_vdw_lst;
        printf("Used pair_coul_vdw_lst\n");

     }
   else
     sim->pair = pair_coul_vdw;

   sim->engVdW = 0.0;
   sim->engEwaldReal = 0.0;

   //initial forces
   clear_force(atoms, field->species, sim, box);
   sim->engEwaldRec = ewald_rec(atoms, field->species, box, sim, ewald);
   all_pair(atoms, field, box, sim);
   sim->engEwaldConst = ewald_const(atoms, field->species, sim, box);
   printf("initial values:\n");
   printf("x0=%4.3f vx=%4.3f fx=%4.3f VdW+Real=%E Ewald=%E\n", atoms->xs[0], atoms->vxs[0], atoms->fxs[0], sim->engVdW+sim->engEwaldReal, /*potE1*/ sim->engEwaldRec);

   //! Сделать функцию вывода ВдВ пар и вызвать её здесь
/*
   printf("Short range intractions:\n");
   for (i = 0; i < field->nSpec; i++)
     for (j = i; j < field->nSpec; j++)
       if (field->vdws[i][j] != NULL)
         printf("%s %s %s r2cut=%f p1=%f p2=%f p3=%f p4=%f p5=%f\n", field->species[i].name, field->species[j].name, vdw_names[field->vdws[i][j]->type], field->vdws[i][j]->r2cut, field->vdws[i][j]->p0, field->vdws[i][j]->p1, field->vdws[i][j]->p2, vdws[i][j]->p3, vdws[i][j]->p4);
       else
         printf("%s %s NULL\n", species[i].name, species[j].name);
*/


   /*
   printf("Bonded intractions:\n");
   for (i = 0; i < sim->nBtypes; i++)
    printf("Bonds[%d]: %d %d %d %d\n", i+1, bTypes[i].spec1, bTypes[i].spec2, bTypes[i].spec1br, bTypes[i].spec2br);
   */

/*
   printf("Species that can bond:\n");
   for (i = 0; i < sim->nSpec; i++)
    if (species[i].canBond)
      printf("%s[%d,%d,%d]\n", species[i].name, species[i].bondKeys[0], species[i].bondKeys[1], species[i].bondKeys[2]);
*/



   if (sim->useClist)
     {
        integrator1 = integrate1_lst;
        forcefield = cell_list_coul_vdw;

        /*
        if (sim->eJump)
          forcefield = cell_list_coul_vdw;
        else
          forcefield = cell_list_coul_vdw;
        */
     }
   else
     {
        integrator1 = integrate1;
        forcefield = all_pair;
     }


   // MAIN MD LOOP:
   //! разбить на 2 цикла - эквилибровочный и настоящий
   while (iSt < sim->nSt)
   {
     iSt++;

     // reset some quantities
     //sim->engPot = 0.0;
     sim->engVdW = 0.0;
     sim->engEwaldReal = 0.0;
     sim->engElecField = 0.0;
     sim->engBond = 0.0;

     // THE FIRST STAGE OF THE MOITION EQUATION INTEGRATION :  X += V*dt + F/m *dt^2/2;   V += F/m dt/2
     integrator1(atoms, field->species, sim, box, chit, conint);

     //if (iSt % 500 == 0)
     //  printf("dev by Raskovalov A.A. Stat: t=%d x=%f  K=%f\n", iSt, atoms->xs[0], sim->engKin);


     // FORCE CALCULATION :  F
     clear_force(atoms, field->species, sim, box);
     //fill_clist(atoms, nAt, head, list, &clist);  //! если вызываются все пары, то в этом алгоримте нет надобности
     sim->engEwaldRec = ewald_rec(atoms, field->species, box, sim, ewald);
     forcefield(atoms, field, box, sim);
     if (field->nBonds)
        exec_bondlist(atoms, field, sim, box);
     if (field->nAngles)
        exec_anglelist(atoms, field, sim, box);

     // if negative ejump do ejump every -Nth step
     if (sim->eJump < 0)
     {
        if (iSt % (-sim->eJump) == 0)
           ejumper(atoms, field->species, field->vdws, sim, box);
     }
     else if (sim->eJump) // if ejump > 0 do N times ejump
     {
        //if ejump has a result repeat sim->eJump times:
        i = 0;
        for (j = 0; j < sim->eJump; j++)
          {
             if (!ejumper(atoms, field->species, field->vdws, sim, box))
               break;

          }
        //! counter (temp)
        if (j > mxJump)
            mxJump = j;

        tJump += j;
        //if (iSt % 100 == 0) printf("nJ=%d\n", i);
     }

     // THE SECOND STAGE OF THE MOITION EQUATION INTEGRATION :   // V += F/m dt/2
     //printf("before: engKin=%f\n", sim->engKin);
     if (iSt > sim->nEq) // end equilibration, no temperature scaling
       {
          integrate2(atoms, field->species, sim, 0, chit, conint);
          if (iSt % sim->frRDF == 0)
            get_rdf(atoms, rdf, sim, field, box, nRDFout);
       }
     else             // equilibration period, T-Scale
       {
         if ((iSt % sim->freqEq) == 0)
           integrate2(atoms, field->species, sim, 1, chit, conint);
         else
           integrate2(atoms, field->species, sim, 0, chit, conint);

         //update initial coordinates (for MSD calculation)
         if (iSt == sim->nEq)
           for (i = 0; i < atoms->nAt; i++)
             {
                atoms->x0s[i] = atoms->xs[i];
                atoms->y0s[i] = atoms->ys[i];
                atoms->z0s[i] = atoms->zs[i];
             }
       }
     sim->Temp = 2.0 * sim->engKin * sim->revDegFree * rkB;

     // STATISTICS :
     tSim += sim->tSt;
     sim->engTot = sim->engElecField + sim->engVdW + sim->engEwaldReal + sim->engEwaldRec + sim->engEwaldConst + sim->engKin + sim->engBond;

     if (iSt % sim->stat == 0)
     {
         out_stat(stat_file, tSim, iSt, sim, field->species);
         out_msd(msd_file, atoms, atoms->nAt, field->species, field->nSpec, box, tSim, iSt);
         out_info(info_file, tSim, iSt, atoms, field->species);
     }
     // hist file output
     if (iSt % sim->hist == 0)
     {
       fprintf(hf, "%f %d %f %f %f %f %f %f %f %f %f %f %f\n", tSim, iSt, sim->engTot, sim->Temp, atoms->xs[0], atoms->ys[0], field->species[atoms->types[0]].charge, box->momXn, box->momXp, box->momYn, box->momYp, box->momZn, box->momZp);
       if (sim->eJump)
         {
            fprintf(jf, "%f\t%d\t%d\t%d\t%d\t%d\t%d", tSim, iSt, sim->nJump, sim->elJump.px, sim->elJump.nx, sim->pEjump, sim->nEjump);
            for (i = 0; i < field->nSpec; i++)
              if ((field->species[i].donacc >> 0) & 1)
                for (j = 0; j < field->nSpec; j++)
                  if ((field->species[j].donacc >> 1) & 1)
                    fprintf(jf, "\t%d", sim->jumps[i][j]);

            fprintf(jf, "\n");
         }

     } // end history output

     if (sim->revcon)
       if (iSt % sim->revcon == 0)
         {
             sprintf(revfname, "revcon%d.xyz\0", iSt);
             out_atoms(atoms, atoms->nAt, field->species, box, revfname);
         }
   }
   // end main MD loop

   //if (sim->eJump)
     //printf("number of eJump procedure per timestep: max=%d, mean=%f\n", mxJump, (double)(tJump / sim->nSt));


   if (field->nBonds)
     save_bondlist("revbonds.txt", field);
   if (field->nAngles)
     save_anglelist("revangles.txt", field);

   out_rdf(rdf, field, box, sim, "rdf.dat", nRDFout);
   out_atoms(atoms, atoms->nAt, field->species, box, "revcon.xyz");

   if (hf != NULL) fclose(hf);
   if (stat_file != NULL) fclose(stat_file);
   if (msd_file != NULL) fclose(msd_file);
   if (info_file != NULL) fclose(info_file);


   free_rdf(&rdf, field->nPair);
   free_ewald(ewald, atoms);

   if (sim->useClist)
     free_clist(sim);

   if (sim->eJump)
   {
     fclose(jf);
     free_ejump(sim, field);
   }

   if (sim->nVarSpec)
     delete[] sim->varSpecs;

   free_neighbors(atoms, sim);
   //free_md(&atoms, field);
   free_field(field);
   delete field;
   delete ewald;
   delete sim;
   free_atoms(&atoms);

   // Spended time calculation
   int end_time = time(NULL);
   int spend_time = end_time - start_time;
   int hours = spend_time / 3600;
   int minutes = spend_time - 3600 * hours;
   int secs = minutes % 60;
   minutes = minutes / 60;

   printf("The program's just finished correctly, the running time: %d s (%d h, %d min, %d sec)\n", spend_time, hours, minutes, secs);
   scanf("%c", &c);
}
