#include <stdlib.h>  // malloc, alloc, rand, NULL
#include <math.h>   // log, sqrt
#include <stdio.h>   // printf (temp)

#include "dataStruct.h"  // Sim, Box, Atoms ....
#include "utils.h"  // int_size, pointer_size, etc...
#include "const.h"  // sphera consts, kB
#include "box.h"  // rect_periodic
#include "vdw.h"  // vdw_iter
#include "ejump.h"

// move electron [ind] from iat to jat,
void electron_move(int ind, int iat, int jat, int ti1, int ti2, int tj1, int tj2, int px, Spec *sp, Atoms *atm, Sim *sim, Box *bx)
{
    sim->electrons[ind] = jat;        // move electron to j-atom
    // if new electron localization center can't donor electron, we need to decrease free electrons
    if (!((sp[tj2].donacc >> 0) & 1))
    {
        sim->electrons[ind] = sim->electrons[sim->nFreeEl-1];
        //! ������ � ���� ������ ��������� �������� ��������� ���� ���
        sim->nFreeEl--;
        printf("The free electron number is decreased!\n");
    }

    atm->types[iat] = ti2;
    atm->types[jat] = tj2;
    sim->jumps[ti1][tj1]++;

    // change number of particles:
    sp[ti1].number--;
    sp[ti2].number++;
    sp[tj1].number--;
    sp[tj2].number++;

    if (px > 0)
      sim->elJump.px++;
    else if (px < 0)
      sim->elJump.nx++;
    //! �������� ������ ���������

    // counters for middle box
    if (atm->xs[iat] <= bx->ha)
    {
       if (atm->xs[jat] > bx->ha)
          if (atm->xs[iat] > (bx->ha - sim->rElec))
              sim->pEjump++;
    }
    else // atm->xs[iat] > box->ha
    {
         if (atm->xs[jat] <= bx->ha)
            if (atm->xs[iat] <= (bx->ha + sim->rElec))
                sim->nEjump++;
    }
}


// initialize arrays and vars for e-jump
void init_ejump(Atoms *atm, Field *field, Sim *sim)
//int **Ns, int **tProbs, int ***acs, double ***dists, int ***probs, , ThrowBox **eThrow, int ***jmps
// Ns - ��� ��� nPairs

 /*  int *nPairs;  // number of neightbours for jumping near atom [i]
   double **dist; // distances
   int **acceptors; // neightbours
   int **probs;  // probability of jump
   int *totProbs; // total probability
   int nJump;
   ThrowBox *elJump;
   int **jumps; // array for keeping nJump between different donor-acceptor pairs
*/

{
   //double **ad;
   //int **ai;
   //int *a, *a1;
   int i, j, k;
   Spec *sp = field->species;

   //sim->nPairs = (int*)malloc(sim->nAt * int_size);
   //sim->totProbs = (int*)malloc(sim->nAt * int_size);
   //sim->acceptors = (int**)malloc(sim->nAt * pointer_size);
   //sim->dist = (double**)malloc(sim->nAt * pointer_size);
   //sim->probs = (int**)malloc(sim->nAt * pointer_size);
   sim->jumps = (int**)malloc(field->nSpec * pointer_size);
   sim->elJump.nx = 0;
   sim->elJump.px = 0;
   sim->nJump = 0;
   //sim->nJumpVar = 10; //! temp
   //sim->EngDifs = (double*)malloc(sim->nJumpVar * double_size);


   // define electrons array and their positions
   k = 0;
   sim->electrons = (int*)malloc(sim->nFreeEl * int_size);
   for (i = 0; i < atm->nAt; i++)
     for (j = 0; j < sp[atm->types[i]].nFreeEl; j++)
       {
          sim->electrons[k] = i;
          k++;
       }

   //a1 = *tProbs;
   /*
   for (i = 0; i < sim->nAt; i++)
     {
        sim->nPairs[i] = 0;
        //a1[i] = 0;

        sim->acceptors[i] = (int*)malloc(sim->maxAcc * int_size);
        sim->dist[i] = (double*)malloc(sim->maxAcc * double_size);
     }
   */

   for (i = 0; i < field->nSpec; i++)
     {
       sim->jumps[i] = (int*)malloc(field->nSpec * int_size);
       for (j = 0; j < field->nSpec; j++)
        sim->jumps[i][j] = 0;
     }

}
// end 'init_ejump' function


int ejump(Atoms *atm, Spec *sp, VdW ***vdws, Sim *sim, Box *bx)
// do electron jump (return the number of electron jumps)
//    nJump is a counter for total number of Jump
//    px, py - counters for ejumps across box border
{
   int i, j, k/*, r, sum*/;
   int iat, jat, kat, ktype;
   int tai1, taj1, tai2, taj2; //type of atom i/j  before jump(1) after jump(2)
   double U1, U2, dU; // energies before(1) and after jump(2), delta U
   VdW *vdw;
   double dx, dy, dz, r2, rad;
   int px, py, pz;
   //double xi, xj;
   //double p;
   double kcharge;

   //int nAcc = 0;    // number of acceptor
   //int nNoRad = 0;


   int result = 0;  // the number of electron jumps that we need to return
   for (i = 0; i < sim->nFreeEl; i++)
   {
      iat = sim->electrons[i]; //! �������� ��� ���� ��������, ��� ��� �� �������� �� ���� �����? ���� ��� ����� ��������� ����? ����, ��������, �� ��������� � ���� bond breaking ��� formation
      tai1 = atm->types[iat];
      tai2 = sp[tai1].oxForm - 1;

      for (j = 0; j < sim->nNbors[iat]; j++)
        if ((sim->tnbors[iat][j] >> 2) & 1) //!temp  // verify neighbors distance type: 0 - Coul, 1 - vdw, 2 - eJump ...
        {
           jat = sim->nbors[iat][j];
           taj1 = atm->types[jat];
           if ((sp[taj1].donacc >> 1) & 1)  // verify, can the neighbor be an acceptor
           {
              taj2 = sp[taj1].redForm - 1;

              //ENERGY CALCULATION
              U1 = 0.0; U2 = 0.0;

              // loop by iat neighbors:
              //! ��������, ����� ���������� � ���� ���� �� ������ iat � �� ������� jat? � ����� ������ �������
              for (k = 0; k < sim->nNbors[iat]; k++)
              {
                 //nAcc++;
                 kat = sim->nbors[iat][k];
                 ktype = atm->types[kat];
                 kcharge = sp[ktype].charge;

                 //! ��������� ejump ���������� ����� forcefield ���� ������ ������������� �������� �� ������
                 //!  ��� ��� �� ���? ���� ��� ������� �������� � ���������� ������ ������ ��������� �� ������, �� ���
                 rad = sim->distances[iat][k];
                 if (rad < 0.3)
                 {
                    //printf("ejump: rad(%f) < 0.3\n", rad);
                    //nNoRad++;
                 }
                 r2 = rad * rad;

                 // VdW energy before jump:
                 vdw = vdws[tai1][ktype];
                 if (vdw != NULL)
                   if (r2 <= vdw->r2cut)  //! ����������, ����� �������� ��� ��������, ���� �� �� ����� � ������, �� ������ ������������ �� ���
                     //vdw_iter(r2, vdw, U1);
                     U1 += vdw->eng_r(r2, rad, vdw);

                 // Columbic energy before jump:
                 //! �������� �� �������?
                 U1 += Fcoul_scale * kcharge * sp[tai1].charge / rad;

                 //vdw energy after jump:
                 vdw = vdws[tai2][ktype];
                   if (vdw != NULL)
                     if (r2 <= vdw->r2cut)
                       //vdw_iter(r2, vdw, U2);
                       U2 += vdw->eng_r(r2, rad, vdw);
                 // Columbic energy after jump
                 //! �������� �� �������?
                 U2 += Fcoul_scale * kcharge * sp[tai2].charge / rad;
              } // end loop by iat neigbors

              // loop by jat neighbors:
              for (k = 0; k < sim->nNbors[jat]; k++)
              {
                 //nAcc++;
                 kat = sim->nbors[jat][k];
                 if (kat == iat)  // this interaction was calculated in previous cycle, skip
                   continue;

                 ktype = atm->types[kat];
                 kcharge = sp[ktype].charge;

                 //! ��������� ejump ���������� ����� forcefield ���� ������ ������������� �������� �� ������
                 //!  ��� ��� �� ���? ���� ��� ������� �������� � ���������� ������ ������ ��������� �� ������, �� ���
                 rad = sim->distances[jat][k];
                 if (rad < 1e-6)
                 {
                    //! ��������� ������ ����������. ���� ����� �� ������ ���������� ����� ��������
                 }
                 r2 = rad * rad;

                 // VdW energy before jump:
                 vdw = vdws[taj1][ktype];
                 if (vdw != NULL)
                   if (r2 <= vdw->r2cut)  //! ����������, ����� �������� ��� ��������, ���� �� �� ����� � ������, �� ������ ������������ �� ���
                     //vdw_iter(r2, vdw, U1);
                     U1 += vdw->eng_r(r2, rad, vdw);

                 // Columbic energy before jump:
                 //! �������� �� �������?
                 U1 += Fcoul_scale * kcharge * sp[taj1].charge / rad;

                 //vdw energy after jump:
                 vdw = vdws[taj2][ktype];
                   if (vdw != NULL)
                     if (r2 <= vdw->r2cut)
                       //vdw_iter(r2, vdw, U2);
                       U2 += vdw->eng_r(r2, rad, vdw);
                 // Columbic energy after jump
                 //! �������� �� �������?
                 U2 += Fcoul_scale * kcharge * sp[taj2].charge / rad;
              }  // end loop by jat neigbors

              dU = U2 - U1;
              //if (isnan(dU))
              //  printf("eq_min: dU is NAN\n");
              dU += (sp[tai2].energy + sp[taj2].energy - sp[tai1].energy - sp[taj1].energy);

              dx = atm->xs[jat] - atm->xs[iat];
              dy = atm->ys[jat] - atm->ys[iat];
              dz = atm->zs[jat] - atm->zs[iat];
              pass_rect_box(dx, dy, dz, bx, px, py, pz);
              //r2 = dx*dx + dy*dy + dz*dz;


              //external electric field addition:
              //don't forget about periodic conditions!
              //! ����� ���� ��� � �� ���������... ��. ������� ��� ���� �����

              //! ��. ����������� � ������� how to change energy via jump away box
              dU += sim->Ux * bx->ra * (atm->xs[iat] * (sp[tai2].charge - sp[tai1].charge) + (atm->xs[jat] + px * bx->la) * (sp[taj2].charge - sp[taj1].charge));


              if (dU < sim->dEjump/*dU < 0*/) //! ��������, ��� �� ��� ���������, �� ������ ��� �� ������ ������� ������������ ���0-��
                if (dU > -sim->dEjump)  // energy equality (the Frank-Condon principle), doing e-jump
                  {
                     result++;
                     electron_move(i, iat, jat, tai1, tai2, taj1, taj2, px, sp, atm, sim, bx);
                     break; // very small chance of the energy equality with several neighbors
                  }


           }  // end if jat is acceptor
        }  // end loop by neighbors of electron
   }  // end loop by electrons

   sim->nJump += result;
   //printf("%d iterations (without distance: %d)\n", nAcc, nNoRad);
   return result;
}
// end 'ejump' function

int ejump_min(Atoms *atm, Spec *sp, VdW ***vdws, Sim *sim, Box *bx)
// do electron jump (return the number of electron jumps), chose jump with minimal energy
//    nJump is a counter for total number of Jump
//    px, py - counters for ejumps across box border
{
   int i, j, k/*, r, sum*/;
   int iat, jat, kat, ktype;
   int tai1, taj1, tai2, taj2; //type of atom i/j  before jump(1) after jump(2)
   double U1, U2, dU; // energies before(1) and after jump(2), delta U
   VdW *vdw;
   double dx, dy, dz, r2, rad;
   int px, py, pz;
   //double xi, xj;
   //double p;
   double kcharge;
   int indMin;
   int pxMin;
   double minE;

   //vars for verification: (delete in final version)
   //int nJat = 0;    // number of neighbors
   //int nAcc = 0;    // number of acceptor
   //int nNeg = 0;    // number of negative dU
   //double dUmin, dUmax, dUmean; // energy diferences...
   //int first;   // flag first dU saving

   int result = 0;  // the number of electron jumps that we need to return
   for (i = 0; i < sim->nFreeEl; i++)
   {
      iat = sim->electrons[i];
      tai1 = atm->types[iat];
      tai2 = sp[tai1].oxForm - 1;

      //nVar = 0;
      minE = 0.0;
      indMin = 0;
      //first = 1;
      for (j = 0; j < sim->nNbors[iat]; j++)
        if ((sim->tnbors[iat][j] >> 2) & 1) //!temp  // verify neighbors distance type: 0 - Coul, 1 - vdw, 2 - eJump ...
        {
           //nJat++;
           jat = sim->nbors[iat][j];
           taj1 = atm->types[jat];
           if ((sp[taj1].donacc >> 1) & 1)  // verify, can the neighbor be an acceptor
           {
              taj2 = sp[taj1].redForm - 1;

              //ENERGY CALCULATION
              U1 = 0.0; U2 = 0.0;
              dU = 0.0;

              // loop by iat neighbors:
              //! ��������, ����� ���������� � ���� ���� �� ������ iat � �� ������� jat? � ����� ������ �������
              for (k = 0; k < sim->nNbors[iat]; k++)
              {
                 //nAcc++;
                 kat = sim->nbors[iat][k];
                 ktype = atm->types[kat];
                 kcharge = sp[ktype].charge;

                 //! ��������� ejump ���������� ����� forcefield ���� ������ ������������� �������� �� ������
                 //!  ��� ��� �� ���? ���� ��� ������� �������� � ���������� ������ ������ ��������� �� ������, �� ���
                 rad = sim->distances[iat][k];
                 r2 = rad * rad;
                 //if (rad < 0.3)
                   // printf("ejump: rad<0.3\n");

                 // VdW energy before jump:
                 vdw = vdws[tai1][ktype];
                 if (vdw != NULL)
                   if (r2 <= vdw->r2cut)  //! ����������, ����� �������� ��� ��������, ���� �� �� ����� � ������, �� ������ ������������ �� ���
                     //vdw_iter(r2, vdw, U1);
                     U1 += vdw->eng_r(r2, rad, vdw);

                 //vdw energy after jump:
                 vdw = vdws[tai2][ktype];
                   if (vdw != NULL)
                     if (r2 <= vdw->r2cut)
                       //vdw_iter(r2, vdw, U2);
                       U2 += vdw->eng_r(r2, rad, vdw);

                 // Columbic energy difference:
                 //! �������� �� �������?
                 dU += Fcoul_scale * kcharge * (sp[tai2].charge - sp[tai1].charge) / rad;
              } // end loop by iat neigbors

              // loop by jat neighbors:
              for (k = 0; k < sim->nNbors[jat]; k++)
              {
                 //nAcc++;
                 kat = sim->nbors[jat][k];
                 if (kat == iat)  // this interaction was calculated in previous cycle, skip
                   continue;

                 ktype = atm->types[kat];
                 kcharge = sp[ktype].charge;

                 //! ��������� ejump ���������� ����� forcefield ���� ������ ������������� �������� �� ������
                 //!  ��� ��� �� ���? ���� ��� ������� �������� � ���������� ������ ������ ��������� �� ������, �� ���
                 rad = sim->distances[jat][k];
                 r2 = rad * rad;
                 //if (rad < 0.3)
                   // printf("ejump: rad<0.3\n");

                 // VdW energy before jump:
                 vdw = vdws[taj1][ktype];
                 if (vdw != NULL)
                   if (r2 <= vdw->r2cut)  //! ����������, ����� �������� ��� ��������, ���� �� �� ����� � ������, �� ������ ������������ �� ���
                     //vdw_iter(r2, vdw, U1);
                     U1 += vdw->eng_r(r2, rad, vdw);

                 //vdw energy after jump:
                 vdw = vdws[taj2][ktype];
                   if (vdw != NULL)
                     if (r2 <= vdw->r2cut)
                       //vdw_iter(r2, vdw, U2);
                       U2 += vdw->eng_r(r2, rad, vdw);

                 // Columbic energy difference
                 //! �������� �� �������?
                 dU += Fcoul_scale * kcharge * (sp[taj2].charge - sp[taj1].charge) / rad;
              }  // end loop by jat neigbors

              dU = dU + U2 - U1;
              // own energy
              dU += (sp[tai2].energy + sp[taj2].energy - sp[tai1].energy - sp[taj1].energy);

              dx = atm->xs[jat] - atm->xs[iat];
              dy = atm->ys[jat] - atm->ys[iat];
              dz = atm->zs[jat] - atm->zs[iat];
              pass_rect_box(dx, dy, dz, bx, px, py, pz);
              //r2 = dx*dx + dy*dy + dz*dz;


              //external electric field addition:
              //don't forget about periodic conditions!
              //! ����� ���� ��� � �� ���������... ��. ������� ��� ���� �����
              //dU += sim->Ux * bx->ra * (atm->xs[iat] * (sp[tai2].charge - sp[tai1].charge) + (atm->xs[jat] + px * bx->la) * (sp[taj2].charge - sp[taj1].charge));

              //! � ������ ��� ������, ��� �������� ����� � ��������� ������ ����� +/- �������
              dU += sim->Ux * bx->ra * (atm->xs[iat] - atm->xs[jat] - px * bx->la);

              /*
              // the first dU calculation for this electron (iat)
              if (first)
              {
                 dUmin = dU;
                 dUmax = dU;
                 dUmean = dU;
                 first = 0;
              }
              else
              {
                 if (dU < dUmin)
                    dUmin = dU;
                 else if (dU > dUmax)
                    dUmax = dU;
                 dUmean += dU;
              }

              if (dU < 0.0)
                nNeg++;
              */

              if (dU < minE)
              {
                 // save variant
                 minE = dU;
                 indMin = jat + 1; // 0 - empty
                 pxMin = px;

                 //! �������� ��� ����� ��������� break,
                 //!  ��������� �������� ����������, ��� � ����������� �������, ������������� ��������� ������ � ����� ������ �� ����
              }


           }  // end if jat is acceptor
        }
      // end loop by neighbors of electron

      // select variant with minimal energy
      if (indMin)
      {
           result++;

           jat = indMin - 1;
           taj1 = atm->types[jat];
           taj2 = sp[taj1].redForm - 1;
           electron_move(i, iat, jat, tai1, tai2, taj1, taj2, px, sp, atm, sim, bx);
      }

   }  // end loop by electrons

   sim->nJump += result;
   //printf("nJ=%d\n", result);
   //printf("%d iterations\n", nAcc);
   return result;
}
// end 'ejump_min' function

//(�� ������ ����� ������� ejump_min)
int ejump_metr(Atoms *atm, Spec *sp, VdW ***vdws, Sim *sim, Box *bx)
// do electron jump (return the number of electron jumps), chose jump with metropolis scheme
//    nJump is a counter for total number of Jump
//    px, py - counters for ejumps across box border
{
   int i, j, k/*, r, sum*/;
   int iat, jat, kat, ktype;
   int tai1, taj1, tai2, taj2; //type of atom i/j  before jump(1) after jump(2)
   double U1, U2, dU; // energies before(1) and after jump(2), delta U
   VdW *vdw;
   double dx, dy, dz, r2, rad;
   int px, py, pz;
   //double xi, xj;
   //double p;
   double kcharge;
   int indMin;
   int pxMin;
   double minE;
   double r;
   int succ; // flag of success move

   //vars for verification: (delete in final version)
   //int nJat = 0;    // number of neighbors
   //int nAcc = 0;    // number of acceptor
   //int nNeg = 0;    // number of negative dU
   //double dUmin, dUmax, dUmean; // energy diferences...
   //int first;   // flag first dU saving

   int result = 0;  // the number of electron jumps that we need to return
   for (i = 0; i < sim->nFreeEl; i++)
   {
      iat = sim->electrons[i];
      tai1 = atm->types[iat];
      tai2 = sp[tai1].oxForm - 1;

      //nVar = 0;
      minE = 0.0;
      indMin = 0;
      //first = 1;
      for (j = 0; j < sim->nNbors[iat]; j++)
        if ((sim->tnbors[iat][j] >> 2) & 1) //!temp  // verify neighbors distance type: 0 - Coul, 1 - vdw, 2 - eJump ...
        {
           //nJat++;
           jat = sim->nbors[iat][j];
           taj1 = atm->types[jat];
           if ((sp[taj1].donacc >> 1) & 1)  // verify, can the neighbor be an acceptor
           {
              //nAcc++;
              taj2 = sp[taj1].redForm - 1;

              //ENERGY CALCULATION
              U1 = 0.0; U2 = 0.0;
              dU = 0.0;

              // loop by iat neighbors:
              //! ��������, ����� ���������� � ���� ���� �� ������ iat � �� ������� jat? � ����� ������ �������
              for (k = 0; k < sim->nNbors[iat]; k++)
              {
                 kat = sim->nbors[iat][k];
                 ktype = atm->types[kat];
                 kcharge = sp[ktype].charge;

                 //! ��������� ejump ���������� ����� forcefield ���� ������ ������������� �������� �� ������
                 //!  ��� ��� �� ���? ���� ��� ������� �������� � ���������� ������ ������ ��������� �� ������, �� ���
                 rad = sim->distances[iat][k];
                 r2 = rad * rad;
                 //if (rad < 0.3)
                   // printf("ejump: rad<0.3\n");

                 // VdW energy before jump:
                 vdw = vdws[tai1][ktype];
                 if (vdw != NULL)
                   if (r2 <= vdw->r2cut)  //! ����������, ����� �������� ��� ��������, ���� �� �� ����� � ������, �� ������ ������������ �� ���
                     //vdw_iter(r2, vdw, U1);
                     U1 += vdw->eng_r(r2, rad, vdw);

                 //vdw energy after jump:
                 vdw = vdws[tai2][ktype];
                   if (vdw != NULL)
                     if (r2 <= vdw->r2cut)
                       //vdw_iter(r2, vdw, U2);
                       U2 += vdw->eng_r(r2, rad, vdw);

                 // Columbic energy difference:
                 //! �������� �� �������?
                 dU += Fcoul_scale * kcharge * (sp[tai2].charge - sp[tai1].charge) / rad;
              } // end loop by iat neigbors

              // loop by jat neighbors:
              for (k = 0; k < sim->nNbors[jat]; k++)
              {
                 kat = sim->nbors[jat][k];
                 if (kat == iat)  // this interaction was calculated in previous cycle, skip
                   continue;

                 ktype = atm->types[kat];
                 kcharge = sp[ktype].charge;

                 //! ��������� ejump ���������� ����� forcefield ���� ������ ������������� �������� �� ������
                 //!  ��� ��� �� ���? ���� ��� ������� �������� � ���������� ������ ������ ��������� �� ������, �� ���
                 rad = sim->distances[jat][k];
                 //if (rad < 0.3)
                   // printf("ejump: rad<0.3\n");
                 r2 = rad * rad;

                 // VdW energy before jump:
                 vdw = vdws[taj1][ktype];
                 if (vdw != NULL)
                   if (r2 <= vdw->r2cut)  //! ����������, ����� �������� ��� ��������, ���� �� �� ����� � ������, �� ������ ������������ �� ���
                     //vdw_iter(r2, vdw, U1);
                     U1 += vdw->eng_r(r2, rad, vdw);

                 //vdw energy after jump:
                 vdw = vdws[taj2][ktype];
                   if (vdw != NULL)
                     if (r2 <= vdw->r2cut)
                       //vdw_iter(r2, vdw, U2);
                       U2 += vdw->eng_r(r2, rad, vdw);

                 // Columbic energy difference
                 //! �������� �� �������?
                 dU += Fcoul_scale * kcharge * (sp[taj2].charge - sp[taj1].charge) / rad;
              }  // end loop by jat neigbors

              dU = dU + U2 - U1;
              //if (isnan(dU))
                //printf("dU is NAN\n");
              // own energy
              dU += (sp[tai2].energy + sp[taj2].energy - sp[tai1].energy - sp[taj1].energy);

              dx = atm->xs[jat] - atm->xs[iat];
              dy = atm->ys[jat] - atm->ys[iat];
              dz = atm->zs[jat] - atm->zs[iat];
              pass_rect_box(dx, dy, dz, bx, px, py, pz);
              //r2 = dx*dx + dy*dy + dz*dz;


              //external electric field addition:
              //don't forget about periodic conditions!
              //! ����� ���� ��� � �� ���������... ��. ������� ��� ���� �����
              //dU += sim->Ux * bx->ra * (atm->xs[iat] * (sp[tai2].charge - sp[tai1].charge) + (atm->xs[jat] + px * bx->la) * (sp[taj2].charge - sp[taj1].charge));

              //! � ������ ��� ������, ��� �������� ����� � ��������� ������ ����� +/- �������
              dU += sim->Ux * bx->ra * (atm->xs[iat] - atm->xs[jat] - px * bx->la);

              /*
              // the first dU calculation for this electron (iat)
              if (first)
              {
                 dUmin = dU;
                 dUmax = dU;
                 dUmean = dU;
                 first = 0;
              }
              else
              {
                 if (dU < dUmin)
                    dUmin = dU;
                 else if (dU > dUmax)
                    dUmax = dU;
                 dUmean += dU;
              }

              if (dU < 0.0)
                nNeg++;
              */
              if (indMin)
              {
                if (dU < minE) // rewrite old variant
                {
                   // save variant
                   indMin = jat + 1; // 0 - empty
                   minE = dU;
                   pxMin = px;

                   //! �������� ��� ����� ��������� break,
                   //!  ��������� �������� ����������, ��� � ����������� �������, ������������� ��������� ������ � ����� ������ �� ����
                }
              }
              else // the first variant, save it in any case
              {
                  indMin = jat + 1;
                  minE = dU;
                  pxMin = px;
              }


           }  // end if jat is acceptor
        }
      // end loop by neighbors of electron

      // accept electron move according to Metropolis scheme:
      succ = 0;
      if (indMin)   // only if there are acceptors near electrons[i]
      {
        if (minE < 0.0)
          succ = 1;
        else
        {
           //printf("minE = %f  metr= %f\n", minE, exp(-rkB*minE/sim->tTemp));
           r = rand01();
           if (r < exp(-rkB*minE/sim->tTemp)) //! ����� ����� �������� �����������
             succ = 1;
        }
      }

      if (succ)
      {
           result++;

           jat = indMin - 1;
           taj1 = atm->types[jat];
           taj2 = sp[taj1].redForm - 1;
           electron_move(i, iat, jat, tai1, tai2, taj1, taj2, px, sp, atm, sim, bx);
      }

   }  // end loop by electrons

   sim->nJump += result;
   return result;
}
// end 'ejump_metr' function

int free_ejump(Sim *sim, Field *field)
// free arrays for e-jump
{

   int i;

/*
   for (i = 0; i < sim->nAt; i++)
     {
        delete[] sim->acceptors[i];
        delete[] sim->dist[i];
     }
   //delete[] sim->nPairs;
   //sim->totProbs = (int*)malloc(sim->nAt * int_size);
   delete[] sim->acceptors;
   delete[] sim->dist;
   //sim->probs = (int**)malloc(sim->nAt * pointer_size);
   delete[] sim->jumps;
*/

   for (i = 0; i < field->nSpec; i++)
     {
       delete[] sim->jumps[i];
     }

   delete[] sim->electrons;
   //delete[] sim->EngDifs;

   return 1;
}
// end 'free_ejump' function
