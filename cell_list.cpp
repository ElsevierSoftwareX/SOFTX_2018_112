//#include "box.h"  // box
#include <stdlib.h>  // malloc, alloc, rand, NULL
#include <math.h>  // floor

//! temp
#include <stdio.h>   // printf

#include "utils.h"  // int_size, pointer_size, etc...
#include "dataStruct.h"  // Sim, Box, Atoms ....
#include "cell_list.h"


int cell_index_sim(double x, double y, double z, Sim *sim)
// return the number of a cell by real{x;y;z} coordintates and sim record
{
   int i = x / sim->clX;
   int j = y / sim->clY;
   int k = z / sim->clZ;

   return i * sim->cnYZ + j * sim->cnZ + k;
}

int cell_index(int x, int y, int z, int mx, int my, int mz)
// return the number of a cell by integer{x;y;z} coordintates
//  this function is for internal usage only
{
   int i = x;
   int j = y;
   int k = z;
   if (i >= mx)
     i = 0;
   else
     if (i < 0)
       i = mx - 1;
   if (j >= my)
     j = 0;
   else
     if (j < 0)
       j = my - 1;
   if (k >= mz)
     k = 0;
   else
     if (k < 0)
       k = mz - 1;

   return i * my * mz + j * mz + k;
}

int init_clist(Atoms* atm, Sim *sim, Box *bx, double rCut)
// create cell-list arrays (clist and hlist)
//   rCut - maximal cutoff radius
{
    int i, j, k, ci;
    int dx, dy, dz;
    int *fd, *sd, *td; // first, second and third dimensions (pointers)
    int *d;         // dimension (pointer)
    int nX, nY, nZ, nYZ;
    int nN;     // the number of neighbors
    //int nZneig;

    //define size of the cell
    nX = floor(bx->la / rCut);
    nY = floor(bx->lb / rCut);
    nZ = floor(bx->lc / rCut);
    if (nX == 0) nX = 1;
    if (nY == 0) nY = 1;
    if (nZ == 0) nZ = 1;

    sim->cnX = nX;
    sim->cnY = nY;
    sim->cnZ = nZ;

    if ((nX < 4)&&(nY < 4)&&(nZ < 4))  // cell list 3x3x3 and smaller has no meaning, use all_pair instead
      return 0;

    nYZ = nY * nZ;
    sim->cnYZ = nYZ;


    sim->clX = bx->la / sim->cnX;
    sim->clY = bx->lb / sim->cnY;
    sim->clZ = bx->lc / sim->cnZ;

    // create list and head arrays:
    sim->clist = (int*)malloc(atm->nAt * int_size);
    sim->nHead = sim->cnX * sim->cnY * sim->cnZ;
    sim->chead = (int*)malloc(sim->nHead * int_size);

    // create array of neighbors for all cells?:
    sim->nHeadNeig = (int*)malloc(sim->nHead * int_size);
    sim->lstHNeig = (int**)malloc(sim->nHead * pointer_size);

/*
    if (nX > 3)
      i = 3;
    else
      i = nX;

    if (nY > 3)
      j = 3;
    else
      j = nY;
    nZneig = i * j;  // the number of the upper neigbors
*/


    if ((nX > 2)&&(nY > 2)&&(nZ > 2))   // all dimensions are full, every cell has 13 neighbours
      for (i = 0; i < nX; i++)
        for (j = 0; j < nY; j++)
          for (k = 0; k < nZ; k++)
            {
               ci = i * nYZ + j * nZ + k;   // cell index
               sim->nHeadNeig[ci] = 13;
               sim->lstHNeig[ci] = (int*)malloc(13 * int_size);

               // upper neighbors
               sim->lstHNeig[ci][0] = cell_index(i+1, j+1, k+1, nX, nY, nZ);
               sim->lstHNeig[ci][1] = cell_index(i+1, j, k+1, nX, nY, nZ);
               sim->lstHNeig[ci][2] = cell_index(i+1, j-1, k+1, nX, nY, nZ);
               sim->lstHNeig[ci][3] = cell_index(i, j+1, k+1, nX, nY, nZ);
               sim->lstHNeig[ci][4] = cell_index(i, j, k+1, nX, nY, nZ);
               sim->lstHNeig[ci][5] = cell_index(i, j-1, k+1, nX, nY, nZ);
               sim->lstHNeig[ci][6] = cell_index(i-1, j+1, k+1, nX, nY, nZ);
               sim->lstHNeig[ci][7] = cell_index(i-1, j, k+1, nX, nY, nZ);
               sim->lstHNeig[ci][8] = cell_index(i-1, j-1, k+1, nX, nY, nZ);

               // neighbors in the same plane
               sim->lstHNeig[ci][9] = cell_index(i+1, j+1, k, nX, nY, nZ);
               sim->lstHNeig[ci][10] = cell_index(i+1, j, k, nX, nY, nZ);
               sim->lstHNeig[ci][11] = cell_index(i+1, j-1, k, nX, nY, nZ);
               sim->lstHNeig[ci][12] = cell_index(i, j+1, k, nX, nY, nZ);

            }


    // one direction has only one cell, other are full (slice geometry):
    if (((nX == 1)&&(nY > 2)&&(nZ > 2)) || ((nX > 2)&&(nY == 1)&&(nZ > 2)) || ((nX > 2)&&(nY > 2)&&(nZ == 1)))
    {
      if (nX == 1)
        {
           dx = 0;
           fd = &dy;
           sd = &dz;
        }
      if (nY == 1)
        {
           dy = 0;
           fd = &dx;
           sd = &dz;
        }
      if (nZ == 1)
        {
           dz = 0;
           fd = &dx;
           sd = &dy;
        }
      for (i = 0; i < nX; i++)
        for (j = 0; j < nY; j++)
          for (k = 0; k < nZ; k++)
            {
               ci = i * nYZ + j * nZ + k;   // cell index
               sim->nHeadNeig[ci] = 4;
               sim->lstHNeig[ci] = (int*)malloc(4 * int_size);

               *fd = 1;
               *sd = 1;
               sim->lstHNeig[ci][0] = cell_index(i+dx, j+dy, k+dz, nX, nY, nZ);
               *sd = 0;
               sim->lstHNeig[ci][1] = cell_index(i+dx, j+dy, k+dz, nX, nY, nZ);
               *sd = -1;
               sim->lstHNeig[ci][2] = cell_index(i+dx, j+dy, k+dz, nX, nY, nZ);
               *fd = 0; *sd = 1;
               sim->lstHNeig[ci][3] = cell_index(i+dx, j+dy, k+dz, nX, nY, nZ);
            }

    } // end 'slice' geometry


    // one direction has only two cells, other are full (stack geometry):
    if (((nX == 2)&&(nY > 2)&&(nZ > 2)) || ((nX > 2)&&(nY == 2)&&(nZ > 2)) || ((nX > 2)&&(nY > 2)&&(nZ == 2)))
    {
      if (nX == 2)
        {
           dx = 0;
           fd = &dy;
           sd = &dz;
           td = &dx;
           d = &i;
        }
      if (nY == 2)
        {
           dy = 0;
           fd = &dx;
           sd = &dz;
           td = &dy;
           d = &j;
        }
      if (nZ == 2)
        {
           dz = 0;
           fd = &dx;
           sd = &dy;
           td = &dz;
           d = &k;
        }
      for (i = 0; i < nX; i++)
        for (j = 0; j < nY; j++)
          for (k = 0; k < nZ; k++)
            {
               ci = i * nYZ + j * nZ + k;   // cell index
               if (*d == 1)
                 nN = 4;
               else
                 nN = 13;

               sim->nHeadNeig[ci] = nN;
               sim->lstHNeig[ci] = (int*)malloc(nN * int_size);

               *td = 0;
               *fd = 1;
               *sd = 1;
               sim->lstHNeig[ci][0] = cell_index(i+dx, j+dy, k+dz, nX, nY, nZ);
               *sd = 0;
               sim->lstHNeig[ci][1] = cell_index(i+dx, j+dy, k+dz, nX, nY, nZ);
               *sd = -1;
               sim->lstHNeig[ci][2] = cell_index(i+dx, j+dy, k+dz, nX, nY, nZ);
               *fd = 0; *sd = 1;
               sim->lstHNeig[ci][3] = cell_index(i+dx, j+dy, k+dz, nX, nY, nZ);

               if (*d == 0) // in this case the next layer is in neighbors
               {
                  *td = 1;

                  *fd = 1;
                  *sd = 1;
                  sim->lstHNeig[ci][4] = cell_index(i+dx, j+dy, k+dz, nX, nY, nZ);
                  *sd = 0;
                  sim->lstHNeig[ci][5] = cell_index(i+dx, j+dy, k+dz, nX, nY, nZ);
                  *sd = -1;
                  sim->lstHNeig[ci][6] = cell_index(i+dx, j+dy, k+dz, nX, nY, nZ);

                  *fd = 0;
                  *sd = 1;
                  sim->lstHNeig[ci][7] = cell_index(i+dx, j+dy, k+dz, nX, nY, nZ);
                  *sd = 0;
                  sim->lstHNeig[ci][8] = cell_index(i+dx, j+dy, k+dz, nX, nY, nZ);
                  *sd = -1;
                  sim->lstHNeig[ci][9] = cell_index(i+dx, j+dy, k+dz, nX, nY, nZ);

                  *fd = -1;
                  *sd = 1;
                  sim->lstHNeig[ci][10] = cell_index(i+dx, j+dy, k+dz, nX, nY, nZ);
                  *sd = 0;
                  sim->lstHNeig[ci][11] = cell_index(i+dx, j+dy, k+dz, nX, nY, nZ);
                  *sd = -1;
                  sim->lstHNeig[ci][12] = cell_index(i+dx, j+dy, k+dz, nX, nY, nZ);

               }
            }

    } // end 'stack' geometry

    //! добавить верификацию, что список создан и не избыточен

/*
    for (i = 0; i < sim->nHead; i++)
      {
         printf("neighobrs[%d](%d): ", i, sim->nHeadNeig[i]);
         for (j = 0; j < sim->nHeadNeig[i]; j++)
           printf("%d ", sim->lstHNeig[i][j]);
         printf("\n");
      }
*/
    //! ¬ќ“ Ё“ј ¬≈–»‘» ј÷»я
    int c1, c2, l, m;
    int is[26] = {1, 1, 1, 1,  1,  1,  1,  1,  1, 0, 0, 0,  0,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    int js[26] = {1, 1, 0, 0, -1,  1, -1,  0, -1, 1, 1, 0, -1,  1, -1,  0, -1,  1,  1,  0,  0, -1,  1, -1,  0, -1};
    int ks[26] = {1, 0, 1, 0,  1, -1,  0, -1, -1, 1, 0, 1,  1, -1,  0, -1, -1,  1,  0,  1,  0,  1, -1,  0, -1, -1};
    for (i = 0; i < nX; i++)
      for (j = 0; j < nY; j++)
        for (k = 0; k < nZ; k++)
          {
             c1 = cell_index(i, j, k, nX, nY, nZ);

             for (m = 0; m < 26; m++)
             {
                c2 = cell_index(i+is[m], j+js[m], k+ks[m], nX, nY, nZ);
                if (c1 != c2)
                   {
                      nN = 0;
                      for (l = 0; l < sim->nHeadNeig[c1]; l++)
                        if (sim->lstHNeig[c1][l] == c2)
                          nN++;

                      for (l = 0; l < sim->nHeadNeig[c2]; l++)
                        if (sim->lstHNeig[c2][l] == c1)
                          nN++;

                      if (nN == 0)
                        printf("WARNING![134] pair %d-%d is absent!\n", c1, c2);
                      if (nN > 1)
                        printf("WARNING![135] pair %d-%d is excess!, %d\n", c1, c2, nN);
                   }
               }


          }


    return nX * nY * nZ;
}
// end init_clist function

void free_clist(Sim *sim)
// free memory from cell-list method arrays
{
    int i;

    delete[] sim->clist;
    delete[] sim->chead;

    for (i = 0; i < sim->nHead; i++)
      delete[] sim->lstHNeig[i];

    delete[] sim->nHeadNeig;
    delete[] sim->lstHNeig;
}
// end free_clist function


/*

double force_clist(Atoms *atm, int ***head, int *list, CList *cl, VdW ***vdw, Box *bx)
// conventional cell-list calculation of forces
//   return energy
{
    int i, j, k, in;
    int i1, j1, k1;
    int ia, ja;
    double eng; // energy
    //double r2, f;
    //double dx, dy, dz;
    double shX, shY, shZ;
    int N = atm->nAt;

    // reset forces and energy
    eng = 0.0;
    for (i = 0; i < N; i++)
      {
        atm->fxs[i] = 0.0;
        atm->fys[i] = 0.0;
        atm->fzs[i] = 0.0;
      }

    for (i = 0; i < cl->nx; i++)
      for (j = 0; j < cl->ny; j++)
        for (k = 0; k < cl->nz; k++)
          {
              ia = head[i][j][k];
              //printf("head[%d][%d][%d]\n", i, j, k);
              while (ia)
                {

                   //! взаимодействи€ внтури €чейки
                   //! а может тут сделать ja = list[ia-1]; ??? надо подумать, тогда это исключит саму необходиомсть проверки ia <> ja
                   ja = head[i][j][k];  // the same cell!
                   while (ja)
                     {
                        // skip i-i iteration
                        if (ja == ia)
                          break;

                        //! REMEMBER! in HEAD and LIST indexes shifted to +1
                        eng += pair1_nobox(atm, ia-1, ja-1, vdw);
                        ja = list[ja - 1]; // next particle (-1) as we keep (index+1) to indicate empty as 0
                     }

                   //по сосед€м (13 соседей)
                   //  тут не надо будет делать проверочку (ia == ja)
                   //   но надо будет добавл€ть сдвиги, если вышли за границу
                   //! важно! в отличие от традиц алгоритмов делаем только по половине соседей
                   //! Ќќ! это работает только дл€ периодических границ!
                   //printf("begin neig\n");
                   for (in = 0; in < nNei; in++)
                     {
                        i1 = i + neis[in][0];
                        j1 = j + neis[in][1];
                        k1 = k + neis[in][2];

                        //! вообще по хорошему надо наверное с else сделать
                        shX = 0.0; shY = 0.0; shZ = 0.0;
                        if (i1 >= cl->nx)
                          {i1 = 0; shX = bx->ax;}
                        if (i1 < 0)
                          {i1 = cl->nx - 1; shX = -bx->ax;}

                        if (j1 >= cl->ny)
                          {j1 = 0; shY = bx->by;}
                        if (j1 < 0)
                          {j1 = cl->ny - 1; shY = -bx->by;}

                        //! а тут не надо провер€ть что меньше 0, k мы всегда не уменьшаем
                        if (k1 >= cl->nz)
                          {k1 = 0; shZ = bx->cz;}

                        //printf("nei[%d][%d][%d]\n", i1, j1, k1);
                        ja = head[i1][j1][k1];
                        while (ja)
                          {
                             //printf("pair: ia=%d, ja=%d\n", ia, ja);
                             eng += pair1_nobox_sh(atm, ia-1, ja-1, vdw, shX, shY, shZ);
                             ja = list[ja-1];
                          }


                     }
                   // end cycle by neigbours cells

                   ia = list[ia - 1]; // next particle
                }
              // end cycle by ia

          }
          // end cycle by head[i][j][k]
    //printf("end force\n");

    return eng;
}
//end fill_clist function
*/

/*
//neighbourgs for cell-list method
const int nNei = 13;
const int neis[nNei][3] = {{1,1,1},{0,1,1},{1,0,1},{0,0,1},{-1,1,1},{1,-1,1},{-1,0,1},{0,-1,1},{-1,-1,1},{1,1,0},{1,0,0},{0,1,0},{1,-1,0}};
*/

/*
void fill_clist(Atoms *atm, int ***head, int *list, CList *cl)
// fill cell list arrays by actual values
{
    int i, cx, cy, cz;
    int nat = atm->nAt;

    // clear list
    for (cx = 0; cx < cl->nx; cx++)
      for (cy = 0; cy < cl->ny; cy++)
        for (cz = 0; cz < cl->nz; cz++)
          head[cx][cy][cz] = 0;


    for (i = 0; i < nat; i++)
      {
         cx = (int)(atm->xs[i] / cl->Lx);
         cy = (int)(atm->ys[i] / cl->Ly);
         cz = (int)(atm->zs[i] / cl->Lz);
         list[i] = head[cx][cy][cz];
         head[cx][cy][cz] = i + 1;
         //printf("iter: %d, head[%d][%d][%d]=%d list=%d\n", i, cx, cy, cz, head[cx][cy][cz], list[i]);
      }
}
*/

