// ћќƒ”Ћ№ ƒЋя –ј—„≈“ј » ¬џ¬ќƒј –јƒ»јЋ№Ќќ… ‘”Ќ ÷»» –ј—ѕ–≈ƒ≈Ћ≈Ќ»я

#include <stdlib.h>  // malloc, alloc, rand, NULL
#include <stdio.h>   // FILE
#include <math.h>   // log, sqrt

#include "dataStruct.h"  // Sim, Box, Atoms ....
#include "utils.h"  // int_size, pointer_size, etc...
#include "const.h"  // sphera consts
#include "box.h"  // rect_periodic
#include "rdf.h"

int init_rdf(double ***rdf, Sim *sim, Field *field, Box *box, int &nRDF)
// create arrays and values for RDF calculation and print
{
   int i, j;
   double minR = sim->maxRDF;
   double **gr;

   nRDF = 0;

   //seek min of sim.maxR or box.maxLength
   if (minR > box->maxLength)
     minR = box->maxLength;

   sim->nRDF = minR / sim->dRDF;
   if (!sim->nRDF)
     return 0;  // dR > R - no RDF points, go away

   // create arrays of nPair arrays
   *rdf = (double**)malloc(field->nPair * pointer_size);
   gr = *rdf;

   if (gr == NULL)
     return 0;

   for (i = 0; i < field->nPair; i++)
     {
        gr[i] = (double*)malloc(sim->nRDF * double_size);
        for (j = 0; j < sim->nRDF; j++)
          gr[i][j] = 0.0;
     }

   //printf("parameters: nRDF=%d, dRDF=%f, maxRDF=%f, box.Length=%f\n", sim->nRDF, sim->dRDF, sim->maxRDF, box->maxLength);
   return 1;
}
// end 'init_rdf' function

void clear_rdf(double **rdf, Sim *sim, Field *field, int &nRDF)
// set RDF values to 0.0
{
   int i, j;

   for (i = 0; i < field->nPair; i++)
     for (j = 0; j < sim->nRDF; j++)
       rdf[i][j] = 0.0;

   nRDF = 0;
}
// end 'clear_rdf' function

int get_rdf(Atoms *atm, double **rdf, Sim *sim, Field *field, Box *box, int &nRDF)
// calculate radial distribution function
{
   int i, j, iR;
   int iMin, iMax, iPair;
   double dx, dy, dz, r2, r;
   int m = field->nSpec - 1;
   int N = atm->nAt;

   //printf("get_rdf: nRDF=%d  dRDF=%f   maxRDF=%f\n", sim->nRDF, sim->dRDF, sim->maxRDF);
   for(i = 0; i < N-1; i++)
     for(j = i + 1; j < N; j++)
       {
          dx = atm->xs[i] - atm->xs[j];
          dy = atm->ys[i] - atm->ys[j];
          dz = atm->zs[i] - atm->zs[j];

          //periodic boundary:
          rect_periodic(dx, dy, dz, box);

          r2 = dx * dx + dy * dy + dz * dz;
          r = sqrt(r2);
          iR = r / sim->dRDF;
          if (iR < sim->nRDF)
            {
               iMin = atm->types[i];
               iMax = atm->types[j];
               if (iMin > iMax)
                 {
                    iMin = atm->types[j];
                    iMax = atm->types[i];
                 }
               iPair = iMin * m + iMin * (1 - iMin) / 2 + iMax;
               rdf[iPair][iR] += 1.0;
            }
       }

   nRDF++;
}
// end 'get_rdf' function

int out_rdf(double **rdf, Field *field, Box *box, Sim *sim, char *fname, int nRDF)
// write radial distribution function
{
   int i, j, k;
   int nP = field->nPair;
   double *nAnB = (double*)malloc(nP * double_size); // array for keeping production Na * Nb for every pair
   FILE *of;
   double dr3 = sim->dRDF * sim->dRDF * sim->dRDF;

   //! почему rho не используетс€ нигде?
   // плотность системы:
   //double rho = sim->nAt * box->rvol;  //! перенести в sim или box const? среднее число атомов в единице объЄма

   of = fopen(fname, "w");
   if (of == NULL)
     return 0; // error, can't open file, exit from function


   k = 0;
   fprintf(of, "r ");
   for (i = 0; i < field->nSpec; i++)
     for (j = i; j < field->nSpec; j++)
       {
         fprintf(of, "%s-%s ", field->species[i].name, field->species[j].name);
         nAnB[k] = field->species[i].number * field->species[j].number;
         if (i == j)
           nAnB[k] *= 0.5; // дл€ одинаковой пары мы должны домножить результат на 2 (или в два раза уменьшить знаментаель)
         //printf("n[%s]n[%s]=%d\n", spec[i].name, spec[j].name, nAnB[k]);
         k++;
       }
   fprintf(of, "\n");

   for (i = 0; i < sim->nRDF; i++)
     {
       fprintf(of, "%f ", (i + 0.5)*sim->dRDF);
       for (j = 0; j < nP; j++)
         {
            //! нормализуем на плотность частиц ид газа. т.е. найти, сколько в среднем частиц приходитс€ на этот шаровой слой и нормализовать.
            //rdf[j][i] /= (sphera * (3 * i * (i + 1) + 1) * dr3 * rho * sim->nAt * nRDF);

            //! у мен€ получилась формула g(AB) = npair(AB) * V / v(сло€) / Na / Nb
            //! может быть она неправильна€. ну и конечно нужно разделить на число выборок
            if (nAnB[j])
              rdf[j][i] = rdf[j][i] * box->vol / (sphera * (3 * i * (i + 1) + 1) * dr3) / nAnB[j] / nRDF;
            fprintf(of, "%4.2E ", rdf[j][i]);
         }
       fprintf(of, "\n");
     }

   fclose(of);
   delete[] nAnB;
   return 1; // success
}
// end 'out_rdf' function

int free_rdf(double ***rdf, int nPair)
// clear RDF arrays
{
   int i;
   double **gr = *rdf;

   for (i = 0; i < nPair; i++)
     delete[] gr[i];
     //gr[i][0] = 0.0;

   delete[] gr;
   //sim->nRDF = 0;

   return 1;
}
// end 'free_rdf' function
