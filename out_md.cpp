#include <stdio.h>   // FILE, fprintf, scanf

#include "dataStruct.h"  // Sim, Box, Atoms ....
#include "box.h"
#include "out_md.h"

void info_header(FILE *f, Atoms *at, Spec *sp)
{
   char *s0 = sp[at->types[0]].name;
   char *s1 = sp[at->types[1]].name;
   char *s2 = sp[at->types[2]].name;
   char *s3 = sp[at->types[3]].name;
   char *s18 = sp[at->types[18]].name;
   char *s20 = sp[at->types[20]].name;
   fprintf(f, "time iStep x0 y%s0 x1 y%s1 x2 y%s2 x3 y%s3 x18 y%s18 x20 y%s20\n", s0, s1, s2, s3, s18, s20);
   fprintf(f, "time,ps iStep x,A y,A x,A y,A x,A y,A x,A y,A x,A y,A x,A y,A\n");
}

void out_info(FILE *f, double tm, int step, Atoms *at, Spec *sp)
{
   fprintf(f, "%f %d %f %f %f %f %f %f %f %f %f %f\n", tm, step, at->xs[0], at->ys[0], at->xs[1], at->ys[1], at->xs[2], at->ys[2], at->xs[3], at->ys[3], at->xs[18], at->ys[18], at->xs[20], at->ys[20]);
}

void history_header(FILE *f)
{
   //f = fopen("hist.dat", "w");
   fprintf(f, "time iStep totEn temp atm1x atm1y atm1ch momXn momXp momYn momYp momZn momZp\n");
   fprintf(f, "time,ps iStep totEn,eV temp,K atm[1].x,A atm[1].y,A atm1ch,e momXn momXp momYn momYp momZn momZp\n");
}

void msd_header(FILE *f, Sim *sim, Field *field)
{
   int i;

   fprintf(f, "Time\tStep");
   for (i = 0; i < field->nSpec; i++)
     {
        fprintf(f, "\t%s%s\t%s%s\t%s%s", field->species[i].name, "-msd", field->species[i].name, "-nOyz", field->species[i].name, "-pOyz");
     }
   fprintf(f, "\n");
}

void stat_header(FILE *f, Sim *sm, Spec *sp)
{
   int i;

   fprintf(f, "Time\tStep\tTemp\tpotE\tpotE1\tkinE\ttotE");
   // the number of species which can change
   for (i = 0; i < sm->nVarSpec; i++)
     fprintf(f, "\t%s", sp[sm->varSpecs[i]].name);
   fprintf(f, "\n");

   fprintf(f, "Time,ps\tStep\tTemp,K\tpotE,eV\tpotE1,eV\tkinE,eV\ttotE,eV");
   for (i = 0; i < sm->nVarSpec; i++)
     fprintf(f, "\t%s", sp[sm->varSpecs[i]].name);
   fprintf(f, "\n");
}

void out_stat(FILE *f, double tm, int step, Sim *sm, Spec *sp)
{
   int i;

   fprintf(f, "%f\t%d\t%f\t%f\t%f\t%f\t%f", tm, step, sm->Temp, sm->engVdW + sm->engEwaldReal, sm->engEwaldRec, sm->engKin, sm->engTot);
   if (sm->nVarSpec)
     for (i = 0; i < sm->nVarSpec; i++)
       fprintf(f, "\t%d", sp[sm->varSpecs[i]].number);
   fprintf(f, "\n");
}

/*
void stat_out(FILE *f, )
{
   fprintf(f, "%f\t%d\t%f\t%f\t%f\t%f\t%f\n", tSim, iSt, Temp, potE, potE1, kinE, totE);
}
*/

int out_atoms(Atoms *atm, int N, Spec *spec, Box *bx, char *fname)
// write coordinates of atoms (.XYZ specification)
{
   int i;
   FILE *of;

   of = fopen(fname, "w");
   if (of == NULL)
     return 0; // error

   fprintf(of, "%d\n", N);

   //box parameters saving:
   if (bx -> type == 1)
     fprintf(of, "%d %f %f %f\n", bx->type, bx->ax, bx->by, bx->cz);

   for (i = 0; i < N; i++)
     {
        fprintf(of, "%s %f %f %f\n", spec[atm->types[i]].name, atm->xs[i], atm->ys[i], atm->zs[i]);
     }
   fclose(of);
   return 1; // success
}
// end 'out_atoms' function

int out_msd(FILE *f, Atoms *atm, int N, Spec *spec, int NSp, Box *bx, double tm, int tst)
// write msd in open file f
{
   int i, j;
   double dx, dy, dz;


   for (i = 0; i < NSp; i++)
     {
        spec[i].number = 0; //! вообще нужно убрать это отсюда, а менять эти количества всегда, когда до этого доходит дело
        spec[i].displ = 0.0;
     }

   fprintf(f, "%f\t%d", tm, tst);
   for (i = 0; i < N; i++)
     {
        dx = atm -> xs[i] - atm -> x0s[i];
        dy = atm -> ys[i] - atm -> y0s[i];
        dz = atm -> zs[i] - atm -> z0s[i];
        rect_periodic(dx, dy, dz, bx);


        j = atm->types[i];

        spec[j].number++;
        spec[j].displ += dx * dx + dy * dy + dz * dz;

     }

   for (i = 0; i < NSp; i++)
     {
        fprintf(f, "\t%f\t%d\t%d", spec[i].displ / spec[i].number, spec[i].nOyz, spec[i].pOyz);
     }
   fprintf(f, "\n");
}
