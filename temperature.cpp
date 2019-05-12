#include <math.h>   // log, sqrt, fabs

#include "temperature.h"
#include "dataStruct.h"
#include "md_utils.h"  // gauss(), rand01

void gauss_temp(Atoms *atm, Spec *spec, Sim *sim)
{
   int i;
   double /*VX, VY, VZ,*/ k, m;

   // centre of mass
   double cmx = 0.0;
   double cmy = 0.0;
   double cmz = 0.0;

   // centre of momentum
   double cpx = 0.0;
   double cpy = 0.0;
   double cpz = 0.0;

   double totMass = 0.0;

   double mean = 0.0;
   double stdev = 0.5;
   double kE = 0.0;

   //VX = 0.0; VY = 0.0; VZ = 0.0; // velocity summs
   for (i = 0; i < atm->nAt; i++)
     {
        atm->vxs[i] = gauss(stdev, mean);
        atm->vys[i] = gauss(stdev, mean);
        atm->vzs[i] = gauss(stdev, mean);

        m = spec[atm->types[i]].mass;

        cmx += atm->xs[i] * m;
        cmy += atm->ys[i] * m;
        cmz += atm->zs[i] * m;
        cpx += atm->vxs[i] * m;
        cpy += atm->vys[i] * m;
        cpz += atm->vzs[i] * m;
        totMass += m;


        //VX += atm->vxs[i];
        //VY += atm->vys[i];
        //VZ += atm->vzs[i];

     }

   // total momentum must be zero:
   cmx /= totMass;
   cmy /= totMass;
   cmz /= totMass;
   cpx /= totMass;
   cpy /= totMass;
   cpz /= totMass;
   for (i = 0; i < atm->nAt; i++)
     {
        //atm->vxs[i] -= VX;
        //atm->vys[i] -= VY;
        //atm->vzs[i] -= VZ;

        atm->vxs[i] -= cpx;
        atm->vys[i] -= cpy;
        atm->vzs[i] -= cpz;

        kE += spec[atm->types[i]].mass * (atm->vxs[i] * atm->vxs[i] + atm->vys[i] * atm->vys[i] + atm->vzs[i] * atm->vzs[i]);
     }

   // normalization:
   kE *= 0.5;
   k = sqrt(sim->tKin / kE);
   for (i = 0; i < atm->nAt; i++)
     {
        atm->vxs[i] *= k;
        atm->vys[i] *= k;
        atm->vzs[i] *= k;
     }

   //! если косяк - переделываем
   if (isnan(atm->vxs[0]))
     gauss_temp;
}

// Nose-Hoover thermostat (as in DL_POLY source), return new kinetic energy
double nvtscale(Atoms *atm, Sim *sim, double &chit, double &conint)
{
   int i;
   double kinE;
   double scale;

   //printf("rQ:%f  kE:%f  chit:%f conint:%f\n", sim->rQmass, kE, *chit, *conint);
   chit += sim->tSt * (sim->engKin - sim->tKin) * sim->rQmass;
   scale = 1 - sim->tSt * chit;
   for (i = 0; i < atm->nAt; i++)
     {
        atm->vxs[i] *= scale;
        atm->vys[i] *= scale;
        atm->vzs[i] *= scale;
     }
   kinE = sim->engKin * scale * scale;
   conint += sim->tSt * chit * sim->qMassTau2;
   chit += sim->tSt * (kinE - sim->tKin) * sim->rQmass; // new kinetic energy (отличие от первого действия этой процедуры)

   //if (isnan(kinE))
     //printf("ERROR nvtsacle is NAN\n");
   //printf("kinE=%f\n", kinE);
   return kinE;
}
// end nvtscale function
