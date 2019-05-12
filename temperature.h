#ifndef TEMPERATURE_H
#define TEMPERATURE_H

#include "dataStruct.h"

void gauss_temp(Atoms *atm, Spec *spec, Sim *sim);
double nvtscale(Atoms *atm, Sim *sim, double &chit, double &conint);

#endif // TEMPERATURE_H
