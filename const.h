#ifndef CONST_H
#define CONST_H
#include <math.h> // sqrt

//physical constant:
const double pi = 3.14159265359;
const double twopi = 2.0 * pi;         // 2pi
const double twopi2 = twopi * twopi;   // (2pi)^2
const double sqrtpi = sqrt(pi);
const double sphera = 4.0 * pi / 3.0;  // factor before volume of sphere

//internal units in CI
const double r_SI  = 1.0E-10;  // internal units of length in SI (m) //
const double t_SI  = 1.0E-12;  // internal untis of time (pikosecond) in SI (s)
const double E_SI = 1.60217733E-19;  // internal value of engergy SI(J) - electronvolt
const double q_SI = 1.60217657E-19; // internal units of charge in SI(C) - charge of proton

const double kB_SI = 1.3806488E-23;  // Boltzman constant in J/K
const double e0_SI = 8.854187817E-12;// value of electric constant in F/m

//derived untis in SI
const double m_SI = E_SI * t_SI * t_SI / r_SI / r_SI;  // E = m*r^2/t^2 - internal units of mass in SI(kg)
const double F_SI = E_SI / r_SI;  // F = m*r/t^2 - internal units of force in SI (N)
const double Fcoul_SI = 0.25 / pi / e0_SI * q_SI * q_SI / r_SI / r_SI; // force electron affected on another electron at unit length in SI(N)

//traditional units in SI
const double eV_SI = 1.60217733E-19; // electronovolt in SI
const double ang_SI = 1.0E-10;       //  angstrom in SI
const double ps_SI = 1.0E-12;        // pikosecond in SI
const double amu_SI = 1.6605402E-27; // value of atomic mass unit in SI (kg)
const double echarge_SI = 1.60217657E-19; // charge of proton

//factors for translating traditional units (inputs) to internal units:
const double r_scale = ang_SI / r_SI;  // must be ang_SI / r_SI
const double t_scale = 1.0;  //! must be ps_SI / t_SI
const double E_scale = 1.0;  //! must be eV_SI / E_SI
const double q_scale = 1.0; //!  must be echarge_SI / q_SI
const double m_scale = amu_SI / m_SI;
const double Fcoul_scale = Fcoul_SI / F_SI;

//fundamental constants in internal units
const double kB = kB_SI / (E_scale * eV_SI);    // Boltzman konstant - energy of 1K of temperature  (E/K)
//const double r3kB = 1.0 / kB / 3.0;             //  1/3kB,  kB - Boltzman const  // 3 - dimensionity of the system
const double rkB =  1.0 / kB;              // invert Boltzman const

#endif  /* CONST_H */
