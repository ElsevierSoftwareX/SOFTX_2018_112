#ifndef DATASTR_H
#define DATASTR_H

// electron jump types
const int ejt_eq = 0;       // equality (Frank-Condon principle)
const int ejt_min = 1;      // minimization
const int ejt_metr = 2;     // Metropolis scheme (comparision with kB*T)

// pre definitions...
struct Atoms;
struct VdW;
struct Spec;
struct Box;
struct Bond;
struct BondArrays;  //! устареет, будет в sim
struct Field;

// structure for keeping number of particles flow through the box:
struct ThrowBox
{
  int nx, px; // negative x, positive x
  int ny, py; //...
  int nz, pz; //...
};

// structure for simulation parameters
struct Sim
{
   double tSt;    // (d-p) value of timestep (in MD units)      (d-p) directive - parameter
   //! не обязательно его хранить
   double tSim;   // (d-p) length of simulation (in MD time units)
   int nSt;       // (d-p) number of timestep

   //! не обязательно его хранить
   double tEq;    // (d-p) length of equilibration period (in MD time units)
   int nEq;       // (d-p) number of timestep in equilibration period
   int freqEq;    // (d-p) freqence of equilibration (equilibration will be called every freqEq step)

   //SYSTEM PROPERTIES:
   //double engPot;           // (pr) potential energy
   double engTot;           // (pr) total energy
   double engKin;           // (pr) kinetic energy
   double engVdW;           // (pr) Van der Waals energy
   double engEwaldConst;    // (pr) 'Constant' part of Ewald Sum (depends on volume and charge of system)
   double engEwaldRec;      // (pr) Recipropal part of Ewald Sum
   double engEwaldReal;     // (pr) Real part of Ewald Sum
   double engElecField;     // (pr) Energy of electric field
   double engBond;          // (pr) Energy of bonded interaction
   double engAngle;         //  (pr) Energy of valent angles
   double Temp;             // (pr) Instanteous temperature


   double tTemp;   // (d-p) target temperature (in K)
   double tKin;    // (dp) target kinetic energy, according to temperature (sigma in DL_POLY source)  (dp) - derived parameter
   int degFree;     // (p) number of degree freedom (3*N) in the simplest case
   double revDegFree;       // (p)  1/degFree   p - parameter

   //int nAt;         // (n)  number of atoms                     n - the number
   //int nSpec;       // (n)  number of species (atom types)
   //int nBtypes;     // (n)  number of bond types

   //vars for cell-list method
   int *clist;        // (a)    atom list [NAtom] for cell-list method  (a) - array
   int *chead;     //   (a) head list [Ncell] for cell-list method
   int nHead;         // (n)    the number of cells (cnX*cnY*cnZ)
   int *nHeadNeig;      // (a)  the numbers of cells neighbors
   int **lstHNeig;      // (a)  lists of the cells neighbors
   int cnX, cnY, cnZ; // (n)    numbers of cells in X, Y and Z-directions
   int cnYZ;            // (n)  cnY * cnZ
   double clX, clY, clZ;  // (dp)    cell length in X,Y and Z-directions
   int useClist;        // (d-f)    flag

   //vars for neighbors list
   int *nNbors;      // numbers[nAt] of neighbors of atoms
   int **nbors;        // indexes of neighbors of atom
   int **tnbors;        // types of neighbors (VdW, Coulomb, bonding, e-jump, bond formation)
   double **distances;     // distances to neighbors
   int maxNbors;        // maximal number of neighbors in list


   //max VdW range (square):
   double mxRvdw2;      // maximal cutoff distance for VdW interactions

   double r2Max;        // maximal cutoff distance (among VdW, real Ewald, etc)


   //vars for Ewald
   int no_elec;     // (d-f) flag to ignore Coulombic forces
   double eps;   // permittivity,  e in 1/(4pi * e * e0)
   int kx;
   int ky;
   int kz;
   double alpha;
   double rkcut2;
   double r2Real;   //r^2 of cutoff radii for real part of Ewald summ

   //vars for thermostates
   int gaussian;       // (d-f) use or not initial velocities
   int nvt;            // (d-f) nvt ensemble
   double tauT;        // termostat relaxation time
   double qMass;       // qMass = 2*tKin*tauT^2;
   double rQmass;      // 1/Qmass = 1/(2*tKin*tauT^2)
   double qMassTau2;   // qmass/tauT^2 = 2*tKin

   //vars for RDF output
   double dRDF; // width of RDF step
   double maxRDF; // maximum of RDF //! attention: if max(box.length) > maxR, RDF is only to box.length
   int nRDF; // number of RDF points (maxR / dRDF)
   int frRDF; // frequency of RDF statistics

   //vars for E-Jump
   int eJump;      //    the number of eJumps routine at every timestep (if negative - every Nth step)
   double rElec;    // (d-p)    maximal length of electron jump
   double r2Elec; //    (dp)    square of rElec
   double dEjump; // criteria then dE is assumed to zero
   int ejtype;     // (d-f) type of electron jump (see constants above)
   //double *EngDifs;     // for saving energy difference between jumps
   //int nJumpVar;        // the number for jumps variants

   //int varSpecNumb; // (f)  the variable number of atoms of one specie
   int *varSpecs;   // indexes of species with variable quantity
   int nVarSpec;    // the number of species with variable quantity

   int nFreeEl;     // number of free electrons in the system (number of electrons avaliable for jumps)
   int* electrons;      // array of electron number sites
   //int *nPairs;  // number of neightbours for jumping near atom [i]
   //double **dist; // distances
   //int **acceptors; // neightbours
   //int **probs;  // probability of jump
   //int *totProbs; // total probability
   int nJump;       // counter of electron jumps
   ThrowBox elJump;    //! заменить тупо на 2 переменные, убрать эту структуру
   int **jumps; // array for keeping nJump between different donor-acceptor pairs
   int pEjump, nEjump;  // (c)  the number of electron jumped in positive and negative direction through the center of box

   //vars for bonding
   int Bonding;  // (d-f) flag for bonding during MD run
   double rBond, r2Bond; // r and r^2 for bonding distance
   int only_twoatomic;      // flag for two-atomic bond creating optimization
   int *bondAtoms;      // array for keeping indexes of bonded atom. -1 value if no bonding (for only_twoatomic)
   int *bondTypes;      // array for keeping bond types (for only_twoatomic)
   int nBonds;          // number of bonds in the system
   int nBndForm, nBndBr;    // (c) the number of formed and breaking bonds (counters)
   int maxCanBondAtm;  // (d-m) the maximal number of the atoms which can form bond
   int nCanBondAtm;    // the number of atoms which can form bond
   int *canBondAtms;    // indexes of atoms which can form bond

   // (d) - directive
   // (d-f) - directive flag
   // (d-m) - directive for memory allocating



   //var for External Field
   //int eField;      // use or not external field
   double Ux;       // (d-p)    delta U in x-direction


   //vars for output
   int hist;      // as freq will be outputed history file (every 'hist' step)
   int revcon;    // frequency of REVCON
   int stat;     // frequency of stat file output

   //functions for operating!
     //enery of i-j pair calculations (используется для обхода пар в функциях cell_list или all_pair)
   void (*pair)(int i, int j, Atoms *atm, Field *field, Box *bx, Sim *sim);

};

// box structure
struct Box
{
   int type;

   // a,b,c vectors
   double ax, ay, az;
   double bx, by, bz;
   double cx, cy, cz;

   // invert vectors:
   double iax, iay, iaz;
   double ibx, iby, ibz;
   double icx, icy, icz;

   //cell perpendicular widths for the invert matrix:
   double ip1, ip2, ip3;

   double det; // determinant

   double cosA, cosB, cosC; // cos of angles between a,b,c.  A - between b and c; B - a and c; C - a and b

   //for rectangular geometry simple la = xx, lb = yy, lc = zz;
   double la, lb, lc; // length of box vectors
   double maxLength;  // maximum of la, lb, lc
   double ra, rb, rc; //  reversible length of vectors a, b, c
   double ha, hb, hc; // half of length a, b, c
   double nha, nhb, nhc;  // negative ha, hb, hc = -0.5*la, ...

   double vol; // volume
   double rvol; // 1/vol

   double ra2, rb2, rc2; // ra^2, rb^2, ...

   // Impulses from particles for pressure determination
   // momentums in X,Y,Z coord in positive(p) or negative (n) direction
   double momXp, momXn, momYp, momYn, momZp, momZn;
};

/*
// structure for cell-list
struct CList
{
    int nx, ny, nz;         // number of cells in corresponding directions
    double Lx, Ly, Lz;     // length of cell in corresponding direction
};
*/

// structure for Atom type (Species)
//struct Spec;  // predeclaration
struct Spec
{
  int number; // количество частиц данного сорта

  char name[8];
  double mass;
  double rmass; // 1/m
  double charge;
  double energy; // own energy for dE calculation during jump

  int charged; // 0 - neitral
  int donacc; // donor/acceptor  binary flags:  01 - donor 10 - acceptor 00 - no donor, no acceptor, 11 - both
  //Spec *oxForm;    // link to Spec after e-donoring
  //Spec *redForm;   //  ....              e-accepting
  int oxForm;
  int redForm;  // index of spec - redForm of this Spec
  int varNumber;    // (f)  number of particles is variable

  int nFreeEl;  // number of electron available for donoring

  int canBond;  // flag: can create bond with some spec or not
  int *bondKeys; // array of keys+1 for bonds that can be assigned to species with corresponding index
  int angleType;    // (f) possibility to form angle (=0 - can not form angle, >0 - angle with id=angleType


  //counters for moving through edges of the box:
  int pOyz; // icrease, then particle away from box throw the Oyz plane in positive direction
  int nOyz; //                                                              negative

  int pOxz;
  int nOxz;
  int pOxy;
  int nOxy;

  double displ; // displacement (for MSD calculation)
};

//structure for short-range pair iteraction (Van-der-Waals)
struct VdW  // pair potential
{
   int type;
   double p0, p1, p2, p3, p4;
   double r2cut; // cuttoff^2
   double (*eng)(double r2, VdW *vdw); // function to calculate energy
   double (*feng)(double r2, VdW *vdw, double &eng); // function to calculate force (return) & energy (save in eng)
   double (*feng_r)(double r2, double &r, VdW *vdw, double &eng); // function to calculate force (return) & energy (save in eng) if r may be knonw
   double (*eng_r)(double r2, double r, VdW *vdw); // return energy by r and r^2
};

struct Atoms
{
   int nAt;
   int* types;  //! заменить на Spec*
   double *xs, *ys, *zs;
   double *vxs, *vys, *vzs;
   double *fxs, *fys, *fzs;
   int *nBonds;  // the number of bonds (provided by bond_list)
   int *parents;    //! index of atom connected with this, однозначно определен только для атомов с одной связью

   //initial coordinates (for MSD calculation)
   double *x0s, *y0s, *z0s;
};

// structure for bond type
struct Bond
{
  int type;  // type of potential
  int breakable;  // flag breaking bond or not
  int spec1, spec2; // type of atoms that connected by this bond type
  int spec1br, spec2br; // type of atoms after bond breaking
  double p0, p1, p2, p3;    // potential parameters
  double rMax2;          // square of distance of bond breaking (if applicapble)
  double energy;        // energy difference between and after bond breaking
};

struct Angle
{
  int type; // potential type (now = 1 harmonic cosinus)
  int central;  // of central atom species type
  double p0, p1, p2;    // parameters
};

struct Ewald
{
    double alpha;
    double scale;
    double scale2;
    double mr4a2;
    double rkcut, rkcut2;

    //arrays:
    double **elc, **els, **emc, **ems, **enc, **ens;
    double *lmc, *lms, *ckc, *cks;
};

// for force field keeping
struct Field
{
    int nSpec;
    int nPair;       // (n)  number of pairs between Spec, including self pair: 0.5*Nsp*(Nsp-1) + Nsp;
    Spec *species;

    VdW   *pairpots;    // array for keeping all possble pair potentials
    VdW ***vdws;    // matrix of pointer to pairpots: vdws[iSpec][iSpec] = &pairpot

    int nBdata;     // number of bond types
    Bond *bdata;    // array of bond types

    // valent bonds:
    int *nbonds;
    int **bonds;
    int **neighs;

    int nBonds;
    int mxBonds;  // maximal number of bonds (memory allocation)
    int *at1;
    int *at2;
    int *bTypes;

    int** bond_matrix;

    int nAdata;     // the number of angle types
    Angle *adata;   // array of angle types

    // valent angles part
    int nAngles;    // the number of valency angles
    int mxAngles;   // maximal number for angles allocation
    int *centrs;    // indexes of central atom
    int *lig1;      // indexes of the first ligand
    int *lig2;      // indexes of the second ligand
    int *angTypes;  // angles types

};

const int spec_size = sizeof(Spec);
const int vdw_size  = sizeof(VdW);
//const int atom_size = sizeof(Atom);

#endif  /* DATASTR_H */
