// ћќƒ”Ћ№ ƒЋя »Ќ»÷»јЋ»«ј÷»» —»—“≈ћџ » ≈® ¬џ—¬ќЅќ∆ƒ≈Ќ»я
#include <stdlib.h>  // malloc, alloc, rand, NULL
#include <stdio.h>   // FILE
#include <math.h>   // log, sqrt, fabs
#include <string.h>  //strcmp

#include "utils.h"
#include "const.h"
#include "dataStruct.h"  // Sim, Box, Atoms ....
#include "vdw.h"
#include "sys_init.h"
#include "box.h"
#include "cell_list.h"
#include "bonds.h"
#include "angles.h"

//default all - ok  (result = 1), but if something wrong we must set res as 0
//! тут хитра€ конструкци€, € еЄ сам придумал, но не уверен, как она работает,похоже в i кладЄтс€ сравнение fscanf <= 0
//! поэтому если по окончанию цикал i = 1, т.е. тру, то цикл не нашЄл того, что искал
int find_int(FILE *f, const char *templ, int &value)
{
  int i, x;

  rewind(f);
  while (!feof(f) && (i = fscanf(f, templ, &x) <= 0))
    fscanf(f, "%s");

  if (i)
    return 0;
  else
    {
        value = x;
        return 1;
    }
}

int find_int_def(FILE *f, const char *templ, int &value, int def_val)
// with default value
{
  int i, x;

  rewind(f);
  while (!feof(f) && (i = fscanf(f, templ, &x) <= 0))
    fscanf(f, "%s");

  if (i)
    {
        value = def_val;
        return 0;
    }
  else
    {
        value = x;
        return 1;
    }
}

int find_double(FILE *f, const char *templ, double &value)
{
  int i;

  rewind(f);
  while (!feof(f) && (i = fscanf(f, templ, &value) <= 0))
    fscanf(f, "%s");

  if (i)
    return 0;
  else
    return 1;
}

// seek spec id by name and save it in id-variable. return success or not
int spec_by_name(Field *field, char *name, int &id)
{
    int i;

    for (i = 0; i < field->nSpec; i++)
      if (strcmp(field->species[i].name, name) == 0)
      {
         id = i;
         return 1;
      }

    return 0;
}

//create atoms struct
Atoms* init_atoms_and_box(Field *field, Sim *sim, Box *&box)
{
   int i, j, n, stop;
   int res = 1;
   char aname[8];
   double x, y, z;

   FILE  *f = fopen("atoms.xyz", "r");
   if (f == NULL)
     {
         printf("ERROR[007] Can't open file 'atoms.xyz'\n");
         return NULL;
     }

   Atoms *atm = (Atoms*)malloc(sizeof(Atoms));
   fscanf(f, "%d", &n);
   atm->nAt = n;

   box = (Box*)malloc(sizeof(Box));
   // box reading (the second line)
   fscanf(f, "%d", &box->type);
   if (box->type == 1)
     {
        fscanf(f, "%lf %lf %lf", &box->ax, &box->by, &box->cz);
     }
   else
     {
        printf("ERROR[008] Unknown box type!\n");
        fclose(f);
        delete atm;
        delete box;
        return NULL;
     }
   box_prop(box); // calculate all box properties

   atm->types = (int*)malloc(int_size * n);
   atm->xs = (double*)malloc(double_size * n);
   atm->ys = (double*)malloc(double_size * n);
   atm->zs = (double*)malloc(double_size * n);
   atm->vxs = (double*)malloc(double_size * n);
   atm->vys = (double*)malloc(double_size * n);
   atm->vzs = (double*)malloc(double_size * n);
   atm->fxs = (double*)malloc(double_size * n);
   atm->fys = (double*)malloc(double_size * n);
   atm->fzs = (double*)malloc(double_size * n);
   atm->x0s = (double*)malloc(double_size * n);
   atm->y0s = (double*)malloc(double_size * n);
   atm->z0s = (double*)malloc(double_size * n);
   atm->nBonds = (int*)malloc(int_size * n);
   atm->parents = (int*)malloc(int_size * n);

   //Reading atoms:
   sim->nFreeEl = 0;
   for (i = 0; i < n; i++)
   {
        fscanf(f, "%s %lf %lf %lf", aname, &x, &y, &z);
        // seek spec that corresponds to atom
        stop = 1;
        for (j = 0; j < field->nSpec; j++)
          if (strcmp(aname, field->species[j].name) == 0)
            {
               atm->types[i] = j;
               field->species[j].number++;
               sim->nFreeEl += field->species[j].nFreeEl;
               stop = 0;
               break;
            }

        if (stop)
          {
              printf("ERROR[009]! unknown atom[%d] type=%s in atoms.xyz file\n", i+1, aname);
              res = 0;
          }
        atm->xs[i] = x;
        atm->ys[i] = y;
        atm->zs[i] = z;
        atm->x0s[i] = x;
        atm->y0s[i] = y;
        atm->z0s[i] = z;
        atm->vxs[i] = 0;
        atm->vys[i] = 0;
        atm->vzs[i] = 0;
        atm->fxs[i] = 0;
        atm->fys[i] = 0;
        atm->fzs[i] = 0;
        atm->nBonds[i] = 0;
        atm->parents[i] = -1;
   }  // cycle by read atoms
   fclose(f);
   if (res)
     return atm;
   else
   {
      delete atm;
      delete box;
      return NULL;
   }
}

void free_atoms(Atoms **atm)
{
   delete[] (*atm)->types;
   delete[] (*atm)->xs;
   delete[] (*atm)->ys;
   delete[] (*atm)->zs;
   delete[] (*atm)->vxs;
   delete[] (*atm)->vys;
   delete[] (*atm)->vzs;
   delete[] (*atm)->fxs;
   delete[] (*atm)->fys;
   delete[] (*atm)->fzs;
   delete[] (*atm)->x0s;
   delete[] (*atm)->y0s;
   delete[] (*atm)->z0s;
   delete[] (*atm)->nBonds;
   delete[] (*atm)->parents;

   delete (*atm);
}

Field* init_field(Sim *sim)
{
  int i, j, k, n, stop, at1, at2;
  VdW pp;
  Field *field;
  FILE *f;
  char aname[8], bname[8], cname[8];

  f = fopen("field.txt", "r");
  if (f == NULL)
    {
        printf("ERROR[001]! Fatal Error. Can't open Field.txt file\n");
        return NULL;
    }
  field = (Field*)malloc(sizeof(Field));

  //RESET COUNTERS and ETC..
  field->nBonds = 0;
  field->nAngles = 0;

  //READ SPECIES:
  if (find_int(f, " spec %d", n))
  {
     field->nSpec = n;
     field->nPair = (field->nSpec)*(field->nSpec - 1) / 2 + field->nSpec;
     field->species = (Spec*)malloc(field->nSpec * spec_size);
     for (i = 0; i < field->nSpec; i++)
     {
        fscanf(f, "%s %lf %lf %d %d %lf", field->species[i].name, &field->species[i].mass, &field->species[i].charge, &n, &j, &field->species[i].energy);
        field->species[i].number = 0;
        field->species[i].varNumber = 0;
        //! recalculate mass in OUR system
        field->species[i].mass *= m_scale;
        field->species[i].rmass = 1.0 / field->species[i].mass;
        field->species[i].charge *= q_scale;
        if (fabs(field->species[i].charge) < 1.0E-10)
          field->species[i].charged = 0;
        else
          field->species[i].charged = 1;
        field->species[i].donacc = n + 2 * j;
        //s[i].redForm = NULL;
        //s[i].oxForm = NULL;
        field->species[i].redForm = 0;
        field->species[i].oxForm = 0;
        field->species[i].pOyz = 0;
        field->species[i].nOyz = 0;
        field->species[i].canBond = 0;
        field->species[i].bondKeys = NULL;
        field->species[i].angleType = 0;    // ability to form angles

        printf("spec[%d]:'%s': m=%f, q=%f, E=%f\n", i, field->species[i].name, field->species[i].mass, field->species[i].charge, field->species[i].energy);
     }
  }
  else
  {
     printf("ERROR[004]! There is no 'spec' section in the Field.txt file\n");
     return NULL;
  }

  //READ DONOR-ACCEPTORS
  if (find_int(f, " red-ox %d ", n))
  {
    for (i = 0; i < n; i++)
    {
       fscanf(f, " %8s %8s %8s ", aname, bname, cname);
       stop = 1;
       for (j = 0; j < field->nSpec; j++)
         if (strcmp(field->species[j].name, aname) == 0)
           {
              stop = 0;
              break;
           }
       if (stop)
         {
             printf("ERROR[002]: Unknown species in red-ox section: %s %s %s", aname, bname, cname);
             return NULL;
         }

       // define Ox Form (after e-donating)
       if (strcmp(bname, "null") == 0)
         field->species[j].oxForm = 0;
       else
         {
            stop = 1;
            for (k = 0; k < field->nSpec; k++)
              if (strcmp(field->species[k].name, bname) == 0)
                {
                   stop = 0;
                   //s[j].oxForm = &s[k];
                   field->species[j].oxForm = k + 1;
                   field->species[j].varNumber = 1;
                   field->species[k].varNumber = 1;
                   break;
                }
            if (stop)
              {
                 printf("ERROR[003]: Unknown OxForm in red-ox section: %s %s %s", aname, bname, cname);
                 return NULL;
              }
         }

       // define RedForm (after e-accepting)
       if (strcmp(cname, "null") == 0)
         field->species[j].redForm = 0;
       else
         {
            stop = 1;
            for (k = 0; k < field->nSpec; k++)
              if (strcmp(field->species[k].name, cname) == 0)
                {
                   stop = 0;
                   //s[j].redForm = &s[k];
                   field->species[j].redForm = k + 1;
                   field->species[j].varNumber = 1;
                   field->species[k].varNumber = 1;
                   break;
                }
            if (stop)
              {
                 printf("ERROR[004]: Unknown RedForm in red-ox section: %s %s %s", aname, bname, cname);
                 return NULL;
              }
         }

       printf("Red-ox[%d]: %s %s %s\n", i, aname, bname, cname);
    }

     //define nFreeElectrons for all species
     for (i = 0; i < field->nSpec; i++)
     {
       field->species[i].nFreeEl = 0;
       j = i;
       while (field->species[j].oxForm && ((field->species[j].donacc >> 0) & 1))  // while oxForm exists (species can be donor)
        {
           field->species[i].nFreeEl++;
           j = field->species[j].oxForm - 1;
        }
      }

  } // end red-ox Section

  //READ VAN DER WAALS INTERACTIONS
  sim->mxRvdw2 = 0.0;  //! перенести это в ‘иелд?
  if (find_int(f, " vdw %d", n))
  {
     //vdws array allocation
     field->vdws = (VdW***)malloc(field->nSpec * pointer_size);
     for (i = 0; i < field->nSpec; i++)
     {
       //printf("i=%d\n", i);
       field->vdws[i] = (VdW**)malloc(field->nSpec * pointer_size);
       for (j = 0; j < field->nSpec; j++)
         field->vdws[i][j] = NULL;
     }

     field->pairpots = (VdW*)malloc(n * sizeof(VdW));
     for (i = 0; i < n; i++)
     {
            fscanf(f, " %8s %8s %8s %lf %lf %lf %lf", aname, bname, cname, &pp.r2cut, &pp.p0, &pp.p1, &pp.p2);

            //save cuttoff radii for cell_list method:
            if (pp.r2cut > sim->mxRvdw2)
              sim->mxRvdw2 = pp.r2cut;

            if (strcmp(cname, BHM_name) == 0)
              fscanf(f, " %lf %lf", &pp.p3, &pp.p4);
            //printf("a=%s; b=%s; v=%s;\n", aname, bname, vdwnm);
            at1 = 0; at2 = 0;
            for (j = 0; j < field->nSpec; j++)
              {
                 if (strcmp(aname, field->species[j].name) == 0)
                   at1 = j + 1;
                 if (strcmp(bname, field->species[j].name) == 0)
                   at2 = j + 1;
                 if (at1 && at2)
                   break;
              }
            if (!(at1 && at2))
              {
                  printf("ERROR[005]! Unknown atom type in vdw-line: %s   %s   %s\n", aname, bname, cname);
                  return NULL;
              }

            if (!prepare_vdw(cname, pp))
              {
                  printf("ERROR[006]! Unknown potential type in vdw-line: %s   %s   %s\n", aname, bname, cname);
                  return NULL;
              }

            field->pairpots[i] = pp;

            if (field->vdws[at1-1][at2-1] != NULL)
                printf("WARNING[002]: Pair potential between %s and %s redeclarated\n", aname, bname);

            field->vdws[at1-1][at2-1] = &field->pairpots[i];
            field->vdws[at2-1][at1-1] = &field->pairpots[i];
     } // end loop by pair pots number
  }
  else
    {
       printf("WARING[001] no Van-der-Waals iteractions!\n");
    }
  sim->mxRvdw2 *= sim->mxRvdw2;

  //READ BOND TYPES
  field->nBdata = 0;
  if (find_int(f, " bonds %d ", n))
  {
     field->nBdata = n;
     field->bdata = (Bond*)malloc(n * sizeof(Bond));
     //bt = *bTypes;

     // array for default bonds
     field->bond_matrix = (int**)malloc(field->nSpec * pointer_size);
     for (i = 0; i < field->nSpec; i++)
     {
        field->bond_matrix[i] = (int*)malloc(field->nSpec * int_size);
     }
     // read bonds from file:
     for (i = 0; i < n; i++)
       read_bond(i, f, field);
  }

  //READ ANGLE TYPES
  if (find_int(f, " angles %d ", n))
  {
     field->nAdata = n + 1; //[0] reserved for empty angle
     field->adata = (Angle*)malloc(field->nAdata * sizeof(Angle));

     // read angles from file (angle[0] reserved for no angle):
     field->adata[0].type = 0;
     for (i = 1; i < field->nAdata; i++)
       read_angle(i, f, field);

     //printf("adata[0]=%d [1]=%d\n", field->adata[0].central, field->adata[1].central);
  }

  //READ AUTO FORMING ANGLES
  if (find_int(f, " angle_forming %d ", n))
  {
     for (i = 0; i < n; i++)
     {
        fscanf(f, "%s %d", aname, &k);
        if (spec_by_name(field, aname, at1))
          field->species[at1].angleType = k;
        else
          printf("ERROR[017] wrong species(%s) in angle_formin section(%d)\n", aname, i+1);
     }
  }

  fclose(f);
  return field;
}

int init_md(Atoms *atm, Field *field, Sim *sim, Box *bx)
// set up starting MD parameters
{
  int i, j, n;
  double x;
  char s[5];
  double rCutMax = 0.0;
  FILE *f;
  int res = 1; // result of function 0 - wrong 1 - true

  // reset all counters and .... calculated parameters?
  sim->nBonds = 0;
  sim->pEjump = 0;
  sim->nEjump = 0;
  sim->engKin = 0.0;
  sim->Temp = 0.0;

  f = fopen("field.txt", "r");
  if (f == NULL)
    {
        printf("ERROR[001]! Can't open file 'field.txt'\n");
        return 0;
    }

  if (find_int(f, " bond_list %d", n))
  {
      if (!read_bondlist(atm, field))
        res = 0;
      printf("the list of bonds is used(N=%d)!\n", field->nBonds);
  }

  if (find_int(f, " angle_list %d", n))
  {
      read_anglelist(atm, field);
      printf("the list of angles is used(N=%d)!\n", field->nAngles);
  }

  // automatic initial bond connection
  if (find_int(f, " autobonding %d", n))
  {
       //printf("autobond n=%d\n", n);
       fscanf(f, "%lf", &x); // cutoff
       if (!autobonding(atm, field, sim, bx, n, x, f, field->nbonds, field->bonds, field->neighs, 6))
         return 0;
  }

  // define ability of bond creation
  if (find_int(f, " linkage %d", n))
  {
     if (!read_linkage(f, n, field, field->nBdata))
       return 0;
  }
  fclose(f);

  // define array of atoms, which can form bond
  //! дл€ чего это?
  /*
  j = 0;
  for (i = 0; i < atm->nAt; i++)
    if (field->species[atm->types[i]].canBond)
      {
         if (j < sim->maxCanBondAtm)
           {
              sim->canBondAtms[j] = i;
              j++;
           }
         else
           printf("WARNING[115] the maximal number of atoms which can form bond is reached!\n");
      }
  sim->nCanBondAtm = j;
  */


  //4. READ 'CONTROL.TXT' file (analog CONTROL in DLPOLY)
  f = fopen("control.txt", "r");
  if (f == NULL)
    {
        printf("ERROR! Can't open file 'control.txt'\n");
        return 0;
    }

  //! непон€тно, что это за параметр, раньше он где-то использовалс€ в том числе в процедуре only_twoatomic
  /*
  if (!find_int(f, " maxCanBondAtm %d ", sim->maxCanBondAtm))
   {
      printf("ERROR[114]: maxCanBondAtm is not specified. Use 'maxCanBondAtm' directive\n");
      res = 0;
   }
  */

  if (!find_double(f, " timestep %lf ", sim->tSt))
    {
       printf("ERROR: timestep must be declared in control.txt file!\n");
       res = 0;
    }

  // seeking 'timesim' or 'nstep'
  //! нужно сделать проверку на противоречие директив timesim и nstep
  if (!find_double(f, " timesim %lf ", sim->tSim))
    {
       rewind(f);
       while (!feof(f) && (j = fscanf(f, " nstep %d ", &sim->nSt) <= 0))
        fscanf(f, "%s");

       if (!find_int(f, " nstep %d", sim->nSt))
         {
            printf("ERROR! no 'nstep' or 'timesim' directives in control.txt file!\n");
            res = 0;
         }
       else
          sim->tSim = double(sim->nSt * sim->tSt);

    }
  else
    sim->nSt = sim->tSim / sim->tSt;

  // seeking 'timeequil' or 'nequil'
  //! добавить проверку на противоречивость
  sim->nEq = 0;
  if (!find_double(f, " timeequil %lf ", sim->tEq))
    {
       if (!find_int(f, " nequil %d ", sim->nEq))
         {
            printf("WARNING: no 'nequil' or 'timeequil' directives in control.txt file  - there is no equilibration period!\n");
         }
        else
          sim->tEq = double(sim->nEq * sim->tSt);

    }
  else
    sim->nEq = sim->tEq / sim->tSt;

  // seeking 'eqfreq'
  if (sim->nEq)
    if (sim->nEq)
      {
         if (!find_int(f, " eqfreq %d ", sim->freqEq))
           {
              printf("ERROR: no eqfreq directive!\n");
              res = 0;
           }
      }


  // seeking 'temperature'
  if (!find_double(f, " temperature %lf ", sim->tTemp))
    {
       printf("WARNING: temperature was not defined in control.txt file!\n");
    }

  //EWALDs VARIABLES:
  find_int_def(f, " no_elec %d ", sim->no_elec, 0);  // do not use electrostatic flag
  if (!sim->no_elec)
  {
    // seeking 'alpha'
    if (!find_double(f, " alpha %lf ", sim->alpha))
      {
         printf("WARNING: alpha was not defined in control.txt file!\n");
      }

    if (!find_double(f, " rReal %lf ", sim->r2Real))
      {
         printf("WARNING[112]: cutoff radii for real part of Ewald was not defined in control.txt file! Default 9.0 values was used\n");
         sim->r2Real = 9.0;
      }
    if (sim->r2Real > rCutMax)
      rCutMax = sim->r2Real;

    sim->r2Real *= r_scale;
    sim->r2Real = sim->r2Real * sim->r2Real; // square found

    // seeking 'kx'
    if (!find_int(f, " kx %d ", sim->kx))
      {
         printf("WARNING: kx was not defined in control.txt file!\n");
      }

    // seeking 'ky'
    if (!find_int(f, " ky %d ", sim->ky))
      {
         printf("WARNING: ky was not defined in control.txt file!\n");
      }

    // seeking 'kz'
    if (!find_int(f, " kz %d ", sim->kz))
      {
         printf("WARNING: kz was not defined in control.txt file!\n");
      }
    //printf("sim kx = %d, ky = %d, kz = %d\n", sim->kx, sim->ky, sim->kz);


    //! почему-то этот блок вызывает сбой в значении kx, если подставить sim->eps вместо x ???
    sim->eps = 1.0;
    // seeking 'permittivity'
    if (!find_double(f, " permittivity %lf ", x))
      {
         printf("WARNING[131]: permittivity was not defined in control.txt file. Set as 1.0!\n");
      }
    else
      sim->eps = x;

    // seeking 'rkcut2'
    //! вроде это вычисл€ема€ величина, убираем
    /*
    if (!find_double(f, " rkcut2 %lf ", sim->rkcut2))
      {
         printf("WARNING: rkcut2 was not defined in control.txt file!\n");
      }
    */
  }  // end if (!no_elec)


  sim->nvt = 0;
  // seeking 'nvt'
  find_int(f, " nvt %d ", sim->nvt);
  if (sim->nvt)
    {
       //! наверное это нужно в конец
       sim->degFree -= 1; // themperature is constant
       //printf("ATTENTION: NVT simulation defined!\n");
    }

  find_int_def(f, " gaussian %d ", sim->gaussian, 0);

  if (!find_double(f, " tautemp %lf ", sim->tauT))
    {
       if (sim->nvt)
         {
            printf("ERROR[132]: thermostat relaxation time (tautemp) was not defined in control.txt file!\n");
            res = 0;
         }
    }
  else
    {
    }

  //RDF OUTPUT
  // seeking 'maxRDF'
  if (!find_double(f, " maxRDF %lf ", sim->maxRDF))
    {
       printf("WARNING[133]: maxRDF was not defined in control.txt file!\n");
    }

  // seeking 'dRDF' - width of RDF step
  if (!find_double(f, " dRDF %lf ", x))
    {
       printf("ERROR[134]: dRDF was not defined in control.txt file!\n");
       res = 0;
    }
  else
    sim->dRDF = x;

  //seeking frRDF  - frequency of RDF statictics
  if (!find_int(f, " frRDF %d ", n))
    {
       printf("ERROR[135]: frRDF was not defined in control.txt file!\n");
       res = 0;
    }
  else
    sim->frRDF = n;

  // seeking 'eJump' - directive for do ejump procedure
  find_int(f, " eJump %d ", sim->eJump);
  if (sim->eJump != 0)
    {
        fscanf(f, "%lf %s ", &sim->rElec, s);
        if (strcmp(s, "eq") == 0)
        {
          sim->ejtype = ejt_eq;
          // read admissible energy deviation from zero
          fscanf(f, "%lf", &sim->dEjump);
        }
        else if (strcmp(s, "min") == 0)
          sim->ejtype = ejt_min;
        else if (strcmp(s, "metr") == 0)
          sim->ejtype = ejt_metr;
        else
          {
                printf("ERROR[121]: unknown electron jump type in control file!\n");
                return 0;
          }

        sim->r2Elec = sim->rElec * sim->rElec;
    } // if (eJump)
  else
    {}//printf("PARAMETER: eJump swith off!\n");



  // seeking 'Ux' - delta U of external field in x-direction
  if (!find_double(f, " Ux %lf ", x))
    {
       //printf("INFO: no external field\n");
       sim->Ux = 0.0;
    }
  else
    sim->Ux = x;


  // seeking 'hist' - frequency steps for output
  //find_int(f, " hist %d ", sim->hist);
  if (!find_int(f, " hist %d ", n))
    {
    }
  else
    sim->hist = n;

  if (!find_int(f, " stat %d ", n))
    {
        printf("WARNING[133]: stat directive is not specified, default value of 1000 is used\n");
        sim->stat = 1000;
    }
  else
    sim->stat = n;

  // seeking 'revcon' - frequency steps for revcon
  //find_int(f, " hist %d ", sim->revcon);
  sim->revcon = 0;
  if (!find_int(f, " revcon %d ", n))
    {
    }
  else
    sim->revcon = n;

  find_int_def(f, " only_twoatomic %d ", sim->only_twoatomic, 0);
  if (sim->only_twoatomic)
    {
    }

  if (find_double(f, " bonding %lf ", sim -> rBond))
    {
       sim -> Bonding = 1;
       //printf("Bonding with Rcut = %f\n", sim -> rBond);
       sim -> r2Bond = sim->rBond * sim -> rBond;
    }
  else
      sim -> Bonding = 0;

  // cell list:
  sim->useClist = 0;
  n = 0;
  find_int(f, " cell_list %d ", n);
  if (n || sim->eJump)
    {
       sim->useClist = init_clist(atm, sim, bx, rCutMax);    //! перенести это в другое место! хот€ может и не надо
    }

  if (!find_int_def(f, " max_neigh %d ", sim->maxNbors, 50))
   {
      printf("WARNING[113]: max_neigh is not specified. Set to default value, 50!\n");
   }

  //printf("end control reading\n");
  fclose(f);

 return res;
}
// end 'init_md' function

void prepare_md(Atoms *atm, Sim *sim, Field *field)
// calculate derivated parameters
{
  int i, k;

  sim->r2Max = sim->mxRvdw2;

  if (!sim->no_elec)
    if (sim->r2Max < sim->r2Real)
      sim->r2Max = sim->r2Real;

  //! дописать сюда и проверку на остальные рассто€ни€ (eJump, bonded, bondcreation..)

  sim->nVarSpec = 0;
  for (i = 0; i < field->nSpec; i++)
    if (field->species[i].varNumber)
      sim->nVarSpec++;

  sim->varSpecs = (int*)malloc(sim->nVarSpec * int_size);
  k = 0;
  for (i = 0; i < field->nSpec; i++)
    if (field->species[i].varNumber)
      {sim->varSpecs[k] = i; k++;};


  //FINAL CALCULATION FROM INPUT PARAMETERS

  sim->degFree = 3 * atm->nAt - sim->nBonds;
  sim->revDegFree = (double)(1.0 / sim->degFree);
  sim->tKin = 0.5 * sim->tTemp * kB * sim->degFree;
  if (sim->nvt)
    {
       sim->qMass = 2 * sim->tKin * sim->tauT * sim->tauT;
       sim->rQmass = 0.5 / sim->tKin / sim->tauT / sim->tauT; // 1/Qmass = 1/(2*sigma*tauT^2)
       sim->qMassTau2 = 2 * sim->tKin; // qMass / tauT^2
    }
  //set all counters to zero
  sim->nBndBr = 0;
  sim->nBndForm = 0;
}
// end 'prepare_md' function

const char jtypes[3][30] = {"equality(Frank-Condon)", "minimal", "Metropolis"};

void info_md(Sim *sim)
{

   //printf("box[%d](%f x %f x %f)\n", box->type, box->ax, box->by, box->cz);
   printf("MD long %d timesteps of %f ps\n", sim->nSt, sim->tSt);
   if (sim->nvt)
     printf("ensemble: NVT(T = %f K) thermostat: Nose-Hoover(rt = %f ps)\n", sim->tTemp, sim->tauT);
   else
     printf("ensemble: NVE\n");

   if (sim->gaussian)
     printf("initital velocities are scaled with gaussian\n");

   if (sim->useClist)
     printf("used cell list: %d cells\n", sim->useClist);

   printf("ewald parameters: %f[%d %d %d]\n", sim->alpha, sim->kx, sim->ky, sim->kz);

   if (sim->eJump != 0)
   {
       printf("electron jumps (max r: %f A, type: %s) ", sim->rElec, jtypes[sim->ejtype]);
       if (sim->eJump < 0)
         printf("every %d-th timestep", sim->eJump);
       else
         printf("%d times per timestep", sim->eJump);

       if (sim->ejtype == ejt_eq)
         printf(" dE=%f\n", sim->dEjump);
       else
         printf("\n", sim->dEjump);


   }
   printf("nFreeEl: %d\n", sim->nFreeEl);

   //printf("initial : degFree=%d revDegFree=%f\n", ssim->degFree, sim->revDegFree);

}

void free_md(Atoms *atm, Field *field, Sim *sim)
// free memory after MD
{
}
// end 'free_md' function

int init_neighbors(Atoms *atm, Sim *sim)
// create neighbors arrays for different purposes
{
   int i, j;
   int res = 1;

   //sim->maxNbors = 52; //! TEMP (перенести в контрол)
   sim->nNbors = (int*)malloc(atm->nAt * int_size);                 // the number of neighbors
   sim->nbors = (int**)malloc(atm->nAt * pointer_size);             // indexes of neighbors
   sim->distances = (double**)malloc(atm->nAt * pointer_size);      // r to neighbors
   sim->tnbors = (int**)malloc(atm->nAt * pointer_size);            // type of neighbors
   for (i = 0; i < atm->nAt; i++)
   {
      sim->nNbors[i] = 0;
      sim->nbors[i] = (int*)malloc(sim->maxNbors * int_size);
      sim->distances[i] = (double*)malloc(sim->maxNbors * double_size);
      sim->tnbors[i] = (int*)malloc(sim->maxNbors * int_size);
      for (j = 0; j < sim->maxNbors; j++)
        sim->tnbors[i][j] = 0;
   }

   return res;
}
// end 'init_neighbors' function

void free_neighbors(Atoms *atm, Sim *sim)
// free neighbors arrays
{
   int i;
   //printf("nAt=%d\n", atm->nAt);

   for (i = 0; i < atm->nAt; i++)
   {
      delete[] sim->nbors[i];
      delete[] sim->distances[i];
      delete[] sim->tnbors[i];
   }
   delete[] sim->nNbors;
   delete[] sim->nbors;
   delete[] sim->distances;
   delete[] sim->tnbors;
}
// end 'free_neighbors' function

void free_field(Field *field)
{
    int i;

    //! что это?
    //free_bonds(sim, field->nbonds, field->bonds, field->neighs, sim->nAt);

    if (field->nBonds)
    {
          delete[] field->at1;
          delete[] field->at2;
          delete[] field->bTypes;
    }

    if (field->nAngles)
      free_angles(field);

    if (field->nBdata)
    {
       for (i = 0; i < field->nSpec; i++)
         delete[] field->bond_matrix[i];
       delete[] field->bond_matrix;
    }

    for (i = 0; i < field->nSpec; i++)
      delete[] field->vdws[i];
    delete[] field->vdws;

    delete[] field->pairpots;
    delete[] field->bdata;
    delete[] field->adata;
}

