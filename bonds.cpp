//MODULE bonds.cpp
//with HEADER bonds.h CONTAINS CONST FOR VDW READING AND APPLYING
#include <stdlib.h>  // malloc, alloc, rand, NULL
#include <stdio.h>   // *FILE
#include <string.h>  //strcmp

#include "utils.h"  // int_size, pointer_size
#include "dataStruct.h" // struct Bond, Atoms
#include "const.h" // r_scale
#include "box.h" // rect, box
#include "vdw.h"  // vdw_iter
#include "angles.h"

// init arrays for bonds:
//      nbonds - array of bond numbers;
//      bonds - array[][] of bonds keys
//      neigs - array[][] of neighbours

/*
int init_bonds(Sim *sim, int **nbonds, int ***bonds, int ***neigs, int def_nbonds)
{
   int i;
   int *nb, **b, **n;
   int nAtm = sim -> nAt;


   *nbonds = (int*)malloc(nAtm * int_size);
   nb = *nbonds;
   *bonds = (int**)malloc(nAtm * pointer_size);
   b = *bonds;
   *neigs = (int**)malloc(nAtm * pointer_size);
   n = *neigs;

   // default values and allocate subarrays:
   for (i = 0; i < nAtm; i++)
     {
        nb[i] = 0;
        b[i] = (int*)malloc(def_nbonds * int_size);
        n[i] = (int*)malloc(def_nbonds * int_size);
     }

   //! это нужно, если св€зи образуютс€ и создаютс€
   sim->canBondAtms = (int*)malloc(sim->maxCanBondAtm * int_size);
   //printf("memory for %d canBondAtoms was allocated\n", sim->maxCanBondAtm);

   if (sim->only_twoatomic)
   {
       sim->bondAtoms = (int*)malloc(sim->nAt * int_size);
       for (i = 0; i < sim->nAt; i++)
         {
            sim->bondAtoms[i] = -1;

         }


    }

    return 1;
}
*/

int read_bondlist(Atoms *atm, Field *field)
{
   int i, j, at1, at2, n, k;
   int res = 1;     // result
   FILE *f = fopen("bonds.txt", "r");
   Bond *btypes = field->bdata;

   fscanf(f, "%d", &field->nBonds);
   field->mxBonds = 1.1 * field->nBonds;
   field->at1 = (int*)malloc(field->mxBonds * int_size);
   field->at2 = (int*)malloc(field->mxBonds * int_size);
   field->bTypes = (int*)malloc(field->mxBonds * int_size);

   for (i = 0; i < field -> nBonds; i++)
   {
      fscanf(f, "%d %d %d", &at1, &at2, &k);
      field->bTypes[i] = k-1;

      // сразу развернем св€зь как надо и проверим, правильно ли св€зана
      if (btypes[k-1].spec1 == atm->types[at1])
      {
          if (btypes[k-1].spec2 != atm->types[at2])
            printf("WARNING [121] wrong bond type in bond list, line: %d\n", i);
      }
      else
        if (btypes[k-1].spec1 == atm->types[at2])
        {
           if (btypes[k-1].spec2 == atm->types[at1])
           {
             n = at1;
             at1 = at2;
             at2 = n;
           }
           else
            printf("WARNING [121] wrong bond type in bond list, line: %d\n", i);
        }
        else
        {
            printf("ERROR [122] wrong bond type in bond list, line: %d\n", i);
            res = 0;
        }


      field->at1[i] = at1;
      field->at2[i] = at2;

      atm->nBonds[at1]++;
      atm->nBonds[at2]++;
      atm->parents[at1] = at2;
      atm->parents[at2] = at1;
      //printf("%d %d %d\n", barrs->at1[i], barrs->at2[i], barrs->types[i]);
   }

/*
   //autoparents verification:
   for (i = 0; i < atm->nAt; i++)
     if (atm->parents[i] == i)
       printf("autoparent %d\n", i);
*/
   fclose(f);
   return res;
}

void save_bondlist(char *fname, Field *field)
{
   int i;
   FILE *f = fopen(fname, "w");

   fprintf(f, "%d\n", field->nBonds);
   for (i = 0; i < field->nBonds; i++)
   {
      fprintf(f, "%d %d %d\n", field->at1[i], field->at2[i], field->bTypes[i] + 1);
   }
   fclose(f);
}

// read bond parameters from file and put them into bond structure
//   (id) - bond number
int read_bond(int id, FILE *f, Field *field)
{
   int n, i, br;
   double r, p0, p1, p2, p3;
   char s1[8], s2[8], key[8];
   int ind1, ind2;
   Bond *bond = &field->bdata[id];

   fscanf(f, "%d %8s %8s %d", &n, s1, s2, &br);
   //printf("reading: %d %s %s %d\n", n, s1, s2, br);

   // define spec indexes
   ind1 = 0; ind2 = 0;
   //printf("nspec = %d\n", nsp);
   for (i = 0; i < field->nSpec; i++)
    {
       if (!ind1)
         if (strcmp(field->species[i].name, s1) == 0)
         {
            ind1 = i + 1;
         }
       if (!ind2)
         if (strcmp(field->species[i].name, s2) == 0)
         {
            ind2 = i + 1;
         }
       if (ind1 && ind2)
         {
            break;
         }
    }
   if (!(ind1 && ind2))
   {
      printf("ERROR[124]: Unknown species in bonds declaration: %d %s %s", n, s1, s2);
      return 0;
   }

   //printf("start bond reading %d %d\n", ind1, ind2);
   bond -> spec1 = ind1 - 1;
   bond -> spec2 = ind2 - 1;
   bond -> breakable = br;

   //! сохраним св€зь как дефолтную между двум€ типами атомов, позже добавим флаг дефолтный в файле spec.txt
   field->bond_matrix[ind1 - 1][ind2 - 1] = id + 1;
   if (ind1 != ind2)
     field->bond_matrix[ind2 - 1][ind1 - 1] = -(id + 1); // negative index means that we need change atom places
   else
     field->bond_matrix[ind2 - 1][ind1 - 1] = id + 1;

   // read breaking bond parameters
   if (br)
   {
      fscanf(f, "%lf %8s %8s", &r, s1, s2);
      ind1 = 0; ind2 = 0;
      for (i = 0; i < field->nSpec; i++)
        {
          if (!ind1)
           if (strcmp(field->species[i].name, s1) == 0)
           {
              ind1 = i + 1;
           }
         if (!ind2)
           if(strcmp(field->species[i].name, s2) == 0)
           {
              ind2 = i + 1;
           }
         if (ind1 && ind2)
           {
              break;
           }
        }
      if (!(ind1 && ind2))
        {
          printf("ERROR[125]: Unknown species in bonds declaration: %d %s %s", n, s1, s2);
          return 0;
        }

      bond -> rMax2 = r * r;
      bond -> spec1br = ind1 - 1;
      bond -> spec2br = ind2 - 1;
      bond -> energy = field->species[bond->spec1].energy + field->species[bond->spec2].energy - field->species[bond->spec1br].energy - field->species[bond->spec2br].energy;
   }

   //read potential parameter
   fscanf(f, "%8s", key);
   if (strcmp(key, "harm") == 0)  // U = 1/2 k (r-r0)^2
     {
        bond -> type = 1;
        fscanf(f, "%lf %lf", &p0, &p1);
        bond -> p0 = p0;  // k (Eng / r^2)
        bond -> p0 *= E_scale;
        bond -> p0 /= (r_scale * r_scale);
        bond -> p1 = p1;  // r0 (r)
        bond -> p1 *= r_scale;
     }
   // Morse potential:
   else if (strcmp(key, "mors") == 0)  // U = D[1 - exp(-a(r-r0))]^2 - C
     {
        bond -> type = 2;
        fscanf(f, "%lf %lf %lf %lf", &p0, &p1, &p2, &p3);
        bond -> p0 = p0;            // D (Eng)
        bond -> p0 *= E_scale;

        p1 /= (r_scale * r_scale);  // a (1/r^2)
        bond -> p1 = p1;

        p2 *= r_scale;              // r0 (r)
        bond -> p2 = p2;

        p3 *= E_scale;              // C (ENg)
        bond -> p3 = p3;
     }
   else
     {
        printf("ERROR[126]: Unknown potential type in bonds declaration: %d %s", n, key);
        return 0;
     }

   return 1;
}

// try to create new bond in bond_list and return success or not
int create_bond(int at1, int at2, int type, Atoms *atm, Field *field)
{
   int j, n, k, t1aft, t2aft, t1bef, t2bef, /*old_type,*/ new_type, bt1, bt2, left;
   int i = field->nBonds;

   //printf("bond[%d]: (%s %s  -> %s %s for %s and %s\n", type, sp[btypes[type].spec1].name, sp[btypes[type].spec2].name, sp[btypes[type].spec1br].name, sp[btypes[type].spec2br].name, sp[atm->types[at1]].name, sp[atm->types[at2]].name);

   // сразу развернем св€зь как надо и проверим, правильно ли св€зана
   if (field->bdata[type].spec1br == atm->types[at1])
      {
          if (field->bdata[type].spec2br != atm->types[at2])
            printf("(A)create wrong bond type[%d] between %s and %s\n", type+1, field->species[atm->types[at1]].name, field->species[atm->types[at2]].name);

          //printf("ok\n");
      }
   else
        //! тут ещЄ проверку правильный ли первый атом

        if (field->bdata[type].spec1br == atm->types[at2])
        {
           n = at1;
           at1 = at2;
           at2 = n;
           //printf("at1=%d, at2=%d n=%d\n", at1, at2, n);
        }
        else
          printf("(B)create wrong bond type[%d] between %s and %s\n", type+1, field->species[atm->types[at1]].name, field->species[atm->types[at2]].name);

   if (i < field->mxBonds)
   {
      field->at1[i] = at1;
      field->at2[i] = at2;
      field->bTypes[i] = type;
      t1bef = atm->types[at1];
      t2bef = atm->types[at2];
      //printf("1\n");

      //change atoms and old bonds to him (+add valent angle! - this is in external place)
      t1aft = field->bdata[type].spec1;
      t2aft = field->bdata[type].spec2;

      // change types of existing bond:
      n = atm->nBonds[at1] + atm->nBonds[at2];
      k = 0;
      for (j = 0; j < field->nBonds; j++)
        if ((field->at1[j] == at1) || (field->at2[j] == at1))
        {
            bt1 = atm->types[field->at1[j]];
            bt2 = atm->types[field->at2[j]];

            if (field->at1[j] == at1)
            {
              new_type = field->bond_matrix[t1aft][bt2];
            }
            else
            {
              new_type = field->bond_matrix[bt1][t1aft];
            }

            if (new_type)
            {
               if (new_type < 0)
               {
                  field->bTypes[j] = - new_type + 1; // as index=0 reserved
                  //change atom places
                  left = field->at1[j];
                  field->at1[j] = field->at2[j];
                  field->at2[j] = left;
               }
               else
               field->bTypes[j] = new_type - 1; // as index=0 reserved
            }

            k++;
            if (k == n)
              break;
        }
        else
        if ((field->at1[j] == at2) || (field->at2[j] == at2))
        {
            bt1 = atm->types[field->at1[j]];
            bt2 = atm->types[field->at2[j]];

            if (field->at1[j] == at2)
            {
              new_type = field->bond_matrix[t2aft][bt2];
            }
            else
            {
              new_type = field->bond_matrix[bt1][t2aft];
            }

            if (new_type)
            {
               if (new_type < 0)
               {
                  field->bTypes[j] = -new_type + 1; // as index=0 reserved
                  //change atom places
                  left = field->at1[j];
                  field->at1[j] = field->at2[j];
                  field->at2[j] = left;
               }
               else
               field->bTypes[j] = new_type - 1; // as index=0 reserved
            }

            k++;
            if (k == n)
              break;
        }
      //printf("3\n");

      //! temp: add angle!
      if (field->species[t1aft].angleType) // O(c)
      {
         //printf("new angle: %s[%d]: %s[%d] %s[%d]\n", field->species[t1aft].name, at1, field->species[t2aft].name, at2, field->species[atm->types[atm->parents[at1]]].name, atm->parents[at1]);
         create_angle(at1, at2, atm->parents[at1], field->species[t1aft].angleType, field);
      }
      else if (field->species[t2aft].angleType)
      {
         //printf("new angle2: %s[%d]\n", field->species[t2aft].name, at2);
         create_angle(at2, at1, atm->parents[at2], field->species[t2aft].angleType, field);
      }
      //printf("end add angle\n");

      atm->types[at1] = t1aft;
      atm->types[at2] = t2aft;
      field->species[t1bef].number--;
      field->species[t2bef].number--;
      field->species[t1aft].number++;
      field->species[t2aft].number++;
      atm->nBonds[at1]++;
      atm->nBonds[at2]++;
      if (atm->parents[at1] < 0)
        atm->parents[at1] = at2;
      if (atm->parents[at2] < 0)
        atm->parents[at2] = at1;
      field->nBonds++;
      //printf("4\n");
      //printf("Bond[%d] %s[%d]-%s[%d] between %s and %s is created!\n", type, sp[t1aft].name, at1, sp[t2aft].name, at2, sp[t1bef].name, sp[t2bef].name);

      return 1;
   }
   else
   {
      printf("WARNING[115] maximal number of bonds is reached!\n");
      return 0;
   }
}

// destroy bond with bnd index
void destroy_bond(int bnd, Atoms *atm, Field *field, Bond *bond)
{
   int a1, a2, p;
   int j, n, k, t1aft, t2aft, t1bef, t2bef, /*old_type, */left, new_type;

   a1 = field->at1[bnd];
   a2 = field->at2[bnd];

   t1bef = atm->types[a1];
   t2bef = atm->types[a2];

   t1aft = bond->spec1br;
   t2aft = bond->spec2br;

   //change types of resting bonds
   k = 0;
   n = atm->nBonds[a1] + atm->nBonds[a2] - 1;
   for (j = 0; j < field->nBonds; j++)
     if (j != bnd)
     {
        //! temp (:???)
        if ((field->at1[j] == a1) || (field->at2[j] == a1))
        {

            //rewrite parent if its necessary:
            p = field->at1[j];
            if (p == a1)
              p = field->at2[j];

            // замен€ем в случае если родитель указан a2
            if (atm->parents[a1] == a2)
               atm->parents[a1] = p;


            if (field->at1[j] == a1)
            {
              new_type = field->bond_matrix[t1aft][atm->types[field->at2[j]]];
            }
            else
            {
              new_type = field->bond_matrix[atm->types[field->at1[j]]][t1aft];
            }

            if (new_type)
            {
               if (new_type < 0)
               {
                  field->bTypes[j] = -new_type + 1; // as index=0 reserved
                  //change atom places
                  left = field->at1[j];
                  field->at1[j] = field->at2[j];
                  field->at2[j] = left;
               }
               else
               field->bTypes[j] = new_type - 1; // as index=0 reserved
            }

            k++;
            if (k == n)
              break;
        }
        else
        if ((field->at1[j] == a2) || (field->at2[j] == a2))
        {
            //rewrite parent if its necessary:
            p = field->at1[j];
            if (p == a2)
              p = field->at2[j];

            if (atm->parents[a2] == a1)
               atm->parents[a2] = p;

            if (field->at1[j] == a2)
            {
              new_type = field->bond_matrix[t2aft][atm->types[field->at2[j]]];
            }
            else
            {
              new_type = field->bond_matrix[atm->types[field->at1[j]]][t2aft];
            }

            if (new_type)
            {
               if (new_type < 0)
               {
                  field->bTypes[j] = -new_type + 1; // as index=0 reserved
                  //change atom places
                  left = field->at1[j];
                  field->at1[j] = field->at2[j];
                  field->at2[j] = left;
               }
               else
               field->bTypes[j] = new_type - 1; // as index=0 reserved
            }

            k++;
            if (k == n)
              break;
        }

     }

   atm->types[a1] = t1aft;
   atm->types[a2] = t2aft;
   field->species[t1bef].number--;
   field->species[t2bef].number--;
   field->species[t1aft].number++;
   field->species[t2aft].number++;
   atm->nBonds[a1]--;
   atm->nBonds[a2]--;
   field->at1[bnd] = -1;   // flag of destroying
   //printf("Bond[%d] between %s and %s is destroyed! Now they are: %s and %s\n", bnd, sp[t1bef].name, sp[t2bef].name, sp[t1aft].name, sp[t2aft].name);
}

// procudure of building bonds in initial configuration
//   rc - cuttoff radius
//      nbonds - array of bond numbers;
//      bonds - array[][] of bonds keys
//      neigs - array[][] of neighbours
//    mxbond - maximal number of bonds
int autobonding(Atoms *atm, Field *field, Sim *sim, Box *bx, int Naut, double rc, FILE *f, int *nbonds, int **bonds, int **neigs, int mxbond)
{
   int i, j, k;
   int nb1, nb2; // number of bonds
   //int ivar;  // index of variant and number of variant
   int **keys, *nvar; // keys and number of variants for autobonding paramenters
   double **rads;  // radiuses for autobonding paramenters
   int N = atm -> nAt;
   int **vars; //variant numbers for [i][j] spec
   char s1[8], s2[8];
   int ind1, ind2;
   double r, dx, dy, dz;

   printf("autobonding: %d %f maxBond=%d\n", Naut, rc, mxbond);
   // create array for a number of variants of bonding pairs
   vars = (int**)malloc(field -> nSpec * pointer_size);
   nvar = (int*)malloc(Naut * int_size);
   keys = (int**)malloc(Naut * pointer_size);
   rads = (double**)malloc(Naut * pointer_size);
   for (i = 0; i < field -> nSpec; i++)
   {
      vars[i] = (int*)malloc(field -> nSpec * int_size);
      for (j = 0; j < field -> nSpec; j++)
        vars[i][j] = 0;
   }

   // read autobonding parameters from FILE f and fill vars array
   for (i = 0; i < Naut; i++)
   {
       fscanf(f, "%8s %8s %d", s1, s2, &k);

       // seeking species numbers
       ind1 = 0; ind2 = 0;
       for (j = 0; j < field -> nSpec; j++)
       {
         if (!ind1)
           if (strcmp(field->species[j].name, s1) == 0)
             {
                ind1 = j + 1;
             }
         if (!ind2)
           if(strcmp(field->species[j].name, s2) == 0)
             {
                ind2 = j + 1;
             }
       if (ind1 && ind2)
         {
            break;
         }
       }

       if (!(ind1 && ind2))
        {
          printf("ERROR[127]: Unknown species in autobonding declaration: %s %s\n", s1, s2);
          return 0;
        }

       nvar[i] = k;
       vars[ind1 - 1][ind2 - 1] = i + 1;
       if (ind1 != ind2)
         vars[ind1 - 1][ind2 - 1] = i + 1;

       //fill keys and rads arrays:
       keys[i] = (int*)malloc(k * int_size);
       rads[i] = (double*)malloc(k * double_size);

       for (j = 0; j < k; j++)
         {
            fscanf(f, "%lf %d", &r, &ind1);
            keys[i][j] = ind1;
            rads[i][j] = r;
         }

   }

   // autobonding
   ind2 = 0;  // will keep number of bonds here...
   for (i = 0; i < N - 1; i++)
     for (j = i + 1; j < N; j++)
       {
          // verify, can this pair be bonded
          ind1 = vars[atm -> types[i]][atm -> types[j]];
          if (ind1)
          {
             ind1 -= 1;
             dx = atm->xs[i] - atm->xs[j];
             dy = atm->ys[i] - atm->ys[j];
             dz = atm->zs[i] - atm->zs[j];

             rect_periodic(dx, dy, dz, bx);  //! need be replaced for NOT rectangular geomery

             r = dx*dx + dy*dy + dz*dz;
             if (r < (rc*rc))
             {
                r = sqrt(r);
                // select correct bond variant
                for (k = 0; k < nvar[ind1]; k++)
                  if (r < rads[ind1][k])
                    {
                       //bond creating
                       nb1 = nbonds[i];
                       nb2 = nbonds[j];

                       if ((nb1 == mxbond)||(nb2 == mxbond))
                         {
                            printf("WARNING[128]: a maximal number of bonds exceeded\n");
                            break;
                         }
                       nbonds[i]++;
                       nbonds[j]++;
                       bonds[i][nb1] = keys[ind1][k] - 1; // because keys indexed from 1, but bonds from 0
                       bonds[j][nb2] = keys[ind1][k] - 1;
                       // save neighbours:
                       neigs[i][nb1] = j;
                       neigs[j][nb2] = i;
                       ind2++;
                       break;
                    }
             }
          }
       }
   printf("%d bonds were created\n", ind2);

   //free arrays
   for (i = 0; i < Naut; i++)
   {
      delete[] keys[i];
      delete[] rads[i];
   }
   delete[] nvar;
   delete[] keys;
   delete[] rads;

   for (i = 0; i < field -> nSpec; i++)
   {
      delete[] vars[i];
   }
   delete[] vars;

   return 1;
}
// end autobonding function

// read ability of bond creating
//   nBonds = number of bond types
int read_linkage(FILE *f, int Nlnk, Field *field, int nBonds)
{
   int i, j, k, n, bond;
   int sp1, sp2; // indexeses of species
   char ion[8], ion2[8];
   int nSp = field -> nSpec;


   //read lines in format: "ion N ion1 k1 ion2 k2 ..."
   //! надо ещЄ где-то определить рассто€ние дл€ автосв€зывани€
   for (i = 0; i < Nlnk; i++)
   {
      fscanf(f, "%s %d", ion, &n);
      sp1 = -1;
      for (j = 0; j < nSp; j++)
        if (strcmp(field->species[j].name, ion) == 0)
        {
           sp1 = j;
           break;
        }

      if (sp1 < 0)
      {
         printf("ERROR[128]: Unknown species %s in linkage definition!\n", ion);
         return 0;
      }

      field->species[sp1].canBond = 1;
      field->species[sp1].bondKeys = (int*)malloc(nSp * int_size);
      for (j = 0; j < nSp; j++)
        field->species[sp1].bondKeys[j] = 0;

      // n - the number of possible species to form bond
      for (j = 0; j < n; j++)
      {
         fscanf(f, "%s %d", ion2, &bond);
         if (bond > nBonds + 1)
         {
            printf("ERROR[129]: Unknown bond number %d in linkage definition!\n", bond);
            return 0;
         }
         sp2 = -1;
         for (k = 0; k < nSp; k++)
           if (strcmp(field->species[k].name, ion2) == 0)
             {
               sp2 = k;
               break;
             }

         if (sp2 < 0)
         {
            printf("ERROR[130]: Unknown species %s in linkage definition!\n", ion2);
            return 0;
         }

         field->species[sp1].bondKeys[sp2] = bond;
      }

   }
   return 1;
}

double bond_iter(double r2, Bond *bnd, double &eng)
// calculate energy and return force of bonded(intramolecular) iteraction  (   Fx/dx = -(1/r)*dU(r)/dr   )
//  r2 - square of distance
{
   double r, x, y;

   switch (bnd->type)
    {
       case 1:   // U = 1/2 k (r-r0)^2   (k = p0, r0 = p1)
         r = sqrt(r2);
         x = r - bnd->p1; // r - r0

         eng += 0.5 * bnd->p0 * x * x;
         return -bnd->p0 / r * x;       // -dU/dr * (1/r)

       // Morse:
       case 2:   // U = D[1 - exp(-a(r-r0))]^2 - C   (D = p0, a = p1, r0 = p2, C = p3)
         r = sqrt(r2);
         x = r - bnd->p2; // r - r0
         x = exp(-bnd->p1 * x); // exp(-a(r-r0))
         y = 1 - x;

         eng += bnd->p0 * y * y - bnd->p3;
         return -2.0 * bnd->p0 * bnd->p1 * x * y / r;
    }

}

double bond_iter_r(double r2, double &r, Bond *bnd, double &eng)
// calculate energy and return force of bonded(intramolecular) iteraction  (   Fx/dx = -(1/r)*dU(r)/dr   )
//  with r version. r - distance (if known, else 0) r2 - square of distance
{
   double x;

   switch (bnd->type)
    {
       case 1:   // U = 1/2 k (r-r0)^2   (k = p0, r0 = p1)

         if (r == 0.0)    // calculate if unkonwn (zero) and use otherwise
           r = sqrt(r2);

         x = r - bnd->p1; // r - r0

         eng += 0.5 * bnd->p0 * x * x;
         return -bnd->p0 / r * x;

    }

}

void exec_bondlist(Atoms *atm, Field *field, Sim *sim, Box *bx)
{
   int i, j, ia, ja, tp, k;
   int nbr = 0;  // the number of breaking bond
   int *bInds;   // indexes of breaking bond
   int mxBr = 20;   // maximum for breaking bond keeping
   double dx, dy, dz, r2, f;
   double eng = 0.0;
   Bond *btypes = field->bdata;
   Spec *sp = field->species;

   k = 0;
   bInds = (int*)malloc(mxBr * int_size);

   i = field->nBonds-1;
   //printf("exec_bondlist. Last bond: %d-%d %d\n", barrs->at1[i], barrs->at2[i], barrs->types[i]);

   for (i = 0; i < field->nBonds; i++)
   {
      ia = field->at1[i];
      if (ia < 0)   // flag of breaking bond
        continue;
      ja = field->at2[i];

/*
      //! temp! bond type correction:
      t1 = atm->types[ia];
      t2 = atm->types[ja];
      barrs->types[i] = barrs->bond_matrix[t1][t2];
*/
      dx = atm->xs[ia] - atm->xs[ja];
      dy = atm->ys[ia] - atm->ys[ja];
      dz = atm->zs[ia] - atm->zs[ja];
      rect_periodic(dx, dy, dz, bx);
      r2 = dx*dx + dy*dy + dz*dz;
      tp = field->bTypes[i];

      if ((btypes[tp].spec1 != atm->types[ia])||(btypes[tp].spec2 != atm->types[ja]))
        printf("!!!unkn bond type[%d]=%d between %s[%d] and %s[%d]\n", i, tp, sp[atm->types[ia]].name, ia, sp[atm->types[ja]].name, ja);

      if (btypes[tp].breakable)
        if (r2 > btypes[tp].rMax2)
        {
           nbr++;
           destroy_bond(i, atm, field, &btypes[tp]);
           destroy_angles(ia, ja, field);
           if (k < mxBr)
             {bInds[k] = i; k++;}
           continue;

        }

      f = bond_iter(r2, &btypes[tp], eng);
      //if (isnan(f))
        //{printf("bond[%d] f is nan\n", i);break;}
      //printf("f=%f\n", f);

      atm->fxs[ia] += f * dx;
      atm->fxs[ja] -= f * dx;
      atm->fys[ia] += f * dy;
      atm->fys[ja] -= f * dy;
      atm->fzs[ia] += f * dz;
      atm->fzs[ja] -= f * dz;

      //if (isnan(atm->fxs[ia]))
        //printf("fx[%d] f is nan\n", ia);
   }

   //remove breaking bonds
   if (nbr)
   {
      j = field->nBonds - 1;
      for (i = 0; i < k; i++)
      {
         while (field->at1[j] < 0)
           j--;
         tp = bInds[i];
         if (j <= tp)
            break;
         field->at1[tp] = field->at1[j];
         field->at2[tp] = field->at2[j];
         field->bTypes[tp] = field->bTypes[j];
         j--;
      }
      field->nBonds -= nbr;
   }

   delete[] bInds;
   //if (isnan(eng))
     // printf("bonds eng: NAN\n");
   sim->engBond = eng;
   //printf("eng=%f f=%f r2=%f  r0^2=%f\n", eng, f, r2, btypes[barrs->types[i]].p1*btypes[barrs->types[i]].p1);
}

/*
void bonding(Atoms *atm, Spec *sp, VdW ***vdws, Sim *sim, Box *bx, Bond *btps, Field *field)
// create new bonds during MD run
{
   int i, j, k;
   int ni, nj; // number of bonds
   int ktype;
   int iBond;
   int tai1, taj1, tai2, taj2; //type of atom i/j  before jump(1) after jump(2)
   double U1, U2, dU; // energies before(1) and after jump(2), delta U
   VdW *vdw;
   double dx, dy, dz, r2;
   int mx = 6; //! temp maximum number of bonds

   // brute force algorithm!
   //! вообще тут идЄт проверка каждой пары, а надо бы заранее запомнить соседей
   for (i = 0; i < atm->nAt; i++)
     if (sp[atm->types[i]].canBond)
     {
        tai1 = atm->types[i];
        for (j = 0; j < atm->nAt; j++)
        {
           // type[iatm] can change, so we need verify, that bonding ability is kept
           if (!sp[atm->types[i]].canBond)
             break;

           if (i == j)
             continue;

           if (!sp[tai1].bondKeys[atm->types[j]])
             continue;

           dx = atm->xs[i] - atm->xs[j];
           dy = atm->ys[i] - atm->ys[j];
           dz = atm->zs[i] - atm->zs[j];

           rect_periodic(dx, dy, dz, bx);
           r2 = dx * dx + dy * dy + dz * dz;
           if (r2 > sim->r2Bond)
             continue;


           //define specie's types:
           //    tai1 - before (non-bonding) i-atom
           //    atm->types[j] (non-bonding j-atom
           iBond = sp[tai1].bondKeys[atm->types[j]] - 1; // shifted, because 0 reserved for no-bonding
           if (tai1 == btps[iBond].spec1br)
             {
                tai2 = btps[iBond].spec1;
                taj2 = btps[iBond].spec2;
             }
           else
             {
                tai2 = btps[iBond].spec2;
                taj2 = btps[iBond].spec1;
             }
           taj1 = atm->types[j];
           //printf("tai1 2 j1 j2: %d %d %d %d bond[%d].p0=%f\n", tai1, tai2, taj1, taj2, iBond, btps[iBond].p0);


           // Energy difference calculation
           //! coulombic energy difference omitted
           U1 = 0.0; U2 = 0.0;
           // before bonding vdw energy
           vdw = vdws[tai1][taj1];
           if (vdw != NULL)
             if (r2 <= vdw->r2cut)
               vdw_iter(r2, vdw, U1);
           //after bonding bond energy
           bond_iter(r2, &btps[iBond], U2);
           for (k = 0; k < sim->nAt; k++)
             if (k != i)     // exclude i-i and i-j iteration (as i-j the same do and after jump)
               if (k != j)  // exclude j-j iteration
                 {
                    ktype = atm->types[k];

                    // I - K interaction:
                    dx = atm->xs[i] - atm->xs[k];
                    dy = atm->ys[i] - atm->ys[k];
                    dz = atm->zs[i] - atm->zs[k];
                    rect_periodic(dx, dy, dz, bx);
                    r2 = dx * dx + dy * dy + dz * dz;
                    //rad = sqrt(r2);

                    //vdw energy from i atom before bonding
                    vdw = vdws[tai1][ktype];
                    if (vdw != NULL)
                      if (r2 <= vdw->r2cut)
                        vdw_iter(r2, vdw, U1);

                    //vdw energy from i atom after bonding
                    vdw = vdws[tai2][ktype];
                      if (vdw != NULL)
                        if (r2 <= vdw->r2cut)
                          vdw_iter(r2, vdw, U2);

                    // J - K interaction:
                    dx = atm->xs[j] - atm->xs[k];
                    dy = atm->ys[j] - atm->ys[k];
                    dz = atm->zs[j] - atm->zs[k];
                    rect_periodic(dx, dy, dz, bx);
                    r2 = dx * dx + dy * dy + dz * dz;
                    //rad = sqrt(r2);

                    //vdw energy from j atom before jump
                    vdw = vdws[taj1][ktype];
                      if (vdw != NULL)
                        if (r2 <= vdw->r2cut)
                          vdw_iter(r2, vdw, U1);

                                 //vdw energy from j atom after jump
                    vdw = vdws[taj2][ktype];
                      if (vdw != NULL)
                        if (r2 <= vdw->r2cut)
                          vdw_iter(r2, vdw, U2);

                 } // end k-loop

                 dU = U2 - U1;
                 // own energy addition
                 dU += (sp[tai2].energy + sp[taj2].energy - sp[tai1].energy - sp[taj1].energy);

                 //BONDING!
                 //! here I add not only equality condition, but this energy decreasing too
                 if (dU <= 1e-6) //! temp parameter
                 {
                    ni = field -> nbonds[i];
                    nj = field -> nbonds[j];

                    if (ni == mx)
                    {
                       printf("WARNING: maximal numbers of bonds for atom[%d](%s) reached\n", i, sp[tai1].name);
                       continue; //! на самом деле тут нужно выйти не из j-цикла, а из i-
                    }

                    if (nj == mx)
                    {
                       printf("WARNING: maximal numbers of bonds for atom[%d](%s) reached\n", j, sp[taj1].name);
                       continue;
                    }


                    field -> bonds[i][ni] = iBond;
                    field -> bonds[j][nj] = iBond;
                    field -> neighs[i][ni] = j;
                    field -> neighs[j][nj] = i;
                    field -> nbonds[i]++;
                    field -> nbonds[j]++;
                    atm -> types[i] = tai2;
                    atm -> types[j] = taj2;

                    //printf("new bond created! %d(%s -> %s)-%d(%s -> %s)\n", i, sp[tai1].name, sp[tai2].name, j, sp[taj1].name, sp[taj2].name);

                 } // end bonding


        }    // end j-loop

     } // end main loop
}
*/

/*
void only_twoatomic_bonding(Atoms *atm, Spec *sp, VdW ***vdws, Sim *sim, Box *bx, Bond *btps, BondArrays *barr)
// create new bonds during MD run
{
   int i, j, k;
   int iat, jat, kat;
   int ktype;
   int iBond;
   int tai1, taj1, tai2, taj2; //type of atom i/j  before jump(1) after jump(2)
   double U1, U2, dU; // energies before(1) and after jump(2), delta U
   VdW *vdw;
   double dx, dy, dz, r2, rad;
   double kcharge;

   for (i = 0; i < sim->nCanBondAtm; i++)
   {
      iat = sim->canBondAtms[i];

      //verify that current atom are still can form bond:
      while ((!sp[atm->types[iat]].canBond) && (i < sim->nCanBondAtm))
      {
         sim->canBondAtms[i] = sim->canBondAtms[sim->nCanBondAtm - 1];
         iat = sim->canBondAtms[i];
         sim->nCanBondAtm--;
      }
      if (i >= sim->nCanBondAtm)
        break;

      tai1 = atm->types[iat];
      //cycle by neighbors:
      for (j = 0; j < sim->nNbors[iat]; j++)
        if ((sim->tnbors[iat][j] >> 3) & 1)  // neighbor in bonding distance
          {
             jat = sim->nbors[iat][j];
             // verify that neighbor j can form bond with iAtom:
             iBond = sp[tai1].bondKeys[atm->types[jat]];
             if (iBond)
             {
                //define indexes
                iBond--;
                taj1 = atm->types[jat];
                if (tai1 == btps[iBond].spec1br)
                {
                    tai2 = btps[iBond].spec1;
                    taj2 = btps[iBond].spec2;
                }
                else
                {
                    tai2 = btps[iBond].spec2;
                    taj2 = btps[iBond].spec1;
                }

                // energy difference calculation
                U1 = 0.0; U2 = 0.0;
                for (k = 0; k < sim->nNbors[iat]; k++)
                {
                    if (k == jat)
                      continue;

                    kat = sim->nbors[iat][k];
                    ktype = atm->types[kat];
                    kcharge = sp[ktype].charge;

                    //! поскольку ejump вызываетс€ после forcefield этот массив гарантировано заполнен не нул€ми
                    //!  или это не так? если все частицы зар€жены и рассто€ние прыжка меньше обрезани€ по кулону, то так
                    rad = sim->distances[iat][k];
                    r2 = rad * rad;

                    // VdW energy before bond formation:
                    vdw = vdws[tai1][ktype];
                    if (vdw != NULL)
                        if (r2 <= vdw->r2cut)  //! посмотреть, может опустить эту проверку, если он уж попал в соседи, то должен дот€гиватьс€ по ¬д¬
                            U1 += vdw->eng_r(r2, rad, vdw);

                    // Columbic energy before bond formation:
                    //! заменить на Ёвальда?
                    if (!sim->no_elec)
                      U1 += Fcoul_scale * kcharge * sp[tai1].charge / rad;

                    //vdw energy after bond formation:
                    vdw = vdws[tai2][ktype];
                    if (vdw != NULL)
                        if (r2 <= vdw->r2cut)
                            U2 += vdw->eng_r(r2, rad, vdw);
                    // Columbic energy after jump
                    //! заменить на Ёвальда?
                    if (!sim->no_elec)
                      U2 += Fcoul_scale * kcharge * sp[tai2].charge / rad;
                } // end loop by iat neigbors

                // loop by jat neighbors:
                for (k = 0; k < sim->nNbors[jat]; k++)
                {
                    if (k == iat)
                        continue;

                    kat = sim->nbors[jat][k];
                    ktype = atm->types[kat];
                    kcharge = sp[ktype].charge;

                    //! поскольку ejump вызываетс€ после forcefield этот массив гарантировано заполнен не нул€ми
                    //!  или это не так? если все частицы зар€жены и рассто€ние прыжка меньше обрезани€ по кулону, то так
                    rad = sim->distances[jat][k];
                    r2 = rad * rad;

                    // VdW energy before bond formation:
                    vdw = vdws[taj1][ktype];
                    if (vdw != NULL)
                        if (r2 <= vdw->r2cut)  //! посмотреть, может опустить эту проверку, если он уж попал в соседи, то должен дот€гиватьс€ по ¬д¬
                            U1 += vdw->eng_r(r2, rad, vdw);

                    // Columbic energy before bond formation:
                    //!     заменить на Ёвальда?
                    if (!sim->no_elec)
                      U1 += Fcoul_scale * kcharge * sp[taj1].charge / rad;

                    //vdw energy after bond formation:
                    vdw = vdws[taj2][ktype];
                        if (vdw != NULL)
                            if (r2 <= vdw->r2cut)
                                U2 += vdw->eng_r(r2, rad, vdw);
                    // Columbic energy after bond formation
                    //! заменить на Ёвальда?
                    if (!sim->no_elec)
                      U2 += Fcoul_scale * kcharge * sp[taj2].charge / rad;
                }  // end loop by jat neigbors

                dx = atm->xs[iat] - atm->xs[jat];
                dy = atm->ys[iat] - atm->ys[jat];
                dz = atm->zs[iat] - atm->zs[jat];
                rect_periodic(dx, dy, dz, bx);
                r2 = dx * dx + dy * dy + dz * dz;
                rad = sim->distances[iat][j];

                // VdW to Bond
                vdw = vdws[taj1][taj2];
                if (vdw != NULL)
                  U1 += vdw->eng_r(r2, rad, vdw);
                bond_iter_r(r2, rad, &btps[iBond], U2);

                dU = U2 - U1;
                dU += btps[iBond].energy;

                // bond forming condition
                if (dU <= 1e-6) //! temp
                {
                    sim->nBndForm++;    // inc the number of crated bonds
                    sim->nBonds++;
                    sim->degFree--;
                    sim->revDegFree = 1.0 / sim->degFree; //! надо бы это вынести за цикл, ну какой смысл это каждую разрушенную св€зь перевычисл€ть

                    sim->bondAtoms[iat] = jat;
                    sim->bondAtoms[jat] = iat;
                    sim->bondTypes[iat] = iBond;
                    sim->bondTypes[jat] = iBond;

                    atm->types[iat] = tai2;
                    atm->types[jat] = taj2;


                    //! наверное можно пропустить изменение массива canBondAtmBond поскольку здесь выполн€ютс€ проверки и чистки

                }

             }
             //end verify that pair iat-jat can form bond
          }
      //end cycle by neighbors
   } // end main loop
}
*/

// free arrays for bond procedures
int free_bonds(Sim *sim, int **nbonds, int ***bonds, int ***neigs/*, Box *box*/, int nAtm)
{
   int i;
   int **b, **n;

   b = *bonds;
   n = *neigs;
   for (i = 0; i < nAtm; i++)
     {
        delete b[i];
        delete n[i];
     }

   delete[] *nbonds;
   delete[] **bonds;
   delete[] **neigs;

   delete[] sim->canBondAtms;
   if (sim->only_twoatomic)
   {
       delete[] sim->bondAtoms;
   }



   return 1;
}

