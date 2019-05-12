#include <math.h>   // log, sqrt
#include <stdio.h>   //!временно, для вывода тестов

#include "dataStruct.h"  // Sim, Box, Atoms ....
#include "box.h"

void box_prop(Box *bx)
// refresh box properties according to xx, xy, xz, ...
{
   double axb1, axb2, axb3;
   double bxc1, bxc2, bxc3;
   double cxa1, cxa2, cxa3;

   // the same for invert matrix
   double iaxb1, iaxb2, iaxb3;
   double ibxc1, ibxc2, ibxc3;
   double icxa1, icxa2, icxa3;

   //double b1, b2, b3, b4, b5, b6, b7, b8, b9;
   double rdet;

   //! temp (for rectangular geometry only:)
   bx->ay = 0.0; bx->az = 0.0; bx->bx = 0.0; bx->bz = 0.0; bx->cx = 0.0; bx->cy = 0.0;

    //bbb(1)-(3) in DLPOLY  a,b,c ? length of box vectors:
   bx->la = sqrt(bx->ax*bx->ax + bx->ay*bx->ay + bx->az*bx->az);
   bx->lb = sqrt(bx->bx*bx->bx + bx->by*bx->by + bx->bz*bx->bz);
   bx->lc = sqrt(bx->cx*bx->cx + bx->cy*bx->cy + bx->cz*bx->cz);
   //printf("box prop, la lb lc: %f %f %f\n", bx->la, bx->lb, bx->lc);


   bx->maxLength = bx->la;
   if (bx->maxLength < bx->lb)
     bx->maxLength = bx->lb;
   if (bx->maxLength < bx->lc)
     bx->maxLength = bx->lc;

   // cosines of cell angles
   bx->cosC = (bx->ax * bx->bx + bx->ay * bx->by + bx->az * bx->bz)/(bx->la * bx->lb); // bbb(4)
   bx->cosB = (bx->ax * bx->cx + bx->ay * bx->cy + bx->az * bx->cz)/(bx->la * bx->lc); // bbb(5)
   bx->cosA = (bx->bx * bx->cx + bx->by * bx->cy + bx->bz * bx->cz)/(bx->lb * bx->lc); // bbb(6)

   //   calculate vector products of cell vectors
   axb1 = bx->ay * bx->bz - bx->az * bx->by;
   axb2 = bx->az * bx->bx - bx->ax * bx->bz;
   axb3 = bx->ax * bx->by - bx->ay * bx->bx;

   bxc1 = bx->by * bx->cz - bx->bz * bx->cy;
   bxc2 = bx->bz * bx->cx - bx->bx * bx->cz;
   bxc3 = bx->bx * bx->cy - bx->by * bx->cx;

   cxa1 = bx->az * bx->cy - bx->ay * bx->cz;
   cxa2 = bx->ax * bx->cz - bx->az * bx->cx;
   cxa3 = bx->ay * bx->cx - bx->ax * bx->cy;

   // invert matrix calculation:
   //    1. calculate adjoint matrix
   /*
   b1 = box->by * box->cz - box->bz * box->cy;  // bxc1
   b2 = box->az * box->cy - box->ay * box->cz;  // cxa1
   b3 = ay * bz - az * by;                      // axb1

   b4 = bz * cx - bx * cz;                      // bxc2
   b5 = ax * cz - az * cx;                      // cxa2
   b6 = az * bx - ax * bz;                      // axb2

   b7 = bx * cy - by * cx;                      // bxc3
   b8 = ay * cx - ax * cy;                      // cxa3
   b9 = ax * by - ay * bx;                      // axb3
   */

   //   2. calculate determinant:
   //! надо посмотреть, а нужно ли его хранить, det. я подозреваю, что det это и есть объём матрицы
   bx->det = bx->ax * bxc1 + bx->bx * cxa1 + bx->cx * axb1;
   rdet = 0.0;
   if (fabs(bx->det) > 0.0)
     rdet = 1.0 / bx->det;

   //   3. complete inverse matrix
   bx->iax = rdet * bxc1;
   bx->iay = rdet * cxa1;
   bx->iaz = rdet * axb1;

   bx->ibx = rdet * bxc2;
   bx->iby = rdet * cxa2;
   bx->ibz = rdet * axb2;

   bx->icx = rdet * bxc3;
   bx->icy = rdet * cxa3;
   bx->icz = rdet * axb3;

   //   calculate vector products of cell vectors for the invert matrix
   iaxb1 = bx->iay * bx->ibz - bx->iaz * bx->iby;
   iaxb2 = bx->iaz * bx->ibx - bx->iax * bx->ibz;
   iaxb3 = bx->iax * bx->iby - bx->iay * bx->ibx;

   ibxc1 = bx->iby * bx->icz - bx->ibz * bx->icy;
   ibxc2 = bx->ibz * bx->icx - bx->ibx * bx->icz;
   ibxc3 = bx->ibx * bx->icy - bx->iby * bx->icx;

   icxa1 = bx->iaz * bx->icy - bx->iay * bx->icz;
   icxa2 = bx->iax * bx->icz - bx->iaz * bx->icx;
   icxa3 = bx->iay * bx->icx - bx->iax * bx->icy;

   // volume and invert volume:
   bx->vol = fabs(bx->det);
   bx->rvol = 1.0 / bx->vol;
   //printf("box prop, vol rvol: %f %f\n", bx->vol, bx->rvol);


   // нужно эту хрень выполнить для обратной матрицы
   bx->ip1 = bx->rvol / sqrt(ibxc1 * ibxc1 + ibxc2 * ibxc2 + ibxc3 * ibxc3);
   bx->ip2 = bx->rvol / sqrt(icxa1 * icxa1 + icxa2 * icxa2 + icxa3 * icxa3);
   bx->ip3 = bx->rvol / sqrt(iaxb1 * iaxb1 + iaxb2 * iaxb2 + iaxb3 * iaxb3);

   // OTHER properties:
   //   half of vector length:
   bx->ha = bx->la * 0.5;  // half a
   bx->hb = bx->lb * 0.5;
   bx->hc = bx->lc * 0.5;

   //   negative of half vector length:
   bx->nha = -bx->ha;  // -half a
   bx->nhb = -bx->hb;
   bx->nhc = -bx->hc;

   //  invert vector length
   bx->ra = 1.0 / bx->la;  // reversible a
   bx->rb = 1.0 / bx->lb;
   bx->rc = 1.0 / bx->lc;
   bx->ra2 = bx->ra * bx->ra;  // 1/a^2
   bx->rb2 = bx->rb * bx->rb;
   bx->rc2 = bx->rc * bx->rc;

   bx->momXp = 0.0;
   bx->momXn = 0.0;
   bx->momYp = 0.0;
   bx->momYn = 0.0;
   bx->momZp = 0.0;
   bx->momZn = 0.0;
   //printf("box prop: ip1, ip2, ip3: %f %f %f\n", bx->ip1, bx->ip2, bx->ip3);
}
// end 'box_prop' function

void rect_periodic(double &dx, double &dy, double &dz, Box *box)
// apply periodic boundary to dx, dy, dz. Only for rectangular geometry!
//! please correct 'rect_periodic_ejump' function if correct this funciton (they are similar)
{
   //double x0 = dx;

   // x
   if (dx > box->ha)
     dx -= box->la;
   else
     if (dx < box->nha)
       dx += box->la;

   //if ((!isnormal(dx))&&(fabs(dx) > 1e6))
   //  printf("rect_periodic: input dx=%f out=%f  box[ha,la]=%d, %d\n", x0, dx, box->ha, box->la);

   // y
   if (dy > box->hb)
     dy -= box->lb;
   else
     if (dy < box->nhb)
       dy += box->lb;

   // z
   if (dz > box->hc)
     dz -= box->lc;
   else
     if (dz < box->nhc)
       dz += box->lc;

   /*
   // x-coordinate:
   if (*dx > box->ha)
     *dx -= box->la;
   else
     if (*dx < -box->ha)
       *dx += box->la;

   // y
   if (*dy > box->hb)
     *dy -= box->lb;
   else
     if (*dy < -box->hb)
       *dy += box->lb;

   // z
   if (*dz > box->hc)
     *dz -= box->lc;
   else
     if (*dz < -box->hc)
       *dz += box->lc;
   */
}
// end 'rect_periodic' function

void pass_rect_box(double dx, double dy, double dz, Box *box, int &px, int &py, int &pz)
// put in 'px', 'py' and 'pz' vars values -1, 0 or 1 depends on how particles with dx dy and dz are situated in box
{
   if (dx > box->ha) // второй атом в отрицательном отображении
     {
       px = -1;
     }
   else
     if (dx < box->nha) //  второй атом в положительном отображении
       {
          px = 1;
       }
     else
       px = 0;

   //! add y- and z- directions
}

// verify that ejump is throw the box
void rect_periodic_ejump(double dx, double dy, double dz, Box *box, ThrowBox *eThrow)
{
   if (dx > box->ha)
     {
        eThrow->nx++;
     }
   else
     if (dx < box->nha)
       {
          eThrow->px++;
       }

   //! add y- and z- dimensions
   // y
/*
   if (dy > box->hb)
     dy -= box->lb;
   else
     if (dy < box->nhb)
       dy += box->lb;

   // z
   if (dz > box->hc)
     dz -= box->lc;
   else
     if (dz < box->nhc)
       dz += box->lc;
   */
}
// end 'rect_periodic_ejump' function

int rect_box(Atoms *atm, int index, Spec *sp, Box *box)
// put atom[index] in periodic box. Only for rectangular geometry!
{
   float x0;
   int t = atm->types[index];

   //! важно, импульс тут без множителя 2, который видимо надо потом поставить


   if (atm->xs[index] < 0)
     {
        x0 = atm->xs[index];


        if (x0 < -box->la)
        {
            printf("x[%d] < -size of box (%f)\n", index, x0);
            return 0;
        }


        //! здесь и далее не учтено, что частица может улететь дальше чем на период бокса
        //atm->xs[index] += box->la;
        //! теперь учтено,но непонятно, откуда берутся такие скорости, которые вышвыривают частицу далеко за бокс
        atm->xs[index] += ((int)(-atm->xs[index]*box->ra) + 1) * box->la;


        //if (atm->xs[index] < 0)
        //  printf("a paricle[%d] has still negative X(%f), prevX=%f Vx=%f\n", index, atm->xs[index], x0, atm->vxs[index]);
        sp[t].nOyz++;
        box->momXn += sp[t].mass * (-atm->vxs[index]); // i suppose that vx in this case is negative
     }
   else
     if (atm->xs[index] > box->la)
       {
          //atm->xs[index] -= box->la;
          //! теперь учтено,но непонятно, откуда берутся такие скорости, которые вышвыривают частицу далеко за бокс
          atm->xs[index] -= ((int)(atm->xs[index]*box->ra)) * box->la;

          sp[t].pOyz++;
          box->momXp += sp[t].mass * atm->vxs[index];
       }


   //! написать аналогичные счетчики для других направлений
   if (atm->ys[index] < 0)
     {
        //atm->ys[index] += box->lb;
        //! теперь учтено,но непонятно, откуда берутся такие скорости, которые вышвыривают частицу далеко за бокс
        atm->ys[index] += ((int)(-atm->ys[index]*box->rb) + 1) * box->lb;

        sp[t].nOxz++;
        box->momYn += sp[t].mass * (-atm->vys[index]); // i suppose that vy in this case is negative
     }
   else
     if (atm->ys[index] > box->lb)
       {
            //atm->ys[index] -= box -> lb;
            //! теперь учтено,но непонятно, откуда берутся такие скорости, которые вышвыривают частицу далеко за бокс
            atm->ys[index] -= ((int)(atm->ys[index]*box->rb)) * box->lb;

            sp[t].pOxz++;
            box->momYp += sp[t].mass * atm->vys[index];
       }

   if (atm->zs[index] < 0)
     {
        //atm->zs[index] += box -> lc;
        //! теперь учтено,но непонятно, откуда берутся такие скорости, которые вышвыривают частицу далеко за бокс
        atm->zs[index] += ((int)(-atm->zs[index]*box->rc) + 1) * box->lc;

        sp[t].nOxy++;
        box->momZn += sp[t].mass * (-atm->vzs[index]); // i suppose that vy in this case is negative
     }
   else
     if (atm->zs[index] > box -> lc)
       {
            //atm->zs[index] -= box -> lc;
            //! теперь учтено,но непонятно, откуда берутся такие скорости, которые вышвыривают частицу далеко за бокс
            atm->zs[index] -= ((int)(atm->zs[index]*box->rc)) * box->lc;

            sp[t].pOxy++;
            box->momZp += sp[t].mass * atm->vzs[index];
       }

   return 1;
}
// end 'rect_box' function
