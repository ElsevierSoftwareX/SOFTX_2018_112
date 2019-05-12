#ifndef BOX_H
#define BOX_H

void box_prop(Box *box);
void rect_periodic(double &dx, double &dy, double &dz, Box *box);
void pass_rect_box(double dx, double dy, double dz, Box *box, int &px, int &py, int &pz);
void rect_periodic_ejump(double dx, double dy, double dz, Box *box, ThrowBox *eThrow);

// put atom[index] in periodic box. Only for rectangular geometry!
int rect_box(Atoms *atm, int index, Spec *sp, Box *box);

#endif  /* BOX_H */
