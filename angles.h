#ifndef ANGLES_H
#define ANGLES_H

void read_anglelist(Atoms *atm, Field *field);
void save_anglelist(char *fname, Field *field);
int read_angle(int id, FILE *f, Field *field);
void destroy_angles(int a1, int a2, Field *field);
int create_angle(int c, int l1, int l2, int type, Field *field);
void exec_anglelist(Atoms *atm, Field *field, Sim *sim, Box *bx);
void free_angles(Field *field);

#endif // ANGLES_H
