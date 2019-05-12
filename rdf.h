// rdf.h
//  заголовочный файл модуля, работающего с радиальной функцией распределения
#ifndef RDF_H
#define RDF_H

int init_rdf(double ***rdf, Sim *sim, Field *field, Box *box, int &nRDF);
void clear_rdf(double **rdf, Sim *sim, Field *field, int &nRDF);
int get_rdf(Atoms *atm, double **rdf, Sim *sim, Field *field, Box *box, int &nRDF);
int out_rdf(double **rdf, Field *field, Box *box, Sim *sim, char *fname, int nRDF);
int free_rdf(double ***rdf, int nPair);

#endif  /* RDF_H */
