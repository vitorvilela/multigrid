//UTILITARIOS


#ifndef utilitariosH
#define utilitariosH


double norma(double **P, int nx, int ny, double h);

void normaliza(double **P, int nx, int ny, double h);

double maximo(double **R, int nx, int ny);

void zero(double ***e, int *nx, int *ny, int niveis);

void igual(double **rhs, double **R, int nx, int ny);

void igual3D(double ***M1, double ***M2, int *nx, int *ny, int niveis);

double escalar(double **M1, double **M2, int nx, int ny);

void corrigeMAIS(double **P, double **e, int nx, int ny);

void corrigeMENOS(double **P, double **e, int nx, int ny);

void interpola(double ***e, int *nx, int *ny, int t);

void restringe(double **R, double ***rhs, int *nx, int *ny, int t);   



#endif
