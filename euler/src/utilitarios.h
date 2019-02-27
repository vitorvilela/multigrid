//utilitarios.h

#ifndef utilitariosH
#define utilitariosH

void erro(double ***Pex, double ***P, double ***E, int nx, int ny, int nz);

void residuo(double ***R, double ***P, double ***rhs, double ***Ap, double ***Ae, double ***Aw, double ***An, double ***As, double ***Au, double ***Ad, int nx, int ny, int nz);

double maximo(double ***R, int nx, int ny, int nz);

void zero(double ****e, int *nx, int *ny, int *nz, int niveis);

void corrige(double ***P, double ***e, int nx, int ny, int nz);

void normaliza(double ***P, int nx, int ny, int nz, double h);

void igual(double ***rhs, double ***R, int nx, int ny, int nz);

double escalar(double ***M1, double ***M2, int nx, int ny, int nz);

void interpola3d(double ****e, int *nx, int *ny, int *nz, int t);

void restringe3d(double ***R, double ****rhs, int *nx, int *ny, int *nz, int t);

void campo(int nx, int ny, int nz, double h, double ***U, double ***V, double ***W, double ***P, int tempo);



#endif


