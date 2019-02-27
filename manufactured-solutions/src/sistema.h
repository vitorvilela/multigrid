//sistema.h

#ifndef sistemaH
#define sistemaH

void ME(double ****Ro, int *nx, int *ny, int *nz, double *h, int K, char imp, int niveis);

void coeficientes(double ****Ap, double ****Aw, double ****Ae, double ****An, double ****As, double ****Au, double ****Ad, double ****Ro, double *h, int *nx, int *ny, int *nz, int niveis);

void f(double ***B, double ***Pex, double ***Ap, double ***Aw, double ***Ae, double ***An, double ***As, double ***Au, double ***Ad, int nx, int ny, int nz);


#endif

