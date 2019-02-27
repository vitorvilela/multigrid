//sistema.h

#ifndef sistemaH
#define sistemaH

void ME1(double ****Ro, int *nx, int *ny, int *nz, double *h, int K, char imp, int niveis);

void ME2(double ****Ro, int *nx, int *ny, int *nz, double *h, int K, char imp, int niveis);

void ME3(double ****Ro, int *nx, int *ny, int *nz, double *h, int K, char imp, int niveis);

void ME4(double ****Ro, int *nx, int *ny, int *nz, double *h, int K, char imp, int niveis);

void ME5(double ****Ro, int *nx, int *ny, int *nz, double *h, int K, char imp, int niveis);

void coeficientes(double ****Ap, double ****Aw, double ****Ae, double ****An, double ****As, double ****Au, double ****Ad, double ****Ro, double *h, int *nx, int *ny, int *nz, int niveis);

void precondicionador(double ****Mp, double ****Mw, double ****Me, double ****Mn, double ****Ms, double ****Mu, double ****Md, double ****Ro, double *h, int *nx, int *ny, int *nz, int niveis, char imp);


#endif

