//SISTEMA


#ifndef sistemaH
#define sistemaH

void coefficients(char **flag, double ***Ap, double ***Ae, double ***Aw, double ***An, double ***As, double **B, double ***Cp, double ***K, double **T0, double **Cp0, double **f, double **f0, int *nx, int *ny, double *h, double ro, double L, double dt, int levels, double Tw, double Te);

void residuo(double **R, double **T, double **rhs, double **Ap, double **Ae, double **Aw, double **An, double **As, int nx, int ny);

#endif


