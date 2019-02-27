#ifndef impressaoH
#define impressaoH

void perfil_cavidade(double **U, double **V, int nx, int ny, double h, FILE *U_y, FILE *V_x);

void campo(int nx, int ny, double h, double **Ro, double **K, double **T, double **f, int t);

#endif

