//solvers.h

#ifndef solversH
#define solversH


void SOR(double ***P, double ***B, double ***Ap, double ***Ae, double ***Aw, double ***An, double ***As, double ***Af, double ***Ab, double X, double Y, double Z, double h, int nx, int ny, int nz, double w, int itr, int sit, char boundaryType);

void MG(double ***Pex, double ***E, double ****e, double ****R, double ****rhs, double ***P, double ***B, double ****Ap, double ****Ae, double ****Aw, double ****An, double ****As, double ****Af, double ****Ab, int *upTimes, int *downTimes, double X, double Y, double Z, double *h, int *nx, int *ny, int *nz, int levels, double *w, double eps, char boundaryType);

void JACOBI(double ***P, double ***B, double ***Ap, int nx, int ny, int nz);

void CONJUGATE_GRADIENT(double ***Pex, double ***P, double ***B, double ***E, double ***z1, double ***z2, double ***R1, double ***R2, double ***p, double ***W, double ****Ap, double ****Ae, double ****Aw, double ****An, double ****As, double ****Af, double ****Ab, double X, double Y, double Z, double *h, int *nx, int *ny, int *nz, double *w, double eps, int levels, char boundaryType, char preconditionType);


#endif

