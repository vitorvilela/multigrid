//exact.h

#ifndef exactH
#define exactH

void EXACT(double ****Pex, double *h, int *nx, int *ny, int *nz, int levels, char imp);

void SET_BOUNDARY(double ***P, double X, double Y, double Z, double h, int nx, int ny, int nz, char boundaryType, int sit);

void DIRICHLET(double ***P, double X, double Y, double Z, double h, int nx, int ny, int nz, int sit);

void NEUMANN(double ***P, double X, double Y, double Z, double h, int nx, int ny, int nz, int sit);

#endif



