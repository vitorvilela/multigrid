//system.h

#ifndef systenH
#define systemH

void DENSITY_FIELD(double ****Ro, double Xo, double Yo, double Zo, double *h, int *nx, int *ny, int *nz, double d, double ro, double W0, double W1, int levels, char imp);

void COEFFICIENTS(double ****Ap, double ****Aw, double ****Ae, double ****An, double ****As, double ****Af, double ****Ab, double ****Ro, double *h, int *nx, int *ny, int *nz, int levels);

void FORCE_TERM(double ****B, double ****Pex, double ****Ap, double ****Aw, double ****Ae, double ****An, double ****As, double ****Au, double ****Ad, int *nx, int *ny, int *nz, int levels);

void RESIDUE(double ***R, double ***P, double ***rhs, double ***Ap, double ***Ae, double ***Aw, double ***An, double ***As, double ***Af, double ***Ab, int nx, int ny, int nz);


#endif

