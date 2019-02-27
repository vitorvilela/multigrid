/*
Domain.cpp
*/

#include "Domain.h"

Domain::Domain(double W, double E)
{
   setWest(W);
   setEast(E);
}

void Domain::setWest(double W) { ( W<EAST ? WEST=W : WEST=EAST=0 ); }
   
void Domain::setEast(double E) { ( E>WEST ? EAST=E : WEST=EAST=0 ); }


