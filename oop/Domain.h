/*
Domain.h
*/

#ifndef Domain_H
#define Domain_H

class Domain {
public:
   Domain(double W = 0.0, double E = 1.0);
   void setWest(double);
   void setEast(double);
   double getWest() const { return WEST; }
   double getEast() const { return EAST; }   
private:
   double WEST;
   double EAST;   
};

#endif

