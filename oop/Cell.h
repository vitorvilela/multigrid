/*
Cell.h
*/

#ifndef Cell_H
#define Cell_H

class Cell {
public:
   Cell(double, double, double, double);
   void setValue(double V) { VALUE = V; }
   void setWface(double WF) { WFACE = WF; }
   void setCenter(double C) { CENTER = C; }
   void setEface(double EF) { EFACE = EF; }
   double getValue() const { return VALUE; }
   double getWface() const { return WFACE; }
   double getCenter() const { return CENTER; }
   double getEface() const { return EFACE; }
private:  
   double WFACE;
   double CENTER;
   double EFACE;
   double DX;
   double VALUE;
   void setDx() { DX = EFACE-WFACE; }
        
   Cell *nextPtr;
   Cell *lastPtr;
};

#endif
