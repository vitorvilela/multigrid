/*
Variable.h
*/

#ifndef Variable_H
#define Variable_H

class Variable {
public:
   Variable(char *VM = "VAR", char *VT = "CENTERED");
   void setName(char *VN);
   void setType(char *VT);
   const char *getName() const { return NAME; }
   const char *getType() const { return TYPE; }
private:
   char *NAME;
   char *TYPE;
};

#endif
