//The class inherits from TObject 

#ifndef HIT_H
#define HIT_H

#include "TObject.h"

class Hit : public TObject{

public:
    
    Hit();//default constructor
    Hit(double x, double y, double z);//custom contructor
    virtual ~Hit();//destructor
    
    //member function
    
    double GetX() const; //getters
    double GetY() const;
    double GetZ() const;
    
private:
    
    //data member
    
    double fX; //values of intersection 
    double fY;
    double fZ;
    
    ClassDef(Hit,1)
    
};

#endif
