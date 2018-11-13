//The class inherits from TObject

#ifndef HIT_H
#define HIT_H

#include "TObject.h"
#include <TRandom3.h>
#include <TH1F.h>

class Hit : public TObject{

public:

    Hit();//default constructor
    Hit(double x, double y, double z);//custom contructor
    Hit(double meanv, double sigmaxy, double sigmaz, TH1F *distr_mult);//event constructor
    virtual ~Hit();//destructor

    //member function

    double GetX() const; //getters
    double GetY() const;
    double GetZ() const;
    int GetMult() const;

private:

    //data member

    double fX; //values of intersection
    double fY;
    double fZ;
    int fMult;

    ClassDef(Hit,2)

};

#endif
