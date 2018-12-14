//The class inherits from TObject

#ifndef HIT_H
#define HIT_H

#include "TObject.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "Layer.h"

class Hit : public TObject{

public:

    Hit(); //default constructor, def1
    Hit(double x, double y, double z); //hit contructor, def1
    Hit(double R, double H); //spurious hit constructor (noise), def1
    Hit(double meanv, double sigmaxy, double sigmaz, TH1F *distr_mult); //event constructor, def2
    Hit(double meanv, double sigmaxy, double z_custom, int mult_custom); //custom event, def2
    virtual ~Hit(); //destructor

    //member function

    double GetX() const;
    double GetY() const;
    double GetZ() const;
    int GetMult() const; //def2
    void SetX(double xnew){fX=xnew;}
    void SetY(double ynew){fY=ynew;}
    void SetZ(double znew){fZ=znew;}
    void SetMult(int multnew){fMult=multnew;} //def2
    void Customize(int custom, double z_custom, int mult_custom); //def2

private:

    //data member

    double fX;
    double fY;
    double fZ;
    int fMult; //def2

    ClassDef(Hit,2)

};

#endif
