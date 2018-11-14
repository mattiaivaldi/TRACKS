//Class inherits from TObject

#ifndef LAYER_H
#define LAYER_H

#include "TObject.h"

class Layer : public TObject{

public:

    Layer();//default constructor
    Layer(double W, double R, double T, double RMS);//custom contructor
    virtual ~Layer();//destructor

    //member function

    double GetWidth() const;
    double GetRadius() const;
    double GetThick() const;
    double GetRMS() const;

private:

    //data member

    double fWidth;
    double fRadius;
    double fThick;
    double fRMS;

    ClassDef(Layer,1)

};

#endif
