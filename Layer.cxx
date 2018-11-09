//Class used in order to define the different layers in our detector system

#include <Riostream.h>
#include "Layer.h"

ClassImp(Layer)

//Initialize all data memebers to 0
Layer::Layer():TObject(),
fWidth(0.),
fRadius(0.),
fThick(0.),
fRMS(0.)
{
    //Default constructor
}

//Initialize data members to values provided in the main program
Layer::Layer(double W, double R, double T, double RMS): TObject(),
fWidth(W),
fRadius(R),
fThick(T),
fRMS(RMS)
{
    //Standard constructor
}

Layer::~Layer() {
    //Default destructor
}

Double_t Layer::GetWidth() const {
    //Returns Width
    return fWidth;
}

Double_t Layer::GetRadius() const {
    //Returns Radius
    return fRadius;
}

Double_t Layer::GetRMS() const {
    //Returns Thickness
    return fRMS;
}
