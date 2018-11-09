//Class used in order to define the hit coorinates of the particle with both the beam pipe and the two detector layers 

#include <Riostream.h>
#include "Hit.h"

ClassImp(Hit)

//Initialize all data members to "0"
Hit::Hit():TObject(),
fX(0.),
fY(0.),
fZ(0.)
{
    //Default constructor
}

//Initialize data members to the value of the intersection point 
Hit::Hit(double x, double y, double z): TObject(),
fX(x),
fY(y),
fZ(z)
{
    //Standard constructor
}

Hit::~Hit() {
    //Default destructor
}

Double_t Hit::GetX() const {
    //Returns X
    return fX;
}

Double_t Hit::GetY() const {
    //Returns Y
    return fY;
}

Double_t Hit::GetZ() const {
    //Returns Z
    return fZ;
}
