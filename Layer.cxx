//Class used in order to define the different layers in our detector system

#include <Riostream.h>
#include "Layer.h"

ClassImp(Layer)

Layer::Layer():TObject(),
fWidth(0.),
fRadius(0.),
fThick(0.),
fRMS(0.)
{
  //default constructor
}

Layer::Layer(double W, double R, double T, double RMS): TObject(),
fWidth(W),
fRadius(R),
fThick(T),
fRMS(RMS)
{
  //standard constructor, W R T RMS are given by the user
}

Layer::~Layer() {
  //default destructor
}

Double_t Layer::GetWidth() const {
  //returns Width
  return fWidth;
}

Double_t Layer::GetRadius() const {
  //returns Radius
  return fRadius;
}

Double_t Layer::GetRMS() const {
  //returns Thickness
  return fRMS;
}
