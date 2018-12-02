//Class used in order to define the different layers in our detector system

#include <Riostream.h>
#include "Layer.h"
#include "TString.h"

ClassImp(Layer)

Layer::Layer():TObject(),
fName("null"),
fWidth(0.),
fRadius(0.),
fThick(0.),
fRMS(0.)
{
  //default constructor
}

Layer::Layer(TString N, double W, double R, double T, double RMS): TObject(),
fName(N),
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

TString Layer::GetLayerName() const {
  //returns name
  return fName;
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
