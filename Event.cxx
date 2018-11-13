//Class used in order to generate the event and its related multiplicity
#include <TRandom2.h>
#include <TRandom3.h>
#include <Riostream.h>
#include <TFile.h>
#include <TH1F.h>
#include "Hit.h"
#include "Event.h"

ClassImp(Event)

//Initialize all data members to 0
Event::Event():Hit(),
fX0(0.),
fY0(0.),
fZ0(0.),
fMult(0.)
{
    //Default constructor
}

//Generate event coordinates and multiplicity
Event::Event(double meanv, double sigmaxy, double sigmaz, TH1F *distr_mult): Hit(),
fX0(gRandom->Gaus(meanv,sigmaxy)), //generates random values according to given data
fY0(gRandom->Gaus(meanv,sigmaxy)),
fZ0(gRandom->Gaus(meanv,sigmaz)),
fMult(fMult=distr_mult->GetRandom())
{
  //Standard constructor
}

Event::~Event() {
    // Default destructor
}

Double_t Event::GetX() const {
    //Returns X0
    return fX0;
}

Double_t Event::GetY() const {
    //Returns Y0
    return fY0;
}

Double_t Event::GetZ() const {
    //Returns Z0
    return fZ0;
}


Float_t Event::GetMult() const {
    //Returns multiplicity
    return fMult;
}
