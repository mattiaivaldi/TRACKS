//Class used in order to generate the event and its related multiplicity 

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
Event::Event(double meanv, double sigmaxy, double sigmaz, const char *distr): Hit(),
fX0(gRandom->Gaus(meanv,sigmaxy)), //generates random values according to given data
fY0(gRandom->Gaus(meanv,sigmaxy)),
fZ0(gRandom->Gaus(meanv,sigmaz)),
fMult(0)
									      
{
  //Standard constructor
  TFile hfile(distr); //opens file
  TH1F *multiplicity = (TH1F*)hfile.Get("hmul"); //assign to TH1F multiplicity the histogram in the .root file 
  fMult=multiplicity->GetRandom(); //get a random value from assigned distribution 
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

