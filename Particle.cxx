//Class used for object of type particle

#include <TRandom3.h>
#include <Riostream.h>
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include "Particle.h"

using namespace TMath;

ClassImp(Particle)

//Initialize all data members to 0
Particle::Particle():TObject(),
  fTheta(0.),
  fPhi(0.),
  fRape(0.)
{
  //Default constructor
}

//Initialize to values given by the user
Particle::Particle(const char *distr): TObject(),
fTheta(0),
fPhi(gRandom->Uniform(2*Pi())),
fRape(0)
{
  //Standard constructor
  TFile hfile(distr); //opens file
  TH1F *pseudorape = (TH1F*)hfile.Get("heta"); //assign to TH1F pseudorape the histogram in the .root file
  fRape = pseudorape->GetRandom(); //get a random value from assigned distribution
  fTheta = 2*ATan(Exp(-(double)fRape)); //get theta from the pseudorapidity distribution
}

Particle::~Particle() {
  // Default destructor
}


Double_t Particle::GetTheta() const {
  //Returns theta
  return fTheta;
}

Double_t Particle::GetPhi() const {
  //Returns phi
  return fPhi;
}

Float_t Particle::GetRaped() const {
  //Returns pseudorapidity
  return fRape;
}
