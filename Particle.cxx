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

void Particle::Rotate(double rms){

  gRandom->SetSeed(0);

  double theta0=rms/Sqrt(2);
  double thetap=gRandom->Gaus(0,theta0);
  double phip=gRandom->Uniform(2*Pi());

  double mr[3][3], pol[3], rot[3], r;

  mr[0][0]=-Sin(fPhi);
  mr[1][0]=Cos(fPhi);
  mr[2][0]=0;
  mr[0][1]=-Cos(fPhi)*Cos(fTheta);
  mr[1][1]=-Cos(fTheta)*Sin(fPhi);
  mr[2][1]=Sin(fTheta);
  mr[0][2]=Sin(fTheta)*Cos(fPhi);
  mr[1][2]=Sin(fTheta)*Sin(fPhi);
  mr[2][2]=Cos(fTheta);

  pol[0]=Sin(thetap)*Cos(phip);
  pol[1]=Sin(thetap)*Sin(phip);
  pol[2]=Cos(thetap);

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      rot[i]+=mr[i][j]*pol[j];
    }
  }

  r=Sqrt(rot[0]*rot[0]+rot[1]*rot[1]+rot[2]*rot[2]);

  fTheta=ACos(rot[2]/r);
  fPhi=phip;

}
