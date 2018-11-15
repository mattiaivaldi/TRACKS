//Class used for object of type particle

#include <TRandom3.h>
#include <Riostream.h>
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include "Particle.h"

using namespace TMath;

ClassImp(Particle)

Particle::Particle():TObject(),
fRap(0.),
fTheta(Pi()/4.),
fPhi(0.)
{
  //default constructor
}

Particle::Particle(TH1F *distr_rap): TObject(),
fRap(distr_rap->GetRandom()),
fTheta(2*ATan(Exp(-(double)fRap))),
fPhi(gRandom->Uniform(2*Pi()))
{
  //standard constructor
}

Particle::~Particle() {
  //default destructor
}


Double_t Particle::GetTheta() const {
  //returns theta
  return fTheta;
}

Double_t Particle::GetPhi() const {
  //returns phi
  return fPhi;
}

Float_t Particle::GetRap() const {
  //returns pseudorapidity
  return fRap;
}

void Particle::Rotate(double rms){

  //rotation matrix for multiple scattering

  gRandom->SetSeed(0);

  double theta0=rms/Sqrt(2);
  double thetap=gRandom->Gaus(0,theta0);//scattering angle theta
  double phip=gRandom->Uniform(2*Pi());//scattering angle phi

  double mr[3][3], pol[3], rot[3], r;//rotation matrixes

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
    rot[i]=0;
    for(int j=0; j<3; j++){
      rot[i]+=mr[i][j]*pol[j];
    }
  }

  fTheta=ACos(rot[2]);

  if(ATan(rot[1]/rot[0]) >= 0) {
    if(double x = gRandom->Rndm() < 0.5){
      fPhi = ATan(rot[1]/rot[0]);
    }
    else {
      fPhi = ATan(rot[1]/rot[0])+Pi();
    }
  }

  else if (ATan(rot[1]/rot[0]) < 0) {
    if(double x = gRandom->Rndm() < 0.5){
      fPhi = ATan(rot[1]/rot[0]) + 2*Pi();
    }
    else {
      fPhi = ATan(rot[1]/rot[0]) + Pi();
    }
  }

}
