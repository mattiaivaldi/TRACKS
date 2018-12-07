//Class used in order to define the hit coorinates of the particle with both the beam pipe and the two detector layers

#include "Riostream.h"
#include "Hit.h"
#include "Layer.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TH1F.h"

//using namespace TMath;

ClassImp(Hit)

Hit::Hit():TObject(),
fX(0.),
fY(0.),
fZ(0.)
{
  //Default constructor
}

Hit::Hit(double x, double y, double z): TObject(),
fX(x),
fY(y),
fZ(z)
{
  //hit constructor
}

Hit::Hit(double R, double H): TObject(),
fX(gRandom->Uniform(-R, R)),
fY(0.),
fZ(gRandom->Uniform(-H/2,H/2))
{
  if(gRandom->Rndm()<0.5){
    fY=TMath::Sqrt(R*R-fX*fX);
  }else{fY=-1*TMath::Sqrt(R*R-fX*fX);}

  //spurious hit constructor
}

Hit::Hit(double meanv, double sigmaxy, double sigmaz, TH1F *distr_mult): TObject(),
fX(gRandom->Gaus(meanv,sigmaxy)),
fY(gRandom->Gaus(meanv,sigmaxy)),
fZ(gRandom->Gaus(meanv,sigmaz)),
fMult((int)distr_mult->GetRandom())
{
  /*double zgen;
  do{zgen=gRandom->Gaus(meanv,sigmaz);}
  while(TMath::Abs(zgen)>13.5);
  fZ=zgen;*/
  //event constructor, def2
}

Hit::Hit(double meanv, double sigmaxy, double z_custom, int mult_custom): TObject(),
fX(gRandom->Gaus(meanv,sigmaxy)),
fY(gRandom->Gaus(meanv,sigmaxy)),
fZ(z_custom),
fMult(mult_custom)
{
  //custom event, def2
}

Hit::~Hit() {
  //default destructor
}

Double_t Hit::GetX() const {
  //returns X
  return fX;
}

Double_t Hit::GetY() const {
  //returns Y
  return fY;
}

Double_t Hit::GetZ() const {
  //returns Z
  return fZ;
}

int Hit::GetMult() const {
  //returns Mult, def2
  return fMult;
}
