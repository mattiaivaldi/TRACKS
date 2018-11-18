#include "vector"
#include "Tools.h"
#include "Layer.h"
#include "Hit.h"
#include "Particle.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TClonesArray.h"
#include "TSystem.h"

//using namespace TMath;

void verbosities(bool b_verbose, bool b_multiscatter, bool b_noise){

  printf("\n\n+++ TRACKS START - tracks generation +++\n\n");

  if (b_verbose==true) {
    printf("Printing vertex and hit coordinates: ON\n\n");
  }else{printf("Printing vertex and hit coordinates: OFF\n\n");}

  if (b_multiscatter==true) {
    printf("Applying multiple scattering: ON\n\n");
  }else{printf("Applying multiple scattering: OFF\n\n");}

  if (b_noise==true) {
    printf("Applying noise: ON\n\nAll distances are in cm, all angles are in rad.\n\n\n");
  }else{printf("Applying noise: OFF\n\nAll distances are in cm, all angles are in rad.\n\n\n");}

}

double *hit_point(double x0, double y0, double z0, double theta, double phi, double R) {

  static double hit[3];
  double c1 = TMath::Sin(theta)*TMath::Cos(phi), c2 = TMath::Sin(theta)*TMath::Sin(phi), c3 = TMath::Cos(theta);//direction cosines
  double delta = 2*x0*y0*c1*c2 - c1*c1*y0*y0 + c1*c1*R*R -c2*c2*x0*x0 + c2*c2*R*R;//delta of II degree equation (>= 0 by construction)
  double t_p = (-(x0*c1 - y0*c2) + TMath::Sqrt(delta))/(c1*c1 + c2*c2);//+ solution
  double t_m = (-(x0*c1 - y0*c2) - TMath::Sqrt(delta))/(c1*c1 + c2*c2);//- solution

  //calculate the values of the intersection points (x,y,z)
  if (t_p >= 0) {
    hit[0] = x0 + c1*t_p;
    hit[1] = y0 + c2*t_p;
    hit[2] = z0 + c3*t_p;
  }

  else {
    hit[0] = x0 + c1*t_m;
    hit[1] = y0 + c2*t_m;
    hit[2] = z0 + c3*t_m;
  }

  return hit;
}

void detect(int &index, Hit* vtx, Layer* L, Particle &part, TClonesArray &cross, bool b_verbose, bool b_multiscatter, char const *detector, int &counter){

  double *hit_buffer;
  bool b_cross=false;

  hit_buffer=hit_point(vtx->GetX(),vtx->GetY(),vtx->GetZ(),part.GetTheta(),part.GetPhi(),L->GetRadius());//evaluate hit point coordinates

  if(*(hit_buffer+2) >= -(L->GetWidth()/2.) && *(hit_buffer+2) <= (L->GetWidth()/2.)) {

    b_cross = true;//yes we have detection
    gSystem->Beep(440,1);
    counter++;

    new(cross[index])Hit(*(hit_buffer+0),*(hit_buffer+1),*(hit_buffer+2));//fill with hit

    if (b_cross == true && b_multiscatter == true) {
      part.Rotate(L->GetRMS());//if multiscattering ON set new angles
    }

    if (b_verbose==true) {
      printf("Hit with %s at (%f, %f, %f)\nAngles after: theta %f - phi %f\n\n",detector,((Hit*)cross[index])->GetX(),((Hit*)cross[index])->GetY(),((Hit*)cross[index])->GetZ(),part.GetTheta(),part.GetPhi());
    }

    index++;

  }//else{new(cross[index])Hit();}//otherwise fill with 0

}

void noise(bool b_verbose, int Noise, int Mult, TClonesArray cross, Layer* L, char const *detector){

  int index_noise=Mult;

  for(int i=0; i<Noise; i++){
    new(cross[index_noise])Hit(L->GetRadius(),L->GetWidth());//random spurious hit
    if(b_verbose==true){
        printf("> Noise hit with %s at (%f, %f,%f) <\n\n",detector,((Hit*)cross[index_noise])->GetX(),((Hit*)cross[index_noise])->GetY(),((Hit*)cross[index_noise])->GetZ());
      }
    index_noise++;
  }

}

void smeagol(int index, double sigmaz, double sigmarf, double R, TClonesArray &cross){

  Hit *buffer=(Hit*)cross.At(index);
  double x=buffer->GetX();
  double y=buffer->GetY();
  double z=buffer->GetZ();
  double phi=0.;
  double theta=TMath::ACos(z/TMath::Sqrt(x*x+y*y+z*z));
  //double theta=TMath::ACos(z);

  /*printf("Dentro smeagol usando i getter su buffer (%f %f %f)\n\n",buffer->GetX(), buffer->GetY(), buffer->GetZ());
  printf("Dentro smeagol usando x y z (%f %f %f)\n\n",x, y, z);

  printf("Dentro smeagol theta %f\n\n",theta);*/

  if(x>0&&y>0){
    phi=TMath::ATan(y/x);
  }else if(x>0&&y<0){
    phi=TMath::ATan(y/x)+2*TMath::Pi();
  }
  else if(x<0&&y<0){
    phi=TMath::ATan(y/x)+TMath::Pi();
  }
  else if(x<0&&y>0){
    phi=TMath::ATan(y/x)+TMath::Pi();
  }

  if(gRandom->Rndm()<0.5){
    z+=gRandom->Gaus(0,sigmaz);
    theta=TMath::ACos(z);
    phi+=gRandom->Gaus(0,sigmarf)/R;
    x=TMath::Sin(theta)*TMath::Cos(phi);
    y=TMath::Sin(theta)*TMath::Sin(phi);
  }else{
    z-=gRandom->Gaus(0,sigmaz);
    theta=TMath::ACos(z);
    phi-=gRandom->Gaus(0,sigmarf)/R;
    x=TMath::Sin(theta)*TMath::Cos(phi);
    y=TMath::Sin(theta)*TMath::Sin(phi);
  }

  //theta=TMath::ACos(z/TMath::Sqrt(x*x+y*y+z*z));

  /*printf("dentro smeagol dopo lo smeagol theta %f phi %f\n\n",theta,phi);

  printf("dentro smeagol dopo lo smeagol direttamente getter (%f %f %f)\n\n",buffer->GetX(), buffer->GetY(), buffer->GetZ());
  printf("dentro smeagol dopo lo smeagol x y z (%f %f %f)\n\n",x, y, z);*/

  buffer->SetX(x);
  buffer->SetY(y);
  buffer->SetZ(z);

  //printf("dentro smeagol dopo lo smeagol dopo i setter (%f %f %f)\n\n",buffer->GetX(), buffer->GetY(), buffer->GetZ());
}
