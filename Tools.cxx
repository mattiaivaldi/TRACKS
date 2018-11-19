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

  if (b_verbose) {
    printf("Printing vertex and hit coordinates: ON\n\n");
  }else{printf("Printing vertex and hit coordinates: OFF\n\n");}

  if (b_multiscatter) {
    printf("Applying multiple scattering: ON\n\n");
  }else{printf("Applying multiple scattering: OFF\n\n");}

  if (b_noise) {
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

    if (b_cross&&b_multiscatter) {
      part.Rotate(L->GetRMS());//if multiscattering ON set new angles
    }

    if (b_verbose) {
      printf("Hit with %s at (%f, %f, %f)\nAngles after: theta %f - phi %f\n\n",detector,((Hit*)cross[index])->GetX(),((Hit*)cross[index])->GetY(),((Hit*)cross[index])->GetZ(),part.GetTheta(),part.GetPhi());
    }

    index++;

  }else{
    if(b_verbose){
      printf("Does not hit %s\n\n",detector);
    }
  }

}

void noise(bool b_verbose, int Noise, int Mult, TClonesArray cross, Layer* L, char const *detector){

  int index_noise=Mult;

  for(int i=0; i<Noise; i++){
    new(cross[index_noise])Hit(L->GetRadius(),L->GetWidth());//random spurious hit
    if(b_verbose){
        printf("> Noise hit with %s at (%f, %f,%f) <\n\n",detector,((Hit*)cross[index_noise])->GetX(),((Hit*)cross[index_noise])->GetY(),((Hit*)cross[index_noise])->GetZ());
      }
    index_noise++;
  }

}

void smeagol(int index, double sigmaz, double sigmarf, double R, TClonesArray &cross){

    double dphi=0;
    double dz=0;

    Hit *hit_buffer=(Hit*)cross.At(index);

    double x=hit_buffer->GetX();
    double y=hit_buffer->GetY();
    double z=hit_buffer->GetZ();

    double phi=TMath::ACos(x/R);
    double theta=TMath::ACos(z/TMath::Sqrt(x*x+y*y+z*z));

    printf("[Entro in smeagol]\n\nGetter su buffer (%f %f %f)\n\nx, y, z (%f %f %f)\n\ntheta %f phi %f\n\n",hit_buffer->GetX(), hit_buffer->GetY(), hit_buffer->GetZ(),x,y,z,theta,phi);

    if (gRandom->Rndm()<0.5) {
        dphi=(gRandom->Gaus(0,sigmarf))/R;
        x=TMath::Cos(phi+dphi);
        y=TMath::Sin(phi+dphi);
        dz=gRandom->Gaus(0,sigmaz);
        z+=dz;
        theta=TMath::ACos(z/TMath::Sqrt(x*x+y*y+z*z));
    } else {
        dphi=(gRandom->Gaus(0,sigmarf))/R;
        x=TMath::Cos(phi-dphi);
        y=TMath::Sin(phi-dphi);
        dz=gRandom->Gaus(0,sigmaz);
        z-=dz;
        theta=TMath::ACos(z/TMath::Sqrt(x*x+y*y+z*z));
    }

    printf("[dopo i conti]\n\nGetter su buffer (%f %f %f)\n\nx, y, z (%f %f %f)\n\ntheta %f phi %f\n\n",hit_buffer->GetX(), hit_buffer->GetY(), hit_buffer->GetZ(),x,y,z,theta,phi);

    hit_buffer->SetX(x);
    hit_buffer->SetY(y);
    hit_buffer->SetZ(z);

    printf("[dopo i setter]\n\nGetter su buffer (%f %f %f)\n\nx, y, z (%f %f %f)\n\ntheta %f phi %f\n\n",hit_buffer->GetX(), hit_buffer->GetY(), hit_buffer->GetZ(),x,y,z,theta,phi);

}
