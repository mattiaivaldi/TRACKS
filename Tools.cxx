#include "vector"
#include "Event.h"
#include "Tools.h"
#include "Layer.h"
#include "Hit.h"
#include "Event.h"
#include "Particle.h"
#include "TMath.h"

//using namespace TMath;

double *hit_point(double x0, double y0, double z0, double theta, double phi, double R) {

    static double hit[3];
    double c1 = TMath::Sin(theta)*TMath::Cos(phi), c2 = TMath::Sin(theta)*TMath::Sin(phi), c3 = TMath::Cos(theta); //direction cosines
    double delta = 2*x0*y0*c1*c2 - c1*c1*y0*y0 + c1*c1*R*R -c2*c2*x0*x0 + c2*c2*R*R; //delta of II degree equation ( >= 0 by construction)
    double t_p = (-(x0*c1 - y0*c2) + TMath::Sqrt(delta))/(c1*c1 + c2*c2); //solution with the "+" sign
    double t_m = (-(x0*c1 - y0*c2) - TMath::Sqrt(delta))/(c1*c1 + c2*c2); //solution with the "-" sign

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

    //returns the pointer to the first element of the array hit which contains the values x, y and z
    return hit;
}

void detect(Event* vtx, Layer* L, double &theta, double &phi, Particle* &part, vector<Hit*> &cross, bool b_verbose, bool b_multiscatter, char const *detector){

  double *hit_buffer;
  bool b_cross=false;

  hit_buffer=hit_point(vtx->GetX(),vtx->GetY(),vtx->GetZ(),theta,phi,L->GetRadius());

  if(*(hit_buffer+2) >= -(L->GetWidth()/2.) && *(hit_buffer+2) <= (L->GetWidth()/2.)) {

    b_cross = true;

    Hit *hit = new Hit(*(hit_buffer+0),*(hit_buffer+1),*(hit_buffer+2));

    cross.push_back(hit);

    if (b_verbose==true) {
      printf("Hit with %s at (%f, %f, %f)\n",detector,*(hit_buffer+0),*(hit_buffer+1),*(hit_buffer+2));
    }

    if (b_cross == true && b_multiscatter == true) {
      part->Rotate(0.001);
      theta=part->GetTheta();
      phi=part->GetPhi();
    }

    if (b_multiscatter == true || b_multiscatter==false) {
      printf("Angles after: theta %f - phi %f\n\n",theta, phi);
    }

  }

}
