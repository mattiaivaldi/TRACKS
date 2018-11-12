//
//TRACKS - generation and recontruction of particle tracks in a detector
//with multiple scattering and noise
//developed by Luca Quaglia and Mattia Ivaldi, 2018
//
//START

#include "TCanvas.h"
#include "TH2F.h"
#include "TGraphPolar.h"
#include "Riostream.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "vector"
#include "Event.h"
#include "Tools.h"
#include "Layer.h"
#include "Hit.h"
#include "Event.h"
#include "Particle.h"

using namespace TMath;

//PrintParticles is a variable used in order to decide wheter or not to print all the info about a particle (verbose)
//multiscatman is used in order to toogle on or off the multiscattering
void test(bool PrintParticles, bool multiscatman) {

  printf("\n\n+++ TRACKS START - tracks generation +++\n\n");

  TStopwatch timer; //declared a timer to perform cpu efficiency measurements

  timer.Start(true); //start the timer

  //3 objects of the class layer which represent the beam pipe, first layer and second layer
  Layer *BP = new Layer(27.,3.,0.8,0.001); //length, radius, thickness and multiscattering RMS
  Layer *L1 = new Layer(27.,4.,0.2,0.001);
  Layer *L2 = new Layer(27.,7.,0.2,0.001);

  //object of class event
  Event *vgen = new Event(0, 0.001, 5.3, "kinem.root");
  //vectors of the class hit used in order to memorize all the coordinates of all the intersections
  vector <Hit*> cross_BP, cross_L1, cross_L2;
  //j, k and l are counters used later on and mult is the cast of the multiplicity (float) to an int value
  int mult = (int)vgen->GetMult();

  if (PrintParticles==true) {
    printf("Printing vertex and hit coordinates: ON\n\n");
  }else{printf("Printing vertex and hit coordinates: OFF\n\n");}

  if (multiscatman==true) {
    printf("Applying multiple scattering: ON\n\n");
  }else{printf("Applying multiple scattering: OFF\n\n");}

  printf("All distances are in cm, all angles are in rad.\n\nGenerated vertex with coordinates (%f, %f, %f) and multiplicity %d\n\n",vgen->GetX(),vgen->GetY(),vgen->GetZ(),mult);

  //cycle over all the particles in current event
  for (int i = 0; i < mult; i++) {

    //create an object of the class particle
    Particle *part = new Particle("kinem.root");

    //print them out only if verbose is on
    if (PrintParticles==true) {
      printf(">>> Particle %i: theta %f - phi %f <<<\n\n",i+1,part->GetTheta(),part->GetPhi());
    }

    detect(*vgen, *BP, *part, cross_BP, PrintParticles, multiscatman, "BP");

    printf(">>> Particle %i: theta %f - phi %f <<<\n\n",i+1,part->GetTheta(),part->GetPhi());

    detect(*vgen, *L1, *part, cross_L1, PrintParticles, multiscatman, "L1");

    detect(*vgen, *L2, *part, cross_L2, PrintParticles, multiscatman, "L2");

    delete part;

  }

  //pronted out how many particles have crossed which layer
  printf("Out of %d generated particles:\n\n%lu crossed BP\n%lu crossed L1\n%lu crossed L2\n\n+++ END generation +++",mult,cross_BP.size(),cross_L1.size(),cross_L2.size());

  TCanvas * CPol = new TCanvas("CPol","TGraphPolar Example",500,500);

  Double_t theta[8];
  Double_t radius[8];
  Double_t etheta[8];
  Double_t eradius[8];

  for (int i=0; i<8; i++) {
     theta[i]   = (i+1)*(TMath::Pi()/4.);
     radius[i]  = (i+1)*0.05;
     etheta[i]  = TMath::Pi()/8.;
     eradius[i] = 0.05;
  }

  TGraphPolar * grP1 = new TGraphPolar(8, theta, radius, etheta, eradius);
  grP1->SetTitle("#theta - MC data");

  grP1->SetMarkerStyle(20);
  grP1->SetMarkerSize(2.);
  grP1->SetMarkerColor(4);
  grP1->SetLineColor(2);
  grP1->SetLineWidth(3);
  grP1->Draw("PE");

  TCanvas *c4 = new TCanvas("c4","c4",600,400);
  TH2F *hscc = new TH2F("hscc","Cylindrical coordinates",20,-4,4,20,-20,20);
  Float_t px, py;
  for (Int_t i = 0; i < 25000; i++) {
      gRandom->Rannor(px,py);
      hscc->Fill(px-1,5*py);
      hscc->Fill(2+0.5*px,2*py-10.,0.1);
   }
   hscc->Draw("SURF1 PSR");
   hscc->SetTitle("PseudoRapidity/Phi coordinates");

  //random info on the cpu usage
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;
  printf("\n\nCPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\n",cpu_time,real_time, cpu_efficiency);

} //END
