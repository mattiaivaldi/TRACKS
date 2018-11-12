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

void test(bool PrintParticles, bool multiscatman) {

  //PrintParticles activates verbose mode
  //multiscatman activates multiple scattering

  printf("\n\n+++ TRACKS START - tracks generation +++\n\n");

  TStopwatch timer;
  timer.Start(true);

  ////length, radius, thickness, multiscattering RMS
  Layer *BP = new Layer(27.,3.,0.8,0.001);
  Layer *L1 = new Layer(27.,4.,0.2,0.001);
  Layer *L2 = new Layer(27.,7.,0.2,0.001);

  //vertex mean, sigmaxy, sigmaz, kinematics file
  Event *vgen = new Event(0, 0.001, 5.3, "kinem.root");
  int mult = (int)vgen->GetMult();

  double *hit_buffer_BP, *hit_buffer_L1, *hit_buffer_L2;
  vector <Hit*> cross_BP, cross_L1, cross_L2;

  int j = 0, k = 0, l = 0;
  double theta[mult], phi[mul];

  //verbosities

  if (PrintParticles==true) {
    printf("Printing vertex and hit coordinates: ON\n\n");
  }else{printf("Printing vertex and hit coordinates: OFF\n\n");}

  if (multiscatman==true) {
    printf("Applying multiple scattering: ON\n\n");
  }else{printf("Applying multiple scattering: OFF\n\n");}

  printf("All distances are in cm, all angles are in rad.\n\nGenerated vertex with coordinates (%f, %f, %f) and multiplicity %d\n\n",vgen->GetX(),vgen->GetY(),vgen->GetZ(),mult);

  //start tracks generation

  /*for (int i = 0; i < mult; i++) {

    Particle *part = new Particle("kinem.root");

    if (PrintParticles==true) {
      printf(">>> Particle %i: theta %f - phi %f <<<\n\n",i+1,part->GetTheta(),part->GetPhi());
    }

    detect(vgen, BP, *part, cross_BP, PrintParticles, multiscatman, "BP");

    printf("Check after detect: theta %f - phi %f\n\n",part->GetTheta(),part->GetPhi());

    detect(vgen, L1, *part, cross_L1, PrintParticles, multiscatman, "L1");

    printf("Check after detect: theta %f - phi %f\n\n",part->GetTheta(),part->GetPhi());

    detect(vgen, L2, *part, cross_L2, PrintParticles, multiscatman, "L2");

    delete part;

    cout<<endl;

  }*/

  bool bBP = false, bL1 = false, bL2 = false;

  for (int i = 0; i < mult; i++) {

    Particle *part = new Particle("kinem.root");

    phi[i] = part->GetPhi();
    theta[i] = part->GetTheta();

    if (PrintParticles==true) {
      printf("Particle %i: phi %f - theta %f\n",i,phi[i],theta[i]);
    }

    //intersection with beam pipe
    hit_buffer_BP=hit_point(vgen->GetX(),vgen->GetY(),vgen->GetZ(),theta[i],phi[i],BP->GetRadius());

    if(*(hit_buffer_BP+2) >= -(BP->GetWidth()/2.) && *(hit_buffer_BP+2) <= (BP->GetWidth()/2.)) {
      bBP = true;
      Hit *hit_BP = new Hit(*(hit_buffer_BP+0),*(hit_buffer_BP+1),*(hit_buffer_BP+2));
      cross_BP.push_back(hit_BP);
      if (PrintParticles==true) {
	       printf("Hit with BP at (%f, %f, %f)\n",cross_BP[j]->GetX(),cross_BP[j]->GetY(),cross_BP[j]->GetZ());
      }
      cout << "Angoli prima " << phi[i] << " " << theta[i] << endl;
      if (bBP == true && multiscatman == true) {
        part->Rotate(BP->GetRMS());
        phi[i] = part->GetPhi();
        theta[i] = part->GetTheta();
      }
      j++;
      cout << "Angoli dopo " <<phi[i] << " " << theta[i] << endl;
    }else{bBP=false;}

    //intersection with L1
    hit_buffer_L1=hit_point(vgen->GetX(),vgen->GetY(),vgen->GetZ(),theta[i],phi[i],L1->GetRadius());

    if(*(hit_buffer_L1+2) >= -(L1->GetWidth()/2.) && *(hit_buffer_L1+2) <= (L1->GetWidth()/2.)) {
      bL1 = true;
      Hit *hit_L1 = new Hit(*(hit_buffer_L1+0),*(hit_buffer_L1+1),*(hit_buffer_L1+2));
      cross_L1.push_back(hit_L1);
      if (PrintParticles==true) {
	       printf("Hit with L1 at (%f, %f, %f)\n",cross_L1[k]->GetX(),cross_L1[k]->GetY(),cross_L1[k]->GetZ());
      }
      cout << "Angoli prima " << phi[i] << " " << theta[i] << endl;
      if (bL1 == true && multiscatman == true) {
        part->Rotate(L1->GetRMS());
        phi[i] = part->GetPhi();
        theta[i] = part->GetTheta();
      }
      k++;
      cout << "Angoli dopo " << phi[i] << " " << theta[i] << endl;
    }else{bL1=false;}

    //intersection with L2
    hit_buffer_L2=hit_point(vgen->GetX(),vgen->GetY(),vgen->GetZ(),theta[i],phi[i],L2->GetRadius());

    if(*(hit_buffer_L2+2) >= -(L2->GetWidth()/2.) && *(hit_buffer_L2+2) <= (L2->GetWidth()/2.)) {
      Hit *hit_L2 = new Hit(*(hit_buffer_L2+0),*(hit_buffer_L2+1),*(hit_buffer_L2+2));
      cross_L2.push_back(hit_L2);
      if (PrintParticles==true) {
	printf("Hit with L2 at (%f, %f, %f)\n",cross_L2[l]->GetX(),cross_L2[l]->GetY(),cross_L2[l]->GetZ());
      }
      l++;
    }else{bL2=false;}

    if (PrintParticles==true){printf("\n");}

    delete part; //deleting the object at the end of the for cycle

  }

  printf("Out of %d generated particles:\n\n%lu crossed BP\n%lu crossed L1\n%lu crossed L2\n\n+++ END generation +++",mult,cross_BP.size(),cross_L1.size(),cross_L2.size());

  /*TCanvas * CPol = new TCanvas("CPol","TGraphPolar Example",500,500);

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
   hscc->SetTitle("PseudoRapidity/Phi coordinates");*/

  //cpu info
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;
  printf("\n\nCPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\n",cpu_time,real_time, cpu_efficiency);

} //END
