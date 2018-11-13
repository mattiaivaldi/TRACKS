//
//TRACKS - generation and recontruction of particle tracks in a detector
//with multiple scattering and noise
//developed by Luca Quaglia and Mattia Ivaldi, 2018
//
//START

#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
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

void test(bool PrintParticles, bool multiscatman, int kExp) {

  //PrintParticles activates verbose mode
  //multiscatman activates multiple scattering

  printf("\n\n+++ TRACKS START - tracks generation +++\n\n");

  if (PrintParticles==true) {
    printf("Printing vertex and hit coordinates: ON\n\n");
  }else{printf("Printing vertex and hit coordinates: OFF\n\n");}

  if (multiscatman==true) {
    printf("Applying multiple scattering: ON\n\nAll distances are in cm, all angles are in rad.\n\n\n");
  }else{printf("Applying multiple scattering: OFF\n\nAll distances are in cm, all angles are in rad.\n\n\n");}

  TStopwatch timer;
  timer.Start(true);

  TFile h_gen("gen.root","RECREATE");
  TTree *tree_gen=new TTree("TG","PORCODIO");

  ////length, radius, thickness, multiscattering RMS
  Layer *BP = new Layer(27.,3.,0.8,0.001);
  Layer *L1 = new Layer(27.,4.,0.2,0.001);
  Layer *L2 = new Layer(27.,7.,0.2,0.001);

  TString distr="kinem.root";
  TFile hfile(distr);
  TH1F *pseudorap = (TH1F*)hfile.Get("heta");
  TH1F *multiplicity = (TH1F*)hfile.Get("hmul");

  TClonesArray *cross_BP=new TClonesArray("Hit",50);
  TClonesArray *cross_L1=new TClonesArray("Hit",50);
  TClonesArray *cross_L2=new TClonesArray("Hit",50);

  TClonesArray &hits_BP=*cross_BP;
  TClonesArray &hits_L1=*cross_L1;
  TClonesArray &hits_L2=*cross_L2;

  tree_gen->Branch("HIT_BP",&cross_BP);
  tree_gen->Branch("HIT_L1",&cross_L1);
  tree_gen->Branch("HIT_L2",&cross_L2);

  for(int i=0; i<kExp; i++){

    printf("> RUN %d <\n\n",i+1);

    //vertex mean, sigmaxy, sigmaz, kinematics file
    Event *vgen = new Event(0, 0.001, 5.3, multiplicity);
    int mult = (int)vgen->GetMult();

    //vector <Hit*> cross_BP, cross_L1, cross_L2;

    printf("Generated vertex with coordinates (%f, %f, %f)\nand multiplicity %d\n\n",vgen->GetX(),vgen->GetY(),vgen->GetZ(),mult);

    //start tracks generation

    for (int j=0; j<mult; j++) {

      Particle *part = new Particle(pseudorap);

      if (PrintParticles==true) {
        printf(">>> Particle %i: theta %f - phi %f <<<\n\n",j+1,part->GetTheta(),part->GetPhi());
      }

      detect(j, vgen, BP, *part, *cross_BP, PrintParticles, multiscatman, "BP");
      detect(j, vgen, L1, *part, *cross_L1, PrintParticles, multiscatman, "L1");
      detect(j, vgen, L2, *part, *cross_L2, PrintParticles, multiscatman, "L2");

      delete part;

    }

    delete vgen;

    //printf("Out of %d generated particles:\n\n%lu crossed BP\n%lu crossed L1\n%lu crossed L2\n\n",mult,cross_BP.size(),cross_L1.size(),cross_L2.size());

    tree_gen->Fill();

    cross_BP->Clear();
    cross_L1->Clear();
    cross_L2->Clear();

  }

  h_gen.Write();
  h_gen.Close();

  printf("+++ END generation +++\n\n");

  //cpu info
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;
  printf("CPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\nScroll up for info and verbosities. Thanks for using TRACKS!\n\n-> DONATE <-\n\n",cpu_time,real_time, cpu_efficiency);

} //END
