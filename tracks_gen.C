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
#include "Tools.h"
#include "Layer.h"
#include "Hit.h"
#include "Particle.h"

bool war=true;//declaring war

using namespace TMath;

void test(bool PrintParticles, bool multiscatman, bool paolonoise, int kExp) {

  //PrintParticles activates verbose mode
  //multiscatman activates multiple scattering
  //paolonoise activates the noise
  //kExp is the number of collisions

  TStopwatch timer;
  timer.Start(true);//start cpu monitor

  gRandom->SetSeed(0);

  //verbosities

  printf("\n\n+++ TRACKS START - tracks generation +++\n\n");

  if (PrintParticles==true) {
    printf("Printing vertex and hit coordinates: ON\n\n");
  }else{printf("Printing vertex and hit coordinates: OFF\n\n");}

  if (multiscatman==true) {
    printf("Applying multiple scattering: ON\n\nAll distances are in cm, all angles are in rad.\n\n\n");
  }else{printf("Applying multiple scattering: OFF\n\nAll distances are in cm, all angles are in rad.\n\n\n");}

  TFile h_gen("gen.root","RECREATE");
  TTree *tree_gen=new TTree("TG","PORCODIO");

  //length, radius, thickness, multiscattering RMS
  Layer *BP = new Layer(27.,3.,0.8,0.001);
  Layer *L1 = new Layer(27.,4.,0.2,0.001);
  Layer *L2 = new Layer(27.,7.,0.2,0.001);

  TString distr="kinem.root";//get starting kinematics
  TFile hfile(distr);
  TH1F *pseudorap = (TH1F*)hfile.Get("heta");//get pseudorapidity distribution
  TH1F *multiplicity = (TH1F*)hfile.Get("hmul");//get multiplicity distribution

  int kNoise=0;
  if(paolonoise==true){
    kNoise=(int)gRandom->Integer(5);//number of spurious hits
  }

  int size=multiplicity->FindLastBinAbove(0,1)+kNoise;

  TClonesArray *cross_BP=new TClonesArray("Hit",size);
  TClonesArray *cross_L1=new TClonesArray("Hit",size);
  TClonesArray *cross_L2=new TClonesArray("Hit",size);//TCA booking

  printf("%d\n\n",size);

  TClonesArray& hits_BP=*cross_BP;
  TClonesArray& hits_L1=*cross_L1;
  TClonesArray& hits_L2=*cross_L2;//TCA aliases to be filled

  tree_gen->Branch("HIT_BP",&cross_BP);
  tree_gen->Branch("HIT_L1",&cross_L1);
  tree_gen->Branch("HIT_L2",&cross_L2);//branch booking

  //start experiment

  for(int i=0; i<kExp; i++){

    //vertex mean, sigmaxy, sigmaz, kinematics file
    Hit *vgen = new Hit(0, 0.001, 5.3, multiplicity);
    int mult = vgen->GetMult();
    int index_noise=0;

    printf("> EVENT %d <\n\nGenerated vertex with coordinates (%f, %f, %f)\nand multiplicity %d\n\n",i+1,vgen->GetX(),vgen->GetY(),vgen->GetZ(),mult);

    //start tracks generation

    for (int j=0; j<mult; j++) {

      Particle *part = new Particle(pseudorap);

      if (PrintParticles==true) {
        printf(">>> Particle %i: theta %f - phi %f <<<\n\n",j+1,part->GetTheta(),part->GetPhi());
      }

      //if particle hits layer TCA is filled, otherwhise gets 0
      detect(j, vgen, BP, *part, hits_BP, PrintParticles, multiscatman, "BP");
      detect(j, vgen, L1, *part, hits_L1, PrintParticles, multiscatman, "L1");
      detect(j, vgen, L2, *part, hits_L2, PrintParticles, multiscatman, "L2");

      delete part;

    }

    delete vgen;

    //randomly add or not add noise

    index_noise=0;

    for(int k=0; k<kNoise; k++){
      if(gRandom->Rndm()>0.5){
        new(hits_L1[index_noise])Hit();
        new(hits_L2[index_noise])Hit();
      }else{
        new(hits_L1[index_noise])Hit();
        new(hits_L2[index_noise])Hit();
      }
      index_noise++;
    }

    tree_gen->Fill();

    cross_BP->Clear();
    cross_L1->Clear();
    cross_L2->Clear();

  }//end for up to kExp

  printf("+++ END generation +++\n\n");

  h_gen.Write();
  h_gen.Close();

  //cpu info
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;
  printf("CPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\nScroll up for info and verbosities. Thanks for using TRACKS!\n\n-> DONATE <-\n\n",cpu_time,real_time, cpu_efficiency);

} //END
