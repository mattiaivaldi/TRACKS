//
//TRACKS - MC generation of particle tracks in a detector
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

bool war2=true;//declaring war

using namespace TMath;

void tracks_reco(bool PrintParticles) {

  //PrintParticles activates verbose mode (0,1)
  printf("\n\n+++ START reconstruction +++\n\n");

  TStopwatch timer;
  timer.Start(true);//start cpu monitor

  gRandom->SetSeed(0);

  //verbosities
  //verbosities(PrintParticles, multiscatman, paolonoise);

  TFile h_reco("gen.root","READ");
  TTree *tree_reco=(TTree*)h_reco.Get("TG");

  int kExp=tree_reco->GetEntries();

  int mult_ev=0;
  TClonesArray *cross_L1=new TClonesArray("Hit",kExp);
  TClonesArray *cross_L2=new TClonesArray("Hit",kExp);

  TClonesArray& hits_L1=*cross_L1;
  TClonesArray& hits_L2=*cross_L2;//TCA aliases to be smeagled

  TBranch *bmult=tree_reco->GetBranch("LAMULTIANI");
  TBranch *b1=tree_reco->GetBranch("HITL1");
  TBranch *b2=tree_reco->GetBranch("HITL2");

  bmult->SetAddress(&mult_ev);
  b1->SetAddress(&cross_L1);
  b2->SetAddress(&cross_L2);

  int mult;

  for(int i=0;i<tree_reco->GetEntries();i++){

    tree_reco->GetEvent(i);
    int multi_ev=cross_L1->GetEntries();

    printf("> EVENT %i - # of hits %i <\n\n",i+1,multi_ev);

    for(int j=0;j<multi_ev;j++){

      Hit *hit_buffer1=(Hit*)cross_L1->At(j);
      //Hit *hit_buffer2=(Hit*)cross_L2->At(j);

      printf("Hit with L1 prima at (%f %f %f)\n\n",hit_buffer1->GetX(),hit_buffer1->GetY(),hit_buffer1->GetZ());

      if(hit_buffer1->GetX() != 0){
        smeagol(j,0.00120,0.0003,4,hits_L1);
      }

      printf("Fuori da smeagol dopo at (%f %f %f)\n\n",hit_buffer1->GetX(),hit_buffer1->GetY(),hit_buffer1->GetZ());

      //printf("Hit with L1 at (%f %f %f)\nHit with L2 at (%f %f %f)\n\n",hit_buffer1->GetX(),hit_buffer1->GetY(),hit_buffer1->GetZ(),hit_buffer2->GetX(),hit_buffer2->GetY(),hit_buffer2->GetZ());

    }
    cout<<endl;
  }

  printf("+++ END reconstruction +++\n\n");

  //cpu info
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;
  printf("CPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\nScroll up for info and verbosities. Thanks for using TRACKS!\n\n-> DONATE <-\n\n",cpu_time,real_time, cpu_efficiency);

} //END