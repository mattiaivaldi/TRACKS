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
  printf("+++ PORCODITO +++\n\n");

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

  TBranch *bmult=tree_reco->GetBranch("LAMULTIANI");
  TBranch *b1=tree_reco->GetBranch("HITL1");
  TBranch *b2=tree_reco->GetBranch("HITL2");

  bmult->SetAddress(&mult_ev);
  b1->SetAddress(&cross_L1);
  b2->SetAddress(&cross_L2);

  for(int i=0;i<tree_reco->GetEntries();i++){
    tree_reco->GetEvent(i);
    cout<<"evento "<<i<<" Mult "<<mult_ev<<endl;
    for(int j=0;j<mult_ev;j++){
      Hit *hit_buffer=(Hit*)cross_L1->At(j);
      cout<<"Zio: "<<hit_buffer->GetX()<<endl;
      //delete hit_buffer;
    }
  }

  printf("+++ END reconstruction +++\n\n");

  //h_reco.Write();
  //h_reco.Close();

  //cpu info
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;
  printf("CPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\nScroll up for info and verbosities. Thanks for using TRACKS!\n\n-> DONATE <-\n\n",cpu_time,real_time, cpu_efficiency);

} //END
