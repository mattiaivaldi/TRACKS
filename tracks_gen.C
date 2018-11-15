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

bool war1=true;//declaring war

using namespace TMath;

void tracks_gen(bool PrintParticles, bool multiscatman, bool paolonoise, int kExp) {

  //PrintParticles activates verbose mode (0,1)
  //multiscatman activates multiple scattering (0,1)
  //paolonoise activates the noise (0,1)
  //kExp is the number of collisions you want to perform

  TStopwatch timer;
  timer.Start(true);//start cpu monitor

  gRandom->SetSeed(0);

  //verbosities
  verbosities(PrintParticles, multiscatman, paolonoise);

  TFile h_gen("gen.root","RECREATE");
  TTree *tree_gen=new TTree("TG","tree_gen");

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

  int mult_ev;
  TClonesArray& hits_BP=*cross_BP;
  TClonesArray& hits_L1=*cross_L1;
  TClonesArray& hits_L2=*cross_L2;//TCA aliases to be filled

  tree_gen->Branch("LAMULTIANI",&mult_ev);
  tree_gen->Branch("HITL1",&cross_L1);
  tree_gen->Branch("HITL2",&cross_L2);//branch booking

  //start experiment

  for(int i=0; i<kExp; i++){

    //mult_ev=0;

    //vertex mean, sigmaxy, sigmaz, kinematics file
    Hit *vgen = new Hit(0, 0.001, 5.3, multiplicity);
    int mult = vgen->GetMult();
    mult_ev=mult+kNoise;

    int counter_BP=0,counter_L1=0,counter_L2=0;

    printf("> EVENT %d <\n\nGenerated vertex with coordinates (%f, %f, %f)\nand multiplicity %d\n\n",i+1,vgen->GetX(),vgen->GetY(),vgen->GetZ(),mult);

    //start tracks generation

    for (int j=0; j<mult; j++) {

      Particle *part = new Particle(pseudorap);

      if (PrintParticles==true) {
        printf(">>> Particle %i: theta %f - phi %f <<<\n\n",j+1,part->GetTheta(),part->GetPhi());
      }

      //if particle hits layer TCA is filled, otherwhise gets 0
      detect(j, vgen, BP, *part, hits_BP, PrintParticles, multiscatman, "BP", counter_BP);
      detect(j, vgen, L1, *part, hits_L1, PrintParticles, multiscatman, "L1", counter_L1);
      detect(j, vgen, L2, *part, hits_L2, PrintParticles, multiscatman, "L2", counter_L2);

      delete part;

    }

    if(PrintParticles==true){
      printf("Out of %d generated particles\n%d crossed BP\n%d crossed L1\n%d crossed L2\n\n",mult,counter_BP, counter_L1, counter_L2);
    }

    //randomly add or not add noise

    if(paolonoise==true){
      noise(PrintParticles,kNoise,mult,hits_L1, L1,"L1");
      noise(PrintParticles,kNoise,mult,hits_L2, L2,"L2");
    }

    tree_gen->Fill();

    cross_BP->Clear();
    cross_L1->Clear();
    cross_L2->Clear();

    delete vgen;

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
