//
//TRACKS - MC generation of particle tracks in a detector
//with multiple scattering and noise
//developed by Luca Quaglia and Mattia Ivaldi, 2018
//
//START

#include <ctime>
#include "TNtuple.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TString.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1D.h"
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

void tracks_reco(bool PrintParticles, double smear_z, double smear_phi){

  //PrintParticles activates verbose mode (0,1)

  printf("\n\n+++ START reconstruction +++\n\n");

  TStopwatch timer;
  timer.Start(true);//start cpu monitor

  gRandom->SetSeed(0);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);

  TH1D *h_zreco=new TH1D("h_zreco","TRACKS reconstruction - Z Vertex",100,-13.5,13.5);
  h_zreco->GetXaxis()->SetTitle("Z_{V} [cm]");
  h_zreco->GetYaxis()->SetTitle("# [a.u.]");
  h_zreco->GetXaxis()->SetTitleOffset(1.1);
  h_zreco->GetYaxis()->SetTitleOffset(0.8);
  h_zreco->GetXaxis()->SetTitleSize(0.04);
  h_zreco->GetYaxis()->SetTitleSize(0.04);

  TFile h_reco("gen.root","READ");
  TTree *tree_reco=(TTree*)h_reco.Get("TG");

  int kExp=tree_reco->GetEntries();

  TClonesArray *cross_L1=new TClonesArray("Hit",kExp);
  TClonesArray *cross_L2=new TClonesArray("Hit",kExp);

  TClonesArray& hits_L1=*cross_L1;
  TClonesArray& hits_L2=*cross_L2;//TCA aliases to be smeagled

  TBranch *b1=tree_reco->GetBranch("HITL1");
  TBranch *b2=tree_reco->GetBranch("HITL2");

  b1->SetAddress(&cross_L1);
  b2->SetAddress(&cross_L2);

  double phi1=0, phi2=0, x1=0,x2=0, z1=0, z2=0, z_reco=0;

  for(int i=0;i<kExp;i++){

    tree_reco->GetEvent(i);
    int mult_ev1=cross_L1->GetEntries();
    int mult_ev2=cross_L2->GetEntries();

    if(PrintParticles){
      printf("> EVENT %i <\n\n%i hits with L1\n\n%i hits with L2\n\n",i+1,mult_ev1, mult_ev2);
    }else if((i+1)%(kExp/10)==0){
      printf("> EVENT %i <\n\n%i hits with L1\n\n%i hits with L2\n\n[running]\n\n",i+1,mult_ev1, mult_ev2);
    }

    for(int l=0;l<mult_ev2;l++){

      Hit *hit_buffer2=(Hit*)cross_L2->At(l);
      smeagol(l,smear_z,smear_phi,7.,hits_L2);
      phi2=ACos(hit_buffer2->GetX()/7.);

      for(int m=0;m<mult_ev1;m++){

        Hit *hit_buffer1=(Hit*)cross_L1->At(m);
        smeagol(m,smear_z,smear_phi,4.,hits_L1);
        phi1=ACos(hit_buffer1->GetX()/4.);

        if(Abs(phi2-phi1)<0.01){

          x1=hit_buffer1->GetX();
          x2=hit_buffer2->GetX();
          z1=hit_buffer1->GetZ();
          z2=hit_buffer2->GetZ();
          z_reco=z1+((z2-z1)/(x1-x2)*x1);

          if(z_reco<13.5&&z_reco>-13.5){h_zreco->Fill(z_reco);}
        }
      }
    }
  }

  printf("+++ END reconstruction +++\n\nSaving files...\n\n");

  TCanvas *c_zreco=new TCanvas("c_zreco","c_zreco",600,400);
  c_zreco->cd();

  TPaveText *pt_1 = new TPaveText(0.15,0.7,0.35,0.85,"NDC");
  pt_1->SetTextAlign(13);
  pt_1->SetFillStyle(0);
  pt_1->SetShadowColor(0);
  pt_1->SetLineColor(0);
  pt_1->SetBorderSize(0);
  pt_1->SetMargin(0);
  pt_1->SetTextSize(.03);
  pt_1->AddText(Form("%d events",kExp));
  pt_1->AddText(Form("#mu = %f cm",h_zreco->GetMean(1)));
  pt_1->AddText(Form("#sigma = %f cm",h_zreco->GetStdDev(1)));

  TFile *f_reco = new TFile("z_vertex.root","RECREATE"); //file to save histogram
  h_zreco->Draw();
  pt_1->Draw();
  h_zreco->Write();
  c_zreco->SaveAs("c_zreco.eps");
  f_reco->Close();

  //cpu info
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;
  printf("CPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\nScroll up for info and verbosities. Thanks for using TRACKS!\n\n-> DONATE <-\n\n",cpu_time,real_time, cpu_efficiency);

} //END