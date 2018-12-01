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
#include "TH1D.h"
#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
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

struct reco_perform{
  double cputime;
  double runtime;
  double reso;
  double eff;
};

bool war2=true;//declaring war

using namespace TMath;

reco_perform tracks_reco(bool PrintParticles, double smear_z, double smear_phi){

  reco_perform perform;

  //PrintParticles activates verbose mode (0,1)

  printf("\n\n+++ START reconstruction +++\n\n");

  TStopwatch timer;
  timer.Start(true);//start cpu monitor

  TString dirplot=TString("/Users/mattiaivaldi/GitHub/TRACKS/")+"tracksplot/";

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

  TH1D *h_ROI=new TH1D("h_ROI","TRACKS reconstruction - ROI",270,-13.45,13.55);
  TH1D *h_tracklet=new TH1D("h_tracklet","TRACKS reconstruction - tracklet",270000,-13.49995,13.50005);

  TH1F *h_reso=new TH1F("h_reso","TRACKS reconstruction - resolution",200,-0.0995,0.1005);

  TFile h_gen("gen.root","READ");
  TTree *tree_gen=(TTree*)h_gen.Get("TG");
  TNtuple *z_gen=(TNtuple*)h_gen.Get("z_gen");
  float zgen;
  z_gen->SetBranchAddress("z_gen",&zgen);

  int kExp=tree_gen->GetEntries();

  TClonesArray *cross_L1=new TClonesArray("Hit",100);
  TClonesArray *cross_L2=new TClonesArray("Hit",100);

  TClonesArray& hits_L1=*cross_L1;
  TClonesArray& hits_L2=*cross_L2;//TCA aliases to be smeagled

  TBranch *b1=tree_gen->GetBranch("HITL1");
  TBranch *b2=tree_gen->GetBranch("HITL2");

  b1->SetAddress(&cross_L1);
  b2->SetAddress(&cross_L2);

  int goodz=0,left_ROI=0,right_ROI=0,counter_tracklet=0,percent=(int)kExp*0.01;
  double phi1=0,phi2=0,x1=0,x2=0,z1=0,z2=0,z_reco=0,z_reco_fin=0,center_ROI=0,z_event=0,total_reco=0,delta = 0.1;

  for(int i=0;i<kExp;i++){

    tree_gen->GetEvent(i);
    int mult_ev1=cross_L1->GetEntries();
    int mult_ev2=cross_L2->GetEntries();

    z_gen->GetEvent(i);

    if(PrintParticles){
      printf("> EVENT %i <\n\n%i hits with L1\n\n%i hits with L2\n\n",i+1,mult_ev1, mult_ev2);
    }else if(kExp>=100&&((i+1)%percent==0)){
      printf("\r[reconstruction running %3d%%]",100*(i+1)/kExp);
      fflush(stdout);
    }else if(kExp<100){
      printf("\r[reconstruction running]");
      fflush(stdout);
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
          z_reco=z1+((z2-z1)/(x1-x2))*x1;

          if(Abs(z_reco)<13.5){
            goodz++;
            h_zreco->Fill(z_reco);
            h_ROI->Fill(z_reco);
            h_tracklet->Fill(z_reco);
          }
        }
      }
    }

    /*new TCanvas();
    h_ROI->Draw();
    new TCanvas();
    h_tracklet->Draw();*/

    total_reco+=(double)goodz;

    if(goodz!=0){
      center_ROI = h_ROI->GetXaxis()->GetBinCenter(h_ROI->GetMaximumBin());//mm
      left_ROI=h_tracklet->FindBin(center_ROI-(delta/2));
      right_ROI=h_tracklet->FindBin(center_ROI+(delta/2));
      for(int k=left_ROI;k<right_ROI;k++){
        if(h_tracklet->GetBinContent(k)!=0){
          z_event+=(h_tracklet->GetXaxis()->GetBinCenter(k));
          counter_tracklet++;
        }
      }
      if(counter_tracklet!=0){
        z_event/=(double)counter_tracklet;
        h_reso->Fill(zgen-(float)z_event);
      }
    }
    h_ROI->Reset();
    h_tracklet->Reset();
    z_event=0;
    goodz=0;
    counter_tracklet=0;
  }

  printf("\n\n+++ END reconstruction +++\n\nSaving files...\n\n");

  TCanvas *c_zreco=new TCanvas("c_zreco","c_zreco",600,400);
  c_zreco->cd();
  TPaveText *pt_reco = new TPaveText(0.15,0.7,0.35,0.85,"NDC");
  pavestyler(*pt_reco,0.03);
  pt_reco->AddText(Form("%d events",kExp));
  pt_reco->AddText(Form("#mu = %f cm",h_zreco->GetMean(1)));
  pt_reco->AddText(Form("#sigma = %f cm",h_zreco->GetStdDev(1)));
  h_zreco->Draw();
  pt_reco->Draw();
  c_zreco->SaveAs(dirplot+"c_zreco.eps");

  TCanvas *c_reso=new TCanvas("c_reso","c_reso",600,400);
  c_reso->cd();
  h_reso->GetXaxis()->SetTitle("z_{gen} - z_{reco} [cm]");
  h_reso->GetYaxis()->SetTitle("# [a.u.]");
  h_reso->SetLineColor(kBlue+1);
  h_reso->Draw();
  TPaveText *pt_reso = new TPaveText(0.15,0.7,0.35,0.85,"NDC");
  pavestyler(*pt_reso,0.03);
  pt_reso->AddText(Form("%d events",kExp));
  pt_reso->AddText(Form("RMS = %f cm",h_reso->GetRMS()));
  pt_reso->Draw();
  c_reso->SaveAs(dirplot+"c_reso.eps");

  //cpu info
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double run_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/run_time)*100;

  perform.cputime=cpu_time;
  perform.runtime=run_time;
  perform.reso=h_reso->GetRMS();
  perform.eff=(double)kExp/total_reco;

  h_gen.Close();

  printf("CPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\nScroll up for info and verbosities. Thanks for using TRACKS!\n\n-> DONATE <-\n\n",cpu_time,run_time, cpu_efficiency);

  return perform;

} //END
