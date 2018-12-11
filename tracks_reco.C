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
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "Tools.h"
#include "Layer.h"
#include "Hit.h"
#include "Particle.h"
#include <vector>
#include <deque>

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

struct reco_perform{
  double cputime;
  double runtime;
  double reso;
  double e_reso;
  double eff;
  double e_eff;
  //double kExp;
};

bool war2=true;//declaring war

using namespace TMath;
using namespace std;

reco_perform tracks_reco(bool printparticles, bool printplot, double smear_z, double smear_phi, double amplitude, int width){

  reco_perform perform;

  //printparticles activates verbose mode (0,1)
  //printplot plot print (0,1)

  printf("\n\n+++ START reconstruction +++\n\n");

  TStopwatch timer;
  timer.Start(true);//start cpu monitor

  char cwd[FILENAME_MAX];
  GetCurrentDir(cwd,FILENAME_MAX);
  TString dirplot=TString(cwd)+"/tracksplot/";

  gRandom->SetSeed(0);
  gStyle->SetOptStat(1111);
  gStyle->SetLegendBorderSize(0);

  TH1D *h_zreco=new TH1D("h_zreco","TRACKS reconstruction - Z Vertex;z_{V} [cm];# [a.u.]",200,-40.,40.);
  //TH1D *h_zreco=new TH1D("h_zreco","TRACKS reconstruction - Z Vertex;z_{V} [cm];# [a.u.]",200,-13.5,13.5);
  histostyler(*h_zreco,2);

  //TH1D *h_ROI=new TH1D("h_ROI","TRACKS reconstruction - ROI",270,-13.5,13.5);
  TH1D *h_ROI=new TH1D("h_ROI","TRACKS reconstruction - ROI",801,-40.,40.);

  //TH1D *h_tracklet=new TH1D("h_tracklet","TRACKS reconstruction - tracklet",800001,-40.,40.);
  //TH1D *h_tracklet=new TH1D("h_tracklet","TRACKS reconstruction - tracklet",270000,-13.5,13.5);

  vector<double> tracklet;

  TH1F *h_reso=new TH1F("h_reso","TRACKS reconstruction - resolution;z_{gen} - z_{reco} [cm];# [a.u.]",201,-0.1005,0.1005);
  histostyler(*h_reso,2);

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

  bool startsum=false;
  int goodz=0;
  //int left_ROI=0,right_ROI=0;
  int percent=(int)kExp*0.01;
  double phi1=0,phi2=0,x1=0,x2=0,z1=0,z2=0,z_reco=0,z_reco_fin=0,center_ROI=0,z_event=0,total_reco=0, total_good=0,delta = 0.05;

  double left_ROI=0, right_ROI=0,counter_tracklet=0,diffz=0;

  for(int i=0;i<kExp;i++){

    tree_gen->GetEvent(i);
    int mult_ev1=cross_L1->GetEntries();
    int mult_ev2=cross_L2->GetEntries();

    z_gen->GetEvent(i);

    if(printparticles){
      printf("> EVENT %i <\n\n%i hits with L1\n\n%i hits with L2\n\n",i+1,mult_ev1, mult_ev2);
    }else if(kExp>=100&&((i+1)%percent==0)){
      printf("\r[reconstruction running %3d%%]",100*(i+1)/kExp);
      fflush(stdout);
    }else if(kExp<100){
      printf("\r[reconstruction running]");
      fflush(stdout);
    }

    printf("\nlist of z reco produced\n");

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

          goodz++;
          h_zreco->Fill(z_reco);
          h_ROI->Fill(z_reco);
          //h_tracklet->Fill(z_reco);
          tracklet.push_back(z_reco);
          printf("\nevento %d %f %f #reco %d",i,zgen,z_reco,goodz);

        }
      }
    }

    sort(tracklet.begin(),tracklet.end());
    cout<<endl;

    printf("\nnow z reco are ordered\n");

    if(goodz!=0&&peakfinder(h_ROI,amplitude,width)){
      total_good++;
      center_ROI=h_ROI->GetXaxis()->GetBinCenter(h_ROI->GetMaximumBin());//cm
      left_ROI=center_ROI-delta;
      right_ROI=center_ROI+delta;
      printf("\nevento %d left %f center %f right %f\n\n",i,left_ROI,center_ROI,right_ROI);
      for(int k=0;k<tracklet.size();k++){
        //printf("\nevento %d %f %f",i,zgen,tracklet[k]);
        if(tracklet[k]>=center_ROI-delta&&tracklet[k]<=center_ROI+delta){
          printf("%d %f\n",k+1,tracklet[k]);
          z_event+=tracklet[k];
          counter_tracklet++;
        }
      }
      z_event/=counter_tracklet;
    }

    printf("\nevento %d %f %f\n",i,zgen,z_event);

    diffz=zgen-z_event;
    h_reso->Fill(diffz);

    /*if(goodz!=0){
      total_reco++;
      if(peakfinder(h_ROI,amplitude,width)){
        total_good++;
        center_ROI=h_ROI->GetXaxis()->GetBinCenter(h_ROI->GetMaximumBin());//mm
        //left_ROI=h_tracklet->FindBin(center_ROI-(delta));
        //right_ROI=h_tracklet->FindBin(center_ROI+(delta));
        for(int k=0;k<tracklet.size();k++){
          printf("\nevento %d %f %f %f",i,zgen,tracklet[k],center_ROI);
          if(tracklet[k]>=center_ROI-delta&&tracklet[k]<=center_ROI+delta){
            z_event+=tracklet[k];
            counter_tracklet++;
          }
        }*/
        /*for(int k=left_ROI;k<right_ROI;k++){
          if(h_tracklet->GetBinContent(k)!=0){
            z_event+=(h_tracklet->GetXaxis()->GetBinCenter(k));
            counter_tracklet++;
          }
        }*/
        /*if(counter_tracklet!=0){
          z_event/=(double)counter_tracklet;
          h_reso->Fill(zgen-(float)z_event);
        }
      }
    }

    for(int p=0;p<tracklet.size();p++){
      printf("\nevento %d %f %f %f %f",i,zgen,tracklet[p],center_ROI,z_event);
    }*/

    h_ROI->Reset();
    //h_tracklet->Reset();
    vector<double>().swap(tracklet);
    z_event=0;
    goodz=0;
    left_ROI=0;
    right_ROI=0;
    counter_tracklet=0;
    diffz=0;
  }

  printf("\n\n+++ END reconstruction +++\n\nSaving files...\n\n");

  if(printplot){

    TCanvas *c_zreco=new TCanvas("c_zreco","c_zreco",1200,400);
    c_zreco->Divide(2,1);
    c_zreco->cd(1);
    h_zreco->DrawCopy();
    TPaveText *pt_reco = new TPaveText(0.15,0.7,0.35,0.85,"NDC");
    pavestyler(*pt_reco,0.03);
    pt_reco->AddText(Form("%d events",kExp));
    pt_reco->AddText(Form("#mu = %f cm",h_zreco->GetMean(1)));
    pt_reco->AddText(Form("#sigma = %f cm",h_zreco->GetStdDev(1)));
    pt_reco->Draw();
    c_zreco->cd(2);
    h_reso->DrawCopy();
    TPaveText *pt_reso = new TPaveText(0.15,0.7,0.35,0.85,"NDC");
    pavestyler(*pt_reso,0.03);
    pt_reso->AddText(Form("%d events",kExp));
    pt_reso->AddText(Form("RMS = %f cm",h_reso->GetRMS()));
    pt_reso->Draw();
    c_zreco->SaveAs(dirplot+"c_zreco.eps");

    TCanvas *c_reco_perform=new TCanvas("c_reco_perform","c_reco_perform",600,400);
    c_reco_perform->cd();
    c_reco_perform->SetLogx();
    c_reco_perform->SetLogy();
    double exp[10]={100,500,1000,5000,10000,50000,100000,500000,1000000,5000000};
    double cput[10]={0.25,0.33,0.37,0.88,1.6,6.01,14,64,124,680};
    double runt[10]={0.26,0.35,0.39,0.92,1.65,6.08,14.1,66,125,694};

    TGraph *p_cpu=new TGraph(10,exp,cput);
    p_cpu->SetLineColor(kBlue);
    p_cpu->SetLineWidth(2);
    TGraph *p_run=new TGraph(10,exp,runt);
    p_run->SetLineColor(kRed);
    p_run->SetLineWidth(2);

    TMultiGraph *p_time=new TMultiGraph();
    p_time->Add(p_cpu);
    p_time->Add(p_run);
    p_time->Draw("AL");
    p_time->SetTitle("TRACKS reconstruction - performances;# collisions;t [s]");
    p_time->GetXaxis()->SetTitleSize(0.04);
    p_time->GetYaxis()->SetTitleSize(0.035);
    p_time->GetXaxis()->SetTitleOffset(1.1);

    auto legt = new TLegend(0.15,0.65,0.3,0.85);
    legt->SetHeader("#varepsilon_{CPU} > 93%","");
    legt->AddEntry(p_cpu,"CPU time","l");
    legt->AddEntry(p_run,"RUN time","l");
    legt->Draw();

    c_reco_perform->SaveAs(dirplot+"c_reco_perform.eps");

  }

  //cpu info
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double run_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/run_time)*100;

  perform.cputime=cpu_time;
  perform.runtime=run_time;
  perform.reso=h_reso->GetRMS();
  perform.e_reso=h_reso->GetRMSError();
  perform.eff=total_good/(double)kExp;
  perform.e_eff=Sqrt(perform.eff*(1-perform.eff)/(double)kExp);

  h_gen.Close();

  delete h_zreco;
  delete h_ROI;
  //delete h_tracklet;
  delete h_reso;

  printf("Reconstruction info:\n\nCPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\nScroll up for info and verbosities. Thanks for using TRACKS!\n\n-> DONATE <-\n\n",cpu_time,run_time, cpu_efficiency);

  //printf("\n%d  %f  %f\n\n",goodz, total_reco, total_good);

  //printf("\n\n\n%f  %f\n\n\n",total_reco/(double)kExp, perform.eff);

  printf("\n\n%f\n\n",total_good);

  return perform;

} //END
