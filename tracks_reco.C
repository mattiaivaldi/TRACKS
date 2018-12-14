//
//TRACKS - reconstruction algorithm
//
//START

#include <vector>
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
#include "TStopwatch.h"
#include "Tools.h"
#include "Layer.h"
#include "Hit.h"
#include "Particle.h"

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif //needed to get current working directory

struct reco_perform{
  double reso;
  double e_reso;
  double eff;
  double e_eff;
};//an issue of reco_perform will contain resolution, efficiency and respective errors

using namespace TMath;

reco_perform tracks_reco(bool printparticles, bool printplot, double smear_z, double smear_phi, double amplitude, int width){

  reco_perform perform;

  //printparticles activates verbose mode (0,1)
  //printplot plot print (0,1)

  printf("\n\n+++ START reconstruction +++\n\n");

  TStopwatch timer;
  timer.Start(true);//start cpu monitor

  char cwd[FILENAME_MAX];
  GetCurrentDir(cwd,FILENAME_MAX);
  TString dirplot=TString(cwd)+"/tracksplot/";//define the path where the plots will be saved

  gRandom->SetSeed(0);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);

  //histogram booking

  TH1D *h_zreco=new TH1D("h_zreco","TRACKS reconstruction - Z Vertex;z_{V} [cm];# [a.u.]",200,-40.,40.);
  histostyler(*h_zreco,2);

  //TH1D *h_ROI=new TH1D("h_ROI","TRACKS reconstruction - ROI",801,-40.55,40.55);//old version with 1mm binning
  TH1D *h_ROI=new TH1D("h_ROI","TRACKS reconstruction - ROI",161,-40.25,40.25);

  TH1F *h_reso=new TH1F("h_reso","TRACKS reconstruction - resolution;z_{gen} - z_{reco} [cm];# [a.u.]",201,-0.1,0.1);
  histostyler(*h_reso,2);

  vector<double> tracklet;//reconstructed vertexes

  TFile f_gen("gen.root","READ");//accessing hit information
  TTree *tree_gen=(TTree*)f_gen.Get("TG");

  TNtuple *z_gen=(TNtuple*)f_gen.Get("z_gen");//accessing MC truth
  float zgen;
  z_gen->SetBranchAddress("z_gen",&zgen);

  int kExp=tree_gen->GetEntries();

  TClonesArray *cross_L1=new TClonesArray("Hit",100);
  TClonesArray *cross_L2=new TClonesArray("Hit",100);//TCA to retrieve hot coordinates

  TClonesArray& hits_L1=*cross_L1;
  TClonesArray& hits_L2=*cross_L2;//TCA aliases to be smeared

  TBranch *b1=tree_gen->GetBranch("HITL1");
  b1->SetAddress(&cross_L1);
  TBranch *b2=tree_gen->GetBranch("HITL2");
  b2->SetAddress(&cross_L2);//getting branches

  int goodz=0,percent=(int)kExp*0.01;
  double phi1=0,phi2=0,x1=0,x2=0,z1=0,z2=0,z_reco=0,z_reco_fin=0;
  double left_ROI=0,right_ROI=0,center_ROI=0,z_event=0,total_good=0,counter_tracklet=0;
  double diffz=0,delta=0.25;//work variables

  printf("Ready to reconstruct %d events!\n\n",kExp);

  for(int i=0;i<kExp;i++){//start loop over total events

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

    for(int l=0;l<mult_ev2;l++){//loop over total hits on layer #2

      Hit *hit_buffer2=(Hit*)cross_L2->At(l);
      smear(l,smear_z,smear_phi,7.,hits_L2);//smearing on hit point
      phi2=ACos(hit_buffer2->GetX()/7.);

      for(int m=0;m<mult_ev1;m++){//loop over total hits on layer #1

        Hit *hit_buffer1=(Hit*)cross_L1->At(m);
        smear(m,smear_z,smear_phi,4.,hits_L1);//smearing on hit point
        phi1=ACos(hit_buffer1->GetX()/4.);

        if(Abs(phi2-phi1)<0.01){//hit compatibility check

          x1=hit_buffer1->GetX();
          x2=hit_buffer2->GetX();
          z1=hit_buffer1->GetZ();
          z2=hit_buffer2->GetZ();
          z_reco=z1+((z2-z1)/(x1-x2))*x1;//tracklet intersection with z-axis

          goodz++;
          h_zreco->Fill(z_reco);
          h_ROI->Fill(z_reco);
          tracklet.push_back(z_reco);

        }
      }//end loop over total hits on layer #1
    }//loop over total hits on layer #2

    //un-comment to TCanvas-gun every roi histogram
    //new TCanvas();
    //h_ROI->DrawCopy();

    sort(tracklet.begin(),tracklet.end());

    if(goodz!=0&&peakfinder(h_ROI,amplitude,width)){//event ambiguity check
      total_good++;
      center_ROI=h_ROI->GetXaxis()->GetBinCenter(h_ROI->GetMaximumBin());
      left_ROI=center_ROI-delta;
      right_ROI=center_ROI+delta;//define the region of interest in cm
      for(int k=0;k<(int)tracklet.size();k++){//loop over tracklets
        if(tracklet[k]>=center_ROI-delta&&tracklet[k]<=center_ROI+delta){
          z_event+=tracklet[k];
          counter_tracklet++;
        }
      }//end loop over tracklets
      z_event/=counter_tracklet;//mean reconstructed z
      diffz=zgen-z_event;
      h_reso->Fill(diffz);//resolution histogram
    }

    h_ROI->Reset();//get the histogram ready for the next event
    vector<double>().swap(tracklet);//freed tracklet
    z_event=0;
    goodz=0;
    center_ROI=0;
    left_ROI=0;
    right_ROI=0;
    counter_tracklet=0;
    diffz=0;
  }

  printf("\n\n+++ END reconstruction +++\n\nThe plots are saved in the reco.root file and also in .eps format in the tracksplot folder.\n\n");

  f_gen.Close();

  TFile f_reco("reco.root","RECREATE");//to store reconstruction plots

  if(printplot){//draw and save plots

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

    h_zreco->Write();
    h_reso->Write();

    delete c_zreco;

    TCanvas *c_reco_perform=new TCanvas("c_reco_perform","c_reco_perform",600,400);
    c_reco_perform->cd();
    c_reco_perform->SetLogx();
    c_reco_perform->SetLogy();
    double exp[10]={100,500,1000,5000,10000,50000,100000,500000,1000000,5000000};
    double cput[10]={0.08,0.09,0.11,0.31,0.44,2.27,3.19,14.1,23.69,142.42};
    double runt[10]={0.15,0.16,0.18,0.46,0.52,2.40,3.31,14.21,23.79,142.59};

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

    delete c_reco_perform;

  }

  //cpu info
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double run_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/run_time)*100;

  perform.reso=h_reso->GetRMS();
  perform.e_reso=h_reso->GetRMSError();
  perform.eff=total_good/(double)kExp;
  perform.e_eff=Sqrt(perform.eff*(1-perform.eff)/(double)kExp);//calculate resolution and robustness

  f_reco.Close();

  delete h_zreco;//delete all objects
  delete h_ROI;
  delete h_reso;

  printf("Reconstruction info:\n\nCPU time = %f s\nRun time = %f s\nCPU efficiency = %f %%",cpu_time,run_time, cpu_efficiency);

  return perform;

} //END
