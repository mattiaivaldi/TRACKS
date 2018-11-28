//
//TRACKS - MC generation of particle tracks in a detector
//with multiple scattering and noise
//developed by Luca Quaglia and Mattia Ivaldi, 2018
//
//START

#include "TCanvas.h"
#include "TPad.h"
#include "TNtuple.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
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

bool war1=true;//declaring war

using namespace TMath;

void tracks_gen(bool PrintParticles, bool multiscatman, bool paolonoise, bool custom, int kExp, double z_custom, int mult_custom) {

  //PrintParticles activates verbose mode (0,1)
  //multiscatman activates multiple scattering (0,1)
  //paolonoise activates the noise (0,1)
  //custom activates custom vertex z and multiplicity (0,1)
  //kExp is the number of collisions you want to perform
  //z_custom is the imposed vertex z
  //mult_custom is the imposed vertex multiplicity

  TStopwatch timer;
  timer.Start(true);//start cpu monitor

  gRandom->SetSeed(0);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);

  //length, radius, thickness, multiscattering RMS
  Layer *BP = new Layer(27.,3.,0.8,0.001);
  Layer *L1 = new Layer(27.,4.,0.2,0.001);
  Layer *L2 = new Layer(27.,7.,0.2,0.001);

  double R_range=L2->GetRadius();

  THStack *s_hit_x = new THStack("s_hit_x","TRACKS generation - X hit");
  THStack *s_hit_y = new THStack("s_hit_y","TRACKS generation - Y hit");
  TH1D *h_BP[3], *h_L1[3], *h_L2[3];

  h_BP[0] = new TH1D("hbp_x","hbp_x",100,-R_range,R_range);
  h_BP[1] = new TH1D("hbp_y","hbp_y",100,-R_range,R_range);
  h_BP[2] = new TH1D("hbp_z","hbp_z",100,-BP->GetWidth()/2,BP->GetWidth()/2);

  h_L1[0] = new TH1D("hL1_x","hL1_x",100,-R_range,R_range);
  h_L1[1] = new TH1D("hL1_y","hL1_y",100,-R_range,R_range);
  h_L1[2] = new TH1D("hL1_z","hL1_z",100,-L1->GetWidth()/2,L1->GetWidth()/2);

  h_L2[0] = new TH1D("hL2_x","TRACKS generation - X hit",100,-R_range,R_range);
  h_L2[1] = new TH1D("hL2_y","TRACKS generation - Y hit",100,-R_range,R_range);
  h_L2[2] = new TH1D("hL2_z","TRACKS generation - Z hit",100,-L2->GetWidth()/2,L2->GetWidth()/2);

  TH1D *h_vgen=new TH1D("h_vgen","TRACKS generation - Z Vertex",100,-30,30);
  h_vgen->GetXaxis()->SetTitle("z [cm]");
  h_vgen->GetYaxis()->SetTitle("# [a.u.]");
  TH1I *h_mult=new TH1I("h_mult","TRACKS generation - event multiplicity",60,0,60);
  h_mult->GetXaxis()->SetTitle("multiplicity [a.u.]");
  h_mult->GetYaxis()->SetTitle("# [a.u.]");

  TH1D *h_Dphi=new TH1D("h_Dphi","TRACKS generation - #Delta#phi",20,-0.015,0.015);
  h_Dphi->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h_Dphi->GetYaxis()->SetTitle("# [a.u.]");

  //verbosities
  verbosities(PrintParticles, multiscatman, paolonoise, kExp);

  TString distr="kinem.root";//get starting kinematics

  TFile hfile(distr);
  TH1F *pseudorap = (TH1F*)hfile.Get("heta");
  TH1F *multiplicity = (TH1F*)hfile.Get("hmul");//get multiplicity distribution

  TFile h_gen("gen.root","RECREATE");
  TTree *tree_gen=new TTree("TG","tree_gen");
  TNtuple *z_gen=new TNtuple("z_gen","z_gen","z_gen");//generation data

  int kNoise1=0, kNoise2=0;
  if(paolonoise){
    kNoise1=(int)gRandom->Integer(5);
    kNoise2=(int)gRandom->Integer(5);
  }//number of spurious hits

  int size0=multiplicity->FindLastBinAbove(0,1);
  int size1=multiplicity->FindLastBinAbove(0,1)+kNoise1;
  int size2=multiplicity->FindLastBinAbove(0,1)+kNoise2;

  TClonesArray *cross_BP=new TClonesArray("Hit",size0);
  TClonesArray *cross_L1=new TClonesArray("Hit",size1);
  TClonesArray *cross_L2=new TClonesArray("Hit",size2);//TCA booking

  TClonesArray& hits_BP=*cross_BP;
  TClonesArray& hits_L1=*cross_L1;
  TClonesArray& hits_L2=*cross_L2;//TCA aliases to be filled

  tree_gen->Branch("HITL1",&cross_L1);
  tree_gen->Branch("HITL2",&cross_L2);//branch booking

  //start experiment

  int percent=(int)kExp*0.01;
  double phib=0,phia=0;

  for(int i=0; i<kExp; i++){

    Hit *vgen;

    if(custom){
      vgen=new Hit(0, 0.001, z_custom, mult_custom);
    }else{vgen=new Hit(0, 0.001, 5.3, multiplicity);}

    int mult = vgen->GetMult();
    z_gen->Fill((float)vgen->GetZ());
    h_vgen->Fill(vgen->GetZ());
    h_mult->Fill(mult);

    bool flag, flag1=false;
    int counter_BP=0,counter_L1=0,counter_L2=0;

    if(PrintParticles){
      printf("> EVENT %d <\n\nGenerated vertex with coordinates (%f, %f, %f)\nand multiplicity %d\n\n",i+1,vgen->GetX(),vgen->GetY(),vgen->GetZ(),mult);
    }else if((i+1)%percent==0){
      printf("\r[generation running %3d%%]",100*(i+1)/kExp);
      fflush(stdout);
    }

    //start tracks generation

    for (int j=0; j<mult; j++) {

      Particle *part = new Particle(pseudorap);

      if (PrintParticles) {
        printf(">>> Particle %i: theta %f - phi %f <<<\n\n",j+1,part->GetTheta(),part->GetPhi());
      }

      //if particle hits layer TCA is filled, otherwhise gets 0
      flag=detect(vgen, BP, *part, hits_BP, PrintParticles, multiscatman, "BP", counter_BP,h_BP);
      phib=part->GetPhi();
      flag=detect(vgen, L1, *part, hits_L1, PrintParticles, multiscatman, "L1", counter_L1,h_L1);
      phia=part->GetPhi();
      if(flag){flag1=true;}
      flag=detect(vgen, L2, *part, hits_L2, PrintParticles, multiscatman, "L2", counter_L2,h_L2);
      if(flag&&flag1){h_Dphi->Fill(phia-phib);}

      delete part;

    }

    if(PrintParticles){
      printf("Out of %d generated particles\n%d crossed BP\n%d crossed L1\n%d crossed L2\n\n",mult,counter_BP, counter_L1, counter_L2);
    }

    //randomly add or not noise

    if(paolonoise){
      noise(PrintParticles,kNoise1,counter_L1,hits_L1, L1,"L1");
      noise(PrintParticles,kNoise2,counter_L2,hits_L2, L2,"L2");
    }

    tree_gen->Fill();

    cross_BP->Clear();
    cross_L1->Clear();
    cross_L2->Clear();

    delete vgen;

  }//end for up to kExp

  double Dphi_MAX=h_Dphi->GetBinCenter(h_Dphi->FindLastBinAbove());

  printf("\n\n+++ END generation +++\n\nSaving files...\n\nYou will find gen.root containing the detection info in the current directory.\n\n");

  TCanvas *c_gen = new TCanvas("c_gen","c_gen",1200,400);
  c_gen->Divide(2,1);
  c_gen->cd(1);
  h_vgen->SetLineColor(kBlue+1);
  h_vgen->Draw();
  c_gen->cd(2);
  h_mult->SetLineColor(kBlue+1);
  h_mult->Draw();
  c_gen->SaveAs("c_gen.eps");

  TCanvas *c_dphi = new TCanvas("c_dphi","c_dphi",600,400);
  c_dphi->cd();
  h_Dphi->SetLineColor(kBlue+1);
  h_Dphi->Draw();
  TPaveText *pt_dphi = new TPaveText(0.15,0.7,0.35,0.85,"NDC");
  pavestyler(*pt_dphi,0.03);
  pt_dphi->AddText(Form("%d events",kExp));
  pt_dphi->AddText(Form("#Delta#phi_{max} = %f rad",Dphi_MAX));
  pt_dphi->AddText("1 bin = 1 mrad");
  pt_dphi->Draw();
  c_dphi->SetLogy();
  c_dphi->SaveAs("c_dphi.eps");

  TCanvas *c_hit = new TCanvas("c_hit","c_hit",600,400);
  c_hit->Divide(2,2);
  for(int i=1;i<=3;i++){
    c_hit->cd(i);
    h_L2[i-1]->SetLineColor(kGreen+3);
    h_L2[i-1]->SetFillColor(kWhite);
    h_L2[i-1]->SetLineWidth(1);
    h_L1[i-1]->SetLineColor(kRed);
    h_L1[i-1]->SetFillColor(kWhite);
    h_L1[i-1]->SetLineWidth(1);
    h_BP[i-1]->SetLineColor(kBlue+1);
    h_BP[i-1]->SetFillColor(kWhite);
    h_BP[i-1]->SetLineWidth(1);
    if(i==1){
      s_hit_x->Add(h_L2[i-1]);
      s_hit_x->Add(h_L1[i-1]);
      s_hit_x->Add(h_BP[i-1]);
      s_hit_x->Draw();
      s_hit_x->GetXaxis()->SetTitle("x [cm]");
      s_hit_x->GetYaxis()->SetTitle("# [a.u.]");
      gPad->SetLogy();
    }
    else if(i==2){
      s_hit_y->Add(h_L2[i-1]);
      s_hit_y->Add(h_L1[i-1]);
      s_hit_y->Add(h_BP[i-1]);
      s_hit_y->Draw();
      s_hit_y->GetXaxis()->SetTitle("y [cm]");
      s_hit_y->GetYaxis()->SetTitle("# [a.u.]");
      gPad->SetLogy();
    }
    else{
      h_L2[i-1]->Draw();
      h_L1[i-1]->Draw("SAME");
      h_BP[i-1]->Draw("SAME");
      h_L2[i-1]->SetMaximum(h_BP[i-1]->GetMaximum()+5000);
    }
  }
  c_hit->cd(4);
  auto leg_hit = new TLegend(0.15,0.5,0.5,0.85);
  leg_hit->SetHeader(Form("%d events",kExp),"");
  leg_hit->AddEntry(h_BP[0],"BP","l");
  leg_hit->AddEntry(h_L1[0],"L1","l");
  leg_hit->AddEntry(h_L2[0],"L2","l");
  leg_hit->Draw();
  c_hit->SaveAs("c_hit.eps");

  h_gen.Write();
  h_gen.Close();

  //cpu info
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;
  printf("CPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\nScroll up for info and verbosities. Thanks for using TRACKS!\n\n-> DONATE <-\n\n",cpu_time,real_time, cpu_efficiency);

} //END
