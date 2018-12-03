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
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "THStack.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TDirectory.h"
#include "vector"
#include "Tools.h"
#include "Layer.h"
#include "Hit.h"
#include "Particle.h"
#include<iostream>
#include<fstream>
#include<string>
#include<stdio.h>

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#define  MAX_LEN 30

bool war1=true;//declaring war

using namespace TMath;

void tracks_gen(bool printparticles, bool printplot, bool multiscatman, bool paolonoise, int custom, int kExp, double z_custom, int mult_custom) {

  //printparticles activates verbose mode (0,1)
  //printplot plot print (0,1)
  //multiscatman activates multiple scattering (0,1)
  //paolonoise activates the noise (0,1)
  //custom activates custom vertex z and/or multiplicity (10,5,15)
  //kExp is the number of collisions you want to perform
  //z_custom is the imposed vertex z
  //mult_custom is the imposed vertex multiplicity

  TStopwatch timer;
  timer.Start(true);//start cpu monitor

  char cwd[FILENAME_MAX];
  GetCurrentDir(cwd,FILENAME_MAX);
  TString dirplot=TString(cwd)+"/tracksplot/";

  gRandom->SetSeed(0);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetFillStyle(0);
  gStyle->SetLabelFont(42);

  char buffer_name[MAX_LEN+1];
  double W_buffer=0,R_buffer=0,T_buffer=0,RMS_buffer=0;
  TString layer_name;
  vector<Layer*> layer;

  FILE *stream_data;
  char isopen;
  stream_data = fopen("detector_info.txt","r");

  while(1){
    isopen = fgetc(stream_data);
    if(isopen==EOF){break;}
    else{
      fscanf(stream_data,"%s %lf %lf %lf %lf\n",&buffer_name [0],&W_buffer,&R_buffer,&T_buffer,&RMS_buffer);
      layer_name=buffer_name;
      layer.push_back(new Layer(layer_name,W_buffer,R_buffer,T_buffer,RMS_buffer));
    }
  }

  fclose (stream_data);

  double R_range=layer[2]->GetRadius(),width=layer[0]->GetWidth();

  TH1D *h_vgen=new TH1D("h_vgen","TRACKS generation - Z Vertex;z [cm];# [a.u.]",100,-width/2,width/2);
  histostyler(*h_vgen,4);
  TH1D *h_theta=new TH1D("h_theta","TRACKS generation - #theta;#theta [rad];# [a.u.]",100,0,Pi());
  histostyler(*h_theta,4);
  TH1D *h_phi=new TH1D("h_phi","TRACKS generation - #varphi;#varphi [rad];# [a.u.]",100,0,2*Pi());
  histostyler(*h_phi,4);
  TH1I *h_mult=new TH1I("h_mult","TRACKS generation - event multiplicity;multiplicity [a.u.];# [a.u.]",60,-0.5,59.5);
  histostyler(*h_mult,2);
  TH1I *h_rap=new TH1I("h_rap","TRACKS generation - pseudorapidity;#eta;# [a.u.]",100,-6,6);
  histostyler(*h_rap,2);

  THStack *s_hit_x = new THStack("s_hit_x","TRACKS generation - X hit;x [cm];# [a.u.]");
  THStack *s_hit_y = new THStack("s_hit_y","TRACKS generation - Y hit;y [cm];# [a.u.]");

  TH1D *h_BP[3], *h_L1[3], *h_L2[3];

  h_BP[0] = new TH1D("hbp_x","hbp_x",100,-R_range,R_range);
  h_BP[1] = new TH1D("hbp_y","hbp_y",100,-R_range,R_range);
  h_BP[2] = new TH1D("hbp_z","hbp_z",100,-width/2,width/2);

  h_L1[0] = new TH1D("hL1_x","hL1_x",100,-R_range,R_range);
  h_L1[1] = new TH1D("hL1_y","hL1_y",100,-R_range,R_range);
  h_L1[2] = new TH1D("hL1_z","hL1_z",100,-width/2,width/2);

  h_L2[0] = new TH1D("hL2_x","hL2_x",100,-R_range,R_range);
  h_L2[1] = new TH1D("hL2_y","hL2_y",100,-R_range,R_range);
  h_L2[2] = new TH1D("hL2_z","TRACKS generation - Z hit;z [cm];# [a.u.]",100,-width/2,width/2);
  histostyler(*h_L2[2],4);

  TH1D *h_Dphi=new TH1D("h_Dphi","TRACKS generation - #Delta#varphi;#Delta#varphi [rad];# [a.u.]",20,-0.015,0.015);
  histostyler(*h_Dphi,1);

  //verbosities
  verbosities(printparticles, multiscatman, paolonoise, kExp);

  TString distr="kinem.root";//get starting kinematics
  TFile hfile(distr);
  TH1F *pseudorap = (TH1F*)hfile.Get("heta");
  pseudorap->SetLineWidth(1);
  pseudorap->SetLineColor(kRed);
  TH1F *multiplicity = (TH1F*)hfile.Get("hmul");
  multiplicity->SetLineWidth(1);
  multiplicity->SetLineColor(kRed);

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

    Hit *vgen=new Hit(0, 0.001, 5.3, multiplicity);

    if(custom==10){
      vgen->SetZ(z_custom);
    }else if(custom==5){
      vgen->SetMult(mult_custom);
    }else if(custom==15){
      vgen->SetZ(z_custom);
      vgen->SetMult(mult_custom);
    }

    int mult = vgen->GetMult();
    z_gen->Fill((float)vgen->GetZ());
    h_vgen->Fill(vgen->GetZ());
    h_mult->Fill(mult);

    bool flag, flag1=false;
    int counter_BP=0,counter_L1=0,counter_L2=0;

    if(printparticles){
      printf("> EVENT %d <\n\nGenerated vertex with coordinates (%f, %f, %f)\nand multiplicity %d\n\n",i+1,vgen->GetX(),vgen->GetY(),vgen->GetZ(),mult);
    }else if(kExp>=100&&((i+1)%percent==0)){
      printf("\r[generation running %3d%%]",100*(i+1)/kExp);
      fflush(stdout);
    }else if(kExp<100){
      printf("\r[generation running]");
      fflush(stdout);
    }

    //start tracks generation

    for (int j=0; j<mult; j++) {

      Particle *part = new Particle(pseudorap);
      h_rap->Fill(part->GetRap());
      h_theta->Fill(part->GetTheta());
      h_phi->Fill(part->GetPhi());

      if (printparticles) {
        printf(">>> Particle %i: theta %f - phi %f <<<\n\n",j+1,part->GetTheta(),part->GetPhi());
      }

      //if particle hits layer TCA is filled, otherwhise gets 0
      flag=detect(vgen, layer[0], *part, hits_BP, printparticles, multiscatman, counter_BP,h_BP);
      phib=part->GetPhi();
      flag=detect(vgen, layer[1], *part, hits_L1, printparticles, multiscatman, counter_L1,h_L1);
      phia=part->GetPhi();
      if(flag){flag1=true;}
      flag=detect(vgen, layer[2], *part, hits_L2, printparticles, multiscatman, counter_L2,h_L2);
      if(flag&&flag1){h_Dphi->Fill(phia-phib);}

      delete part;

    }

    if(printparticles){
      printf("Out of %d generated particles\n%d crossed BP\n%d crossed L1\n%d crossed L2\n\n",mult,counter_BP, counter_L1, counter_L2);
    }

    //randomly add or not noise

    if(paolonoise){
      noise(printparticles,kNoise1,counter_L1,hits_L1, layer[1]);
      noise(printparticles,kNoise2,counter_L2,hits_L2, layer[2]);
    }

    tree_gen->Fill();

    cross_BP->Clear();
    cross_L1->Clear();
    cross_L2->Clear();

    delete vgen;

  }//end for up to kExp

  double Dphi_MAX=h_Dphi->GetBinCenter(h_Dphi->FindLastBinAbove());

  printf("\n\n+++ END generation +++\n\nSaving files...\n\nYou will find gen.root containing the detection info in the current directory.\n\n");

  if(printplot){

    TCanvas *c_kinem = new TCanvas("c_kinem","c_kinem",1200,400);
    c_kinem->Divide(2,1);
    c_kinem->cd(1);
    h_rap->SetLineColor(kBlue+1);
    h_rap->Draw();
    pseudorap->DrawCopy("SAME");
    auto legk = new TLegend(0.73,0.7,0.94,0.86);
    legk->SetHeader("10^{5} events","");
    legk->AddEntry(pseudorap,"theo","l");
    legk->AddEntry(h_rap,"extracted","l");
    legk->Draw();
    c_kinem->cd(2);
    gPad->SetLogy();
    h_mult->SetLineColor(kBlue+1);
    h_mult->Draw();
    multiplicity->DrawCopy("SAME");
    legk->Draw();
    c_kinem->SaveAs(dirplot+"c_kinem.eps");

    TCanvas *c_gen = new TCanvas("c_gen","c_gen",600,400);
    c_gen->Divide(2,2);
    c_gen->cd(1);
    h_vgen->SetLineColor(kBlue+1);
    h_vgen->Draw();
    TPaveText *pt_gen = new TPaveText(0.15,0.7,0.35,0.85,"NDC");
    pavestyler(*pt_gen,0.03);
    pt_gen->AddText(Form("%d events",kExp));
    pt_gen->AddText(Form("#mu = %f cm",h_vgen->GetMean(1)));
    pt_gen->AddText(Form("#sigma = %f cm",h_vgen->GetStdDev(1)));
    pt_gen->Draw();
    h_vgen->GetYaxis()->SetNdivisions(506);
    c_gen->cd(3);
    gPad->SetLogy();
    h_theta->SetLineColor(kBlue+1);
    h_theta->Draw();
    c_gen->cd(4);
    h_phi->SetLineColor(kBlue+1);
    h_phi->Draw();
    h_phi->SetAxisRange(h_phi->GetBinContent(h_phi->GetMaximumBin())+5000, h_phi->GetBinContent(h_phi->GetMaximumBin())-5000,"Y");
    h_phi->GetYaxis()->SetNdivisions(506);
    c_gen->SaveAs(dirplot+"c_gen.eps");

    TCanvas *c_dphi = new TCanvas("c_dphi","c_dphi",600,400);
    c_dphi->cd();
    h_Dphi->SetLineColor(kBlue+1);
    h_Dphi->Draw();
    TPaveText *pt_dphi = new TPaveText(0.15,0.7,0.35,0.85,"NDC");
    pavestyler(*pt_dphi,0.03);
    pt_dphi->AddText(Form("%d events",kExp));
    pt_dphi->AddText(Form("#Delta#varphi_{max} = %f rad",Dphi_MAX));
    pt_dphi->AddText("1 bin = 1 mrad");
    pt_dphi->Draw();
    c_dphi->SetLogy();
    c_dphi->SaveAs(dirplot+"c_dphi.eps");

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
        s_hit_x->GetXaxis()->SetLabelSize(0.045);
        s_hit_x->GetYaxis()->SetLabelSize(0.045);
        s_hit_x->GetXaxis()->SetTitleSize(0.05);
        s_hit_x->GetYaxis()->SetTitleSize(0.05);
        s_hit_x->GetXaxis()->SetTitleOffset(0.9);
        s_hit_x->GetYaxis()->SetTitleOffset(0.7);
        gPad->SetLogy();
      }
      else if(i==2){
        s_hit_y->Add(h_L2[i-1]);
        s_hit_y->Add(h_L1[i-1]);
        s_hit_y->Add(h_BP[i-1]);
        s_hit_y->Draw();
        s_hit_y->GetXaxis()->SetLabelSize(0.045);
        s_hit_y->GetYaxis()->SetLabelSize(0.045);
        s_hit_y->GetXaxis()->SetTitleSize(0.05);
        s_hit_y->GetYaxis()->SetTitleSize(0.05);
        s_hit_y->GetXaxis()->SetTitleOffset(0.9);
        s_hit_y->GetYaxis()->SetTitleOffset(0.7);
        gPad->SetLogy();
      }
      else{
        h_L2[i-1]->Draw();
        h_L1[i-1]->Draw("SAME");
        h_BP[i-1]->Draw("SAME");
        h_L2[i-1]->SetMaximum(h_BP[i-1]->GetMaximum()+5000);
        auto leg_hit = new TLegend(0.15,0.6,0.4,0.85);
        leg_hit->SetHeader(Form("%d events",kExp),"");
        leg_hit->AddEntry(h_BP[0],"BP","l");
        leg_hit->AddEntry(h_L1[0],"L1","l");
        leg_hit->AddEntry(h_L2[0],"L2","l");
        leg_hit->Draw();
      }
    }
    c_hit->cd(4);

    double exp[10]={100,500,1000,5000,10000,50000,100000,500000,1000000,5000000};
    double cput[10]={0.68,0.72,0.78,1.11,1.5,5.2,9.6,47,89,490};
    double runt[10]={0.71,0.76,0.82,1.16,1.6,5.3,9.9,49,92,510};
    double fsize[10]={0.035156,0.169922,0.353516,1.3,2.8,16,28,189,249,1945.6};

    TPad *pad = new TPad("pad","",0,0,1,1);
    pad->SetFillColor(0);
    pad->Draw();
    pad->SetLogx();
    pad->SetLogy();
    pad->cd();

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
    p_time->SetTitle("TRACKS generation - performances;# collisions;t [s]");
    p_time->GetXaxis()->SetLabelSize(0.045);
    p_time->GetYaxis()->SetLabelSize(0.045);
    p_time->GetXaxis()->SetTitleSize(0.05);
    p_time->GetYaxis()->SetTitleSize(0.05);
    p_time->GetYaxis()->SetTitleOffset(0.5);

    auto legt = new TLegend(0.15,0.65,0.3,0.85);
    legt->SetHeader("#varepsilon_{CPU} > 96%","");
    legt->AddEntry(p_cpu,"CPU time","l");
    legt->AddEntry(p_run,"RUN time","l");
    legt->Draw();

    TGraph *p_size=new TGraph(10,exp,fsize);
    p_size->SetLineColor(kGreen+2);
    p_size->SetLineWidth(2);
    p_size->Draw("SAMEL");

    TPad *overlay = new TPad("overlay","",0,0,1,1);
    overlay->SetFillStyle(4000);
    overlay->SetFillColor(0);
    overlay->SetFrameFillStyle(4000);
    overlay->Draw();
    overlay->cd();
    Double_t xmin = pad->GetUxmin();
    Double_t ymin = 0.25;
    Double_t xmax = pad->GetUxmax();
    Double_t ymax = pad->GetUymax();
    TH1F *hframe = overlay->DrawFrame(xmin,10,xmax,ymax);
    hframe->GetXaxis()->SetNdivisions(0);
    hframe->GetYaxis()->SetNdivisions(0);

    TGaxis *axis = new TGaxis(xmax,ymin,xmax,5.8,0.1,2000,505,"G+L");
    axis->SetLineColor(kGreen+2);
    axis->SetLabelColor(kGreen+2);
    axis->SetTitleColor(kGreen+2);
    axis->SetTitle("gen.root size [MB]");
    axis->SetLabelSize(0.045);
    axis->SetTitleSize(0.05);
    axis->Draw();

    c_hit->SaveAs(dirplot+"c_hit.eps");

  }

  h_gen.Write();
  h_gen.Close();

  //cpu info
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;
  printf("Generation info:\n\nCPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\nScroll up for info and verbosities. Thanks for using TRACKS!\n\n-> DONATE <-\n\n",cpu_time,real_time, cpu_efficiency);

} //END
