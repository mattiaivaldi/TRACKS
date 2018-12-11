#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

void cluster_study(){

  gROOT->Reset();

  TStopwatch timer;
  timer.Start(true);//start cpu monitor

  char cwd[FILENAME_MAX];
  GetCurrentDir(cwd,FILENAME_MAX);
  TString dirplot=TString(cwd)+"/tracksplot/";

  const int kTest=6, kCluster=4;

  double ampli_custom[kCluster]={0.7,1,1.3,1.5};
  double width_custom[kCluster]={1,3,5,7};

  double z_custom[kTest]={0,2,4,7,10,13};
  double mult_custom[kTest]={5,15,30,40,50};

  double resoz_ampli[kCluster][kTest], e_resoz_ampli[kCluster][kTest], resoz_width[kCluster][kTest], e_resoz_width[kCluster][kTest];

  TGraphErrors *p_resoz_ampli[kCluster],*p_resoz_width[kCluster];

  double reso_buffer[kTest], e_reso_buffer[kTest];

  TString z, m, ampli, width, exec;

  printf("+++ START resolution performances +++\n\n");

  for(int i=0; i<kCluster;i++){
    for(int j=0; j<kTest;j++){
      z=Form("%f",z_custom[j]);
      exec="tracks_gen(0,0,1,1,15,100000,"+z+",20)";
      gROOT->ProcessLine(exec);
      reco_perform perform=tracks_reco(0,0,0.0012,0.0003,ampli_custom[i],3);
      resoz_ampli[i][j]=perform.reso;
      e_resoz_ampli[i][j]=perform.e_reso;
      gROOT->Reset();
    }
  }

  for(int i=0; i<kCluster;i++){
    for(int j=0; j<kTest;j++){
      z=Form("%f",z_custom[j]);
      exec="tracks_gen(0,0,1,1,15,100000,"+z+",20)";
      gROOT->ProcessLine(exec);
      reco_perform perform=tracks_reco(0,0,0.0012,0.0003,1,width_custom[i]);
      resoz_width[i][j]=perform.reso;
      e_resoz_width[i][j]=perform.e_reso;
      gROOT->Reset();
    }
  }

  TCanvas *c_study=new TCanvas("c_study","c_study",600,400);
  c_study->Divide(2,2);

  c_study->cd(1);
  TMultiGraph *m_resoz_ampli=new TMultiGraph;
  for(int i=0;i<kCluster;i++){//over ampli custom
    for(int j=0;j<kTest;j++){//over z custom
      reso_buffer[j]=resoz_ampli[i][j];
      e_reso_buffer[j]=e_resoz_ampli[i][j];
    }
    p_resoz_ampli[i]=new TGraphErrors(kTest,z_custom,reso_buffer,NULL,e_reso_buffer);
    p_resoz_ampli[i]->SetMarkerStyle(20+i);
    p_resoz_ampli[i]->SetMarkerSize(0.4);
    p_resoz_ampli[i]->SetMarkerColor(kRed+i);
    p_resoz_ampli[i]->SetLineColor(kRed+i);
    m_resoz_ampli->Add(p_resoz_ampli[i]);
  }
  m_resoz_ampli->Draw("AP");
  m_resoz_ampli->SetTitle("TRACKS performances - resolution vs Z_{gen};z_{gen} [cm];RMS [cm]");

  c_study->cd(2);
  TMultiGraph *m_resoz_width=new TMultiGraph;
  for(int i=0;i<kCluster;i++){//over ampli custom
    for(int j=0;j<kTest;j++){//over z custom
      reso_buffer[j]=resoz_width[i][j];
      e_reso_buffer[j]=e_resoz_width[i][j];
    }
    p_resoz_width[i]=new TGraphErrors(kTest,z_custom,reso_buffer,NULL,e_reso_buffer);
    p_resoz_width[i]->SetMarkerStyle(20+i);
    p_resoz_width[i]->SetMarkerSize(0.4);
    p_resoz_width[i]->SetMarkerColor(kRed+i);
    p_resoz_width[i]->SetLineColor(kRed+i);
    m_resoz_width->Add(p_resoz_width[i]);
  }
  m_resoz_width->Draw("AP");
  m_resoz_width->SetTitle("TRACKS performances - resolution vs Z_{gen};z_{gen} [cm];RMS [cm]");

  /*TGraphErrors *p_resoz=new TGraphErrors(kTest,z_custom,resoz,NULL,e_resoz);
  p_resoz->SetTitle("TRACKS performances - resolution vs Z_{gen};z_{gen} [cm];RMS [cm]");
  graphstyler(*p_resoz,4);
  p_resoz->GetYaxis()->SetTitleOffset(1.1);
  p_resoz->SetMarkerStyle(20);
  p_resoz->SetMarkerSize(0.4);
  p_resoz->SetMarkerColor(1);
  p_resoz->Draw("AP");*/

  /*c_study->cd(2);
  TGraphErrors *p_resom=new TGraphErrors(kTest,mult_custom,resom,NULL,e_resom);
  p_resom->SetTitle("TRACKS performances - resolution vs multiplicity;multiplicity;RMS [cm]");
  graphstyler(*p_resom,4);
  p_resom->GetYaxis()->SetTitleOffset(1.1);
  p_resom->SetMarkerStyle(20);
  p_resom->SetMarkerSize(0.4);
  p_resom->SetMarkerColor(1);
  p_resom->Draw("AP");

  c_study->cd(3);
  TGraphErrors *p_effz=new TGraphErrors(kTest,z_custom,effz,NULL,e_effz);
  p_effz->SetTitle("TRACKS performances - robustness vs Z_{gen};z_{gen} [cm];#tilde{#varepsilon}");
  graphstyler(*p_effz,4);
  p_effz->GetYaxis()->SetTitleOffset(0.5);
  p_effz->GetYaxis()->SetTitleSize(0.07);
  p_effz->SetMarkerStyle(20);
  p_effz->SetMarkerSize(0.4);
  p_effz->SetMarkerColor(1);
  p_effz->Draw("AP");

  c_study->cd(4);
  TGraphErrors *p_effm=new TGraphErrors(kTest,mult_custom,effm,NULL,e_effm);
  p_effm->SetMarkerStyle(20);
  p_effm->SetMarkerSize(0.4);
  p_effm->SetMarkerColor(1);
  TGraphErrors *p_effm1s=new TGraphErrors(kTest,mult_custom,effm1s,NULL,e_effm1s);
  p_effm1s->SetMarkerStyle(22);
  p_effm1s->SetMarkerSize(0.35);
  p_effm1s->SetMarkerColor(2);
  TMultiGraph *m_effm=new TMultiGraph();
  m_effm->Add(p_effm);
  m_effm->Add(p_effm1s);
  m_effm->Draw("AP");
  m_effm->SetTitle("TRACKS performances - robustness vs multiplicity;multiplicity;#tilde{#varepsilon}");
  graphstyler(*p_effz,4);
  m_effm->GetYaxis()->SetTitleOffset(0.5);
  m_effm->GetYaxis()->SetTitleSize(0.07);*/

  c_study->SaveAs(dirplot+"c_study.eps");

  /*for(int i=0; i<kTest;i++){
    m=Form("%f",mult_custom[i]);
    exec="tracks_gen(0,0,1,1,5,1000000,0,"+m+")";
    gROOT->ProcessLine(exec);
    reco_perform perform=tracks_reco(0,0,0.0012,0.0003);
    effz[i]=perform.eff;
    e_effz[i]=TMath::Sqrt(effz[i]*(1-effz[i])/perform.kExp);
    gROOT->Reset();
  }*/

  //cpu info
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;
  printf("Performances info:\n\nCPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\nScroll up for info and verbosities. Thanks for using TRACKS!\n\n-> DONATE <-\n\n",cpu_time,real_time, cpu_efficiency);

}
