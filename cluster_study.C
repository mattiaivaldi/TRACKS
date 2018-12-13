#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

void cluster_study(){

  gROOT->Reset();
  gStyle->SetOptStat(0);

  TStopwatch timer;
  timer.Start(true);//start cpu monitor

  char cwd[FILENAME_MAX];
  GetCurrentDir(cwd,FILENAME_MAX);
  TString dirplot=TString(cwd)+"/tracksplot/";

  const int kTest=6, kCluster=4;

  int color[kCluster]={2,4,6,9};

  double ampli_custom[kCluster]={0.7,1,1.3,1.5};
  double width_custom[kCluster]={1,3,5,7};

  double z_custom[kTest]={0,2,4,7,10,13};
  double mult_custom[kTest]={3,5,15,30,40,50};

  double resoz_ampli[kCluster][kTest], e_resoz_ampli[kCluster][kTest], resoz_width[kCluster][kTest], e_resoz_width[kCluster][kTest], resom_ampli[kCluster][kTest], e_resom_ampli[kCluster][kTest], resom_width[kCluster][kTest], e_resom_width[kCluster][kTest];

  TGraphErrors *p_resoz_ampli[kCluster],*p_resoz_width[kCluster], *p_resom_ampli[kCluster],*p_resom_width[kCluster];

  double reso_buffer[kTest], e_reso_buffer[kTest];

  TString z, m, ampli, width, exec;

  printf("+++ START resolution performances +++\n\n");

  for(int i=0; i<kCluster;i++){
    for(int j=0; j<kTest;j++){
      z=Form("%f",z_custom[j]);
      exec="tracks_gen(0,0,1,1,15,1000000,"+z+",20)";
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
      exec="tracks_gen(0,0,1,1,15,1000000,"+z+",20)";
      gROOT->ProcessLine(exec);
      reco_perform perform=tracks_reco(0,0,0.0012,0.0003,1,width_custom[i]);
      resoz_width[i][j]=perform.reso;
      e_resoz_width[i][j]=perform.e_reso;
      gROOT->Reset();
    }
  }

  for(int i=0; i<kCluster;i++){
    for(int j=0; j<kTest;j++){
      m=Form("%f",mult_custom[j]);
      exec="tracks_gen(0,0,1,1,15,1000000,0,"+m+")";
      gROOT->ProcessLine(exec);
      reco_perform perform=tracks_reco(0,0,0.0012,0.0003,ampli_custom[i],3);
      resom_ampli[i][j]=perform.reso;
      e_resom_ampli[i][j]=perform.e_reso;
      gROOT->Reset();
    }
  }//resolution vs multiplicity for different amplitude

  for(int i=0; i<kCluster;i++){
    for(int j=0; j<kTest;j++){
      m=Form("%f",mult_custom[j]);
      exec="tracks_gen(0,0,1,1,15,1000000,0,"+m+")";
      gROOT->ProcessLine(exec);
      reco_perform perform=tracks_reco(0,0,0.0012,0.0003,1,width_custom[i]);
      resom_width[i][j]=perform.reso;
      e_resom_width[i][j]=perform.e_reso;
      gROOT->Reset();
    }
  }//resolution vs multiplicity for different width

  TCanvas *c_study=new TCanvas("c_study","c_study",600,400);
  c_study->Divide(2,2);

  c_study->cd(1);
  TMultiGraph *m_resoz_ampli=new TMultiGraph();
  for(int i=0;i<kCluster;i++){//over ampli custom
    for(int j=0;j<kTest;j++){//over z custom
      reso_buffer[j]=resoz_ampli[i][j];
      e_reso_buffer[j]=e_resoz_ampli[i][j];
    }
    p_resoz_ampli[i]=new TGraphErrors(kTest,z_custom,reso_buffer,NULL,e_reso_buffer);
    p_resoz_ampli[i]->SetMarkerStyle(20+i);
    p_resoz_ampli[i]->SetMarkerSize(0.4);
    p_resoz_ampli[i]->SetMarkerColor(color[i]);
    p_resoz_ampli[i]->SetLineColor(color[i]);
    m_resoz_ampli->Add(p_resoz_ampli[i]);
  }
  m_resoz_ampli->Draw("AP");
  m_resoz_ampli->SetTitle("TRACKS performances - resolution vs Z_{gen};z_{gen} [cm];RMS [cm]");

  c_study->cd(2);
  TMultiGraph *m_resoz_width=new TMultiGraph();
  for(int i=0;i<kCluster;i++){//over ampli custom
    for(int j=0;j<kTest;j++){//over z custom
      reso_buffer[j]=resoz_width[i][j];
      e_reso_buffer[j]=e_resoz_width[i][j];
    }
    p_resoz_width[i]=new TGraphErrors(kTest,z_custom,reso_buffer,NULL,e_reso_buffer);
    p_resoz_width[i]->SetMarkerStyle(20+i);
    p_resoz_width[i]->SetMarkerSize(0.4);
    p_resoz_width[i]->SetMarkerColor(color[i]);
    p_resoz_width[i]->SetLineColor(color[i]);
    m_resoz_width->Add(p_resoz_width[i]);
  }
  m_resoz_width->Draw("AP");
  m_resoz_width->SetTitle("TRACKS performances - resolution vs Z_{gen};z_{gen} [cm];RMS [cm]");

  c_study->cd(3);
  TMultiGraph *m_resom_ampli=new TMultiGraph();
  for(int i=0;i<kCluster;i++){//over ampli custom
    for(int j=0;j<kTest;j++){//over z custom
      reso_buffer[j]=resom_ampli[i][j];
      e_reso_buffer[j]=e_resom_ampli[i][j];
    }
    p_resom_ampli[i]=new TGraphErrors(kTest,mult_custom,reso_buffer,NULL,e_reso_buffer);
    p_resom_ampli[i]->SetMarkerStyle(20+i);
    p_resom_ampli[i]->SetMarkerSize(0.4);
    p_resom_ampli[i]->SetMarkerColor(color[i]);
    p_resom_ampli[i]->SetLineColor(color[i]);
    m_resom_ampli->Add(p_resom_ampli[i]);
  }
  m_resom_ampli->Draw("AP");
  m_resom_ampli->SetTitle("TRACKS performances - resolution vs multiplicity;multiplicity;RMS [cm]");

  c_study->cd(4);
  TMultiGraph *m_resom_width=new TMultiGraph();
  for(int i=0;i<kCluster;i++){//over ampli custom
    for(int j=0;j<kTest;j++){//over z custom
      reso_buffer[j]=resom_width[i][j];
      e_reso_buffer[j]=e_resom_width[i][j];
    }
    p_resom_width[i]=new TGraphErrors(kTest,mult_custom,reso_buffer,NULL,e_reso_buffer);
    p_resom_width[i]->SetMarkerStyle(20+i);
    p_resom_width[i]->SetMarkerSize(0.4);
    p_resom_width[i]->SetMarkerColor(color[i]);
    p_resom_width[i]->SetLineColor(color[i]);
    m_resom_width->Add(p_resom_width[i]);
  }
  m_resom_width->Draw("AP");
  m_resom_width->SetTitle("TRACKS performances - resolution vs multiplicity;multiplicity;RMS [cm]");

  c_study->SaveAs(dirplot+"c_study.eps");

  //cpu info
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;
  printf("Performances info:\n\nCPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\nScroll up for info and verbosities. Thanks for using TRACKS!\n\n-> DONATE <-\n\n",cpu_time,real_time, cpu_efficiency);

}
