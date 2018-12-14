#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif //needed to get current working directory

void cluster_study(){

  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetFillStyle(0);

  TStopwatch timer;
  timer.Start(true);//start cpu monitor

  char cwd[FILENAME_MAX];
  GetCurrentDir(cwd,FILENAME_MAX);
  TString dirplot=TString(cwd)+"/tracksplot/";//define the path where the plots will be saved

  const int kTest=6, kCluster=4;

  int color[kCluster]={2,4,6,9};

  double ampli_custom[kCluster]={0.7,1,1.3,1.5};
  double width_custom[kCluster]={1,3,5,7};

  double z_custom[kTest]={0,2,4,7,10,13};
  double mult_custom[kTest]={3,5,15,30,40,50};

  double resoz_ampli[kCluster][kTest], e_resoz_ampli[kCluster][kTest], resoz_width[kCluster][kTest], e_resoz_width[kCluster][kTest], resom_ampli[kCluster][kTest], e_resom_ampli[kCluster][kTest], resom_width[kCluster][kTest], e_resom_width[kCluster][kTest];

  TGraphErrors *p_resoz_ampli[kCluster],*p_resoz_width[kCluster], *p_resom_ampli[kCluster],*p_resom_width[kCluster];

  double reso_buffer[kTest], e_reso_buffer[kTest];

  auto leg_z_ampli = new TLegend(0.15,0.65,0.41,0.85);
  leg_z_ampli->SetHeader("500000 events - amplitude","");

  auto leg_z_width = new TLegend(0.15,0.65,0.41,0.85);
  leg_z_width->SetHeader("500000 events - width","");

  auto leg_m_ampli = new TLegend(0.68,0.65,0.93,0.86);
  leg_m_ampli->SetHeader("500000 events - amplitude","");

  auto leg_m_width = new TLegend(0.68,0.65,0.93,0.86);
  leg_m_width->SetHeader("500000 events - width","");

  TString z, m, ampli, width, exec;//to vary vertex z, multiplicity, amplitude and width during the study

  printf("\n\nxxx START performances: amplitude xxx\n\n");

  for(int i=0; i<kCluster;i++){//resolution vs vertex z for different amplitude
    for(int j=0; j<kTest;j++){
      z=Form("%f",z_custom[j]);
      exec="tracks_gen(0,0,1,1,15,1,"+z+",20)";
      gROOT->ProcessLine(exec);
      reco_perform perform=tracks_reco(0,0,0.0012,0.0003,ampli_custom[i],3);
      resoz_ampli[i][j]=perform.reso;
      e_resoz_ampli[i][j]=perform.e_reso;
    }
  }

  for(int i=0; i<kCluster;i++){//resolution vs vertex z for different width
    for(int j=0; j<kTest;j++){
      z=Form("%f",z_custom[j]);
      exec="tracks_gen(0,0,1,1,15,1,"+z+",20)";
      gROOT->ProcessLine(exec);
      reco_perform perform=tracks_reco(0,0,0.0012,0.0003,1,width_custom[i]);
      resoz_width[i][j]=perform.reso;
      e_resoz_width[i][j]=perform.e_reso;
    }
  }

  printf("\n\nxxx START performances: width xxx\n\n");

  for(int i=0; i<kCluster;i++){//resolution vs multiplicity for different amplitude
    for(int j=0; j<kTest;j++){
      m=Form("%f",mult_custom[j]);
      exec="tracks_gen(0,0,1,1,15,1,0,"+m+")";
      gROOT->ProcessLine(exec);
      reco_perform perform=tracks_reco(0,0,0.0012,0.0003,ampli_custom[i],3);
      resom_ampli[i][j]=perform.reso;
      e_resom_ampli[i][j]=perform.e_reso;
    }
  }

  for(int i=0; i<kCluster;i++){//resolution vs multiplicity for different width
    for(int j=0; j<kTest;j++){
      m=Form("%f",mult_custom[j]);
      exec="tracks_gen(0,0,1,1,15,1,0,"+m+")";
      gROOT->ProcessLine(exec);
      reco_perform perform=tracks_reco(0,0,0.0012,0.0003,1,width_custom[i]);
      resom_width[i][j]=perform.reso;
      e_resom_width[i][j]=perform.e_reso;
    }
  }

  TFile f_perform("perform.root","UPDATE");//to store performance plots

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
    leg_z_ampli->AddEntry(p_resoz_ampli[i],Form("%f",ampli_custom[i]),"p");
  }
  m_resoz_ampli->Draw("AP");
  m_resoz_ampli->SetTitle("TRACKS performances - resolution vs Z_{gen};z_{gen} [cm];RMS [cm]");
  m_resoz_ampli->GetXaxis()->SetLabelSize(0.045);
  m_resoz_ampli->GetYaxis()->SetLabelSize(0.045);
  m_resoz_ampli->GetXaxis()->SetTitleSize(0.05);
  m_resoz_ampli->GetYaxis()->SetTitleSize(0.05);
  m_resoz_ampli->GetXaxis()->SetTitleOffset(0.9);
  m_resoz_ampli->GetYaxis()->SetTitleOffset(1.1);
  leg_z_ampli->Draw();

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
    leg_z_width->AddEntry(p_resoz_width[i],Form("%f",width_custom[i]),"p");
  }
  m_resoz_width->Draw("AP");
  m_resoz_width->SetTitle("TRACKS performances - resolution vs Z_{gen};z_{gen} [cm];RMS [cm]");
  m_resoz_width->GetXaxis()->SetLabelSize(0.045);
  m_resoz_width->GetYaxis()->SetLabelSize(0.045);
  m_resoz_width->GetXaxis()->SetTitleSize(0.05);
  m_resoz_width->GetYaxis()->SetTitleSize(0.05);
  m_resoz_width->GetXaxis()->SetTitleOffset(0.9);
  m_resoz_width->GetYaxis()->SetTitleOffset(1.1);
  leg_z_width->Draw();

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
    leg_m_ampli->AddEntry(p_resom_ampli[i],Form("%f",ampli_custom[i]),"p");
  }
  m_resom_ampli->Draw("AP");
  m_resom_ampli->SetTitle("TRACKS performances - resolution vs multiplicity;multiplicity;RMS [cm]");
  m_resom_ampli->GetXaxis()->SetLabelSize(0.045);
  m_resom_ampli->GetYaxis()->SetLabelSize(0.045);
  m_resom_ampli->GetXaxis()->SetTitleSize(0.05);
  m_resom_ampli->GetYaxis()->SetTitleSize(0.05);
  m_resom_ampli->GetXaxis()->SetTitleOffset(0.9);
  m_resom_ampli->GetYaxis()->SetTitleOffset(1.1);
  leg_m_ampli->Draw();

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
    leg_m_width->AddEntry(p_resom_width[i],Form("%f",width_custom[i]),"p");
  }
  m_resom_width->Draw("AP");
  m_resom_width->SetTitle("TRACKS performances - resolution vs multiplicity;multiplicity;RMS [cm]");
  m_resom_width->GetXaxis()->SetLabelSize(0.045);
  m_resom_width->GetYaxis()->SetLabelSize(0.045);
  m_resom_width->GetXaxis()->SetTitleSize(0.05);
  m_resom_width->GetYaxis()->SetTitleSize(0.05);
  m_resom_width->GetXaxis()->SetTitleOffset(0.9);
  m_resom_width->GetYaxis()->SetTitleOffset(1.1);
  leg_m_width->Draw();

  c_study->SaveAs(dirplot+"c_study.eps");

  c_study->Write();

  f_perform.Close();

  timer.Stop();//stop cpu monitoring
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;
  printf("Performances info:\n\nCPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\n",cpu_time,real_time, cpu_efficiency);

}
