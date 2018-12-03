#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

void spit_perform(){

  gROOT->Reset();

  TStopwatch timer;
  timer.Start(true);//start cpu monitor

  char cwd[FILENAME_MAX];
  GetCurrentDir(cwd,FILENAME_MAX);
  TString dirplot=TString(cwd)+"/tracksplot/";

  //TString dirplot=TString("/Users/mattiaivaldi/GitHub/TRACKS/")+"tracksplot/";

  const int kTest=11;

  double z_custom[kTest]={-13.,-10.,-7.,-4.,-2.,0.,2.,4.,7.,10.,13.};
  double mult_custom[kTest]={2,5,10,15,20,25,30,35,40,45,50};
  double resoz[kTest], resom[kTest], effm[kTest];

  TString z, m, exec;

  printf("+++ START resolution performances +++\n\n");

  for(int i=0; i<kTest;i++){
    z=Form("%f",z_custom[i]);
    exec="tracks_gen(0,0,1,1,15,100000,"+z+",20)";
    gROOT->ProcessLine(exec);
    reco_perform perform=tracks_reco(0,0,0.0012,0.0003);
    resoz[i]=perform.reso;
  }

  printf("+++ START efficiency performances +++\n\n");

  for(int i=0; i<kTest;i++){
    m=Form("%f",mult_custom[i]);
    exec="tracks_gen(0,0,1,1,15,100000,0,"+m+")";
    gROOT->ProcessLine(exec);
    reco_perform perform=tracks_reco(0,0,0.0012,0.0003);
    resom[i]=perform.reso;
    effm[i]=perform.eff;
  }

  TCanvas *c_perform=new TCanvas("c_perform","c_perform",600,400);
  c_perform->Divide(2,2);
  c_perform->cd(1);
  TGraph *p_resoz=new TGraph(kTest,z_custom,resoz);
  p_resoz->SetTitle("TRACKS performances - resolution vs Z_{gen};z_{gen} [cm];RMS [cm]");
  graphstyler(*p_resoz,4);
  p_resoz->GetYaxis()->SetTitleOffset(1.1);
  p_resoz->SetMarkerStyle(20);
  p_resoz->SetMarkerSize(0.4);
  p_resoz->SetMarkerColor(1);
  p_resoz->Draw("AP");
  c_perform->cd(2);
  TGraph *p_resom=new TGraph(kTest,mult_custom,resom);
  p_resom->SetTitle("TRACKS performances - resolution vs multiplicity;multiplicity;RMS [cm]");
  graphstyler(*p_resom,4);
  p_resom->GetYaxis()->SetTitleOffset(1.1);
  p_resom->SetMarkerStyle(20);
  p_resom->SetMarkerSize(0.4);
  p_resom->SetMarkerColor(1);
  p_resom->Draw("AP");
  c_perform->cd(3);
  TGraph *p_effm=new TGraph(kTest,mult_custom,effm);
  p_effm->SetTitle("TRACKS performances - #varepsilon vs multiplicity;multiplicity;#varepsilon");
  graphstyler(*p_effm,4);
  p_effm->GetYaxis()->SetTitleOffset(0.5);
  p_effm->GetYaxis()->SetTitleSize(0.07);
  p_effm->SetMarkerStyle(20);
  p_effm->SetMarkerSize(0.4);
  p_effm->SetMarkerColor(1);
  p_effm->Draw("AP");
  c_perform->SaveAs(dirplot+"c_perform.eps");

  //cpu info
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;
  printf("Performances info:\n\nCPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\nScroll up for info and verbosities. Thanks for using TRACKS!\n\n-> DONATE <-\n\n",cpu_time,real_time, cpu_efficiency);

}
