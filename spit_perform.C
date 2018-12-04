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

  const int kTest=11;

  double z_custom[kTest]={-13.,-10.,-7.,-4.,-2.,0.,2.,4.,7.,10.,13.};
  double mult_custom[kTest]={3,5,10,15,20,25,30,35,40,45,50};
  double resoz[kTest], e_resoz[kTest], resom[kTest], e_resom[kTest], effm[kTest], e_effm[kTest];

  double effz[kTest]={0.599157, 0.791143, 0.951297, 0.989414, 0.997564, 0.999472, 0.999870, 0.999968, 0.999991, 0.999999, 1.000000};
  double e_effz[kTest]={0.000490, 0.000406, 0.000215, 0.000102, 0.000049, 0.000023, 0.000011, 0.000006, 0.000003, 0.000001, 0.000000};

  TString z, m, exec;

  printf("+++ START resolution performances +++\n\n");

  for(int i=0; i<kTest;i++){
    z=Form("%f",z_custom[i]);
    exec="tracks_gen(0,0,1,1,15,10000,"+z+",20)";
    gROOT->ProcessLine(exec);
    reco_perform perform=tracks_reco(0,0,0.0012,0.0003);
    resoz[i]=perform.reso;
    e_resoz[i]=perform.e_reso;
    gROOT->Reset();
  }

  printf("+++ START efficiency performances +++\n\n");

  for(int i=0; i<kTest;i++){
    m=Form("%f",mult_custom[i]);
    exec="tracks_gen(0,0,1,1,15,10000,0,"+m+")";
    gROOT->ProcessLine(exec);
    reco_perform perform=tracks_reco(0,0,0.0012,0.0003);
    resom[i]=perform.reso;
    e_resom[i]=perform.e_reso;
    effm[i]=perform.eff;
    e_effm[i]=TMath::Sqrt(effm[i]*(1-effm[i])/perform.kExp);
    gROOT->Reset();
  }

  TCanvas *c_perform=new TCanvas("c_perform","c_perform",600,400);
  c_perform->Divide(2,2);
  c_perform->cd(1);
  TGraphErrors *p_resoz=new TGraphErrors(kTest,z_custom,resoz,NULL,e_resoz);
  p_resoz->SetTitle("TRACKS performances - resolution vs Z_{gen};z_{gen} [cm];RMS [cm]");
  graphstyler(*p_resoz,4);
  p_resoz->GetYaxis()->SetTitleOffset(1.1);
  p_resoz->SetMarkerStyle(20);
  p_resoz->SetMarkerSize(0.4);
  p_resoz->SetMarkerColor(1);
  p_resoz->Draw("AP");
  c_perform->cd(2);
  TGraphErrors *p_resom=new TGraphErrors(kTest,mult_custom,resom,NULL,e_resom);
  p_resom->SetTitle("TRACKS performances - resolution vs multiplicity;multiplicity;RMS [cm]");
  graphstyler(*p_resom,4);
  p_resom->GetYaxis()->SetTitleOffset(1.1);
  p_resom->SetMarkerStyle(20);
  p_resom->SetMarkerSize(0.4);
  p_resom->SetMarkerColor(1);
  p_resom->Draw("AP");
  c_perform->cd(3);
  TGraphErrors *p_effm=new TGraphErrors(kTest,mult_custom,effm,NULL,e_effm);
  p_effm->SetTitle("TRACKS performances - #varepsilon vs multiplicity;multiplicity;#varepsilon");
  graphstyler(*p_effm,4);
  p_effm->GetYaxis()->SetTitleOffset(0.5);
  p_effm->GetYaxis()->SetTitleSize(0.07);
  p_effm->SetMarkerStyle(20);
  p_effm->SetMarkerSize(0.4);
  p_effm->SetMarkerColor(1);
  p_effm->Draw("AP");
  c_perform->cd(4);
  TGraphErrors *p_effz=new TGraphErrors(kTest,mult_custom,effz,NULL,e_effz);
  p_effz->SetTitle("TRACKS performances - #varepsilon vs multiplicity for |z|<1#sigma;multiplicity;#varepsilon");
  graphstyler(*p_effz,4);
  p_effz->GetYaxis()->SetTitleOffset(0.5);
  p_effz->GetYaxis()->SetTitleSize(0.07);
  p_effz->SetMarkerStyle(20);
  p_effz->SetMarkerSize(0.4);
  p_effz->SetMarkerColor(1);
  p_effz->Draw("AP");
  c_perform->SaveAs(dirplot+"c_perform.eps");

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
