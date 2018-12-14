#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif //needed to get current working directory

void spit_perform(){

  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetFillStyle(0);

  TStopwatch timer;
  timer.Start(true);//start cpu monitor

  char cwd[FILENAME_MAX];
  GetCurrentDir(cwd,FILENAME_MAX);
  TString dirplot=TString(cwd)+"/tracksplot/";//define the path where the plots will be saved

  const int kTest=11;

  double z_custom[17]={-40,-30,-20,-13,-10,-7,-4,-2,0,2,4,7,10,13,20,30,40};
  double mult_custom[kTest]={3,5,10,15,20,25,30,35,40,45,50};
  double resoz[17], e_resoz[17], resom[kTest], e_resom[kTest], effm[kTest], e_effm[kTest], effz[17], e_effz[17];

  double effm1s[kTest]={0.870975,0.878662,0.875956,0.875345,0.877736,0.875762,0.879214,0.870685,0.878893,0.868891,0.875211};
  double e_effm1s[kTest]={0.000335,0.000327,0.000330,0.000330,0.000328,0.000330,0.000326,0.000336,0.000326,0.000338,0.000330};

  TString z, m, exec;//to vary vertex z and multiplicity during the study

  printf("\n\nxxx START performances: vertex z xxx\n\n");

  for(int i=0; i<17;i++){//performances varying vertex z
    z=Form("%f",z_custom[i]);
    exec="tracks_gen(0,0,1,-1,15,100000,"+z+",20)";
    gROOT->ProcessLine(exec);
    reco_perform perform=tracks_reco(0,0,0.0012,0.0003,1,5);
    resoz[i]=perform.reso;
    e_resoz[i]=perform.e_reso;
    effz[i]=perform.eff;
    e_effz[i]=perform.e_eff;
  }

  printf("\n\nxxx START performances: multiplicity xxx\n\n");

  for(int i=0; i<kTest;i++){//performances varying multiplicity
    m=Form("%f",mult_custom[i]);
    exec="tracks_gen(0,0,1,-1,15,100000,0,"+m+")";
    gROOT->ProcessLine(exec);
    reco_perform perform=tracks_reco(0,0,0.0012,0.0003,1,5);
    resom[i]=perform.reso;
    e_resom[i]=perform.e_reso;
    effm[i]=perform.eff;
    e_effm[i]=perform.e_eff;
  }

  TFile f_perform("perform.root","RECREATE");//to store performance plots

  TPaveText *pt_perform = new TPaveText(0.22,0.11,0.41,0.26,"NDC");
  pavestyler(*pt_perform,0.03);
  pt_perform->AddText(Form("%d events",1000000));

  TCanvas *c_perform=new TCanvas("c_perform","c_perform",600,400);
  c_perform->Divide(2,2);

  c_perform->cd(1);
  TGraphErrors *p_resoz=new TGraphErrors(17,z_custom,resoz,NULL,e_resoz);
  p_resoz->SetTitle("TRACKS performances - resolution vs Z_{gen};z_{gen} [cm];RMS [cm]");
  graphstyler(*p_resoz,4);
  p_resoz->GetYaxis()->SetTitleOffset(1.1);
  p_resoz->SetMarkerStyle(20);
  p_resoz->SetMarkerSize(0.4);
  p_resoz->SetMarkerColor(1);
  p_resoz->Draw("AP");
  pt_perform->Draw();

  c_perform->cd(2);
  TGraphErrors *p_resom=new TGraphErrors(kTest,mult_custom,resom,NULL,e_resom);
  p_resom->SetTitle("TRACKS performances - resolution vs multiplicity;multiplicity;RMS [cm]");
  graphstyler(*p_resom,4);
  p_resom->GetYaxis()->SetTitleOffset(1.1);
  p_resom->SetMarkerStyle(20);
  p_resom->SetMarkerSize(0.4);
  p_resom->SetMarkerColor(1);
  p_resom->Draw("AP");
  pt_perform->Draw();

  c_perform->cd(3);
  TGraphErrors *p_effz=new TGraphErrors(17,z_custom,effz,NULL,e_effz);
  p_effz->SetTitle("TRACKS performances - efficiency vs Z_{gen};z_{gen} [cm];#tilde{#varepsilon}");
  graphstyler(*p_effz,4);
  p_effz->GetYaxis()->SetTitleOffset(0.5);
  p_effz->GetYaxis()->SetTitleSize(0.07);
  p_effz->SetMarkerStyle(20);
  p_effz->SetMarkerSize(0.4);
  p_effz->SetMarkerColor(1);
  p_effz->Draw("AP");
  pt_perform->Draw();

  c_perform->cd(4);
  TGraphErrors *p_effm=new TGraphErrors(kTest,mult_custom,effm,NULL,e_effm);
  p_effm->SetMarkerStyle(20);
  p_effm->SetMarkerSize(0.4);
  p_effm->SetMarkerColor(1);
  TGraphErrors *p_effm1s=new TGraphErrors(kTest,mult_custom,effm1s,NULL,e_effm1s);
  p_effm1s->SetMarkerStyle(22);
  p_effm1s->SetMarkerSize(0.4);
  p_effm1s->SetMarkerColor(2);
  TMultiGraph *m_effm=new TMultiGraph();
  m_effm->Add(p_effm);
  m_effm->Add(p_effm1s);
  m_effm->Draw("AP");
  m_effm->SetTitle("TRACKS performances - efficiency vs multiplicity;multiplicity;#tilde{#varepsilon}");
  m_effm->GetXaxis()->SetLabelSize(0.045);
  m_effm->GetYaxis()->SetLabelSize(0.045);
  m_effm->GetXaxis()->SetTitleSize(0.05);
  m_effm->GetYaxis()->SetTitleSize(0.07);
  m_effm->GetXaxis()->SetTitleOffset(0.9);
  m_effm->GetYaxis()->SetTitleOffset(0.5);

  auto legt = new TLegend(0.63,0.17,0.82,0.43);
  legt->SetHeader("1000000 events","");
  legt->AddEntry(p_effm,"z_{gen} = 0 cm","p");
  legt->AddEntry(p_effm1s,"|Z_{gen}| < 1#sigma","p");
  legt->Draw();

  c_perform->SaveAs(dirplot+"c_perform.eps");

  c_perform->Write();

  //uncomment to study the efficienciency for generated z whithin 1 sigma
  //a change in the Hit event constructor is required
  /*for(int i=0; i<kTest;i++){
    m=Form("%f",mult_custom[i]);
    exec="tracks_gen(0,0,1,-1,10,1000000,0,"+m+")";
    gROOT->ProcessLine(exec);
    reco_perform perform=tracks_reco(0,0,0.0012,0.0003,1,5);
    effz[i]=perform.eff;
    e_effz[i]=perform.e_eff;
  }

  for(int i=0; i<kTest;i++){
    printf("%f,",effz[i]);
  }
  cout<<endl;
  for(int i=0; i<kTest;i++){
    printf("%f,",e_effz[i]);
  }*/

  f_perform.Close();

  timer.Stop();//stop cpu monitoring
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;

  printf("Performances info:\n\nCPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\n",cpu_time,real_time, cpu_efficiency);

}
