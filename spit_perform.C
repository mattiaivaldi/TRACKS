void spit_perform(){

  gROOT->Reset();

  TStopwatch timer;
  timer.Start(true);//start cpu monitor

  double z_custom[9]={-13.,-9.,-4.,-2.,0.,2.,4.,9.,13.};
  double mult_custom[11]={2,5,10,15,20,25,30,35,40,45,50};
  double resoz[9], resom[11], effm[11];

  TString z, m, exec;

  printf("+++ START resolution performances +++\n\n");

  for(int i=0; i<9;i++){
    z=Form("%f",z_custom[i]);
    exec="tracks_gen(0,0,1,1,15,100000,"+z+",20)";
    gROOT->ProcessLine(exec);
    reco_perform perform=tracks_reco(0,0,0.0012,0.0003);
    resoz[i]=perform.reso;
  }

  printf("+++ START efficiency performances +++\n\n");

  for(int i=0; i<11;i++){
    m=Form("%f",mult_custom[i]);
    exec="tracks_gen(0,0,1,1,15,100000,0,"+m+")";
    gROOT->ProcessLine(exec);
    reco_perform perform=tracks_reco(0,0,0.0012,0.0003);
    resom[i]=perform.reso;
    effm[i]=perform.eff;
  }

  TCanvas *c_reso=new TCanvas("c_reso","c_reso",600,400);
  c_reso->Divide(2,1);
  c_reso->cd(1);
  TGraph *p_resoz=new TGraph(9,z_custom,resoz);
  p_resoz->SetTitle("TRACKS performances - resolution vs Z_{gen};z_{gen} [cm];RMS [cm]");
  graphstyler(*p_resoz,2);
  p_resoz->SetMarkerStyle(20);
  p_resoz->SetMarkerSize(0.6);
  p_resoz->SetMarkerColor(1);
  p_resoz->Draw("AP");
  c_reso->cd(2);
  TGraph *p_resom=new TGraph(11,mult_custom,resom);
  p_resom->SetTitle("TRACKS performances - resolution vs multiplicity;multiplicity;RMS [cm]");
  graphstyler(*p_resom,2);
  p_resom->SetMarkerStyle(20);
  p_resom->SetMarkerSize(0.6);
  p_resom->SetMarkerColor(1);
  p_resom->Draw("AP");

  TCanvas *c_eff=new TCanvas("c_eff","c_eff",600,400);
  c_eff->Divide(2,1);
  c_eff->cd(1);
  TGraph *p_effm=new TGraph(11,mult_custom,effm);
  p_effm->SetTitle("TRACKS performances - #varepsilon vs multiplicity;multiplicity;#varepsilon");
  graphstyler(*p_effm,2);
  p_effm->SetMarkerStyle(20);
  p_effm->SetMarkerSize(0.6);
  p_effm->SetMarkerColor(1);
  p_effm->Draw("AP");

  //cpu info
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;
  printf("Performances info:\n\nCPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\nScroll up for info and verbosities. Thanks for using TRACKS!\n\n-> DONATE <-\n\n",cpu_time,real_time, cpu_efficiency);

}
