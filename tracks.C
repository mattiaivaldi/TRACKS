//
//TRACKS - MC generation of particle tracks in a detector
//with multiple scattering and noise
//developed by Luca Quaglia and Mattia Ivaldi, 2018
//
//START

void tracks(TString myopt="fast"){

  TStopwatch timer;
  timer.Start(true);//start cpu monitor

  printf("\n\n+++ TRACKS START +++\n\n");

  TString opt;

  if(myopt.Contains("force")){
    opt = "kfg";
  }
  else {
    opt = "kg";
  }

  gSystem->CompileMacro("Layer.cxx",opt.Data());//compile classes and libraries
  gSystem->CompileMacro("Hit.cxx",opt.Data());
  gSystem->CompileMacro("Particle.cxx",opt.Data());
  gSystem->CompileMacro("Tools.cxx",opt.Data());
  gSystem->CompileMacro("tracks_gen.C",opt.Data());
  gSystem->CompileMacro("tracks_reco.C",opt.Data());
  gROOT->ProcessLine(".x spit_perform.C");//reconstruction performance study
  gROOT->ProcessLine(".x cluster_study.C");//peak finding performance study
  gROOT->ProcessLine("tracks_gen(0,1,1,-1,0,10,0,0)");//perform generation
  gROOT->ProcessLine("tracks_reco(0,1,0.0012,0.0003,1,5)");//perform reconstruction

  timer.Stop();//stop cpu monitoring
  double cpu_time = timer.CpuTime();
  double run_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/run_time)*100;

  printf("\n\n+++ TRACKS STOP +++\n\nGlobal CPU time = %f s\nGlobal Run time = %f s\nGlobal CPU efficiency = %f %% \n\nScroll up for info and verbosities. Thanks for using TRACKS!\n\n",cpu_time,run_time, cpu_efficiency);

}

//STOP
