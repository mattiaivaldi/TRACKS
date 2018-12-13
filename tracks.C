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

  gSystem->CompileMacro("Layer.cxx",opt.Data());
  gSystem->CompileMacro("Hit.cxx",opt.Data());
  gSystem->CompileMacro("Particle.cxx",opt.Data());
  gSystem->CompileMacro("Tools.cxx",opt.Data());
  gSystem->CompileMacro("tracks_gen.C",opt.Data());
  gSystem->CompileMacro("tracks_reco.C",opt.Data());
  gROOT->ProcessLine("tracks_gen(0,1,1,1,0,10,2.65,30)");
  gROOT->ProcessLine("tracks_reco(0,1,0.0012,0.0003,1,3)");
  //gROOT->ProcessLine(".x spit_perform.C");
  //gROOT->ProcessLine(".x cluster_study.C");

  timer.Stop();
  double cpu_time = timer.CpuTime();
  double run_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/run_time)*100;

  printf("\n\n+++ TRACKS STOP +++\n\nCPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\nScroll up for info and verbosities. Thanks for using TRACKS!\n\n",cpu_time,run_time, cpu_efficiency);

}
