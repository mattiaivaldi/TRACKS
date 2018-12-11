void tracks(TString myopt="fast"){
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
    gROOT->ProcessLine("tracks_gen(0,1,1,1,0,1000000,0,30)");
    gROOT->ProcessLine("tracks_reco(0,1,0.0012,0.0003,1,3)");
    gROOT->ProcessLine(".x spit_perform.C");
    //gROOT->ProcessLine(".x cluster_study.C");
    gROOT->Reset();
}
