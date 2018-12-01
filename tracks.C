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
    gROOT->ProcessLine("tracks_gen(0,1,1,0,100000,1.04,10)");
    //gROOT->ProcessLine(".x gen_performance.C");
    gROOT->ProcessLine("tracks_reco(0,0.0012,0.0003)");
}