//Macro used in order to compile all classes. You have to run this macro with .x before executing the main code

void classcompiler(TString myopt="fast"){
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
}
