# TRACKS

Monte Carlo particles' tracks generation and reconstruction with multiple scattering and noise. Developed by Luca Quaglia and Mattia Ivaldi - Torino, 2018.  

TRACKS requires a working version of Root on your device. It has an easy 3-steps usage:
1. open and edit the file _detector_info.txt_ with the format  
`layer_name width radius thickness theta_rms`

2. open and edit the macro tracks.C with your preferences  
`gROOT->ProcessLine(“tracks_gen(a,b,c,d,e,f,g,h)”)` where  
bool a = verbose mode ON/OFF  
bool b = print and save plots ON/OFF  
bool c = multiple scattering ON/OFF    
bool d = noise ON/OFF  
int e = 5 custom z, 10 custom multiplicity, 15 both custom  
int f = # of collisions performed  
double g = custom vertex z  
double h = custom event multiplicity
`gROOT->ProcessLine(“tracks_reco(a,b,c)”)` where  
bool a = verbose mode ON/OFF  
bool b = print and save plots ON/OFF  
double c = smearing parameter on z  
double d = smearing parameter on phi  
double e = ambiguity check amplitude  
int e = ambiguity check width

3. open a ROOT session and interpret root [0] .x tracks.C

The generated data will be stored in the _gen.root_ file, the reconstruction plots will be stored in the 

Caveat! When you use TRACKS do pay attention to time and CPU efficiency, the software operates with the following performances:

<img src="https://github.com/mattiaivaldi/TRACKS/blob/TRACKSinprogress/c_perform.jpg" alt="alt text" width="400" height="250">

-_reconstruction and wiki under development_-
