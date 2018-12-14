# TRACKS

Monte Carlo particles' tracks generation and reconstruction with multiple scattering and noise. Developed by Luca Quaglia and Mattia Ivaldi - Torino, 2018.  

TRACKS requires a working version of Root on your device. It has an easy 3-steps usage:
1. open and edit the file detector_info.txt with the format  
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

3.

root[0] .x classcompiler.C _to load and compile all the libraries_  
root[1] tracks_gen(a,b,c,d)

where  
_a_ = verbose mode ON/OFF  
_b_ = multiple scattering ON/OFF  
_c_ = noise ON/OFF  
must be 1 or 0 and _d_ is the (integer) number of collisions you want to generate. The output is a file _gen.root_ containing the hit data from the detector layers.

Caveat! When you use TRACKS do pay attention to time and CPU efficiency, the software operates with the following performances:

<img src="https://github.com/mattiaivaldi/TRACKS/blob/TRACKSinprogress/c_perform.jpg" alt="alt text" width="400" height="250">

-_reconstruction and wiki under development_-
