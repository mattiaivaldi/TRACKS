# TRACKS

Monte Carlo particles' tracks generation and reconstruction with multiple scattering and noise. Developed by Luca Quaglia and Mattia Ivaldi - Torino, 2018.  

To use TRACKS you must have a working version of Root installed on your device. To generate the tracks just execute:

root[0] .x classcompiler.C _to load and compile all the libraries_  
root[1] tracks_gen(a,b,c,d)

where  
_a_ = verbose mode ON/OFF  
_b_ = multiple scattering ON/OFF  
_c_ = noise ON/OFF  
must be 1 or 0 and _d_ is the (integer) number of collisions you want to generate. The output is a file _

-_reconstruction under development_-
