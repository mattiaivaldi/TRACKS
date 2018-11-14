# TRACKS

Monte Carlo particles' tracks generation and reconstruction with multiple scattering and noise. Developed by Luca Quaglia and Mattia Ivaldi - Torino, 2018.  

To use TRACKS you must have a working version of Root installed on your device. To generate the tracks just execute:

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
