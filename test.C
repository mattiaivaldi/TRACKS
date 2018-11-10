//
//TRACKS - generation and recontruction of particle tracks in a detector
//developed by Luca Quaglia and Mattia Ivaldi, 2018
//
//START

#include "TCanvas.h"
#include "Riostream.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "vector"
#include "Event.h"
#include "Tools.h"
#include "Layer.h"
#include "Hit.h"
#include "Event.h"
#include "Particle.h"

using namespace TMath;

//PrintParticles is a variable used in order to decide wheter or not to print all the info about a particle (verbose)
//multiscatman is used in order to toogle on or off the multiscattering
void test(bool PrintParticles, bool multiscatman) {

  TStopwatch timer; //declared a timer to perform cpu efficiency measurements

  timer.Start(true); //start the timer

  //3 objects of the class layer which represent the beam pipe, first layer and second layer
  Layer *BP = new Layer(27.,3.,0.8,0.001); //length, radius, thickness and multiscattering RMS
  Layer *L1 = new Layer(27.,4.,0.2,0.001);
  Layer *L2 = new Layer(27.,7.,0.2,0.001);

  //object of class event
  Event *vgen = new Event(0, 0.001, 5.3, "kinem.root");
  //pointers to 3 arrays which keep track of the intersection point between the particle and the beam pipe/detectors
  double *hit_buffer_BP, *hit_buffer_L1, *hit_buffer_L2;
  //vectors of the class hit used in order to memorize all the coordinates of all the intersections
  vector <Hit*> cross_BP, cross_L1, cross_L2;
  //j, k and l are counters used later on and mult is the cast of the multiplicity (float) to an int value
  int j = 0, k = 0, l = 0, mult = (int)vgen->GetMult();

  //print out some info for the user
  printf("\n\n+++ TRACKS START - tracks generation +++\n\n");
  printf("Generated vertex with coordinates (%f, %f, %f)\nEvent multiplicity %d\n\n",vgen->GetX(),vgen->GetY(),vgen->GetZ(),mult);

  if (PrintParticles==true) {
    printf("Printing vertex and hit coordinates: ON\n\n");
  }else{printf("Printing vertex and hit coordinates: OFF\n\n");}

  //arrays which contain the values of the randomly generated theta and phi angles
  double theta[mult], phi[mult];

  //cycle over all the particles in current event
  for (int i = 0; i < mult; i++) {

    //create an object of the class particle
    Particle *part = new Particle("kinem.root");

    //get the values of theta and phi
    phi[i] = part->GetPhi();
    theta[i] = part->GetTheta();
    //print them out only if verbose is on
    if (PrintParticles==true) {
      printf("Particle %i: phi %f - theta %f\n",i,phi[i],theta[i]);
    }

    detect(vgen, &theta[i], &phi[i], BP, part, cross_BP, PrintParticles, multiscatman);

    detect(vgen, &theta[i], &phi[i], L1, part, cross_L1, PrintParticles, multiscatman);

    detect(vgen, &theta[i], &phi[i], L2, part, cross_L2, PrintParticles, multiscatman);

    delete part;

  }

  //pronted out how many particles have crossed which layer
  printf("\n\nOut of %d generated particles:\n\n%lu crossed BP\n%lu crossed L1\n%lu crossed L2\n\n+++ END generation +++",mult,cross_BP.size(),cross_L1.size(),cross_L1.size());

  //random info on the cpu usage
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;
  printf("\n\nCPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\n",cpu_time,real_time, cpu_efficiency);

} //END
