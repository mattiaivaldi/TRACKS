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
  vector <Hit*> ciccioBP, ciccioL1, ciccioL2;
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

  //variables used in order to see if there is an effective intersection between a particle and a piece of the detector
  bool bBP = false, bL1 = false, bL2 = false; //they are all set to false and changed to true if we have an intersection

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

    //---intersection with the beam pipe---//

    //hit buffer is used in order to save the values of x, y and z of the current particle intersection in the event
    hit_buffer_BP=hit_point(vgen->GetX(),vgen->GetY(),vgen->GetZ(),theta[i],phi[i],BP->GetRadius());

    //check to see if the z of the intersection is in the length of the beam pipe (27 cm)
    if(*(hit_buffer_BP+2) >= -(BP->GetWidth()/2.) && *(hit_buffer_BP+2) <= (BP->GetWidth()/2.)) {
      //there has been an interesction with the beam pipe
      bBP = true;
      //object of the class hit is created here, the values calculated for the current particle are used
      Hit *hit_BP = new Hit(*(hit_buffer_BP+0),*(hit_buffer_BP+1),*(hit_buffer_BP+2));
      //the object is pushed into a vector
      ciccioBP.push_back(hit_BP);
      //print the coordinates (x,y,z) of intersection with beam pipe usng the getter function of the hit class
      if (PrintParticles==true) {
	       printf("Hit with BP at (%f, %f, %f)\n",ciccioBP[j]->GetX(),ciccioBP[j]->GetY(),ciccioBP[j]->GetZ());
      }
      //direction before multiple scattering in the beam pipe
      cout << "Before beam pipe, phi = " << phi[i] << " theta = " << theta[i] << endl;
      //if there is an intersection and verbose is set to true
      if (bBP == true && multiscatman == true) {
        //using the member function GetRMS of the class layer.cxx the rotation is calculated
        part->Rotate(BP->GetRMS());
        //the angles are updated using the getter function of the class
        phi[i] = part->GetPhi();
        theta[i] = part->GetTheta();
      }
      //the counter is incremented in order to be able to read and print out the next particle intersection coordinates
      j++;
      //Print to video
      cout << "After beam pipe, phi = " << phi[i] << " theta = " << theta[i] << endl;
    }

    //---intersection with L1 is the same as beam pipe (only some name changes here and there)---//
    hit_buffer_L1=hit_point(vgen->GetX(),vgen->GetY(),vgen->GetZ(),theta[i],phi[i],L1->GetRadius());

    if(*(hit_buffer_L1+2) >= -(L1->GetWidth()/2.) && *(hit_buffer_L1+2) <= (L1->GetWidth()/2.)) {
      bL1 = true;
      Hit *hit_L1 = new Hit(*(hit_buffer_L1+0),*(hit_buffer_L1+1),*(hit_buffer_L1+2));
      ciccioL1.push_back(hit_L1);
      if (PrintParticles==true) {
	       printf("Hit with L1 at (%f, %f, %f)\n",ciccioL1[k]->GetX(),ciccioL1[k]->GetY(),ciccioL1[k]->GetZ());
      }
      cout << "Angoli prima " << phi[i] << " " << theta[i] << endl;
      if (bL1 == true && multiscatman == true) {
        part->Rotate(L1->GetRMS());
        phi[i] = part->GetPhi();
        theta[i] = part->GetTheta();
      }
      k++;
      cout << "Angoli dopo " << phi[i] << " " << theta[i] << endl;
    }

    //---intersection with L2, here there is no multiple scattering for now (we don't care about what happens after layer 2)---//
    hit_buffer_L2=hit_point(vgen->GetX(),vgen->GetY(),vgen->GetZ(),theta[i],phi[i],L2->GetRadius());

    if(*(hit_buffer_L2+2) >= -(L2->GetWidth()/2.) && *(hit_buffer_L2+2) <= (L2->GetWidth()/2.)) {
      Hit *hit_L2 = new Hit(*(hit_buffer_L2+0),*(hit_buffer_L2+1),*(hit_buffer_L2+2));
      ciccioL2.push_back(hit_L2);
      if (PrintParticles==true) {
	printf("Hit with L2 at (%f, %f, %f)\n",ciccioL2[l]->GetX(),ciccioL2[l]->GetY(),ciccioL2[l]->GetZ());
      }
      l++;
    }

    if (PrintParticles==true){printf("\n");}

    //deleting the object of the class particle at the end of the for cycle
    delete part;

  }//end of for cycle over all the particles in the current event

  //pronted out how many particles have crossed which layer
  printf("\n\nOut of %d generated particles:\n\n%lu crossed BP\n%lu crossed L1\n%lu crossed L2\n\n+++ END generation +++",mult,ciccioBP.size(),ciccioL1.size(),ciccioL1.size());

  //random info on the cpu usage
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;
  printf("\n\nCPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\n",cpu_time,real_time, cpu_efficiency);

} //END
