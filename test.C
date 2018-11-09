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

void test(bool PrintParticles, bool multiscatman) {
//gRandom->SetSeed(0);
  //
  //TRACK - generation and recontruction of particle tracks in a detector
  //developed by Luca Quaglia and Mattia Ivaldi, 2018
  //

  TStopwatch timer;

  timer.Start(true);

  Layer *BP = new Layer(27.,3.,0.8,0.001);
  Layer *L1 = new Layer(27.,4.,0.2,0.001);
  Layer *L2 = new Layer(27.,7.,0.2,0.001);

  Event *vgen = new Event(0, 0.001, 5.3, "kinem.root"); //dichiarato oggetto della struct event
  double *hit_buffer_BP, *hit_buffer_L1, *hit_buffer_L2;
  vector <Hit*> ciccioBP, ciccioL1, ciccioL2;
  int j = 0, k = 0, l = 0, mult = (int)vgen->GetMult();

  //stampa a video
  printf("\n\n+++ TRACKS START - tracks generation +++\n\n");
  printf("Generated vertex with coordinates (%f, %f, %f)\nEvent multiplicity %d\n\n",vgen->GetX(),vgen->GetY(),vgen->GetZ(),mult);

  if (PrintParticles==true) {
    printf("Printing vertex and hit coordinates: ON\n\n");
  }else{printf("Printing vertex and hit coordinates: OFF\n\n");}

  //array che contengono i valori di angolo theta e phi per le particelle dell'evento
  double theta[mult], phi[mult];

  //calcolo i valori casauali per theta e phi, phi è uniforme tra 0 e 2 pi mentre theta segue la distribuzione ricavata nel caso del flusso
  //isootropo
  for (int i = 0; i < mult; i++) {

    Particle *part = new Particle("kinem.root");

    phi[i] = part->GetPhi();
    theta[i] = part->GetTheta(); //moltiplico per 2 per avere angoli tra 0 e 2 pi anche per theta
    if (PrintParticles==true) {
      printf("Particle %i: phi %f - theta %f\n",i,phi[i],theta[i]);
    }

    //calcolo intersezioni tracce con la beam pipe
    hit_buffer_BP=hit_point(vgen->GetX(),vgen->GetY(),vgen->GetZ(),theta[i],phi[i],BP->GetRadius());

    if(*(hit_buffer_BP+2) >= -(BP->GetWidth()/2.) && *(hit_buffer_BP+2) <= (BP->GetWidth()/2.)) {
      Hit *hit_BP = new Hit(*(hit_buffer_BP+0),*(hit_buffer_BP+1),*(hit_buffer_BP+2));
      ciccioBP.push_back(hit_BP);
      if (PrintParticles==true) {
	printf("Hit with BP at (%f, %f, %f)\n",ciccioBP[j]->GetX(),ciccioBP[j]->GetY(),ciccioBP[j]->GetZ());
      }
      j++;
    }

    cout << "Angoli prima " << theta[i] << " " << phi[i] << endl;

    if (multiscatman == false) {;} //if multiscattering is off -> do nothing

    else {
      part->Rotate(BP->GetRMS()); //using the member function GetRMS of the class layer.cxx
      phi[i] = part->GetPhi();
      theta[i] = part->GetTheta();
    }

    cout << "Angoli dopo " << theta[i] << " " << phi[i] << endl;

    //intersection with L1
    hit_buffer_L1=hit_point(vgen->GetX(),vgen->GetY(),vgen->GetZ(),theta[i],phi[i],L1->GetRadius());

    if(*(hit_buffer_L1+2) >= -(L1->GetWidth()/2.) && *(hit_buffer_L1+2) <= (L1->GetWidth()/2.)) {
      Hit *hit_L1 = new Hit(*(hit_buffer_L1+0),*(hit_buffer_L1+1),*(hit_buffer_L1+2));
      ciccioL1.push_back(hit_L1);
      if (PrintParticles==true) {
	printf("Hit with L1 at (%f, %f, %f)\n",ciccioL1[k]->GetX(),ciccioL1[k]->GetY(),ciccioL1[k]->GetZ());
      }
      k++;
    }

    if (multiscatman == false) {;} //if multiscattering is off -> do nothing

    else {
      part->Rotate(L1->GetRMS()); //using the member function GetRMS of the class layer.cxx
      phi[i] = part->GetPhi();
      theta[i] = part->GetTheta();
    }

    //intersection with L2
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

    delete part; //deleting the object at the end of the for cycle

  }//fine for

  printf("\n\nOut of %d generated particles:\n\n%lu crossed BP\n%lu crossed L1\n%lu crossed L2\n\n+++ END generation +++",mult,ciccioBP.size(),ciccioL1.size(),ciccioL1.size());

  //cose a caso sulla cpu
  timer.Stop();
  double cpu_time = timer.CpuTime();
  double real_time = timer.RealTime();
  double cpu_efficiency = (cpu_time/real_time)*100;
  printf("\n\nCPU time = %f s\nRun time = %f s\nCPU efficiency = %f %% \n\n",cpu_time,real_time, cpu_efficiency);

}


//END
