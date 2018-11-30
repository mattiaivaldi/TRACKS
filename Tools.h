//This "class" doesn't inherits from any class since it's not a real class but it only contains the function used in order to find the intersection with the detector system

#ifndef TOOLS_H
#define TOOLS_H

#include "Tools.h"
#include "Layer.h"
#include "Hit.h"
#include "Particle.h"
#include "TClonesArray.h"
#include "TPaveText.h"
#include "vector"
#include "TGraph.h"

void verbosities(bool b_verbose, bool b_multiscatter, bool b_noise, int kExp);

void graphstyler(TGraph &graph, int divide);

void histostyler(TH1 &histo, int divide);

void pavestyler(TPaveText &pave, double textsize);

//x0, y0, z0 = generated vertex coordinates; theta, phi = generated angles; R = radius of the detector
double *hit_point(double x0, double y0, double z0, double theta, double phi, double R);

bool detect(Hit* vtx, Layer* L, Particle &part, TClonesArray &cross, bool b_verbose, bool b_multiscatter, int &counter, TH1D** histo);

void noise(bool b_verbose, int Noise, int Mult, TClonesArray &cross, Layer* L);

void smeagol(int index, double sigmaz, double sigmarf, double R, TClonesArray &cross);

double mode(vector<double> v);

#endif
