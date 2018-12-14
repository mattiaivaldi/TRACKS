#ifndef TOOLS_H
#define TOOLS_H

#include "Tools.h"
#include "Layer.h"
#include "Hit.h"
#include "Particle.h"
#include "TClonesArray.h"
#include "TPaveText.h"
#include "THStack.h"
#include "vector"
#include "TGraph.h"
#include "TString.h"

void verbosities(bool b_verbose, bool b_multiscatter, bool b_noise, int kExp);//general information about the simulation

void graphstyler(TGraph &graph, int divide);//TGraph makeup

void histostyler(TH1 &histo, int divide);//TH1 makeup

void stackstyler(THStack &stack);//THStack makeup

void pavestyler(TPaveText &pave, double textsize);//TPave makeup

double *hit_point(double x0, double y0, double z0, double theta, double phi, double R);//evaluate hit coordinates

bool detect(Hit* vtx, Layer* L, Particle &part, TClonesArray &cross, bool b_verbose, bool b_multiscatter, int &counter, TH1D** histo);//hit check

void noise(bool b_verbose, int Noise, int Mult, TClonesArray &cross, Layer* L);//noise hit generator

void smear(int index, double sigmaz, double sigmarf, double R, TClonesArray &cross);//gaussian smearing

bool peakfinder(TH1D* histo, double ampli, int width);//tracklet ambiguity check

#endif
