//This "class" doesn't inherits from any class since it's not a real class but it only contains the function used in order to find the intersection with the detector system

#ifndef TOOLS_H
#define TOOLS_H

#include "Event.h"
#include "Tools.h"
#include "Layer.h"
#include "Hit.h"
#include "Event.h"
#include "Particle.h"

//x0, y0, z0 = generated vertex coordinates; theta, phi = generated angles; R = radius of the detector
double *hit_point(double x0, double y0, double z0, double theta, double phi, double R);

void detect(Event vtx, Layer L, Particle part, vector<Hit*> &cross, bool b_verbose, bool b_multiscatter, char const *detector);

#endif
