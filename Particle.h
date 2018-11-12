//Class inherits from TObject

#ifndef PARTICLE_H
#define PARTICLE_H

#include "TObject.h"

class Particle : public TObject {

public:

	Particle(); //default constructor
	Particle(const char *distr); //custom constructor
	virtual ~Particle(); //destructor

//member functions

	double GetTheta() const; //getters
	double GetPhi() const;
	float GetRap() const;
	void Rotate(double rms);//multiple scattering setter

private:

//data members

	double fTheta; //theta
	double fPhi; //phi
	float fRap; //pseudorapidity

	ClassDef(Particle,1)

};
#endif
