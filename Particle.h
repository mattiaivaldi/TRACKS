//Class inherits from TObject

#ifndef PARTICLE_H
#define PARTICLE_H

#include "TObject.h"
#include "TH1F.h"

class Particle : public TObject {

public:
	
	Particle(); //default constructor
	Particle(TH1F *distr_rap); //custom constructors
	virtual ~Particle(); //destructor

	//member functions

	double GetTheta() const; //getters
	double GetPhi() const;
	float GetRap() const;
	void Rotate(double rms);//multiple scattering setter

private:

	//data members

	float fRap; //pseudorapidity
	double fTheta; //theta
	double fPhi; //phi

	ClassDef(Particle,1)

};
#endif
