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
	float GetRaped() const; 
	void Rotate(double rms);

private:

//data members 	

	double fTheta; //theta
	double fPhi; //phi 
	float fRape; //pseudorapidity	
	
	ClassDef(Particle,1)

};
#endif
