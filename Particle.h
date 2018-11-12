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
	void SetTheta(double newtheta){fTheta=newtheta;}
	void SetPhi(double newphi){fPhi=newphi;}
	void Rotate(double rms);
	void Cazzone(double rms);

private:

//data members

	double fTheta; //theta
	double fPhi; //phi
	float fRape; //pseudorapidity

	ClassDef(Particle,1)

};
#endif
