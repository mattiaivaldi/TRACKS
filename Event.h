//The class inherits from the class Hit

#ifndef EVENT_H
#define EVENT_H

#include "TObject.h"
#include "Hit.h"

class Event : public Hit {

public:

	Event(); //default constructor
	Event(double meanv, double sigmaxy, double sigmaz, const char *distr); //custom constructor
	//meanv is the mean of the vertex coordinates which is 0 for x, y and z, sigmaxy and sigmaz are the RMS of x, y and z
	//distr is the name of the file in which the pseudorapidity distribution is located  
	virtual ~Event(); //destructor 

//member functions

	double GetX() const; //getters 
	double GetY() const;
	double GetZ() const;
	float GetMult() const; 

private:

//data members 	

	double fX0; //primary vertex cooridnates
	double fY0;
	double fZ0;	
	float fMult; //multiplicity 
	
	ClassDef(Event,1)

};
#endif
