#ifndef myParticleInput_h
#define myParticleInput_h 1

#include <string>
#include <iostream>
#include <fstream>
#include "myParticle.hh"

using namespace std;

class myParticleInput
{

  public: 
    myParticleInput(string filename);
    ~myParticleInput();
    
    myParticle getParticle();
	
  private:
    // input file variables
    
    long    fNinput;

    G4int* fEvent = NULL;
    G4String* fParticle = NULL;
    G4double* fEnergy = NULL;     // eV
    G4double* fArrivalTime = NULL; // ms
    double* fPosX = NULL;  // km
    double* fPosY = NULL;  // km
    double* fPosZ = NULL;  // km
		double* fDirX = NULL;        
		double* fDirY = NULL;        
		double* fDirZ = NULL;        

		double* fECD = NULL;   // empirical cumulative distribution
		double* fWeight = NULL;       

		long findx;            // pointer into last fECD bin
		long fpcount = 0l;     // current number of particles chosen

		ofstream fout_ECDlist; 
};

#endif
