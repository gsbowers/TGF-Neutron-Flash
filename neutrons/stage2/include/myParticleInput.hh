#ifndef myParticleInput_h
#define myParticleInput_h 1

#include <string>
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

    G4int* fEvent;
    G4double* fEnergy;     // eV
    G4double* fArrivalTime; // ms
    double* fPosX;         // km
    double* fPosY;         // km
    double* fPosZ;         // km
		double* fDirX;        
		double* fDirY;        
		double* fDirZ;        
		long findx;
};

#endif
