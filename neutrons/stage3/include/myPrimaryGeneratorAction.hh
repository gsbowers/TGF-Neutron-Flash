#ifndef myPrimaryGeneratorAction_h
#define myPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"
#include "myParticleInput.hh"
#include "G4Orb.hh"
#include <fstream>

using namespace std;

class G4GeneralParticleSource;
class G4Event;

class myPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    myPrimaryGeneratorAction(string particleInputFilename);
    virtual ~myPrimaryGeneratorAction();
 
    // method from the base class
    virtual void GeneratePrimaries(G4Event*);
    
    // method to access particle gun
    //const G4GeneralParticleSource* GetParticleGun() 
    const G4ParticleGun* GetParticleGun() 
      const { return fParticleGun; }

  private:
    //G4GeneralParticleSource* fParticleGun;  
    G4ParticleGun* fParticleGun;  
    G4Orb* fEnvelopeOrb;
		myParticleInput* fStage2Output;

    // particle input data
    //ofstream fin;
    

};

#endif
