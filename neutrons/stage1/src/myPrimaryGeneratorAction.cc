#include "myPrimaryGeneratorAction.hh"

#include "G4Event.hh"

myPrimaryGeneratorAction::myPrimaryGeneratorAction()
{
	fParticleGun = new G4GeneralParticleSource();
}

myPrimaryGeneratorAction::~myPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void myPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
