#include "myPrimaryGeneratorAction.hh"
#include "myParticleInput.hh"
#include "myParticle.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"

myPrimaryGeneratorAction::myPrimaryGeneratorAction(string particleInputFilename)
{
	// set number of particles to fire each time
	G4int n_particle = 1e3;
	G4cout << G4endl << "in myPrimaryGeneratorAction: n_particle = " << n_particle << G4endl;
	fParticleGun = new G4ParticleGun(n_particle);
	//fParticleGun = new G4GeneralParticleSource(n_particle);

	// default particle kinematic
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	G4ParticleDefinition* particle
		= particleTable->FindParticle(particleName="neutron");
	fParticleGun->SetParticleDefinition(particle);

	// define particleInput
	cout << "reading particle input: " << particleInputFilename << G4endl;
	fStage1Output = new myParticleInput(particleInputFilename);
	
}

myPrimaryGeneratorAction::~myPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void myPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{

	myParticle particle = fStage1Output->getParticle();
	G4double Energy, ArrivalTime, *StartXYZ, *StartDir;

	Energy = particle.get_Energy();	
	ArrivalTime = particle.get_ArrivalTime();
	StartXYZ = particle.get_StartXYZ();
	StartDir = particle.get_StartDir();

	// set initial particle energy, direction & position
	fParticleGun->SetParticleEnergy(Energy);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(StartDir[0], StartDir[1], StartDir[2]));
  fParticleGun->SetParticlePosition(G4ThreeVector(StartXYZ[0], StartXYZ[1], StartXYZ[2]));
	fParticleGun->SetParticleTime(ArrivalTime);

  fParticleGun->GeneratePrimaryVertex(anEvent);
}
