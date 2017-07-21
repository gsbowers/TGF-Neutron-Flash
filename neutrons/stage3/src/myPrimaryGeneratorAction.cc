#include "myPrimaryGeneratorAction.hh"
#include "myParticleInput.hh"
#include "myParticle.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"

#include "G4Event.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"

myPrimaryGeneratorAction::myPrimaryGeneratorAction(string particleInputFilename)
: G4VUserPrimaryGeneratorAction(), fParticleGun(0), fEnvelopeOrb(0)
{
	// set number of particles to fire each time
	G4int n_particle = 1;
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
	fStage2Output = new myParticleInput(particleInputFilename);
	
}

myPrimaryGeneratorAction::~myPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void myPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  G4double R = 0.6*m; 
  
	/*
  if ( !fEnvelopeOrb ) {
		G4LogicalVolume* envLV
			= G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope");
    if ( envLV ) fEnvelopeOrb = dynamic_cast<G4Orb*>(envLV->GetSolid());
  }

  if ( fEnvelopeOrb ) {
    R = fEnvelopeOrb->GetRadius();
  }
  else {
    G4ExceptionDescription msg;
    msg << "Envelope volume of orb shape not found.\n"; 
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be place at the center.";
    G4Exception("g4PrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002",JustWarning,msg);
	}
	*/

	myParticle particle = fStage2Output->getParticle();
	G4double Energy, ArrivalTime, *StartXYZ, *StartDir;

	Energy = particle.get_Energy();	
	ArrivalTime = particle.get_ArrivalTime();
	StartXYZ = particle.get_StartXYZ();
	StartDir = particle.get_StartDir();

	// set Particle definition 
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* G4particle
		= particleTable->FindParticle(particle.get_ParticleDefinition());
	fParticleGun->SetParticleDefinition(G4particle);


	// set initial particle energy, direction & position
	fParticleGun->SetParticleEnergy(Energy);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(StartDir[0], StartDir[1], StartDir[2]));
  fParticleGun->SetParticlePosition(G4ThreeVector(StartXYZ[0]*R-9.0*2.54*cm, StartXYZ[1]*R, StartXYZ[2]*R));
	fParticleGun->SetParticleTime(ArrivalTime);

  fParticleGun->GeneratePrimaryVertex(anEvent);
}
