#include "myPhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4StepLimiterPhysics.hh"

#include "myNeutronHPphysics.hh"
#include "myGammaPhysics.hh"

// particles 
//
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4ShortLivedConstructor.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics.hh"
#include "G4IonPhysics.hh"

// \file eventgenerator/exgps/src/g4PhysicsList.cc
// \brief Implementation of the myPhysicsList class

myPhysicsList::myPhysicsList()
:G4VModularPhysicsList()
{
	SetVerboseLevel(1);
	
	//add new units
	new G4UnitDefinition(" millielectronVolt", "meV", "Energy", 1.e-3*eV);	

	//EM Physics 
	RegisterPhysics( new G4EmStandardPhysics_option4(1));

	// Ion Physics
	//RegisterPhysics( new G4IonPhysics(1));

	//Gamma Physics
	RegisterPhysics( new myGammaPhysics("gamma"));

	//Neutron Physics
	RegisterPhysics( new myNeutronHPphysics("myNeutronHP"));

}

myPhysicsList::~myPhysicsList()
{ }

void myPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 
	G4BosonConstructor pBosonConstructor;
	pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor  pBaryonConstructor;
  pBaryonConstructor.ConstructParticle(); 

	G4IonConstructor pIonConstructor;
	pIonConstructor.ConstructParticle();

	G4ShortLivedConstructor pShortLivedConstructor;
	pShortLivedConstructor.ConstructParticle();
}

void myPhysicsList::SetCuts()
{
	SetCutValue(0*mm, "proton");
}
