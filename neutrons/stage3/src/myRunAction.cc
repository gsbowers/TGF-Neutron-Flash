#include "myRunAction.hh"
#include "myPrimaryGeneratorAction.hh"
#include "myDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4ParameterManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"

myRunAction::myRunAction()
{
	fgInstance = this;
}

myRunAction::~myRunAction()
{
}

myRunAction* myRunAction::fgInstance = 0;

myRunAction* myRunAction::Instance()
{
	return fgInstance;
}

void myRunAction::BeginOfRunAction(const G4Run* run)
{
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
  // inform the runManager to save random number seed
	nbEventInRun = run->GetNumberOfEventToBeProcessed();
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset parameters to their initial values
  G4ParameterManager* parameterManager = G4ParameterManager::Instance();
  parameterManager->Reset();
}

void myRunAction::EndOfRunAction(const G4Run* run) 
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  G4ParameterManager* parameterManager = G4ParameterManager::Instance();
  parameterManager->Merge();

  // Print
  //  
  if (IsMaster()){
    G4cout
      << G4endl
      << "------------------End of Global Run---------------------"
      << G4endl;
  }
  else {
    G4cout
      << G4endl
      << "------------------End of Local Run----------------------"
      << G4endl;
  }
}
