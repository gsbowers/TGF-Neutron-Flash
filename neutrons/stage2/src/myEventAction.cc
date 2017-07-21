#include "myEventAction.hh"
#include "myRunAction.hh"
#include "myDetectorHit.hh"
#include "myDetectorSD.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"

myEventAction::myEventAction(myRunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction)
{
}

myEventAction::~myEventAction()
{}

void myEventAction::BeginOfEventAction(const G4Event* evt)
{
  G4int evtNb = evt->GetEventID();
  myDetectorSD::Instance()->SetCurrentEventID(evtNb);
  if (evtNb % 10000 == 0){
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
  } 
}

void myEventAction::EndOfEventAction(const G4Event* event)
{
	G4HCofThisEvent* hce = event->GetHCofThisEvent();
	if (!hce)	
	{
		G4ExceptionDescription msg;
		msg << "No hits collection of this event found." << G4endl;
		G4Exception("myEventAction::EndofEventAction()", 
								"myCode001", JustWarning, msg);
		return;
	}
}
