#ifndef myRunAction_h
#define myRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Parameter.hh"
#include "globals.hh"

class G4Run;

class myRunAction : public G4UserRunAction
{
  public:
    myRunAction();
    virtual ~myRunAction();
    
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
		static myRunAction* Instance();
		G4int nbEventInRun;	
	private:
		static myRunAction* fgInstance;
    
};

#endif
