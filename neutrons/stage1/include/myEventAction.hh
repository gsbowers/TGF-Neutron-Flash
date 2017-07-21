#ifndef myEventAction_h
#define myEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class myRunAction;

class myEventAction : public G4UserEventAction
{
  public:
    myEventAction(myRunAction* runAction);
    virtual ~myEventAction();
   
    virtual void BeginOfEventAction(const G4Event *event);
    virtual void EndOfEventAction(const G4Event *event);

    private:
      myRunAction* fRunAction;
};

#endif
