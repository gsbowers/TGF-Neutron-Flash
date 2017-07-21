#include "myActionInitialization.hh"
#include "myPrimaryGeneratorAction.hh"
#include "myRunAction.hh"
#include "myEventAction.hh"

myActionInitialization::myActionInitialization()
 :  G4VUserActionInitialization()
{}

myActionInitialization::~myActionInitialization()
{}

void myActionInitialization::BuildForMaster() const
{
  myRunAction* runAction = new myRunAction;
  SetUserAction(runAction);
}

void myActionInitialization::Build() const
{
  SetUserAction(new myPrimaryGeneratorAction);
  
  myRunAction* runAction = new myRunAction;
  SetUserAction(runAction);
  
  myEventAction* eventAction = new myEventAction(runAction);
  SetUserAction(eventAction);
}
