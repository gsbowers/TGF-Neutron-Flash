#ifndef myPhysicsList_h
#define myPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class myPhysicsList: public G4VModularPhysicsList
{
public:
  myPhysicsList();
 ~myPhysicsList();

public:
  virtual void ConstructParticle();
  virtual void SetCuts();

};

#endif
