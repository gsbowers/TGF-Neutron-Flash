#ifndef myNeutronHPphysics_h
#define myNeutronHPphysics_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"

class myNeutronHPMessenger;

class myNeutronHPphysics : public G4VPhysicsConstructor
{
  public:
    myNeutronHPphysics(const G4String& name="neutron");
   ~myNeutronHPphysics();
  
  public:
    virtual void ConstructParticle() { };
    virtual void ConstructProcess();

	public:
		void SetThermalPhysics(G4bool flag) {fThermal = flag;};
		
	private:
		G4bool fThermal;
		myNeutronHPMessenger* fNeutronMessenger;
};

#endif
