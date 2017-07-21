#ifndef myNeutronHPMessenger_h
#define myNeutronHPMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class myNeutronHPphysics;
class G4UIdirectory;
class G4UIcmdWithABool;

class myNeutronHPMessenger: public G4UImessenger
{
	public:
		myNeutronHPMessenger(myNeutronHPphysics*);
	 ~myNeutronHPMessenger();
	
		virtual void SetNewValue(G4UIcommand*, G4String);

	private:
		myNeutronHPphysics* fNeutronPhysics;

		G4UIdirectory* fPhysDir;
		G4UIcmdWithABool* fThermalCmd;

};

#endif
