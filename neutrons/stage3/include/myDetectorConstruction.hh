#ifndef myDetectorConstruction_h
#define myDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4Element.hh"
#include <string.h>
using namespace std;

class G4VPhysicalVolume;
class G4LogicalVolume;

class myDetectorConstruction : public G4VUserDetectorConstruction 
{
	public:
		myDetectorConstruction(string particleInputPrefix);
		virtual ~myDetectorConstruction();
		
		virtual G4VPhysicalVolume* Construct();
    void ConstructSlab(G4int, G4LogicalVolume**);
		void ConstructSDandField();
		
		G4LogicalVolume* GetScoringVolume() const {return fScoringVolume; }
		
	protected:
		G4LogicalVolume* fScoringVolume;
	  G4LogicalVolume *fLogicSmPl, *fLogicNaI, *fLogicLgPl;	
		G4UserLimits *fStepLimits;

	  G4Element* fTS_H; // Thermal hydrogen for neutron scatter

    G4int fnSlabs;
		G4double *fslabs_z;
		G4double *fslabs_pDz;
		G4double *fslabs_rho;
		G4double *fslabs_T;
		G4double *fslabs_P;

		string *fParticleInputPrefix;
};

#endif
