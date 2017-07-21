#include "myDetectorConstruction.hh"
#include "myDetectorSD.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"

#include "G4ios.hh"

myDetectorConstruction::myDetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0), fLogicLgPl(0)
{ }

myDetectorConstruction::~myDetectorConstruction()
{ }

G4VPhysicalVolume* myDetectorConstruction::Construct()
{
  // Get nist material manager
  //
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm; 
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");

	// Option to switch on/off checking of volume overlaps
	//
	G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

	G4Box* solidWorld =
    new G4Box("World",                 //its name 
      0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ); //its size

  G4LogicalVolume* logicWorld = 
    new G4LogicalVolume(solidWorld,    //its solid
                        world_mat,     //its material
                        "World");      //its name

  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,               //no rotation
                      G4ThreeVector(), //at (0,0,0)
                      logicWorld,      //its logical volume
                      "World",         //its name
                      0,               //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking

  //
  // Envelope
  //
  G4Box* solidEnv = 
    new G4Box("Envelope",              //its name
      0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size

  G4LogicalVolume* logicEnv = 
    new G4LogicalVolume(solidEnv,      //its solid
                        env_mat,       //its material
                        "Envelope");   //its name

    new G4PVPlacement(0,               //no rotation
                      G4ThreeVector(), //at (0,0,0)
                      logicEnv,        //its logical volume
                      "Envelope",      //its name
                      logicWorld,      //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking

	

	// 
	//	Sensitive Detectors
	//	
	G4double large_pRmin = 0.0;
	G4double large_pRmax = 2.5 * 2.54 * cm;
	G4double large_pDz =   2.5 * 2.54 * cm;
	G4double large_pSPhi =  0.0; 
	G4double large_pDPhi =  2.0*pi;

  // detector positions and orientations
  G4RotationMatrix yRot90deg = G4RotationMatrix();
  yRot90deg.rotateY(90*deg);

  G4ThreeVector LgPl_position = G4ThreeVector(0,0,0);
  G4Transform3D LgPl_transform = G4Transform3D(yRot90deg, LgPl_position);

	//
	//  Large Plastic (5"x5")
	//
	
  //G4VSensitiveDetector* LgPl = new g4CMOS("LgPl","LgPlHitsCollection");
  //G4SDManager::GetSDMpointer()->AddNewDetector( LgPl );
  //G4Material *Plastic = 
  //   nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  G4Material *LgPl_mat = nist->FindOrBuildMaterial("G4_Pb");

  G4Tubs* solidLgPl
    = new G4Tubs("LgPl",                    //its name
                large_pRmin,                //its size 
                large_pRmax,
                large_pDz,
                large_pSPhi,
                large_pDPhi);

  //G4LogicalVolume* logicLgPl
  fLogicLgPl
    = new G4LogicalVolume(solidLgPl,       //its solid
                          LgPl_mat,        //its material
                          "LgPl");         //its name 

  //logicLgPl->SetSensitiveDetector( LgPl );
      new G4PVPlacement(LgPl_transform, 
                        fLogicLgPl,           //its logical volume
                        "LgPl",              //its name
                        logicEnv,            //its mother (logical) volume
                        false,               //no boolean operations
                        0,                   //its copy number 
                        checkOverlaps);      //checkOverlap
	
  //
  // sets a max step length in the tracker region
  // 	
   
  //G4double maxStep = large_pDz/16;
  //logicLgPl->SetUserLimits(fStepLimit);

  //
  // Visualization
  //
  logicEnv->SetVisAttributes(new G4VisAttributes(G4Colour(0.,1.,1.)));  
  logicWorld->SetVisAttributes(new G4VisAttributes(G4Colour(1.,1.,0.)));  
  fLogicLgPl->SetVisAttributes(new G4VisAttributes(G4Colour(1.,0.0,1.)));  

	return physWorld;

}

void myDetectorConstruction::ConstructSDandField()
{
	//	Sensitive Detectors
	
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	G4String SDname;	

  G4VSensitiveDetector* LgPl 
		= new myDetectorSD(SDname="/LgPl");
	SDman->AddNewDetector(LgPl);
  fLogicLgPl->SetSensitiveDetector(LgPl);
	
}
