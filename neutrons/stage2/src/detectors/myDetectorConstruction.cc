#include "myDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"

myDetectorConstruction::myDetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
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
  // Visualization
  //
  logicEnv->SetVisAttributes(new G4VisAttributes(G4Colour(0.,1.,1.)));  
  logicWorld->SetVisAttributes(new G4VisAttributes(G4Colour(1.,1.,0.)));  
  

	return physWorld;

}
