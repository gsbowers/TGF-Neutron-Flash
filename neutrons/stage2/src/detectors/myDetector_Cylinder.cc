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

	// Option to switch on/off checking of volume overlaps
	//
	G4bool checkOverlaps = true;

  //
  // World
  //
 
  // TUBS
	//G4double world_pRmin = 0.0;
	//G4double world_pRmax = 0.5 * m;
	//G4double world_pDz =   0.5 * m;
	//G4double world_pSPhi =  0.0; 
	//G4double world_pDPhi =  2.0*pi;
  
  // BOX
  G4double world_hx = 0.5 * m;
  G4double world_hy = 0.5 * m;
  G4double world_hz = 0.5 * m;

  //G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");

  // TUBS
  /*G4Tubs* solidWorld =
    new G4Tubs("World",                 // its name 
               world_pRmin,             // its size
               world_pRmax,
               world_pDz,
               world_pSPhi,
               world_pDPhi);
	*/

  // BOX
  G4Box* solidWorld =
    new G4Box("World",                 // its name 
               world_hx,             // its size
               world_hy,
               world_hz);

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
  // Atmosphere
  //

	/*
  G4double atm_pRmin = world_pRmin;
  G4double atm_pRmax = world_pRmax;
  G4double atm_pDz =  1000 * m - 1 *cm;
  G4double atm_pSPhi = world_pSPhi;   
  G4double atm_pDPhi = world_pDPhi; 

  G4Material* atm_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4Tubs* solidAtm = 
    new G4Tubs("Atmosphere",              //its name
               atm_pRmin,                 // its size 
               atm_pRmax,                 
               atm_pDz,                 
               atm_pSPhi,
               atm_pDPhi);

  G4LogicalVolume* logicAtm = 
    new G4LogicalVolume(solidAtm,      //its solid
                        atm_mat,       //its material
                        "Atmosphere");   //its name

    new G4PVPlacement(0,               //no rotation
                      G4ThreeVector(0, 0, 20 *m + 1*cm), //at (0,0,0)
                      logicAtm,        //its logical volume
                      "Atmosphere",      //its name
                      logicWorld,      //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking

	
	*/

	// 
	//	Sensitive Detectors
	//	
	G4double large_pRmin = 0.0;
	G4double large_pRmax = 2.5 * 2.54 * cm; // radius
	G4double large_pDz =   2.5 * 2.54 * cm; // half length
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
	
  //G4Material *LgPl_mat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  //
  G4int ncomponents, natoms;
	G4double fractionmass, density;

  G4Element* H  = new G4Element("TS_H_of_Water" ,"H" , 1., 1.0079*g/mole);
  G4Element* O  = new G4Element("Oxygen"        ,"O" , 8., 16.00*g/mole);
  G4Element* N  = new G4Element("Nitrogen"      ,"N" , 7., 14.01*g/mole);
  
	// Thermal Water 
  G4Material* H2O = 
  new G4Material("Water_ts", 1.000*g/cm3, ncomponents=2,
                         kStateLiquid, 593*kelvin, 150*bar);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

	// Air
	G4Material* Air = new G4Material("Air ", density=1.290*mg/cm3, ncomponents=2);
	Air->AddElement(N, fractionmass=70*perCent);
	Air->AddElement(O, fractionmass=30*perCent);
	

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
                          Air, //H2O,        //its material
                          "LgPl");         //its name 

      new G4PVPlacement(LgPl_transform, 
                        fLogicLgPl,           //its logical volume
                        "LgPl",              //its name
                        logicWorld,            //its mother (logical) volume
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
  //logicEnv->SetVisAttributes(new G4VisAttributes(G4Colour(0.,1.,1.)));  
  logicWorld->SetVisAttributes(new G4VisAttributes(G4Colour(0.,1.,1.)));  
  fLogicLgPl->SetVisAttributes(new G4VisAttributes(G4Colour(1.,0.0,1.)));  

	G4cout << H2O << G4endl;
	G4cout << Air << G4endl;

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
