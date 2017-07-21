#include "myDetectorConstruction.hh"
#include "myDetectorSD.hh"

#include "G4UserLimits.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"

#include "G4SDManager.hh"
#include "G4GeometryManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4Material.hh"

#include "G4ios.hh"
#include <string>

void ConstructAlEnclosure(G4LogicalVolume *, G4LogicalVolume *, G4double, G4double, G4double);
void ConstructPelicanCase(G4LogicalVolume*);

myDetectorConstruction::myDetectorConstruction(string particleInputPrefix)
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ 

 fParticleInputPrefix = new string(particleInputPrefix);
 // global definitions of atmosphere geometry to be used 
 // by ConstructSlab in myDetectorConstruction::Construct 
 // altitude profile is calculated using makeslabs.pro. 
 // if you wish to change profile, copy and past output 
 // of modified makeslabs.pro into section below
 
 //########## begin output of makeslabs.pro

 // specify number of slabs in atmosphere
 fnSlabs = 16;

 // slab mid points
 G4double slabs_z[] = 
 {-1917.7*m, -1751.7*m, -1582.9*m, -1411.2*m, -1236.5*m, -1058.7*m, -877.7*m, -693.4*m, -505.5*m, -314.1*m, -119.0*m, 80.0*m, 283.1*m, 490.4*m, 702.1*m, 918.3*m};

 // slab half-heights
 G4double slabs_pDz[] = 
 {82.3*m, 83.7*m, 85.1*m, 86.6*m, 88.1*m, 89.7*m, 91.3*m, 93.0*m, 94.8*m, 96.6*m, 98.5*m, 100.5*m, 102.6*m, 104.7*m, 107.0*m, 109.2*m};

 // slab densities
 G4double slabs_rho[] = 
 {1.2154*mg/cm3, 1.1961*mg/cm3, 1.1767*mg/cm3, 1.1573*mg/cm3, 1.1377*mg/cm3, 1.1181*mg/cm3, 1.0984*mg/cm3, 1.0786*mg/cm3, 1.0587*mg/cm3, 1.0387*mg/cm3, 1.0187*mg/cm3, 0.9985*mg/cm3, 0.9782*mg/cm3, 0.9579*mg/cm3, 0.9375*mg/cm3, 0.9169*mg/cm3};

// slab pressures
 G4double slabs_P[] = 
 {100341.381582*pascal, 98378.064526*pascal, 96414.553682*pascal, 94449.922664*pascal, 92485.123261*pascal, 90519.033396*pascal, 88552.241873*pascal, 86585.120224*pascal, 84617.534693*pascal, 82650.852587*pascal, 80683.279004*pascal, 78715.644643*pascal, 76748.894940*pascal, 74782.709118*pascal, 72817.373905*pascal, 70852.632561*pascal};

 // slab temperatures
 G4double slabs_T[] = 
{287.6155*kelvin, 286.5357*kelvin, 285.4388*kelvin, 284.3222*kelvin, 283.1883*kelvin, 282.0329*kelvin, 280.8560*kelvin, 279.6585*kelvin, 278.4373*kelvin, 277.1946*kelvin, 275.9268*kelvin, 274.6347*kelvin, 273.3147*kelvin, 271.9682*kelvin, 270.5940*kelvin, 269.1895*kelvin};

 //########## end output of makeslabs.pro

  // initialize persistent arrays  
  fslabs_z = new double[fnSlabs];
  std::memcpy(fslabs_z, slabs_z, sizeof(slabs_z));
  fslabs_pDz = new double[fnSlabs];
  std::memcpy(fslabs_pDz, slabs_pDz, sizeof(slabs_pDz));
  fslabs_rho = new double[fnSlabs];
  std::memcpy(fslabs_rho, slabs_rho, sizeof(slabs_rho));
  fslabs_T = new double[fnSlabs];
  std::memcpy(fslabs_T, slabs_T, sizeof(slabs_T));
  fslabs_P = new double[fnSlabs];
  std::memcpy(fslabs_P, slabs_P, sizeof(slabs_P));
}

myDetectorConstruction::~myDetectorConstruction()
{ 
  delete[] fslabs_z;
  delete[] fslabs_pDz;
  delete[] fslabs_rho;
  delete[] fslabs_T;
  delete[] fslabs_P;
}

G4VPhysicalVolume* myDetectorConstruction::Construct()
{

  G4bool checkOverlaps = true;
  
  // Get nist material manager
  //
  G4NistManager* nist = G4NistManager::Instance();

  // SetWorldMaximumExtent for scaling checkOverlap precision
  G4double pRmax = 1.0 * m;
  //G4double pDz = 3.0 * m;
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(pRmax);
  
  //
  // World
  //

  // Material
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");

  // Dimensions - Orb
  G4double world_pRmax = 1.1*pRmax;

  // Dimensions - Tubs
	/*
  G4double world_pRmin = 0;
  G4double world_pRmax = 1.1*pRmax;
  G4double world_pDz =   1.1*pDz; 
	G4double world_pSPhi = 0.0;
	G4double world_pDPhi = 2.0*pi;
	*/

  // Geometry - Orb
  G4Orb* solidWorld =
    new G4Orb("World",                 // its name 
               world_pRmax);           // its size

  // Geometry - Tubs 
	/*
  G4Tubs* solidWorld =
    new G4Tubs("World",                 // its name 
               world_pRmin,               // its size
               world_pRmax,
               world_pDz, 
               world_pSPhi, 
               world_pDPhi);              
	*/

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
                      true);           //overlaps checking

	// 
	// Particle Input sphere
	// 

	/*
  // Material
  G4Material* input_mat = nist->FindOrBuildMaterial("G4_Galactic");

  // Dimensions - Orb
	// from R in myPrimaryGenerator.cc
  G4double input_pRmax = 0.6*m; 

  // Geometry - Tubs 
  G4Orb* solidInput =
    new G4Orb("Input",                 // its name 
               input_pRmax);

  G4LogicalVolume* logicInput = 
    new G4LogicalVolume(solidInput,    //its solid
                        input_mat,     //its material
                        "Input");      //its name

    new G4PVPlacement(0,               //no rotation
                      G4ThreeVector(-9*2.54*cm,0,0), //at (0,0,0)
                      logicInput,      //its logical volume
                      "Input",         //its name
                      logicWorld,      //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      true);           //overlaps checking
	*/

  //
  // Envelope
  //

  // Material
  // Thermal Hydrogen in Water
  // N.B. Contribution of O from H2O is negligible compared to O2 in air
  fTS_H = new G4Element("TS_H_of_Water" ,"H" , 1., 1.0079*g/mole);
  G4cout << fTS_H << G4endl;
  
  G4double density, temperature, pressure;
  G4int slabNo = 0;

  // Get material parameters from atmosphere deffinition
  density = fslabs_rho[slabNo];
  temperature = fslabs_T[slabNo];
  pressure = fslabs_P[slabNo];

  char materialName[256];
  sprintf(materialName, "SlabMat%02i", slabNo);
  G4Material *Air = nist->ConstructNewGasMaterial(materialName, "G4_AIR", temperature, pressure);
  sprintf(materialName, "SlabMat%02i", slabNo);
  G4double rho_TS_H = 2.778e-4 * 2.0/18.0; // mg/cm3
  G4double rho_air = density/(mg/cm3); 
  G4double fraction_Air = rho_air/(rho_air + rho_TS_H);
  G4double fraction_TS_H = rho_TS_H/(rho_air + rho_TS_H);
  G4Material* SlabMat = new G4Material(materialName, density, 2, kStateGas, temperature, pressure);  
  SlabMat->AddMaterial(Air, fraction_Air);
  SlabMat->AddElement(fTS_H, fraction_TS_H);

  G4Material *envelope_mat = SlabMat;
  G4cout << envelope_mat << G4endl;
 
  // Dimensions - Orb
  G4double envelope_pRmax = pRmax;

  // Dimensions - Tubs
	/*
  G4double envelope_pRmax = pRmax;
  G4double envelope_pRmin = 0.0;
  G4double envelope_pDz   = pDz;
  G4double envelope_pSPhi = 0.0;
  G4double envelope_pDPhi = 2.0*pi;
	*/

  // Geometry - Orb
  G4Orb* solidEnvelope =
    new G4Orb("Envelope",                 // its name 
               envelope_pRmax);            // its size

  // Geometry - Tubs
	/*
  G4Tubs* solidEnvelope =
    new G4Tubs("Envelope",                 // its name 
               envelope_pRmin,            // its size
               envelope_pRmax,
               envelope_pDz,              
               envelope_pSPhi,             
               envelope_pDPhi);              
	*/

  G4LogicalVolume* logicEnvelope = 
    new G4LogicalVolume(solidEnvelope,    //its solid
                        envelope_mat,     //its material
                        "Envelope");      //its name

   new G4PVPlacement(0,                  //no rotation
                     G4ThreeVector(),    //at (0,0,0)
                     logicEnvelope,      //its logical volume
                     "Envelope",         //its name
                     logicWorld,         //its mother volume
                     false,              //no boolean operation
                     0,                  //copy number
                     true);              //overlaps checking

  //
  // Ground
  //
	/*

  // Material
  G4Material* ground_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
  G4cout << ground_mat << G4endl;

  // Dimensions - Box
  G4double ground_pRmin = 0;
  G4double ground_pRmax = pRmax;
  G4double ground_pDz =   pDz/2.0-0.6*m/2.0; 
	G4double ground_pSPhi = 0.0;
	G4double ground_pDPhi = 2.0*pi;

  // Geometry - Tubs 
  G4Tubs* solidGround =
    new G4Tubs("Ground",                 // its name 
               ground_pRmin,               // its size
               ground_pRmax,
               ground_pDz, 
               ground_pSPhi, 
               ground_pDPhi);              

  G4LogicalVolume* logicGround = 
    new G4LogicalVolume(solidGround,    //its solid
                        ground_mat,     //its material
                        "Ground");      //its name

    new G4PVPlacement(0,               //no rotation
                      G4ThreeVector(0,0,-pDz+ground_pDz), //at (0,0,0)
                      logicGround,     //its logical volume
                      "Ground",        //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      true);           //overlaps checking

	*/

  //
  // GODOT 
  // ##############################################################

  // 
  //  Polyurethane Block
  //  http://www.skandiaupholsterysupplies.com/hr-poly-foam-hr-150/
  //  https://en.wikipedia.org/wiki/Polyurethane
  // 

	// Material
  G4double z, a;
  G4String name, symbol;
  G4int ncomponents, natoms; 

  a = 1.01*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H", z=1., a);

  a = 12.02*g/mole;
  G4Element* elC = new G4Element(name="Carbon",   symbol="C", z=6., a);

  a = 14.01*g/mole;  
  G4Element * elN = new G4Element(name="Nitrogen", symbol="N", z=7., a);

  a = 16.00*g/mole;  
  G4Element * elO = new G4Element(name="Oxygen"  , symbol="O", z=8., a);

  density = 0.07256 * g/cm3;
  G4Material *polyurethane = new G4Material(name="Polyurethane", density, ncomponents=4);
  polyurethane->AddElement(elH, natoms=16);
  polyurethane->AddElement(elC, natoms=17);
  polyurethane->AddElement(elN, natoms=2);
  polyurethane->AddElement(elO, natoms=4);

  G4Material *foam_mat = polyurethane;
  G4cout << foam_mat << G4endl;

	// Dimensions - Box
  static const G4double inch = 2.54*cm;
  G4double foam_hx = 10.0/2.0 * inch;
  G4double foam_hy = 17.0/2.0 * inch;
  G4double foam_hz = 10.0/2.0 * inch;

	// Geometry - Box
  G4Box* solidFoam = 
    new G4Box("Foam",                //its name
               foam_hx,              // its size 
               foam_hy,                 
               foam_hz);                 

  G4LogicalVolume* logicFoam = 
    new G4LogicalVolume(solidFoam,   //its solid
                        foam_mat,    //its material
                        "Foam");     //its name

    new G4PVPlacement(0,               //no rotation
                      G4ThreeVector(0, 0, 0), //
                      logicFoam,       //its logical volume
                      "Foam",          //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking

  //
  //  Aluminum Front Plate
  // 

  // Material
  G4Material *AlPlate_mat 
    = nist->FindOrBuildMaterial("G4_Al");

  // Dimensions - Box
  G4double AlPlate_hx = 0.25/2.0 * inch; //half-length
  G4double AlPlate_hy = 19.0/2.0 * inch; //half-length 
  G4double AlPlate_hz = 10.5/2.0 * inch; //half-length

  // Position
  G4ThreeVector AlPlate_position = G4ThreeVector(foam_hx+AlPlate_hx, 0,0);

	// Geometry - Box
  G4Box* solidAlPlate = 
    new G4Box("AlPlate",                //its name
               AlPlate_hx,              // its size 
               AlPlate_hy,                 
               AlPlate_hz);                 

  G4LogicalVolume* logicAlPlate = 
    new G4LogicalVolume(solidAlPlate,   //its solid
                        AlPlate_mat,    //its material
                        "AlPlate");     //its name

    new G4PVPlacement(0,               //no rotation
                      AlPlate_position,//it's position
                      logicAlPlate,    //its logical volume
                      "AlPlate",       //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking

  
  //
  // Aluminum Enclosure
  //   
  ConstructAlEnclosure(logicFoam, logicEnvelope, foam_hx, foam_hy, foam_hz);

	//
	// Pelican 1730 Large Case
	// 
	ConstructPelicanCase(logicEnvelope);

  // 
  //  Sensitive Detectors
  //  

  // Matrial definitions
  G4Material *SmPl_mat 
    = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  G4Material *LgPl_mat 
    = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  G4Material *NaI_mat 
    = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");

	
	// Dimensions - Tubs
  // large scintillator dimensions
  G4double large_pRmin = 0.0;
  G4double large_pRmax = 2.5 * 2.54 * cm; // radius
  G4double large_pDz =   2.5 * 2.54 * cm; // half length
  G4double large_pSPhi = 0.0; 
  G4double large_pDPhi = 2.0*pi;

  // small scintillator dimensions
  G4double small_pRmin = 0.0;
  G4double small_pRmax = 0.5 * 2.54 * cm; // radius
  G4double small_pDz =   0.5 * 2.54 * cm; // half length
  G4double small_pSPhi = 0.0; 
  G4double small_pDPhi = 2.0*pi;

  // detector positions and orientations
  G4RotationMatrix yRot90deg = G4RotationMatrix();
  yRot90deg.rotateY(90*deg);

  G4ThreeVector SmPl_position = G4ThreeVector(5*cm+0.5*inch, 3*cm,8*cm);
  G4Transform3D SmPl_transform = G4Transform3D(yRot90deg, SmPl_position);

  G4ThreeVector LgPl_position = G4ThreeVector(0.5*inch, 10*cm,-1*cm);
  G4Transform3D LgPl_transform = G4Transform3D(yRot90deg, LgPl_position);

  G4ThreeVector NaI_position = G4ThreeVector(0.5*inch, -10*cm,-1*cm);
  G4Transform3D NaI_transform = G4Transform3D(yRot90deg, NaI_position);

	// Geometries - Tubs

  //
  //  Small Plastic (1"x1")
  //
  
  G4Tubs* solidSmPl
    = new G4Tubs("SmPl",                  //its name
                small_pRmin,              //its size 
                small_pRmax,
                small_pDz,
                small_pSPhi,
                small_pDPhi);

  //G4LogicalVolume* logicSmPl
  fLogicSmPl
    = new G4LogicalVolume(solidSmPl,      //its solid
                          SmPl_mat,       //its material
                          "SmPl");        //its name 

      new G4PVPlacement(SmPl_transform, 
                        fLogicSmPl,       //its logical volume
                        "SmPl",           //its name
                        logicFoam,        //its mother (logical) volume
                        false,            //no boolean operations
                        0,                //its copy number 
                        checkOverlaps);   //checkOverlap

  //
  //  Large Plastic (5"x5")
  //
  
  G4Tubs* solidLgPl
    = new G4Tubs("LgPl",                  //its name
                large_pRmin,              //its size 
                large_pRmax,
                large_pDz,
                large_pSPhi,
                large_pDPhi);

  //G4LogicalVolume* logicLgPl
  fLogicLgPl
    = new G4LogicalVolume(solidLgPl,      //its solid
                          LgPl_mat,       //its material
                          "LgPl");        //its name 

      new G4PVPlacement(LgPl_transform, 
                        fLogicLgPl,       //its logical volume
                        "LgPl",           //its name
                        logicFoam,        //its mother (logical) volume
                        false,            //no boolean operations
                        0,                //its copy number 
                        checkOverlaps);   //checkOverlap

  //
  //  Large NaI(Tl) (5"x5")
  //
  
  G4Tubs* solidNaI
    = new G4Tubs("NaI",                   //its name
                large_pRmin,              //its size 
                large_pRmax,
                large_pDz,
                large_pSPhi,
                large_pDPhi);

  //G4LogicalVolume* logicNaI
  fLogicNaI
    = new G4LogicalVolume(solidNaI,       //its solid
                          NaI_mat,        //its material
                          "NaI");         //its name 

      new G4PVPlacement(NaI_transform, 
                        fLogicNaI,        //its logical volume
                        "NaI",            //its name
                        logicFoam,        //its mother (logical) volume
                        false,            //no boolean operations
                        0,                //its copy number 
                        checkOverlaps);   //checkOverlap


  //
  // Visualization
  //
  logicWorld->  SetVisAttributes(new G4VisAttributes(G4Colour(0.,1.,1.))); 
  logicEnvelope->SetVisAttributes(new G4VisAttributes(G4Colour(0.,0.,1.))); 
  //logicGround-> SetVisAttributes(new G4VisAttributes(G4Colour(1.,1.,0.))); 
  logicWorld->  SetVisAttributes(G4VisAttributes::GetInvisible());
  //logicEnvelope->  SetVisAttributes(G4VisAttributes::GetInvisible());
  logicFoam->   SetVisAttributes(new G4VisAttributes(G4Colour(1.,0,0)));  
  logicAlPlate->SetVisAttributes(new G4VisAttributes(G4Colour(1.,0.5,0))); 
  fLogicSmPl->  SetVisAttributes(new G4VisAttributes(G4Colour(1.,0.,1.))); 
  fLogicLgPl->  SetVisAttributes(new G4VisAttributes(G4Colour(1.,0.,1.))); 
  fLogicNaI->   SetVisAttributes(new G4VisAttributes(G4Colour(1.,1.,0.)));  

  return physWorld;

}

void myDetectorConstruction::ConstructSDandField()
{
  //  Sensitive Detectors
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  char SDName[128];

  sprintf(SDName, "SmPl_%s", fParticleInputPrefix->c_str());
  G4VSensitiveDetector* SmPl 
    = new myDetectorSD(SDName, "SmPlHitsCollection");
  SDman->AddNewDetector(SmPl);
  fLogicSmPl->SetSensitiveDetector(SmPl);

  sprintf(SDName, "LgPl_%s", fParticleInputPrefix->c_str());
  G4VSensitiveDetector* LgPl 
    = new myDetectorSD(SDName, "LgPlHitsCollection");
  SDman->AddNewDetector(LgPl);
  fLogicLgPl->SetSensitiveDetector(LgPl);

  sprintf(SDName, "NaI_%s", fParticleInputPrefix->c_str());
  G4VSensitiveDetector* NaI 
    = new myDetectorSD(SDName, "NaIHitsCollection");
  SDman->AddNewDetector(NaI);
  fLogicNaI->SetSensitiveDetector(NaI);
}

void ConstructPelicanCase(G4LogicalVolume *logicEnvelope){

	//
	// Pelican 1730 Large Case 
	// http://godot.jp/documents/Pelican-1730-protector-Large%20Case.pdf
  G4NistManager* nist = G4NistManager::Instance();
  static const G4double inch = 2.54*cm;
  bool checkOverlaps = true;

	// Material 
	G4Material* pelican_mat = 
    nist->FindOrBuildMaterial("G4_POLYPROPYLENE");
	G4cout << pelican_mat << G4endl;	
	
	// Dimensions
	G4double box_hx = 34.00/2.0 * inch;
	G4double box_hy = 24.13/2.0 * inch;
	G4double box_hz = 12.50/2.0 * inch;
	G4double thick  = 0.0597 * inch;

	// Geometry
	G4Box *outerCase = 
		new G4Box("OuterCase", 
							box_hx+thick/2.0, 
							box_hy+thick/2.0, 
							box_hz+thick/2.0);

	G4Box *innerCase = 
		new G4Box("InnerCase", 
							box_hx, 
							box_hy, 
							box_hz);

	G4SubtractionSolid *solidPelican = 
		new G4SubtractionSolid("PelicanCase", outerCase, innerCase);

  G4LogicalVolume* logicPelican = 
    new G4LogicalVolume(solidPelican,   //its solid
                        pelican_mat,    //its material
                        "Pelican");     //its name

    new G4PVPlacement(0,               //no rotation
                      G4ThreeVector(-9.0*inch,0,0), //it's position
                      logicPelican,    //its logical volume
                      "Pelican",       //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking
	
  logicPelican->  SetVisAttributes(new G4VisAttributes(G4Colour(0.,1.,0.))); 

}

void ConstructAlEnclosure(G4LogicalVolume *logicFoam, G4LogicalVolume *logicEnvelope, G4double foam_hx, G4double foam_hy, G4double foam_hz)
{
  //
  // Al Frame 
  // Originally 1" OD x 1/8" Wall square tubing w/ rho
  // Approximating as 1" OD solid tube w/ rho'
  // volume of 1" OD x 1/8" Wall = 7/16"^2*l 
  // volume of 1" OD Solid tube  = 1"^2*l
  // rho' = volume/volume' * rho = 7/16*rho

  G4NistManager* nist = G4NistManager::Instance();
  static const G4double inch = 2.54*cm;
  bool checkOverlaps = true;

  // Al Tubing material
  G4double density = 2.700*7.0/16.0 * g/cm3;
  G4double a = 26.98*g/mole;
  G4int z;
  G4String name;
  G4Material *AlTubing_mat = new G4Material(name="AlTubing", z=13, a, density);

  // Al Tubing Dimensions
  G4double AlTube_hx = 1.0/2.0 * inch;  //half-length
  G4double AlTube_hz = 1.0/2.0 * inch;  //half-length
  G4double AlTube_hy = 0 * inch; //half-length 

  // Al Tubing Position
  G4ThreeVector AlTubing_Position = G4ThreeVector();

  //#######################################################################
  //FRAME INSIDE FOAM

  // ++++++++++++++++++++++++++++
  // Al Tubing FORWARD TOP
  AlTube_hy = foam_hy;
  AlTubing_Position = G4ThreeVector(foam_hx-AlTube_hx, 0, foam_hz-AlTube_hz);

  G4Box* solidAlTubing_FORWARD_TOP =  
    new G4Box("AlTubing_FORWARD_TOP",               //its name
               AlTube_hx,                           // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_FORWARD_TOP =
    new G4LogicalVolume(solidAlTubing_FORWARD_TOP,  //its solid
                        AlTubing_mat,               //its material
                        "AlTubing_FORWARD_TOP");    //its name

    new G4PVPlacement(0,                            //no rotation
                      AlTubing_Position,            //it's position
                      logicAlTubing_FORWARD_TOP,    //its logical volume
                      "AlTubing_FORWARD_TOP",       //its name
                      logicFoam,                    //its mother volume
                      false,                        //no boolean operation
                      0,                            //copy number
                      checkOverlaps);               //overlaps checking

  // ++++++++++++++++++++++++++++
  // Al Tubing FORWARD BOTTOM
  AlTube_hy = foam_hy;
  AlTubing_Position = G4ThreeVector(foam_hx-AlTube_hx, 0, -(foam_hz-AlTube_hz));

  G4Box* solidAlTubing_FORWARD_BOTTOM =  
    new G4Box("AlTubing_FORWARD_BOTTOM",            //its name
               AlTube_hx,                           // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_FORWARD_BOTTOM =
    new G4LogicalVolume(solidAlTubing_FORWARD_BOTTOM,//its solid
                        AlTubing_mat,               //its material
                        "AlTubing_FORWARD_BOTTOM"); //its name

    new G4PVPlacement(0,                            //no rotation
                      AlTubing_Position,            //it's position
                      logicAlTubing_FORWARD_BOTTOM, //its logical volume
                      "AlTubing_FORWARD_BOTTOM",    //its name
                      logicFoam,                    //its mother volume
                      false,                        //no boolean operation
                      0,                            //copy number
                      checkOverlaps);               //overlaps checking

  // ++++++++++++++++++++++++++++
  // Al Tubing FORWARD RIGHT
  AlTube_hy = foam_hz-2.0*AlTube_hz;
  AlTubing_Position 
    = G4ThreeVector(foam_hx-AlTube_hx, foam_hy-AlTube_hx,0);
  G4RotationMatrix xRot90deg = G4RotationMatrix();
  xRot90deg.rotateX(90*deg);
  G4Transform3D AlTubing_Transform 
    = G4Transform3D(xRot90deg, AlTubing_Position);

  G4Box* solidAlTubing_FORWARD_RIGHT =  
    new G4Box("AlTubing_FORWARD_RIGHT",             //its name
               AlTube_hx,                           // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_FORWARD_RIGHT =
    new G4LogicalVolume(solidAlTubing_FORWARD_RIGHT,//its solid
                        AlTubing_mat,               //its material
                        "AlTubing_FORWARD_RIGHT");  //its name

    new G4PVPlacement(AlTubing_Transform,           //it's position
                      logicAlTubing_FORWARD_RIGHT,  //its logical volume
                      "AlTubing_FORWARD_RIGHT",     //its name
                      logicFoam,                    //its mother volume
                      false,                        //no boolean operation
                      0,                            //copy number
                      checkOverlaps);               //overlaps checking

  // ++++++++++++++++++++++++++++
  // Al Tubing FORWARD LEFT
  AlTube_hy = foam_hz-2.0*AlTube_hz;
  AlTubing_Position 
    = G4ThreeVector(foam_hx-AlTube_hx,-(foam_hy-AlTube_hx),0);
  AlTubing_Transform = G4Transform3D(xRot90deg, AlTubing_Position);

  G4Box* solidAlTubing_FORWARD_LEFT =  
    new G4Box("AlTubing_FORWARD_LEFT",              //its name
               AlTube_hx,                           // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_FORWARD_LEFT =
    new G4LogicalVolume(solidAlTubing_FORWARD_LEFT, //its solid
                        AlTubing_mat,               //its material
                        "AlTubing_FORWARD_LEFT");   //its name

    new G4PVPlacement(AlTubing_Transform,           //its position
                      logicAlTubing_FORWARD_LEFT,   //its logical volume
                      "AlTubing_FORWARD_LEFT",      //its name
                      logicFoam,                    //its mother volume
                      false,                        //no boolean operation
                      0,                            //copy number
                      checkOverlaps);               //overlaps checking

  // ++++++++++++++++++++++++++++
  // Al Tubing RIGHT TOP
  AlTube_hy = foam_hx-AlTube_hx;
  AlTubing_Position 
    = G4ThreeVector(-AlTube_hx, foam_hy-AlTube_hx,foam_hz-AlTube_hx);
  G4RotationMatrix zRot90deg = G4RotationMatrix();
  zRot90deg.rotateZ(90*deg);
  AlTubing_Transform 
    = G4Transform3D(zRot90deg, AlTubing_Position);

  G4Box* solidAlTubing_RIGHT_TOP =  
    new G4Box("AlTubing_RIGHT_TOP",                 //its name
               AlTube_hx,                           // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_RIGHT_TOP =
    new G4LogicalVolume(solidAlTubing_RIGHT_TOP,    //its solid
                        AlTubing_mat,               //its material
                        "AlTubing_RIGHT_TOP");      //its name

    new G4PVPlacement(AlTubing_Transform,           //it's position
                      logicAlTubing_RIGHT_TOP,      //its logical volume
                      "AlTubing_RIGHT_TOP",         //its name
                      logicFoam,                    //its mother volume
                      false,                        //no boolean operation
                      0,                            //copy number
                      checkOverlaps);               //overlaps checking

  // ++++++++++++++++++++++++++++
  // Al Tubing RIGHT BOTTOM
  AlTube_hy = foam_hx-AlTube_hx;
  AlTubing_Position 
    = G4ThreeVector(-AlTube_hx, foam_hy-AlTube_hx,-(foam_hz-AlTube_hx));
  AlTubing_Transform 
    = G4Transform3D(zRot90deg, AlTubing_Position);

  G4Box* solidAlTubing_RIGHT_BOTTOM =  
    new G4Box("AlTubing_RIGHT_BOTTOM",              //its name
               AlTube_hx,                           // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_RIGHT_BOTTOM =
    new G4LogicalVolume(solidAlTubing_RIGHT_BOTTOM, //its solid
                        AlTubing_mat,               //its material
                        "AlTubing_RIGHT_BOTTOM");   //its name

    new G4PVPlacement(AlTubing_Transform,           //it's position
                      logicAlTubing_RIGHT_BOTTOM,   //its logical volume
                      "AlTubing_RIGHT_BOTTOM",      //its name
                      logicFoam,                    //its mother volume
                      false,                        //no boolean operation
                      0,                            //copy number
                      checkOverlaps);               //overlaps checking

  // ++++++++++++++++++++++++++++
  // Al Tubing LEFT TOP
  AlTube_hy = foam_hx-AlTube_hx;
  AlTubing_Position 
    = G4ThreeVector(-AlTube_hx, -(foam_hy-AlTube_hx), foam_hz-AlTube_hx);
  AlTubing_Transform 
    = G4Transform3D(zRot90deg, AlTubing_Position);

  G4Box* solidAlTubing_LEFT_TOP =  
    new G4Box("AlTubing_LEFT_TOP",                  //its name
               AlTube_hx,                           // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_LEFT_TOP =
    new G4LogicalVolume(solidAlTubing_LEFT_TOP,     //its solid
                        AlTubing_mat,               //its material
                        "AlTubing_LEFT_TOP");       //its name

    new G4PVPlacement(AlTubing_Transform,           //it's position
                      logicAlTubing_LEFT_TOP,       //its logical volume
                      "AlTubing_LEFT_TOP",          //its name
                      logicFoam,                    //its mother volume
                      false,                        //no boolean operation
                      0,                            //copy number
                      checkOverlaps);               //overlaps checking

  // ++++++++++++++++++++++++++++
  // Al Tubing LEFT BOTTOM
  AlTube_hy = foam_hx-AlTube_hx;
  AlTubing_Position 
    = G4ThreeVector(-AlTube_hx, -(foam_hy-AlTube_hx),-(foam_hz-AlTube_hx));
  AlTubing_Transform 
    = G4Transform3D(zRot90deg, AlTubing_Position);

  G4Box* solidAlTubing_LEFT_BOTTOM =  
    new G4Box("AlTubing_LEFT_BOTTOM",               //its name
               AlTube_hx,                           // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_LEFT_BOTTOM =
    new G4LogicalVolume(solidAlTubing_LEFT_BOTTOM,  //its solid
                        AlTubing_mat,               //its material
                        "AlTubing_LEFT_BOTTOM");    //its name

    new G4PVPlacement(AlTubing_Transform,           //it's position
                      logicAlTubing_LEFT_BOTTOM,    //its logical volume
                      "AlTubing_LEFT_BOTTOM",       //its name
                      logicFoam,                    //its mother volume
                      false,                        //no boolean operation
                      0,                            //copy number
                      checkOverlaps);               //overlaps checking

  //#######################################################################
  // Al Side Faces

  G4Material *AlFace_mat 
    = nist->FindOrBuildMaterial("G4_Al");

  // ++++++++++++++++++++++++++++
  // TOP FACE
  G4double AlFace_hx = 28.0/2.0 * inch;
  G4double AlFace_hy = foam_hy;
  G4double AlFace_hz = 0.075 * inch;

  G4ThreeVector AlFace_Position = 
    G4ThreeVector(-(AlFace_hx-foam_hx), 0, foam_hz+AlFace_hz);

  G4Box* solidAlFace_TOP =
    new G4Box("AlFace_TOP",           //its name
              AlFace_hx,              // its size 
              AlFace_hy,                 
              AlFace_hz);                 

  G4LogicalVolume* logicAlFace_TOP =
    new G4LogicalVolume(solidAlFace_TOP,   //its solid
                        AlFace_mat,    //its material
                        "AlFace_TOP"); //its name

    new G4PVPlacement(0,               // no rotation
                      AlFace_Position, 
                      logicAlFace_TOP, //its logical volume
                      "AlFace_TOP",    //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking

  // ++++++++++++++++++++++++++++
  // BOTTOM FACE 
  AlFace_Position = 
    G4ThreeVector(-(AlFace_hx-foam_hx), 0, -(foam_hz+AlFace_hz));

  G4Box* solidAlFace_BOTTOM =
    new G4Box("AlFace_BOTTOM",           //its name
              AlFace_hx,              // its size 
              AlFace_hy,                 
              AlFace_hz);                 

  G4LogicalVolume* logicAlFace_BOTTOM =
    new G4LogicalVolume(solidAlFace_BOTTOM,   //its solid
                        AlFace_mat,    //its material
                        "AlFace_BOTTOM"); //its name

    new G4PVPlacement(0,               // no rotation
                      AlFace_Position, 
                      logicAlFace_BOTTOM, //its logical volume
                      "AlFace_BOTTOM",    //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking

  // ++++++++++++++++++++++++++++
  // RIGHT FACE 
  AlFace_hz = foam_hz;  
  AlFace_hy = 0.075 * inch;
  AlFace_Position = 
    G4ThreeVector(-(AlFace_hx-foam_hx), foam_hy+AlFace_hy, 0);

  G4Box* solidAlFace_RIGHT =
    new G4Box("AlFace_RIGHT",           //its name
              AlFace_hx,              // its size 
              AlFace_hy,                 
              AlFace_hz);                 

  G4LogicalVolume* logicAlFace_RIGHT =
    new G4LogicalVolume(solidAlFace_RIGHT,   //its solid
                        AlFace_mat,    //its material
                        "AlFace_RIGHT"); //its name

    new G4PVPlacement(0,               // no rotation
                      AlFace_Position, 
                      logicAlFace_RIGHT, //its logical volume
                      "AlFace_RIGHT",    //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking

  // ++++++++++++++++++++++++++++
  // LEFT FACE 
  AlFace_Position = 
    G4ThreeVector(-(AlFace_hx-foam_hx), -(foam_hy+AlFace_hy), 0);

  G4Box* solidAlFace_LEFT =
    new G4Box("AlFace_LEFT",           //its name
              AlFace_hx,              // its size 
              AlFace_hy,                 
              AlFace_hz);                 

  G4LogicalVolume* logicAlFace_LEFT =
    new G4LogicalVolume(solidAlFace_LEFT,   //its solid
                        AlFace_mat,    //its material
                        "AlFace_LEFT"); //its name

    new G4PVPlacement(0,               // no rotation
                      AlFace_Position, 
                      logicAlFace_LEFT, //its logical volume
                      "AlFace_LEFT",    //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking

  // ++++++++++++++++++++++++++++
  // REAR FACE 
  AlFace_Position = 
    G4ThreeVector(-(2*AlFace_hx-foam_hx)-AlFace_hy, 0, 0);

  G4Box* solidAlFace_REAR =
    new G4Box("AlFace_REAR",           //its name
              AlFace_hy,              // its size 
              foam_hy,                 
              foam_hz);                 

  G4LogicalVolume* logicAlFace_REAR =
    new G4LogicalVolume(solidAlFace_REAR,   //its solid
                        AlFace_mat,    //its material
                        "AlFace_REAR"); //its name

    new G4PVPlacement(0,               // no rotation
                      AlFace_Position, 
                      logicAlFace_REAR, //its logical volume
                      "AlFace_REAR",    //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking
  
  
  

  //#######################################################################
  //FRAME OUTSIDE FOAM

  // ++++++++++++++++++++++++++++
  // Al Tubing REAR TOP
  AlTube_hy = foam_hy;
  AlTubing_Position = 
    G4ThreeVector(-(2*AlFace_hx-foam_hx)+AlTube_hx, 0, foam_hz-AlTube_hz);

  G4Box* solidAlTubing_REAR_TOP =  
    new G4Box("AlTubing_REAR_TOP",                //its name
               AlTube_hx,              // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_REAR_TOP =
    new G4LogicalVolume(solidAlTubing_REAR_TOP,   //its solid
                        AlTubing_mat,    //its material
                        "AlTubing_REAR_TOP");     //its name

    new G4PVPlacement(0,               //no rotation
                      AlTubing_Position,//it's position
                      logicAlTubing_REAR_TOP,    //its logical volume
                      "AlTubing_REAR_TOP",       //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking

  // ++++++++++++++++++++++++++++
  // Al Tubing REAR BOTTOM
  AlTube_hy = foam_hy;
  AlTubing_Position = G4ThreeVector(-(2*AlFace_hx-foam_hx)+AlTube_hx, 0, -(foam_hz-AlTube_hz));

  G4Box* solidAlTubing_REAR_BOTTOM =  
    new G4Box("AlTubing_REAR_BOTTOM",  //its name
               AlTube_hx,              // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_REAR_BOTTOM =
    new G4LogicalVolume(solidAlTubing_REAR_BOTTOM,   //its solid
                        AlTubing_mat,    //its material
                        "AlTubing_REAR_BOTTOM");     //its name

    new G4PVPlacement(0,               //no rotation
                      AlTubing_Position,//it's position
                      logicAlTubing_REAR_BOTTOM,    //its logical volume
                      "AlTubing_REAR_BOTTOM",       //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking


  // ++++++++++++++++++++++++++++
  // Al Tubing REAR RIGHT
  AlTube_hy = foam_hz-2.0*AlTube_hz;
  AlTubing_Position 
    = G4ThreeVector(-(2*AlFace_hx-foam_hx)+AlTube_hx, foam_hy-AlTube_hx,0);
  AlTubing_Transform 
    = G4Transform3D(xRot90deg, AlTubing_Position);

  G4Box* solidAlTubing_REAR_RIGHT =  
    new G4Box("AlTubing_REAR_RIGHT",                //its name
               AlTube_hx,              // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_REAR_RIGHT =
    new G4LogicalVolume(solidAlTubing_REAR_RIGHT,   //its solid
                        AlTubing_mat,    //its material
                        "AlTubing_REAR_RIGHT");     //its name

    new G4PVPlacement(AlTubing_Transform,//it's position
                      logicAlTubing_REAR_RIGHT,    //its logical volume
                      "AlTubing_REAR_RIGHT",       //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking


  // ++++++++++++++++++++++++++++
  // Al Tubing REAR LEFT
  AlTube_hy = foam_hz-2.0*AlTube_hz;
  AlTubing_Position 
    = G4ThreeVector(-(2*AlFace_hx-foam_hx)+AlTube_hx,-(foam_hy-AlTube_hx),0);
  AlTubing_Transform = G4Transform3D(xRot90deg, AlTubing_Position);

  G4Box* solidAlTubing_REAR_LEFT =  
    new G4Box("AlTubing_REAR_LEFT",                //its name
               AlTube_hx,              // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_REAR_LEFT =
    new G4LogicalVolume(solidAlTubing_REAR_LEFT,   //its solid
                        AlTubing_mat,    //its material
                        "AlTubing_REAR_LEFT");     //its name

    new G4PVPlacement(AlTubing_Transform,//it's position
                      logicAlTubing_REAR_LEFT,    //its logical volume
                      "AlTubing_REAR_LEFT",       //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking


  // ++++++++++++++++++++++++++++
  // Al Tubing TOP RIGHT
  AlTube_hy = AlFace_hx-AlTube_hz-foam_hx;
  AlTubing_Position 
    = G4ThreeVector(-foam_hx-AlTube_hy, foam_hy-AlTube_hx,foam_hz-AlTube_hx);
  AlTubing_Transform 
    = G4Transform3D(zRot90deg, AlTubing_Position);

  G4Box* solidAlTubing_TOP_RIGHT =  
    new G4Box("AlTubing_TOP_RIGHT",                //its name
               AlTube_hx,              // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_TOP_RIGHT =
    new G4LogicalVolume(solidAlTubing_TOP_RIGHT,   //its solid
                        AlTubing_mat,    //its material
                        "AlTubing_TOP_RIGHT");     //its name

    new G4PVPlacement(AlTubing_Transform,//it's position
                      logicAlTubing_TOP_RIGHT,    //its logical volume
                      "AlTubing_TOP_RIGHT",       //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking

  // ++++++++++++++++++++++++++++
  // Al Tubing BOTTOM RIGHT
  AlTube_hy = AlFace_hx-AlTube_hz-foam_hx;
  AlTubing_Position 
    = G4ThreeVector(-foam_hx-AlTube_hy, foam_hy-AlTube_hx,-(foam_hz-AlTube_hx));
  AlTubing_Transform 
    = G4Transform3D(zRot90deg, AlTubing_Position);

  G4Box* solidAlTubing_BOTTOM_RIGHT =  
    new G4Box("AlTubing_BOTTOM_RIGHT",                //its name
               AlTube_hx,              // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_BOTTOM_RIGHT =
    new G4LogicalVolume(solidAlTubing_BOTTOM_RIGHT,   //its solid
                        AlTubing_mat,    //its material
                        "AlTubing_BOTTOM_RIGHT");     //its name

    new G4PVPlacement(AlTubing_Transform,//it's position
                      logicAlTubing_BOTTOM_RIGHT,    //its logical volume
                      "AlTubing_BOTTOM_RIGHT",       //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking

  // ++++++++++++++++++++++++++++
  // Al Tubing TOP LEFT 
  AlTube_hy = AlFace_hx-AlTube_hz-foam_hx;
  AlTubing_Position 
    = G4ThreeVector(-foam_hx-AlTube_hy, -(foam_hy-AlTube_hx),foam_hz-AlTube_hx);
  AlTubing_Transform 
    = G4Transform3D(zRot90deg, AlTubing_Position);

  G4Box* solidAlTubing_TOP_LEFT =  
    new G4Box("AlTubing_TOP_LEFT",                //its name
               AlTube_hx,              // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_TOP_LEFT =
    new G4LogicalVolume(solidAlTubing_TOP_LEFT,   //its solid
                        AlTubing_mat,    //its material
                        "AlTubing_TOP_LEFT");     //its name

    new G4PVPlacement(AlTubing_Transform,//it's position
                      logicAlTubing_TOP_LEFT,    //its logical volume
                      "AlTubing_TOP_LEFT",       //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking

  // ++++++++++++++++++++++++++++
  // Al Tubing BOTTOM LEFT 
  AlTube_hy = AlFace_hx-AlTube_hz-foam_hx;
  AlTubing_Position 
    = G4ThreeVector(-foam_hx-AlTube_hy, -(foam_hy-AlTube_hx),-(foam_hz-AlTube_hx));
  AlTubing_Transform 
    = G4Transform3D(zRot90deg, AlTubing_Position);

  G4Box* solidAlTubing_BOTTOM_LEFT =  
    new G4Box("AlTubing_BOTTOM_LEFT",                //its name
               AlTube_hx,              // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_BOTTOM_LEFT =
    new G4LogicalVolume(solidAlTubing_BOTTOM_LEFT,   //its solid
                        AlTubing_mat,    //its material
                        "AlTubing_BOTTOM_LEFT");     //its name

    new G4PVPlacement(AlTubing_Transform,//it's position
                      logicAlTubing_BOTTOM_LEFT,    //its logical volume
                      "AlTubing_BOTTOM_LEFT",       //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking


  // ++++++++++++++++++++++++++++
  // Al Tubing CENTER RIGHT
  AlTube_hy = foam_hz-2.0*AlTube_hz;
  AlTubing_Position 
    = G4ThreeVector(-14*inch+foam_hx+AlTube_hx, foam_hy-AlTube_hx,0);
  AlTubing_Transform 
    = G4Transform3D(xRot90deg, AlTubing_Position);

  G4Box* solidAlTubing_CENTER_RIGHT =  
    new G4Box("AlTubing_CENTER_RIGHT",                //its name
               AlTube_hx,              // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_CENTER_RIGHT =
    new G4LogicalVolume(solidAlTubing_CENTER_RIGHT,   //its solid
                        AlTubing_mat,    //its material
                        "AlTubing_CENTER_RIGHT");     //its name

    new G4PVPlacement(AlTubing_Transform,//it's position
                      logicAlTubing_CENTER_RIGHT,    //its logical volume
                      "AlTubing_CENTER_RIGHT",       //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking


  // ++++++++++++++++++++++++++++
  // Al Tubing CENTER LEFT
  AlTube_hy = foam_hz-2.0*AlTube_hz;
  AlTubing_Position 
    = G4ThreeVector(-14*inch+foam_hx+AlTube_hx,-(foam_hy-AlTube_hx),0);
  AlTubing_Transform = G4Transform3D(xRot90deg, AlTubing_Position);

  G4Box* solidAlTubing_CENTER_LEFT =  
    new G4Box("AlTubing_CENTER_LEFT",                //its name
               AlTube_hx,              // its size 
               AlTube_hy,                 
               AlTube_hz);                 

  G4LogicalVolume* logicAlTubing_CENTER_LEFT =
    new G4LogicalVolume(solidAlTubing_CENTER_LEFT,   //its solid
                        AlTubing_mat,    //its material
                        "AlTubing_CENTER_LEFT");     //its name

    new G4PVPlacement(AlTubing_Transform,//it's position
                      logicAlTubing_CENTER_LEFT,    //its logical volume
                      "AlTubing_CENTER_LEFT",       //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking

  // ++++++++++++++++++++++++++++
  // Al Electronics Tray
  G4double AlTray_hx = 11.5/2.0*inch;
  G4double AlTray_hy = foam_hy;
  G4double AlTray_hz = 1.0/8.0*inch;
  G4ThreeVector AlTray_Position 
    = G4ThreeVector(-15.0*inch-AlTray_hx+foam_hx,0,2.0*inch);

  G4Material *AlTray_mat 
    = nist->FindOrBuildMaterial("G4_Al");

  G4Box* solidAlTray = 
    new G4Box("AlTray",
               AlTray_hx,              // its size 
               AlTray_hy,                 
               AlTray_hz);                 

  G4LogicalVolume* logicAlTray = 
    new G4LogicalVolume(solidAlTray,   //its solid
                        AlTray_mat,    //its material
                        "AlTray");     //its name

    new G4PVPlacement(0,               //no rotation
                      AlTray_Position, //it's position
                      logicAlTray,     //its logical volume
                      "AlTray",        //its name
                      logicEnvelope,   //its mother volume
                      false,           //no boolean operation
                      0,               //copy number
                      checkOverlaps);  //overlaps checking
  
}
