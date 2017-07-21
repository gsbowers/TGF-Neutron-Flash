#include "myDetectorConstruction.hh"
#include "myDetectorSD.hh"

#include "G4UserLimits.hh"
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
#include "G4Material.hh"

#include "G4ios.hh"

myDetectorConstruction::myDetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0), fLogicSD(0)
{ 

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

	fLogicSlabSDs = static_cast<G4LogicalVolume*>(::operator new(sizeof(G4LogicalVolume)*fnSlabs));
	fStepLimits = static_cast<G4UserLimits*>(::operator new(sizeof(G4UserLimits)*fnSlabs));

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
	// destruct fLogicSlabSDs array
	for(int i = fnSlabs-1; i >=0; --i){
		fLogicSlabSDs[i].~G4LogicalVolume();
	} 
	::operator delete[]( fLogicSlabSDs);

	// destruct fStepLimits array
	for(int i = fnSlabs-1; i >=0; --i){
		fStepLimits[i].~G4UserLimits();
	} 
	::operator delete[]( fStepLimits);

	delete[] fslabs_z;
	delete[] fslabs_pDz;
	delete[] fslabs_rho;
	delete[] fslabs_T;
	delete[] fslabs_P;
}

G4VPhysicalVolume* myDetectorConstruction::Construct()
{
  // Get nist material manager
  //
  G4NistManager* nist = G4NistManager::Instance();

  //
  // World
  //

	G4double XYMax = 6*km;  
	G4double ZMax  = 4*km;  
 
  // TUBS
	G4double world_sizeXY = 1.1*XYMax;
	G4double world_sizeZ = 1.1*ZMax;
  
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");

  // TUBS
  G4Box* solidWorld = 
		new G4Box("World", 
      				0.5*world_sizeXY, 
							0.5*world_sizeXY, 
							0.5*world_sizeZ); // its size

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

	// Thermal Hydrogen in Water
	// N.B. Contribution of O from H2O is negligible compared to O2 in air
  fTS_H = new G4Element("TS_H_of_Water" ,"H" , 1., 1.0079*g/mole);
	G4cout << fTS_H << G4endl;

  //
  // Atmosphere
  //
	for(G4int i=0; i<fnSlabs; i++){
		ConstructSlab(i,&logicWorld);
	} 

  //
  // Ground
  //

	// TUBS
	G4double ground_sizeXY = 5*km;
	G4double ground_sizeZ = 100*m;
	
	G4double ground_z = -0.5*(ZMax + ground_sizeZ);

	// Material
  G4Material* ground_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
	G4cout << ground_mat << G4endl;

  // position and orientation
  G4RotationMatrix yRot90deg = G4RotationMatrix();
  yRot90deg.rotateY(0*deg);

  G4ThreeVector ground_position = G4ThreeVector(0,0,ground_z);
  G4Transform3D ground_transform = G4Transform3D(yRot90deg, ground_position);

  G4Box* solidGround = 
		new G4Box("Ground",  //its name
    					0.5*ground_sizeXY, 
							0.5*ground_sizeXY, 
							0.5*ground_sizeZ);  //its size 
		
  G4LogicalVolume* logicGround = 
    new G4LogicalVolume(solidGround,    //its solid
                        ground_mat,     //its material
                        "Ground");      //its name

  new G4PVPlacement(ground_transform, 
            				logicGround,  			//its logical volume
            				"Ground",     			//its name
            				logicWorld,  			  //its mother (logical) volume
            				false,        			//no boolean operations
            				0,            			//its copy number 
            				true);        			//checkOverlap

  // sets a max step length in the tracker region
  // 	
   
  //G4double maxStep = large_pDz/16;
  //logicLgPl->SetUserLimits(fStepLimit);

  //
  // Visualization
  //
  logicWorld->SetVisAttributes(new G4VisAttributes(G4Colour(0.,1.,1.)));  
  logicGround->SetVisAttributes(new G4VisAttributes(G4Colour(1.,0.,0.)));  

	return physWorld;

}

void myDetectorConstruction::ConstructSlab(G4int slabNo, 
	G4LogicalVolume** logicWorld)
{
	// make slab name
	char slabName[50]; 
	sprintf(slabName, "Slab%02i", slabNo);

	G4double slab_sizeXY = 5*km;
	G4double slab_sizeZ = 2.0*fslabs_pDz[slabNo];
	G4double slab_z = fslabs_z[slabNo]; 

  // detector positions and orientations
  G4RotationMatrix yRot90deg = G4RotationMatrix();
  yRot90deg.rotateY(0*deg);

  G4ThreeVector slab_position = G4ThreeVector(0,0,slab_z);
  G4Transform3D slab_transform = G4Transform3D(yRot90deg, slab_position);

	// Material Definitions 
	
	// Slab Air
	char materialName[50]; 
	sprintf(materialName, "Air%02i", slabNo);

	G4double density;
  G4double temperature;
  G4double pressure;
  G4NistManager* nist = G4NistManager::Instance();
	nist->SetVerbose(1);

	density = fslabs_rho[slabNo];
	temperature = fslabs_T[slabNo];
	pressure = fslabs_P[slabNo];

	G4Material *Air = nist->ConstructNewGasMaterial(materialName, "G4_AIR", temperature, pressure);

	// Slab Material - Mixture of Air and Rain 
	sprintf(materialName, "SlabMat%02i", slabNo);
	G4double rho_TS_H = 2.778e-4 * 2.0/18.0; // mg/cm3
	G4double rho_air = density/(mg/cm3); 
	G4double fraction_Air = rho_air/(rho_air + rho_TS_H);
	G4double fraction_TS_H = rho_TS_H/(rho_air + rho_TS_H);
	G4Material* SlabMat = new G4Material(materialName, density, 2, kStateGas, temperature, pressure);  
  SlabMat->AddMaterial(Air, fraction_Air);
  SlabMat->AddElement(fTS_H, fraction_TS_H);

	SlabMat = Air;

	G4cout << SlabMat << G4endl; // output material properties
	G4cout << "\t\t\t\t\tElmMassFraction of TS_H_of_Water: " << fraction_TS_H / perCent << "%" << G4endl << G4endl;
		
  G4Box* solidSlab = new G4Box(slabName,             //its name
    0.5*slab_sizeXY, 0.5*slab_sizeXY, 0.5*slab_sizeZ);

	new (fLogicSlabSDs+slabNo) G4LogicalVolume(solidSlab, SlabMat, slabName);
	
  new G4PVPlacement(slab_transform, 
                    (fLogicSlabSDs+slabNo), //its logical volume
                    slabName,               //its name
                    *logicWorld,            //its mother (logical) volume
                    false,                  //no boolean operations
                    0,                      //its copy number 
                    true);                  //checkOverlap

	// user limits 
	new (fStepLimits+slabNo) G4UserLimits(); 
	// photonuclear production
	G4double minEkine = 8*MeV;
	(fStepLimits+slabNo)->SetUserMinEkine(minEkine);
	// step limit 
	G4double stepMax = slab_sizeZ/4.0;
	(fStepLimits+slabNo)->SetMaxAllowedStep(stepMax);
	// set userlimits for slab
	(fLogicSlabSDs+slabNo)->SetUserLimits(fStepLimits+slabNo);

	// set color of region
  (fLogicSlabSDs+slabNo)->SetVisAttributes(new G4VisAttributes(G4Colour(0.0,0.0,1.-0.7*slabNo/fnSlabs)));  


}

void myDetectorConstruction::ConstructSDandField()
{

	//	Sensitive Detectors
	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	G4VSensitiveDetector *slabSD;

	char slabName[50];
	char slabCollection[50]; 

	// iterate over sensitive detector slabs and 
	// register to G4SDManager
	for (G4int i=0; i<fnSlabs; i++){

		sprintf(slabName, "Slab%02i", i);
		sprintf(slabCollection, "Slab%02iHitsCollection", i);

		slabSD = new  myDetectorSD(slabName, slabCollection); 
		SDman->AddNewDetector(slabSD);
		(fLogicSlabSDs+i)->SetSensitiveDetector(slabSD);
	}

}
