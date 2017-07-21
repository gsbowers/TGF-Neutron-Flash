#include "myDetectorConstruction.hh"
//#include "myPolyBall.hh"
#include "myActionInitialization.hh"
#include "myPhysicsList.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc, char** argv){

	string* inputPrefix;

	if (argc == 2) 
	{	
		string inputFilename = argv[1];
		// get inputFilename Prefix for name mangling 
		int stop = inputFilename.find_last_of('_');
		inputPrefix = new string(inputFilename.substr(0,stop));
	} else {
		inputPrefix = new string("test"); 
	}

	G4RunManager *runManager = new G4RunManager();

	// set mandatory initialization classes
	runManager->SetUserInitialization(new myDetectorConstruction(*inputPrefix));
  runManager->SetUserInitialization(new myPhysicsList);	
	runManager->SetUserInitialization(new myActionInitialization);	

	// initialize G4 kernal
	//runManager->Initialize();

   // get the pointer to the User Interface manager
   G4UImanager* UImanager = G4UImanager::GetUIpointer();

   if (argc==2)   // batch mode
   {
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UImanager->ApplyCommand(command+fileName);
   }
   else
   {  // interactive mode : define UI session
#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

#ifdef G4UI_USE
     G4UIExecutive* ui = new G4UIExecutive(argc, argv);
     G4String command = "/control/execute vis.mac";
     UImanager->ApplyCommand(command);
     ui->SessionStart();
     delete ui;
#endif

#ifdef G4VIS_USE
	delete visManager;
#endif
   }

	// job termination
	delete runManager;
	return 0;

}
