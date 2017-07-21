#include <string>
#include <fstream>
#include "myDetectorConstruction.hh"
#include "myActionInitialization.hh"
#include "myPhysicsList.hh"

using namespace std;

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

	string particleInputFilename = argv[1];
	// get particleInputFilename Prefix for name mangling 
	//int start = 0;
	int start = particleInputFilename.find_last_of('/')+1;
	int stop = particleInputFilename.find_last_of('_');
	string particleInputPrefix = particleInputFilename.substr(start,stop-start);

	// construct the default run manager
	G4RunManager* runManager = new G4RunManager;

	// set mandatory initialization classes
	runManager->SetUserInitialization(new myDetectorConstruction(particleInputPrefix));
	myPhysicsList* phys = new myPhysicsList;
  runManager->SetUserInitialization(phys);	
	runManager->SetUserInitialization(new myActionInitialization(particleInputFilename));	

	// initialize G4 kernal
	//runManager->Initialize();

   // get the pointer to the User Interface manager
   G4UImanager* UImanager = G4UImanager::GetUIpointer();

   if (argc==2)   // batch mode
   {
      G4String fileName = argv[1];
      //G4String command = "/control/execute ";
      //UImanager->ApplyCommand(command+fileName);
		
		  ifstream f;
			unsigned long Ninput;
		  f.open(fileName.c_str(), ios::in);
		  string line;
		 	for (Ninput=0l; getline(f, line); ++Ninput); 
			f.close();

			char command[256];
			sprintf(command, "/run/beamOn %lu", Ninput-1l);
			
      UImanager->ApplyCommand("/run/verbose 2");
      UImanager->ApplyCommand("/testhadr/phys/thermalScattering true");
      UImanager->ApplyCommand("/run/initialize");
      UImanager->ApplyCommand(command);
      UImanager->ApplyCommand("/control/shell ls -l *.out");
      //UImanager->ApplyCommand("/run/beamOn 100");
			
				
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
		G4cout << "foobar" << G4endl;
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
