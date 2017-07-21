//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr04/src/NeutronHPphysics.cc
/// \brief Implementation of the NeutronHPphysics class
//
// $Id: NeutronHPphysics.cc 66587 2012-12-21 11:06:44Z ihrivnac $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "myNeutronHPphysics.hh"

#include "myNeutronHPMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessTable.hh"

// Processes
#include "G4ChipsElasticModel.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4NeutronElasticXS.hh"

#include "G4CascadeInterface.hh"
#include "G4NeutronRadCapture.hh"
#include "G4LFission.hh"

#include "G4HadronElasticProcess.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPThermalScatteringData.hh"
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPThermalScattering.hh"

#include "G4NeutronInelasticProcess.hh"
#include "G4ParticleHPInelasticData.hh"
#include "G4ParticleHPInelastic.hh"

#include "G4HadronCaptureProcess.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4ParticleHPCapture.hh"

#include "G4HadronFissionProcess.hh"
#include "G4ParticleHPFissionData.hh"
#include "G4ParticleHPFission.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

myNeutronHPphysics::myNeutronHPphysics(const G4String& name)
:  G4VPhysicsConstructor(name), fThermal(true), fNeutronMessenger(0)
{
  fNeutronMessenger = new myNeutronHPMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

myNeutronHPphysics::~myNeutronHPphysics()
{
  delete fNeutronMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void myNeutronHPphysics::ConstructProcess()
{
  G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4ProcessManager* pManager = neutron->GetProcessManager();
   
  // delete all neutron processes if already registered
  //
  G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();
  G4VProcess* process = 0;
  process = processTable->FindProcess("hadElastic", neutron);      
  if (process) pManager->RemoveProcess(process);
  //
  process = processTable->FindProcess("neutronInelastic", neutron);
  if (process) pManager->RemoveProcess(process);
  //
  process = processTable->FindProcess("nCapture", neutron);      
  if (process) pManager->RemoveProcess(process);
  //
  process = processTable->FindProcess("nFission", neutron);      
  if (process) pManager->RemoveProcess(process);      
         
	//
  // ############## Elastic Scattering
  // 	
	G4HadronElasticProcess* the_neutron_elastic_process 
		= new G4HadronElasticProcess();

	// Models
	G4ParticleHPThermalScattering* the_thermal_neutron_elastic_model 
		= new G4ParticleHPThermalScattering(); // < 4 ev
	G4ParticleHPElastic* the_HP_neutron_elastic_model 
		= new G4ParticleHPElastic(); //  4 eV - 20 MeV
	G4ChipsElasticModel* the_high_energy_neutron_elastic_model
		= new G4ChipsElasticModel(); // > 20 MeV
	
	// Energy Ranges
	the_HP_neutron_elastic_model->SetMinEnergy(4*eV);
	the_HP_neutron_elastic_model->SetMaxEnergy(20*MeV);
  the_high_energy_neutron_elastic_model->SetMinEnergy(20*MeV);
	
	// Data Sets
	the_neutron_elastic_process->AddDataSet(new G4ParticleHPThermalScatteringData());
	the_neutron_elastic_process->AddDataSet(new G4ParticleHPElasticData());

	// Register models
	the_neutron_elastic_process->RegisterMe(the_thermal_neutron_elastic_model);
	the_neutron_elastic_process->RegisterMe(the_HP_neutron_elastic_model);
	the_neutron_elastic_process->RegisterMe(the_high_energy_neutron_elastic_model);

	// Register Process
  pManager->AddDiscreteProcess(the_neutron_elastic_process);   

	// 
  // ############## Inelastic Scattering
  //
	G4NeutronInelasticProcess* the_neutron_inelastic_process 
		= new G4NeutronInelasticProcess();
  
	// Models
  G4ParticleHPInelastic* the_HP_neutron_inelastic_model 
		= new G4ParticleHPInelastic();  // 0-20MeV
	G4CascadeInterface* the_high_energy_neutron_inelastic_model
		= new G4CascadeInterface(); // > 20MeV

	// Energy ranges
	the_high_energy_neutron_inelastic_model->SetMinEnergy(20*MeV);
  
  // Data Sets
  the_neutron_inelastic_process->AddDataSet(new G4ParticleHPInelasticData());

	// Register Models
  the_neutron_inelastic_process->RegisterMe(the_HP_neutron_inelastic_model);    
  the_neutron_inelastic_process->RegisterMe(the_high_energy_neutron_inelastic_model);    

	// Register Process
  pManager->AddDiscreteProcess(the_neutron_inelastic_process);   

	//  
  // ############## nCapture   
  //
	G4HadronCaptureProcess* the_nCapture_process
		= new G4HadronCaptureProcess("nCapture");

  // Models
  G4ParticleHPCapture* the_HP_nCapture_model 
		= new G4ParticleHPCapture(); // 0-20MeV
	G4NeutronRadCapture* the_high_energy_nCapture_model
		= new G4NeutronRadCapture(); // > 20MeV

	// Energy ranges
	the_high_energy_nCapture_model->SetMinEnergy(20*MeV);
	
  // Data Sets 
  the_nCapture_process->AddDataSet(new G4ParticleHPCaptureData());

  // Register Models
  the_nCapture_process->RegisterMe(the_HP_nCapture_model);
  the_nCapture_process->RegisterMe(the_high_energy_nCapture_model);

	// Register Process
  pManager->AddDiscreteProcess(the_nCapture_process);    

  // ##############  nFission   
  //
	G4HadronFissionProcess* the_nFission_process
		= new G4HadronFissionProcess("nFission");

  // Models
  G4ParticleHPFission* the_HP_nFission_model 
		= new G4ParticleHPFission(); // 0-20MeV
	G4LFission* the_high_energy_nFission_model
		= new G4LFission(); // > 20MeV

	// Energy ranges
	the_high_energy_nFission_model->SetMinEnergy(20*MeV);
	
  // Data Sets 
  the_nFission_process->AddDataSet(new G4ParticleHPFissionData());

  // Register Models
  the_nFission_process->RegisterMe(the_HP_nFission_model);
  the_nFission_process->RegisterMe(the_high_energy_nFission_model);

	// Register Process
  pManager->AddDiscreteProcess(the_nFission_process);    

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
