#include "myParticleInput.hh"
#include "myParticle.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

myParticleInput::myParticleInput(string filename)
{
	ifstream f;	

  // get number of lines in file
  f.open(filename.c_str(), ios::in);
  string line;
  for (fNinput=0l; getline(f, line); ++fNinput);
  // reset file pointer
  f.clear();
  f.seekg(0, ios::beg);

  // allocate memory
  fEvent       = new G4int[fNinput](); 
  fEnergy      = new G4double[fNinput](); 
  fArrivalTime = new G4double[fNinput](); 
	fPosX        = new double[fNinput]();
	fPosY        = new double[fNinput]();
	fPosZ        = new double[fNinput]();
	fDirX        = new double[fNinput]();
	fDirY        = new double[fNinput]();
	fDirZ        = new double[fNinput]();

  // read input file 
  long i = 0l;
	
	string sPos;
	string sDir;
	G4double data;
	vector<G4double> vect;
	
	getline(f, line);
	while (getline(f, line)){ 

		stringstream lineStream(line);
		
		// read line data
		lineStream >> fEvent[i]
							 >> fEnergy[i]
							 >> fArrivalTime[i]
							 >> sPos >> sDir;

		// process 'sPos' and 'sDir' strings
		stringstream ssPos(sPos.erase(0, 1));			
		while (ssPos >> data)
		{	
			vect.push_back(data);
			if (ssPos.peek() == ',' || ssPos.peek() == ')') ssPos.ignore();
		}
		fPosX[i] = vect.at(0);
		fPosY[i] = vect.at(1);
		fPosZ[i] = vect.at(2);
		vect.clear();

		stringstream ssDir(sDir.erase(0, 1));			
		while (ssDir >> data)
		{	
			vect.push_back(data);
			if (ssDir.peek() == ',' || ssDir.peek() == ')') ssDir.ignore();
		}
		fDirX[i] = vect.at(0);
		fDirY[i] = vect.at(1);
		fDirZ[i] = vect.at(2);
		vect.clear();

		/*
		G4cout << fEvent[i] << " "
					 << fEnergy[i] << " "
           << fArrivalTime[i] << " "
           << fPosX[i] << " "
           << fPosY[i] << " "
           << fPosZ[i] << " "
           << fDirX[i] << " "
           << fDirY[i] << " "
           << fDirZ[i] << " "
           << G4endl;
		*/

		i++;
	}
  f.close();

  cout << fNinput << " particles input from " << filename << endl;
	findx = 0l;

  return;
}

myParticleInput::~myParticleInput()
{
  delete[] fEvent;
  delete[] fArrivalTime; 
  delete[] fEnergy;
  delete[] fPosX;
  delete[] fPosY;
  delete[] fPosZ;
  delete[] fDirX;
  delete[] fDirY;
  delete[] fDirZ;
}

myParticle myParticleInput::getParticle()
{
  // 1.  Pick a particle from the particle array according 
  // to the specified fECD
  // 2.  Specify it's direction
  // 3.  Specify it's location on the input sphere
 
  myParticle particle;
	long i = findx;

  particle.set_Energy(fEnergy[i]*eV);

//c CHOOSE A DIRECTION: the third direction cosine is given, pick
//c the horizontal ones using a random azimuthal angle phi.

  particle.set_StartDir(fDirX[i], fDirY[i], fDirZ[i]);

//c CHOOSE A STARTING POSITION randomly on a disk (which will be  the
//c input sphere's shadow cast by a plane wave coming in the input
//c direction generated above for the current photon:

  particle.set_StartXYZ(fPosX[i]*m, fPosY[i]*m, fPosZ[i]*m);

  particle.set_ArrivalTime(fArrivalTime[i]*ms);

	findx++;
	return particle;

}
