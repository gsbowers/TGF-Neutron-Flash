#include "G4UImanager.hh"

#ifndef particle_h
#define particle_h 1

using namespace std;

class myParticle{

	public:
    myParticle(){};
    ~myParticle(){};

    G4double get_Energy(){ 
      return fEnergy; 
    };

    G4double get_ArrivalTime(){
      return fArrivalTime;
    };

    G4double* get_StartXYZ(){
      return fStartXYZ;
    };

    G4double* get_StartDir(){
      return fStartDir;
    };

    void set_Energy(G4double E){
      fEnergy = E;
    };

    void set_ArrivalTime(G4double time){
      fArrivalTime = time;
    };
 
    void set_StartXYZ(G4double x, G4double y, G4double z){
      fStartXYZ[0] = x;
      fStartXYZ[1] = y;
      fStartXYZ[2] = z;
    };
 
    void set_StartDir(G4double a, G4double b, G4double c){
      fStartDir[0] = a;
      fStartDir[1] = b;
      fStartDir[2] = c;
    };

  private:
    G4double fEnergy; 
		G4double fArrivalTime;
    G4double fStartXYZ[3]; 
    G4double fStartDir[3];
};

#endif
