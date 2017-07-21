#ifndef myDetectorHit_h
#define myDetectorHit_h 1

#include "G4Allocator.hh"
#include "G4ParticleDefinition.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"

using namespace std;

class myDetectorHit : public G4VHit
{
  public:
    myDetectorHit();
    virtual ~myDetectorHit();
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);

    void Print();

    // set methods
    void SetParticleDefinition (G4ParticleDefinition *val) { fPD=val; };
    void SetEdep(G4double val) { fEdep = val; };
    void SetKE(G4double val) { fKE = val; };
		void SetPos(G4ThreeVector val) { fPos = val; };
		void SetDir(G4ThreeVector val) { fDir = val; };
		void SetGlobalTime(G4double val) { fGlobalTime = val; };
		void SetProcessName(G4String &val) { fProcessName = val; };
		void SetCreatorProcessName(G4String &val) { fCreatorProcessName = val; };
		void SetEventID(G4int val) { fEventID = val; };
		void SetParentID(G4int val) { fParentID = val; };
		void SetTrackID(G4int val) { fTrackID = val; };

    // get methods
    G4ParticleDefinition* GetParticleDefinition() const { return fPD; };
    G4double GetEdep() const { return fEdep; };
    G4double GetKE() const { return fKE; };
    G4ThreeVector GetPos() const { return fPos; };
    G4ThreeVector GetDir() const { return fDir; };
    G4double GetGlobalTime() const { return fGlobalTime; };
		G4String GetProcessName() const { return fProcessName;};
		G4String GetCreatorProcessName() const { return fCreatorProcessName;};
		G4int GetEventID() const { return fEventID;};
		G4int GetParentID() const { return fParentID;};
		G4int GetTrackID() const { return fTrackID;};

  private:

    // hit variables
    G4ParticleDefinition* fPD;
    G4double fEdep;
    G4double fKE;
		G4ThreeVector fPos;
		G4ThreeVector fDir;
		G4double fGlobalTime;
		G4String fProcessName;
		G4String fCreatorProcessName;
		G4int fEventID;
		G4int fParentID;
		G4int fTrackID;
};

typedef G4THitsCollection<myDetectorHit> myDetectorHitsCollection;

extern G4ThreadLocal G4Allocator<myDetectorHit>* myDetectorHitAllocator;

inline void* myDetectorHit::operator new(size_t)
{
  if (!myDetectorHitAllocator)
    myDetectorHitAllocator = new G4Allocator<myDetectorHit>;
  return (void*) myDetectorHitAllocator->MallocSingle();
}

inline void myDetectorHit::operator delete(void *aHit)
{
  myDetectorHitAllocator->FreeSingle((myDetectorHit*) aHit);
}

#endif
