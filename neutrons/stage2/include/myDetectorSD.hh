#ifndef myDetectorSD_h
#define myDetectorSD_h 1

#include "G4VSensitiveDetector.hh"
#include "myDetectorHit.hh"

#include <iostream>
#include <fstream>

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

// my Sensitive Detector

class myDetectorSD : public G4VSensitiveDetector
{
  public:
    myDetectorSD(const G4String& name, 
								 const G4String& hitsCollectionName);
    virtual ~myDetectorSD();
  
    virtual void Initialize(G4HCofThisEvent *HCE);
    virtual G4bool ProcessHits(G4Step *aStep, G4TouchableHistory* ROhist);
    virtual void EndOfEvent(G4HCofThisEvent* hitCollection);
    static myDetectorSD* Instance();

    G4int GetCurrentEventID() const { return fCurrentEventID; }
    void SetCurrentEventID(G4int val) {fCurrentEventID = val; }

    void initHitOutfile(char* outfile);
    void outputHitData();

  private:
    static myDetectorSD* fgInstance;
    G4int fCurrentEventID;

    myDetectorHitsCollection* fHitsCollection;    
    G4int fHCID;

		// tracking flags
		G4int fTrackEvents;
		G4int fTrackHits;
		G4bool fNeutronFilter;

    // output files
    ofstream fout_events;  
    ofstream fout_hits; 

    // intermediate array storage for event info
    G4int fcount_events;
    G4double* fEdepEvt;
    G4int* fevtNb;

    // intermediate array storage for hit info
    G4int fcount_hits;
		G4ThreeVector* fPos;
		G4ThreeVector* fDir;
		G4double* fGlobalTime;
    G4double* fEdepHit;
		G4double* fKEHit;
		G4int* fhitEventID;
    G4int* fhitNb;
		G4int* fParentID; 
		G4int* fTrackID; 
		G4String* fParticle;
		G4String* fProcess;
		G4String* fCreatorProcess;
};

#endif
