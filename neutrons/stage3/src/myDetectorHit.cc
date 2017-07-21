#include "myDetectorHit.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

G4ThreadLocal G4Allocator<myDetectorHit>* myDetectorHitAllocator;

myDetectorHit::myDetectorHit()
: G4VHit()
{
}

myDetectorHit::~myDetectorHit()
{
}

void myDetectorHit::Print()
{
  G4cout << fEdep/MeV << " (MeV) " << G4endl; 
}
