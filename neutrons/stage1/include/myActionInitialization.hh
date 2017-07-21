#ifndef myActionInitialization_h
#define myActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class myActionInitialization : public G4VUserActionInitialization
{
  public:
    myActionInitialization();
    virtual ~myActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
};

#endif
