//Action Initialization for T1 project
//21.2.2018 Patric Granholm

#include "T1ActionInitialization.hh"
//#include "T1RunAction.hh"
#include "T1PrimaryGeneratorAction.hh"

T1ActionInitialization::T1ActionInitialization()
:G4VUserActionInitialization()
{}

T1ActionInitialization::~T1ActionInitialization()
{}

void T1ActionInitialization::Build() const
{
  SetUserAction(new T1PrimaryGeneratorAction);

  //  T1RunAction* runAction = new T1RunAction;
  //  SetUserAction(runAction);
}