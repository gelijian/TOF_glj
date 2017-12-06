#ifndef LSPhysicsList_h
#define TOFPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "G4VModularPhysicsList.hh"
#include "globals.hh"

#include <vector>
using namespace std;

class G4VPhysicsConstructor;

class TOFPhysicsList: public G4VUserPhysicsList
{
	public:
		TOFPhysicsList();
		~TOFPhysicsList();
    
		void SetCuts();

		// Construct particle and physics
		void ConstructParticle();
		void ConstructProcess();
	
	private:
		// these methods Construct physics processes and register them
		void ConstructHadrPhys();
	
		G4VPhysicsConstructor*  emPhysList;
		G4VPhysicsConstructor*  fParticleList;
		vector<G4VPhysicsConstructor*>  hadronPhys;

};

#endif
