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


#ifndef TOFscinHit_h
#define TOFscinHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

/// Tracker hit class
///
/// It defines data members to store the trackID, chamberNb, energy deposit,
/// and position of charged particles in a selected volume:
/// - TrackID, fChamberNB, fEdep, fPos

class TOFscinHit : public G4VHit
{
	public:
		TOFscinHit();
		TOFscinHit(const TOFscinHit&);
		virtual ~TOFscinHit();

		// operators
		const TOFscinHit& operator=(const TOFscinHit&);
		G4int operator==(const TOFscinHit&) const;

		inline void* operator new(size_t);
		inline void  operator delete(void*);

		// Set methods
		void SetPname 	(G4String pname)	{ Pname = pname; };
		void SetPVname  (G4String pvname)	{ PVname = pvname; };
		void SetDetID  	(G4int detID)      	{ DetID = detID; };
		void SetTrackID (G4int trackID)      	{ TrackID = trackID; };
		void SetEdep    (G4double edep)      	{ Edep = edep; };
		void SetTime	(G4double time)		{ Time = time; };
		
		// Get methods
		G4String	GetPname() const 	{return Pname;};
		G4String	GetPVname() const 	{return PVname;};
		G4int		GetDetID() const 	{return DetID;};
		G4int       GetTrackID() const 	{return TrackID;};
		G4double    GetEdep() const 	{return Edep;};
		G4double	GetTime() const 	{return Time;};
		
		void DisplayHit();

	private:
		
		G4String	  Pname;
		G4String	  PVname;
		G4int		  DetID;
		G4int         TrackID;
		G4double      Edep;
		G4double	  Time;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<TOFscinHit> TOFscinHitsCollection;

extern G4ThreadLocal G4Allocator<TOFscinHit>* TOFscinHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* TOFscinHit::operator new(size_t)
{
	if(!TOFscinHitAllocator)
		TOFscinHitAllocator = new G4Allocator<TOFscinHit>;
	return (void *) TOFscinHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void TOFscinHit::operator delete(void *hit)
{
	TOFscinHitAllocator->FreeSingle((TOFscinHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
