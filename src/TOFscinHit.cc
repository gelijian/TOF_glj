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

#include "TOFscinHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<TOFscinHit>* TOFscinHitAllocator=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TOFscinHit::TOFscinHit()
 : G4VHit(),
   Pname("undifined"),
   PVname("undifined"),
   DetID(-1),
   TrackID(-1),
   Edep(0.),
   Time(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TOFscinHit::~TOFscinHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TOFscinHit::TOFscinHit(const TOFscinHit& right)
  : G4VHit()
{
	Pname	= 	right.Pname;
	PVname 	= 	right.PVname;
	DetID	= 	right.DetID;
	TrackID	= 	right.TrackID;
	Edep	=	right.Edep;
	Time	=	right.Time;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const TOFscinHit& TOFscinHit::operator=(const TOFscinHit& right)
{
	Pname	= 	right.Pname;
	PVname 	= 	right.PVname;
	DetID	= 	right.DetID;
	TrackID	= 	right.TrackID;
	Edep	=	right.Edep;
	Time	=	right.Time;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int TOFscinHit::operator==(const TOFscinHit& right) const
{
	return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TOFscinHit::DisplayHit()
{
	G4cout << " Pname: " << Pname
		   << " PVname: " << PVname
		   << " DetID: " << DetID
		   << " TrackID: " << TrackID
		   << " Edep: " << Edep
		   << " Time: " << Time
		   << G4endl;
}	

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

