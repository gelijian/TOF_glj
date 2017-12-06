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
// $Id: TOFscinSD.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file TOFscinSD.cc
/// \brief Implementation of the TOFscinSD class

#include <stdio.h>
#include "TOFscinSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TOFscinSD::TOFscinSD(const G4String& name,
                         const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TOFscinSD::~TOFscinSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TOFscinSD::Initialize(G4HCofThisEvent* hce)
{
	// Create hits collection

	fHitsCollection 
		= new TOFscinHitsCollection(SensitiveDetectorName, collectionName[0]); 

	// Add this collection in hce

	G4int hcID 
		= G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
	hce->AddHitsCollection( hcID, fHitsCollection ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TOFscinSD::ProcessHits(G4Step* aStep, 
                                     G4TouchableHistory*)
{  
	// energy deposit
	//
	G4int ParentID = aStep->GetTrack()->GetParentID();
	G4int StepID = aStep->GetTrack()->GetCurrentStepNumber();
	G4String Pname = aStep->GetTrack()->GetDefinition()->GetParticleName();
    G4String PVname = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
	
	if (ParentID != 1 || StepID != 1 || Pname == "neutron" || Pname == "gamma" || Pname == "opticalphoton")
		return false;
        
    G4int DetID = -100;
    G4int avnum = -100;
    char temp[100];
    sscanf(PVname.c_str(), "av_%d_impr_%d_%s", &avnum, &DetID, temp);
    DetID = DetID - 1;
	G4int TrackID = aStep->GetTrack()->GetTrackID();
	G4double Edep  = aStep->GetPreStepPoint()->GetKineticEnergy() / keV;
	G4double Time = aStep->GetPreStepPoint()->GetGlobalTime() / ns;
    
	TOFscinHit* newHit = new TOFscinHit();
	newHit->SetPname(Pname);
	newHit->SetPVname(PVname);
	newHit->SetDetID(DetID);
	newHit->SetTrackID(TrackID);
	newHit->SetEdep(Edep);
	newHit->SetTime(Time);
	fHitsCollection->insert( newHit );

	return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TOFscinSD::EndOfEvent(G4HCofThisEvent*)
{
	/*if ( verboseLevel>=0 ) 
	{ 
		G4int nofHits = fHitsCollection->entries();
		if (nofHits == 0)
		{
			return;
		}
		
		G4cout << G4endl
			   << "-------->" << SensitiveDetectorName << ": in this event they are " << nofHits 
			   << G4endl;
            
    for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->DisplayHit();
	}*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
