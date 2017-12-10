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
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "G4ParticleGun.hh"

#include "Analysis.hh"
#include "ionRESP.hh"
#include "TOFEventAction.hh"
#include "TOFscinHit.hh"
#include "TOFscinSD.hh"
#include "TOFPrimaryGeneratorAction.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


TOFEventAction::TOFEventAction()
: G4UserEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TOFEventAction::~TOFEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TOFEventAction::BeginOfEventAction(const G4Event* event)
{
    G4int eventID = event->GetEventID();
    if (eventID % 10000 == 0)
    {
        G4cout << "eventID: " << eventID << G4endl;
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TOFEventAction::EndOfEventAction(const G4Event* event)
{   
    G4int eventID = event->GetEventID();
    
    // get analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
    // Fill the ntuple for each event
    G4int s1_HCID = G4SDManager::GetSDMpointer()->GetCollectionID("s1_hitscollection");
	G4int s2u_HCID = G4SDManager::GetSDMpointer()->GetCollectionID("s2u_hitscollection");
	G4int s2b_HCID = G4SDManager::GetSDMpointer()->GetCollectionID("s2b_hitscollection");
	
	TOFscinHitsCollection* s1_HC = static_cast<TOFscinHitsCollection*>(event->GetHCofThisEvent()->GetHC(s1_HCID));
	TOFscinHitsCollection* s2u_HC = static_cast<TOFscinHitsCollection*>(event->GetHCofThisEvent()->GetHC(s2u_HCID));
	TOFscinHitsCollection* s2b_HC = static_cast<TOFscinHitsCollection*>(event->GetHCofThisEvent()->GetHC(s2b_HCID));

	G4int Nohits_s1 = s1_HC->entries();
	G4int Nohits_s2u = s2u_HC->entries();
	G4int Nohits_s2b = s2b_HC->entries();
	
	if (Nohits_s1 == 0)
	{
		return;
	}
	else if (Nohits_s2u == 0 && Nohits_s2b == 0)
	{
		return;
	}
    
	G4double s1_time[5];
	G4double s2u_time[40];
	G4double s2b_time[40];
	G4double s1_phlight[5];
	G4double s2u_phlight[40];
	G4double s2b_phlight[40];
    
    //init
    for (G4int i = 0; i < 5; i++)
	{
		s1_time[i] = -1;
		s1_phlight[i] = 0;
	}
	
	for (G4int i = 0; i < 40; i++)
	{
		s2u_time[i] = -1;
		s2b_time[i] = -1;
		s2u_phlight[i] = 0;
		s2b_phlight[i] = 0;	
	}
	//
	
	//
	for (int i = 0; i < Nohits_s1; i++)
	{
		G4int DetID = (*s1_HC)[i]->GetDetID();
		G4String Pname = (*s1_HC)[i]->GetPname();
		G4double Time = (*s1_HC)[i]->GetTime();
		G4double Edep = (*s1_HC)[i]->GetEdep();
		s1_phlight[DetID] += Songresp(Edep, Pname);
		if (s1_time[DetID] < 0 || Time < s1_time[DetID])
		{
			s1_time[DetID] = Time;
		}
	}
	
	for (int i = 0; i < Nohits_s2u; i++)
	{
		G4int DetID = (*s2u_HC)[i]->GetDetID();
		G4String Pname = (*s2u_HC)[i]->GetPname();
		G4double Time = (*s2u_HC)[i]->GetTime();
		G4double Edep = (*s2u_HC)[i]->GetEdep();
		s2u_phlight[DetID] += Songresp(Edep, Pname);
		if (s2u_time[DetID] < 0 || Time < s2u_time[DetID])
		{
			s2u_time[DetID] = Time;
		} 
	}
	
	for (int i = 0; i < Nohits_s2b; i++)
	{
		G4int DetID = (*s2b_HC)[i]->GetDetID();
		G4String Pname = (*s2b_HC)[i]->GetPname();
		G4double Time = (*s2b_HC)[i]->GetTime();
		G4double Edep = (*s2b_HC)[i]->GetEdep();
		
		s2b_phlight[DetID] += Songresp(Edep, Pname);
		if (s2b_time[DetID] < 0 || Time < s2b_time[DetID])
		{
			s2b_time[DetID] = Time;
		} 
	}

    G4double s1PhThreshold = 1;
    G4double s2PhThreshold = 1;
    
    //S2_up(detID < 40) S2_bottom(detID >= 40)
	for(int iS1No = 0; iS1No < 5;iS1No++)
	{
		if ((s1_phlight[iS1No] > s1PhThreshold) && (s1_time[iS1No] > 0))
		{
			for (int iS2No = 0;iS2No < 40;iS2No++)
			{
				if ((s2u_phlight[iS2No] > s2PhThreshold) && (s2u_time[iS2No] > 0))
				{
					G4double tflight = s2u_time[iS2No] - s1_time[iS1No];
                    analysisManager->FillNtupleDColumn(0, tflight);
                    analysisManager->FillNtupleDColumn(1, s1_phlight[iS1No]);
                    analysisManager->FillNtupleDColumn(2, s2u_phlight[iS2No]);
                    analysisManager->FillNtupleIColumn(3, iS1No);
                    analysisManager->FillNtupleIColumn(4, iS2No);
                    analysisManager->FillNtupleIColumn(5, eventID);
                    analysisManager->AddNtupleRow();
				}
			}
			
			for (int iS2No = 0;iS2No < 40;iS2No++)
			{
				if ((s2b_phlight[iS2No] > s2PhThreshold) && (s2b_time[iS2No] > 0))
				{
					G4double tflight = s2b_time[iS2No] - s1_time[iS1No];
                    analysisManager->FillNtupleDColumn(0, tflight);
                    analysisManager->FillNtupleDColumn(1, s1_phlight[iS1No]);
                    analysisManager->FillNtupleDColumn(2, s2b_phlight[iS2No]);
                    analysisManager->FillNtupleIColumn(3, iS1No);
                    analysisManager->FillNtupleIColumn(4, iS2No + 40);
                    analysisManager->FillNtupleIColumn(5, eventID);
                    analysisManager->AddNtupleRow();
				}
			}
		}
	} 
	
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
