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

#include "TOFRunAction.hh"
#include "Analysis.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TOFRunAction::TOFRunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each 100 events
  //G4RunManager::GetRunManager()->SetPrintProgress(1000);
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4cout << "Using " << analysisManager->GetType() << G4endl;

    // Create directories 
    //analysisManager->SetHistoDirectoryName("histograms");
    //analysisManager->SetNtupleDirectoryName("ntuple");
    analysisManager->SetVerboseLevel(0);
    
    // Creating ntuple
    //
    analysisManager->SetFirstNtupleId(1);
    
    analysisManager->CreateNtuple("S1", "S1event");
    analysisManager->CreateNtupleDColumn(1, "time");
    analysisManager->CreateNtupleDColumn(1, "Eee_s1");
    analysisManager->CreateNtupleIColumn(1, "s1ID");
    analysisManager->CreateNtupleIColumn(1, "eventID");
    analysisManager->CreateNtupleIColumn(1, "s1ProtonHits");
    analysisManager->CreateNtupleIColumn(1, "s1CarbonHits");
    analysisManager->FinishNtuple(1);
    
    /*
    analysisManager->CreateNtuple("S2", "S2event");
    analysisManager->CreateNtupleDColumn(2, "time");
    analysisManager->CreateNtupleDColumn(2, "Eee_s2");
    analysisManager->CreateNtupleIColumn(2, "s2ID");
    analysisManager->CreateNtupleIColumn(2, "eventID");
    analysisManager->CreateNtupleIColumn(2, "s2trackID");
    analysisManager->FinishNtuple(2);
    */
    
    analysisManager->CreateNtuple("TOF", "TOFevent");
    analysisManager->CreateNtupleDColumn(2, "TimeOfFlight");
    analysisManager->CreateNtupleDColumn(2, "Eee_s1");
    analysisManager->CreateNtupleDColumn(2, "Eee_s2");
    analysisManager->CreateNtupleIColumn(2, "s1ID");
    analysisManager->CreateNtupleIColumn(2, "s2ID");
    analysisManager->CreateNtupleIColumn(2, "eventID");
    //analysisManager->CreateNtupleIColumn(3, "s2trackID");
    analysisManager->FinishNtuple(2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TOFRunAction::~TOFRunAction()
{
    delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TOFRunAction::BeginOfRunAction(const G4Run*)
{ 

    //inform the runManager to save random number seed
    //G4RunManager::GetRunManager()->SetRandomNumberStore(false);
    
    // Get analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    // Open an output file
    //
    G4String fileName = "TOF";
    analysisManager->OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TOFRunAction::EndOfRunAction(const G4Run* )
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    // save histograms & ntuple
    //
    analysisManager->Write();
    analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
