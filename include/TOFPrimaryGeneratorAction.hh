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

#ifndef TOFPrimaryGeneratorAction_h
#define TOFPrimaryGeneratorAction_h 1

#include <iostream>
#include <fstream>
#include <vector>
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

using namespace std;
class G4ParticleGun;
class G4Event;
class TOFPrimaryGeneratorMessenger;

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the Tracker 
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class 
/// (see the macros provided with this example).

class TOFPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
	public:
		TOFPrimaryGeneratorAction();    
		virtual ~TOFPrimaryGeneratorAction();
		
		virtual void GeneratePrimaries(G4Event* );
		G4ParticleGun* GetParticleGun() {return fParticleGun;}
		void SelfDefineGun();
		void ReadEnergyDist();
		
		void SetGPSFlag(G4String data) { GPSFlag = data;};
		void SetSoureType(G4String data)    { SourceType = data;};
		void SetParticleName(G4String data)       { ParticleName = data;};
		void SetSourcePosition_z(G4double data) { SourcePosition_z = data; };
		void SetSourceRadius(G4double data) { SourceRadius = data; };
		void SetSourceEnergy(G4double data)  { SourceEnergy = data; };
		void SetEnergyConst(G4double data)   { EnergyConst  = data; };
		void SetMomentumDirection(G4double data)	{MomentumThetaRange = data;}
  
	private:
		G4ParticleGun*  fParticleGun;
		TOFPrimaryGeneratorMessenger* gunMessenger;
		G4String GPSFlag;
		G4String SourceType;
		G4String ParticleName;
		G4double SourceRadius;
		G4double SourcePosition_z;
		G4double SourceEnergy;
		G4double EnergyConst;
		G4double MomentumThetaRange;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
