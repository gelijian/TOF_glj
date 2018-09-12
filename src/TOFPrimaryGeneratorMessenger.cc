#include "TOFPrimaryGeneratorMessenger.hh"
#include "TOFPrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

TOFPrimaryGeneratorMessenger::TOFPrimaryGeneratorMessenger(TOFPrimaryGeneratorAction* Gun):Action(Gun)
{
	GunDir = new G4UIdirectory("/TOF/Gun/");
	GunDir->SetGuidance("UI commands of Paticle Gun Action");
	GunDir->SetGuidance("Self-defination or GPS Source");

	GPSFlagCmd = new G4UIcmdWithAString("/TOF/Gun/GPSFlag", this);
	GPSFlagCmd->SetGuidance("Set mode of particle generation");
	GPSFlagCmd->SetGuidance("Choice : off, on, self(default)");
	GPSFlagCmd->SetParameterName("choice",true);
	GPSFlagCmd->SetDefaultValue("self");
	GPSFlagCmd->SetCandidates("on off self");
	GPSFlagCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	ParticleNameCmd = new G4UIcmdWithAString("/TOF/Gun/ParticleName",this);
	ParticleNameCmd->SetGuidance("Set Type of source particle");
	ParticleNameCmd->SetGuidance("ParticleType : gamma, neutron(default)");
	ParticleNameCmd->SetParameterName("ParticleName", true);
	ParticleNameCmd->SetDefaultValue("neutron");
	ParticleNameCmd->SetCandidates("gamma neutron");
	ParticleNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	SourceTypeCmd = new G4UIcmdWithAString("/TOF/Gun/SourceType", this);
	SourceTypeCmd->SetGuidance("Set Type of source energy");
	SourceTypeCmd->SetGuidance("SourceType : mono(default), flat, Cf252, expdecay, gauss, rect");
	SourceTypeCmd->SetParameterName("SourceType",true);
	SourceTypeCmd->SetDefaultValue("mono");
	SourceTypeCmd->SetCandidates("mono flat Cf252 expdecay gauss rect Tabulated");
	SourceTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	SourceEnergyCmd = new G4UIcmdWithADoubleAndUnit("/TOF/Gun/SourceEnergy",this);
	SourceEnergyCmd -> SetGuidance("Set energy of particle");
	SourceEnergyCmd -> SetParameterName("SourceEnergy",false);
	SourceEnergyCmd -> SetRange("SourceEnergy>0.");
	SourceEnergyCmd -> SetUnitCategory("Energy");
	SourceEnergyCmd -> AvailableForStates(G4State_PreInit, G4State_Idle);

	SEconst = new G4UIcmdWithADoubleAndUnit("/TOF/Gun/SEconst",this);
	SEconst -> SetGuidance("Set energy constant of sourse");
	SEconst -> SetParameterName("SEconst",false);
	SEconst -> SetRange("SEconst>=0.");
	SEconst -> SetUnitCategory("Energy");
	SEconst -> AvailableForStates(G4State_PreInit,G4State_Idle);

	SourceRadiusCmd = new G4UIcmdWithADoubleAndUnit("/TOF/Gun/SourceRadius",this);
	SourceRadiusCmd -> SetGuidance("Set radius of source");
	SourceRadiusCmd -> SetParameterName("SourceRadius",false);
	SourceRadiusCmd -> SetRange("SourceRadius>=0.");
	SourceRadiusCmd -> SetUnitCategory("Length");
	SourceRadiusCmd -> AvailableForStates(G4State_PreInit, G4State_Idle);

	SourcePosition_zCmd = new G4UIcmdWithADoubleAndUnit("/TOF/Gun/SourcePosition_z",this);
	SourcePosition_zCmd -> SetGuidance("Set z position of source");
	SourcePosition_zCmd -> SetParameterName("SourcePosition_z",false);
	SourcePosition_zCmd -> SetUnitCategory("Length");
	SourcePosition_zCmd -> AvailableForStates(G4State_PreInit,G4State_Idle);

	MomentumThetaRangeCmd = new G4UIcmdWithADoubleAndUnit("/TOF/Gun/MomentumThetaRange", this);
	MomentumThetaRangeCmd -> SetGuidance("Set the source momentum direction range");
	MomentumThetaRangeCmd -> SetParameterName("MomentumThetaRange", false);
	MomentumThetaRangeCmd -> SetRange("MomentumThetaRange>=0.");
	MomentumThetaRangeCmd -> SetUnitCategory("Angle");
	MomentumThetaRangeCmd -> AvailableForStates(G4State_PreInit,G4State_Idle);
}

TOFPrimaryGeneratorMessenger::~TOFPrimaryGeneratorMessenger()
{
	delete GPSFlagCmd;
	delete ParticleNameCmd;
	delete SourceTypeCmd;
	delete SourceEnergyCmd;
	delete SEconst;
	delete SourceRadiusCmd;
	delete SourcePosition_zCmd;
	delete GunDir;
}

void TOFPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String value)
{
	if(command == GPSFlagCmd)
	{  Action -> SetGPSFlag(value); }

	if(command == ParticleNameCmd)
	{  Action -> SetParticleName(value); }

	if(command == SourceTypeCmd)
	{  Action -> SetSoureType(value); }

	if( command == SourceEnergyCmd)
	{ Action -> SetSourceEnergy(SourceEnergyCmd->GetNewDoubleValue(value)); }

	if( command == SourceRadiusCmd)
	{ Action -> SetSourceRadius(SourceRadiusCmd->GetNewDoubleValue(value)); }

	if( command == SourcePosition_zCmd)
	{ Action -> SetSourcePosition_z(SourcePosition_zCmd->GetNewDoubleValue(value)); }

	if( command == SEconst)
	{ Action -> SetEnergyConst(SEconst -> GetNewDoubleValue(value));}

	if( command == MomentumThetaRangeCmd)
	{ Action -> SetMomentumDirection(MomentumThetaRangeCmd -> GetNewDoubleValue(value));}

}
