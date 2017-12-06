
#ifndef TOFPrimaryGeneratorMessenger_h
#define TOFPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class TOFPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class TOFPrimaryGeneratorMessenger: public G4UImessenger
{
public:
	TOFPrimaryGeneratorMessenger(TOFPrimaryGeneratorAction*);
	~TOFPrimaryGeneratorMessenger();

	virtual void SetNewValue(G4UIcommand*, G4String);

private:
	TOFPrimaryGeneratorAction* Action;
	
	G4UIdirectory* GunDir;
	G4UIcmdWithAString* GPSFlagCmd;
	G4UIcmdWithAString* ParticleNameCmd;
	G4UIcmdWithAString* SourceTypeCmd;
	G4UIcmdWithADoubleAndUnit* SourceEnergyCmd;
	G4UIcmdWithADoubleAndUnit* SEconst;
	G4UIcmdWithADoubleAndUnit* SourceRadiusCmd;
	G4UIcmdWithADoubleAndUnit* SourcePosition_zCmd;
	G4UIcmdWithADoubleAndUnit* MomentumThetaRangeCmd;
};

#endif
