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
//  Liqiuid Scintillator  
 
#include "TOFDetectorConstruction.hh"
#include "TOFscinSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4BooleanSolid.hh"

#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4AssemblyVolume.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
TOFDetectorConstruction::TOFDetectorConstruction()
:G4VUserDetectorConstruction(), 
 expHallPV(NULL),
 expHallLV(NULL),
 WallLV(NULL),
 S1_scinLV(NULL),
 S1_MSLV(NULL),
 S1_LGLV(NULL),
 S2_MSLV(NULL),
 S2u_scinLV(NULL),
 S2u_LGLV(NULL),
 S2b_scinLV(NULL),
 S2b_LGLV(NULL),
 AirMaterial(NULL),
 ConcreteMaterial(NULL), 
 S1Material(NULL),
 S2Material(NULL),
 LGMaterial(NULL),
 MSMaterial(NULL),
 fCheckOverlaps(true)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
TOFDetectorConstruction::~TOFDetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* TOFDetectorConstruction::Construct()
{
  // Define m   aterials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TOFDetectorConstruction::DefineMaterials()
{
    // ---------Material definition---------     

    G4NistManager* nistManager = G4NistManager::Instance();

    // ---------Air defined using NIST Manager--------- 
    AirMaterial = nistManager->FindOrBuildMaterial("G4_AIR");
    
    // ---------Concrete defined using NIST Manager---------
    ConcreteMaterial = nistManager->FindOrBuildMaterial("G4_CONCRETE");
    
    // ---------Element---------    
    //G4double a,   z, n;
    G4double density;
    G4int nel;
    nistManager->SetVerbose(1);
    G4Element* H = nistManager->FindOrBuildElement("H");
    G4Element* C = nistManager->FindOrBuildElement("C");
    G4Element* O  = nistManager->FindOrBuildElement("O");
    G4Element* Fe  = nistManager->FindOrBuildElement("Fe");
    G4Element* Mn  = nistManager->FindOrBuildElement("Mn");

    //--------- Plastic Scintilator EJ200 --------- 
    G4Material* EJ200_Scin = new G4Material("EJ200_Scintilator", density = 1.023*g/cm3, nel = 2);
    EJ200_Scin->AddElement(H,  8.467 * perCent);
    EJ200_Scin->AddElement(C, 91.533 * perCent);

    //--------- Plastic Scintilator EJ228 --------- 
    G4Material* EJ228_Scin = new G4Material("EJ228_Scintilator", density = 1.023*g/cm3, nel = 2);
    EJ228_Scin->AddElement(H,  8.4376 * perCent);
    EJ228_Scin->AddElement(C, 91.5624 * perCent);
    
    //--------- PMMA light guide --------- 
    G4Material* PMMA_EJ = new G4Material("PMMA_lightguide_EJ",density= 1.19*g/cm3, nel=3);
    PMMA_EJ->AddElement(H, 8.0605 * perCent);
    PMMA_EJ->AddElement(C, 60.0096 * perCent);
    PMMA_EJ->AddElement(O, 31.9299 * perCent);
    
    //--------Soft iron---------
    G4Material* SoftIron = new G4Material("soft_iron",density= 7.87*g/cm3, nel=2);
    SoftIron->AddElement(Fe, 99.8*perCent);
    SoftIron->AddElement(Mn, 0.2*perCent);

    S1Material = EJ228_Scin;
    S2Material = EJ200_Scin;
    LGMaterial = PMMA_EJ;
    MSMaterial = SoftIron;
    
    //--------- Print   all the materials   defined---------    
    G4cout <<   G4endl <<   "The materials defined are : " <<   G4endl <<   G4endl;
    G4cout <<   *(G4Material::GetMaterialTable())   << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* TOFDetectorConstruction::DefineVolumes()
{
    
    //---------expHall as world---------
    G4double expHall_x = 180.0 * cm;
    G4double expHall_y = 180.0 * cm;
    G4double expHall_z = 180.0 * cm;
    G4Box* expHallSV = new G4Box("expHallSV", expHall_x, expHall_y, expHall_z);
    
    expHallLV = new G4LogicalVolume(expHallSV,
                                    AirMaterial,
                                    "expHallLV");
                                    
    expHallPV = new G4PVPlacement(0,
                                  G4ThreeVector(),
                                  expHallLV,
                                  "expHallPV",
                                  0,
                                  false,
                                  0,
                                  fCheckOverlaps);
    
    G4double Wall_x = 120 * cm;
    G4double Wall_y = 120 * cm;
    G4double Wall_z = 20 * cm;
    G4ThreeVector pos_Wall = G4ThreeVector(0, 0, 100 * cm);
    
    G4Box * WallSV = new G4Box("WallSV", Wall_x, Wall_y, Wall_z);
    
    WallLV = new G4LogicalVolume(WallSV, ConcreteMaterial, "WallLV");
    
    G4VPhysicalVolume* WallPV = new G4PVPlacement(0,
                                                  pos_Wall,
                                                  WallLV,
                                                  "WallPV",
                                                  expHallLV,
                                                  false,
                                                  0,
                                                  fCheckOverlaps);
                               
    //--------- the principal geometrical components---------
    G4double TOFRadius   =  75.0 * cm;
    G4int NumberOfS1 = 5;
    G4int NumberOfPairs = 40;
    G4double S2u_angle = 2 * 25 * deg;  //flight angle
    G4double S2b_angle = 2 * 35 * deg;  //flight angle

    //---------S1detector---------
    
    //S1_scintillator
    G4double S1_scinRadius = 20 * mm; 
    G4double S1_scinHalfThickness = 3 * mm; 
    G4Tubs* S1_scinSV = new G4Tubs("S1_scinSV", 0 * mm, S1_scinRadius, S1_scinHalfThickness, 0 * deg, 360 * deg);
    S1_scinLV = new G4LogicalVolume(S1_scinSV, S1Material, "S1_scinLV");
    
    //S1 magenitic shielding for each PMT
    G4double S1_MSinnerRadius = 24 * mm;
    G4double S1_MSoutterRadius = 30 * mm;
    G4double S1_MSHalflength = 82.5 * mm;
    G4Tubs* S1_MSSV = new G4Tubs("S1_MSSV", S1_MSinnerRadius, S1_MSoutterRadius, S1_MSHalflength, 0 * deg, 360 * deg);
    S1_MSLV = new G4LogicalVolume(S1_MSSV, MSMaterial, "S1_MSLV");
    
    //S1 fishtail Light guide 
    G4Tubs  *S1_tracSV = new G4Tubs("S1_tracSV", 0. , S1_scinRadius, S1_scinHalfThickness + 1 * mm, 0 * deg, 360 * deg);
    G4Box* S1_LGboxSVtemp = new G4Box("S1_LGboxSVtemp", S1_scinRadius, 30 * mm, S1_scinHalfThickness);
    G4SubtractionSolid* S1_LGboxSV = new G4SubtractionSolid("S1_LGboxSV", S1_LGboxSVtemp, S1_tracSV);
    G4Trd* S1_LGFish_up = new G4Trd("S1_LGfish", 11.5 * mm, 20 * mm, 11.5 * mm, 3 * mm, 46 * mm);
    G4Trd* S1_LGFish_down = new G4Trd("S1_LGfish", 11.5 * mm, 20 * mm, 11.5 * mm, 3 * mm, 46 * mm);

    G4RotationMatrix Rot_up;
    Rot_up.rotateX(90 * deg);
    G4ThreeVector trans_up = G4ThreeVector(0.0, 76 * mm, 0.0);
    G4Transform3D transform_up(Rot_up, trans_up);
    
    G4RotationMatrix Rot_down;
    Rot_down.rotateX(-90 * deg);
    G4ThreeVector trans_down = G4ThreeVector(0.0, -76 * mm, 0.0);
    G4Transform3D transform_down(Rot_down, trans_down);

    G4UnionSolid* S1_LGsingleside = new G4UnionSolid("S1_LGsingleside", S1_LGboxSV, S1_LGFish_down, transform_down);
    G4UnionSolid* S1_LGSV = new G4UnionSolid("S1_LGSV", S1_LGsingleside, S1_LGFish_up, transform_up);
    S1_LGLV = new G4LogicalVolume(S1_LGSV, LGMaterial, "S1_LGLV");
    
    //assemble all the parts of S1
    G4AssemblyVolume* S1_detector = new G4AssemblyVolume();
    G4ThreeVector Ta;
    G4Transform3D Tr;
  
    // fill the S1_detector with S1_MSLV (two side)
    Ta.setY(-105 * mm - 82.5 * mm);
    G4RotationMatrix Ra_down;
    Ra_down.rotateX(-90 * deg);
    Tr = G4Transform3D(Ra_down, Ta);
    S1_detector->AddPlacedVolume(S1_MSLV, Tr);

    G4RotationMatrix Ra_up;
    Ra_up.rotateX(90 * deg);
    Ta.setY(105 * mm + 82.5 * mm);
    Tr = G4Transform3D(Ra_up, Ta);
    S1_detector->AddPlacedVolume(S1_MSLV, Tr);
    
    // fill theS1 with S1_LGLV
    Ta.setX(0 * mm);
    Ta.setY(0 * mm);
    Ta.setZ(0 * mm);
    G4RotationMatrix Ra_LG;
    Tr = G4Transform3D(Ra_LG, Ta);
    S1_detector->AddPlacedVolume(S1_LGLV, Tr);
    
    // fill the S1_detector with S1_sinLV
    G4RotationMatrix Ra_scin;
    Tr = G4Transform3D(Ra_scin, Ta);
    S1_detector->AddPlacedVolume(S1_scinLV, Tr);
    
    for (G4int i = 0; i < NumberOfS1 ; i++)
    {
        G4double phi = i * 180 / NumberOfS1 * deg;
        G4double pos_z =-TOFRadius + 2 * (i - 2) * (S1_scinHalfThickness + 0.5 * mm);
        G4ThreeVector pos = G4ThreeVector(0 * mm, 0 * mm, pos_z);
        G4RotationMatrix rot = G4RotationMatrix();
        rot.rotateZ(phi);     
        G4Transform3D transform = G4Transform3D(rot, pos);
        S1_detector->MakeImprint(expHallLV, transform);
    }
    
    //---------S2 magnetic shielding---------
    G4Tubs* S2_MSSV = new G4Tubs("S2_MSSV", 38.5 * mm, 48.5 * mm, 160 * mm, 0. * deg, 360 * deg);
    S2_MSLV = new G4LogicalVolume(S2_MSSV, MSMaterial, "S2_MSLV");  

    //---------S2updetector---------
    
    //S2u_scintillator
    G4double S2u_scinhy = 8.5 * mm;      // Half thickness  of S2u 
    G4double S2u_scinhx1 = 35.0 * mm;    // Half length of S2 (up) top side in x axis
    G4double S2u_scinhx2 = 50.0 * mm;    // Half length of S2 (up) bottom side in x axis
    G4double S2u_scinhz  = 140.0 * mm;   // Half length of S2 (up) in y axis
    //G4double S2u_dev = 10.0 * mm;     // The offset of the S2 (up) center
    G4Trd* S2u_scinSV = new G4Trd("S2u_scinSV", S2u_scinhx2, S2u_scinhx1, S2u_scinhy, S2u_scinhy, S2u_scinhz);
    S2u_scinLV = new G4LogicalVolume(S2u_scinSV, S2Material, "S2u_scinLV");
    
    //S2u fishtail light guide
    G4Box* S2ulg1 = new G4Box("S2ulg1", 50 * mm, 8.5 * mm, 30 * mm);
    G4Trd* S2ulg2 = new G4Trd("S2ulg2", 21.7 * mm, 50 * mm, 21.7 * mm, 8.5 * mm, 35 * mm);
    G4Tubs* S2ulg4 = new G4Tubs("S2ulg4", 0.0 , 25 * mm, 25 * mm, 0. * deg, 360 * deg);
    G4UnionSolid* S2ulg3 = new G4UnionSolid("s2ulg3", S2ulg1, S2ulg2, 0 , G4ThreeVector(0.0, 0.0, -65*mm)); 
    G4UnionSolid* S2u_LGSV = new G4UnionSolid("s2ulg", S2ulg3, S2ulg4, 0, G4ThreeVector(0.0,0.0,-125*mm)); 
    S2u_LGLV = new G4LogicalVolume(S2u_LGSV, LGMaterial, "S2u_LGLV");
    
    //assemble all the parts of S2u
    G4AssemblyVolume* S2u_detector = new G4AssemblyVolume();
    G4ThreeVector Tau;
    G4Transform3D Tru;
    
    // fill the S2u with S2u_LGLV
    Tau.setZ(-170 * mm);
    G4RotationMatrix Rau_LG;
    Tru = G4Transform3D(Rau_LG, Tau);
    S2u_detector->AddPlacedVolume(S2u_LGLV, Tru);
    
    // fill the S2U with S2_MSLV
    Tau.setZ(-410 * mm);
    G4RotationMatrix Rau_MS;
    Tru = G4Transform3D(Rau_MS, Tau);
    S2u_detector->AddPlacedVolume(S2_MSLV, Tru);
    
    // fill the S2u with S2_sinLV
    Tau.setZ(0 * mm);
    G4RotationMatrix Rau_scin;
    Tru = G4Transform3D(Rau_scin, Tau);
    S2u_detector->AddPlacedVolume(S2u_scinLV, Tru);
    
    //S2udetector placement
    for (G4int i = 0; i < NumberOfPairs; i++)
    {
        G4double theta = S2u_angle;
        G4double phi = i * 360 / NumberOfPairs * deg;
        G4double pos_x = TOFRadius * std::sin(theta) * std::cos(phi);
        G4double pos_y = TOFRadius * std::sin(theta) * std::sin(phi);
        G4double pos_z = TOFRadius * std::cos(theta);
        G4ThreeVector pos = G4ThreeVector(pos_x, pos_y, pos_z);
        
        G4RotationMatrix rotm  = G4RotationMatrix();
        rotm.rotateX(-(90 * deg - theta));
        rotm.rotateZ(phi + 90 * deg);
        G4Transform3D transform = G4Transform3D(rotm, pos);
        S2u_detector->MakeImprint(expHallLV, transform);
    }
    
    //---------S2bottomdetector---------

    //S2b_scintillator
    G4double S2b_scinhy = 8.5 * mm;      // Half thickness  of S2b 
    G4double S2b_scinhx1 = 47.5 * mm;    // Half length of S2 (up) top side in x axis
    G4double S2b_scinhx2 = 55.0 * mm;    // Half length of S2 (up) bottom side in x axis
    G4double S2b_scinhz  = 117.5 * mm;   // Half length of S2 (up) in y axis
    //G4double S2b_dev = 12.5 * mm;     // The offset of the S2 (up) center
    G4Trd* S2b_scinSV = new G4Trd("S2b_scinSV", S2b_scinhx2, S2b_scinhx1, S2b_scinhy, S2b_scinhy, S2b_scinhz);
    S2b_scinLV = new G4LogicalVolume(S2b_scinSV, S2Material, "S2b_scinLV");
    
    
    //S2b fishtail light guide
    G4Trd* S2blg1 = new G4Trd("S2blg1", 21.7 * mm, 55 * mm, 21.7 * mm, 8.5 * mm, 40 * mm);
    G4Tubs* S2blg2  =   new G4Tubs("S2blg2", 0.0, 25 * mm, 25 * mm, 0. * deg, 360 * deg);
    G4UnionSolid* S2b_LGSV = new G4UnionSolid("s2b_LGSV", S2blg1, S2blg2, 0, G4ThreeVector(0.0,0.0,-65*mm)); 
    S2b_LGLV = new G4LogicalVolume(S2b_LGSV, LGMaterial, "S2b_LGLV");   
    
    //assemble all the parts of S2b
    G4AssemblyVolume* S2b_detector = new G4AssemblyVolume();
    G4ThreeVector Tab;
    G4Transform3D Trb;
    
    // fill the S2b with S2b_LGLV
    Tab.setZ(-157.5 * mm);
    G4RotationMatrix Rab_LG;
    Trb = G4Transform3D(Rab_LG, Tab);
    S2b_detector->AddPlacedVolume(S2b_LGLV, Trb);
    
    // fill the S2b with S2_MSLV
    Tab.setZ(-410 * mm);
    G4RotationMatrix Rab_MS;
    Trb = G4Transform3D(Rab_MS, Tab);
    S2b_detector->AddPlacedVolume(S2_MSLV, Trb);
    
    // fill the S2b with S2b_sinLV
    Tab.setZ(0 * mm);
    G4RotationMatrix Rab_scin;
    Trb = G4Transform3D(Rab_scin, Tab);
    S2b_detector->AddPlacedVolume(S2b_scinLV, Trb);
    
    //S2bdetector placement
    for (G4int i = 0; i < NumberOfPairs; i++)
    {
        G4double theta = S2b_angle;
        G4double phi = i * 360 / NumberOfPairs * deg;
        G4double pos_x = TOFRadius * std::sin(theta) * std::cos(phi);
        G4double pos_y = TOFRadius * std::sin(theta) * std::sin(phi);
        G4double pos_z = TOFRadius * std::cos(theta);
        G4ThreeVector pos = G4ThreeVector(pos_x, pos_y, pos_z);
        
        G4RotationMatrix rotm  = G4RotationMatrix();
        rotm.rotateX(-(90 * deg -  theta));
        rotm.rotateZ(phi + 90 * deg);
        G4Transform3D transform = G4Transform3D(rotm, pos);
        S2b_detector->MakeImprint(expHallLV, transform);
    }
    
    //--------- Visualization attributes -------------------------------
    
    //G4Colour
    G4Colour white (1.0,1.0,1.0);
    G4Colour grey (0.5,0.5,0.5);
    //G4Colour black (0.0,0.0,0.0);
    G4Colour red (1.0,0.0,0.0);
    G4Colour blue (0.0,0.0,1.0);
    G4Colour green (0.0,1.0,0.0);
    //G4Colour cyan (0.0,1.0,1.0);
    G4Colour magenta (1.0,0.0,1.0);
    //G4Colour yellow(1.0,1.0,0.0);
    
    G4VisAttributes* expHallVisAtt= new G4VisAttributes(G4VisAttributes::Invisible);
    expHallLV -> SetVisAttributes(expHallVisAtt);
    
    G4VisAttributes* WallVisAtt = new G4VisAttributes(white);
    WallLV -> SetVisAttributes(WallVisAtt);  

    G4VisAttributes* S1_MSVisAtt = new G4VisAttributes(grey);
    S1_MSLV -> SetVisAttributes(S1_MSVisAtt);
    
    G4VisAttributes* S1_scinVisAtt = new G4VisAttributes(red);
    S1_scinLV -> SetVisAttributes(S1_scinVisAtt);
    
    G4VisAttributes* S1_LGVisAtt = new G4VisAttributes(blue);
    S1_LGLV -> SetVisAttributes(S1_LGVisAtt);
    
    G4VisAttributes* S2_MSVisAtt = new G4VisAttributes(grey);
    S2_MSLV -> SetVisAttributes(S2_MSVisAtt);
    
    G4VisAttributes* S2u_scinVisAtt = new G4VisAttributes(red);
    S2u_scinLV -> SetVisAttributes(S2u_scinVisAtt);
    
    G4VisAttributes* S2u_LGVisAtt = new G4VisAttributes(blue);
    S2u_LGLV -> SetVisAttributes(S2u_LGVisAtt);
  
    G4VisAttributes* S2b_scinVisAtt = new G4VisAttributes(red);
    S2b_scinLV -> SetVisAttributes(S2b_scinVisAtt);
    
    G4VisAttributes* S2b_LGVisAtt = new G4VisAttributes(blue);
    S2b_LGLV -> SetVisAttributes(S2b_LGVisAtt);

    //------------------------------------------------------------------
    
    return expHallPV;
}
 
void TOFDetectorConstruction::ConstructSDandField()
{
    // Sensitive detectors
    G4SDManager* SDmanager = G4SDManager::GetSDMpointer();
    G4String SDname;
    G4String HCname;
    
    G4VSensitiveDetector* S1_scinSD = new TOFscinSD(SDname = "S1_scinSD", HCname = "s1_hitscollection");
    SDmanager->AddNewDetector(S1_scinSD);
    S1_scinLV->SetSensitiveDetector(S1_scinSD);
    
    G4VSensitiveDetector* S2u_scinSD = new TOFscinSD(SDname = "S2u_scinSD", HCname = "s2u_hitscollection");
    SDmanager->AddNewDetector(S2u_scinSD);
    S2u_scinLV->SetSensitiveDetector(S2u_scinSD);
    
    G4VSensitiveDetector* S2b_scinSD = new TOFscinSD(SDname = "S2b_scinSD", HCname = "s2b_hitscollection");
    SDmanager->AddNewDetector(S2b_scinSD);
    S2b_scinLV->SetSensitiveDetector(S2b_scinSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
