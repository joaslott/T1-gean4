//
// Detector construction for MGC-20 detectorposition T1 
//
// Production of 18F
// Gadolinia using a proton beam from the MGC-20 cyclotron and water target system (PET
// 1054). An amount of 800 mg [18O]H2O was irradiated for about 30 min using a particle beam
// current of 10 &#956;A. The water target chamber is made of silver and the water volume is
// sealed with a 50 &#956;m thick nickel foil (Goodfellow Cambridge Limited, NI000391,
// purity 99.98%). The foil was cooled
// with a flow of helium gas. The proton energy incident on the nickel foil was about 16.6
// MeV and on the water 15.7 MeV (SRIM-2006.02). The water was purchased
// from Rotem Industries Ltd. (Lot #: 14-1445) and contained 98.3 atom% of 18O, 1.0
// atom% of 17O and 0.7 atom% of 16O. The elemental impurities in the water were
// measured by the manufacturer. The metal impurities were < 0.1 mg/L."
//

#include "T1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"

// Constructor
T1DetectorConstruction::T1DetectorConstruction()
  :G4VUserDetectorConstruction(),
   fScoringVolume(0)
{ }

// Destructor
T1DetectorConstruction::~T1DetectorConstruction()
{ }


G4VPhysicalVolume* T1DetectorConstruction::Construct()
{
  
  //************************** Materials *************************************
  G4double fractionmass, density, pressure, temperature;
  G4double abudance;
  G4String name, symbol;
  G4int nel, natoms;
  G4int ncomponents;
  //G4bool isotopes = true;

  //Get nist material
  G4NistManager* nist = G4NistManager::Instance();

  //Isotopes
  G4Isotope* O18 = new G4Isotope("O18", 8, 18, 18.0*g/mole); //z, N, a  
  G4Isotope* O16 = new G4Isotope("O16", 8, 16, 15.99*g/mole);
  G4Isotope* O17 = new G4Isotope("O17", 8, 17, 17.0*g/mole);
  
  //Elements and pure materials
  G4Material* Ag = nist->FindOrBuildMaterial("G4_Ag");
  G4Material* Ni = nist->FindOrBuildMaterial("G4_Ni");
  //G4Material* Fe = nist->FindOrBuildMaterial("G4_Fe");
  //G4Element* elN = nist->FindOrBuildElement("G4_N");
  G4Element* elC = nist->FindOrBuildElement("C");
  G4Element* elFe = nist->FindOrBuildElement("Fe");
  G4Element* elO = nist->FindOrBuildElement("O");
  G4Element* elSi = nist->FindOrBuildElement("Si");
  G4Element* elCa = nist->FindOrBuildElement("Ca");
  G4Element* elAl = nist->FindOrBuildElement("Al");
  //G4Element* elS = nist->FindOrBuildElement("S");
  //G4Element* elBa = nist->FindOrBuildElement("Ba");
  G4Element* elEu = nist->FindOrBuildElement("Eu");
  G4Element* elNa = nist->FindOrBuildElement("Na");
  G4Element* elCl = nist->FindOrBuildElement("Cl");
  G4Element* elH = nist->FindOrBuildElement("H");

  G4Element* elOw = new G4Element("18Oxygen", "18O", ncomponents=3);
  elOw->AddIsotope(O18, abudance = 98.3*perCent);
  elOw->AddIsotope(O17, abudance = 0.7*perCent);
  elOw->AddIsotope(O16, abudance = 1.0*perCent);

//Materials
  // Air
  G4Material* Air = nist->FindOrBuildMaterial("G4_AIR");
 
  //vacuum data from vacaero.com  
  density     = 1.3e-10*g/m3;
  pressure    = 1.3e-8*pascal;
  temperature = 3.1e-28*kelvin;  // From PhysicalConstants.h 
  G4Material* G4vacuum = new G4Material(name="G4vacuum", density, nel=1,
                                          kStateGas, temperature, pressure);
  G4vacuum->AddMaterial(Air, fractionmass=1.);

  // 18O Water
  density = 1.11*g/cm3;
  G4Material* O18Water = new G4Material("H2-18O", density, nel = 2);
  O18Water -> AddElement(elH, natoms = 2);
  O18Water -> AddElement(elOw, natoms = 1);

  //Salt
  density = 2.16*g/cm3;
  G4Material* NaCl = new G4Material("Salt", density, nel = 2);
  NaCl -> AddElement(elNa, natoms = 1);
  NaCl -> AddElement(elCl, natoms = 1);

  // Water
  G4Material* water = nist->FindOrBuildMaterial("G4_WATER");

  //  G4Material* Concrete = nist->FindOrBuildMaterial("G4_CONCRETE");

  //Concrete
  density = 2.3*g/cm3;
  G4Material* Concrete = new G4Material("Concrete", density, nel = 8);
  Concrete -> AddElement(elO, 24.97*perCent);
  Concrete -> AddElement(elC, 24.97*perCent);
  Concrete -> AddElement(elSi, 18.50*perCent);
  Concrete -> AddElement(elCa, 20.35*perCent);
  Concrete -> AddElement(elAl, 3.70*perCent);
  Concrete -> AddElement(elFe, 3.00*perCent);
  Concrete -> AddElement(elEu, (9.0e-7)*perCent);
  Concrete -> AddMaterial(water, 4.50*perCent);

  // Magnetit concrete Fe 61.74 %, magnetit (Fe2O3) 17.89 %, Concrete 20.37 %
  // A.A. El-Sawy, Arab Journal of Nuclear Science and Applications, vol 50, 3, (151-158) 2017 
  density = 5.2*g/cm3;
  G4Material* Fe2O3 = new G4Material("Magnetit", density, nel = 2);
  Fe2O3 -> AddElement(elFe, natoms = 3);
  Fe2O3 -> AddElement(elO, natoms = 4);
  density = 5.11*g/cm3;
  G4Material* MagnetitConcrete = new G4Material("MagConcrete", density, nel = 3);
  MagnetitConcrete -> AddElement(elFe, 61.74*perCent);
  MagnetitConcrete -> AddMaterial(Fe2O3, 26.19*perCent);
  MagnetitConcrete -> AddMaterial(Concrete, 12.07*perCent);



  // variables
  G4double worldXYZ = 4.0*m;              // half of the world size
  G4double magnetitWallx = 3.19*m;        // half of the wall width
  G4double magnetitWally = 2.10*m;        // half of the wall hight
  G4double magnetitWallz = 0.75*m;        // half of the wall thickness
  G4double magnetitShieldingx = 1.25*m;    // half of the shielding width
  G4double magnetitShieldingy = 2.10*m;   // half of the shielding hight
  G4double magnetitShieldingz = 0.25*m;   // half of the shielding thickness
  G4double concreteWallx = 0.15*m;        // half of the wall width
  G4double concreteWally = 2.10*m;        // half of the wall hight
  G4double concreteWallz = 3.75*m;        // half of the wall thickness
  G4double innerRadius = 0.00*mm;         // tubs inner radius
  G4double targetRadius = 9.50*mm;        // radius of the water target
  G4double outerRadius = 10.00*mm;        // outer radius of the target chamber
  G4double beamlineRadius = 10.10*mm;     // beamline vacuum radius
  G4double targetLength = 1.41*mm;        // half length in Z of the water target (could be calculated here)
  G4double chamberLength = 1.66*mm;       // half length in Z of the target chamber 0.5 mm chamber back thickness
  G4double nickelfoilLength = 0.025*mm;   // half length in Z of the nickel foil
  G4double beamlineLength = 0.5*m;        //half length of the beamline
  G4double startAngle = 0.*deg;           // full cylinder
  G4double spanningAngle = 360.*deg;      // full cylinder

  G4bool checkOverlaps = true;  

  // Rotaton matrix
  //G4RotationMatrix rotm = G4RotationMatrix();
  G4ThreeVector position = G4ThreeVector();
  G4Transform3D transform = G4Transform3D();

  // Air filled world
  G4Box* solidWorld = new G4Box("World", worldXYZ, worldXYZ, worldXYZ);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, // solid 
						    Air,        // material
						    "World",    // name
						    0,          // field manager
						    0,          // sensitive detector
						    0);         // user limits

  G4VPhysicalVolume* physicalWorld = new G4PVPlacement(0,                                  // rotation
						       G4ThreeVector(0, 0, 0),             // translation, at (0,0,0)
						       logicWorld,                         // logical volume
						       "World",                            // name
						       0,                                  // mother physical volume
						       false,                              // boolean operator
						       0,                                  // copy number
						       checkOverlaps);                     // check for overlaps

  // Concrete walls
  G4Box* solidConcreteWall = new G4Box("ConcreteWall",concreteWallx, concreteWally, concreteWallz);
  G4LogicalVolume* logicConcreteWall = new G4LogicalVolume(solidConcreteWall,    //solid
							   Concrete,             // materia change!!
							   "ConcreteWall",       // name
							   0,                    // field manager
							   0,                    // sensitive detector
							   0);                   // user limits

  new G4PVPlacement(0,                            // rotation
		    G4ThreeVector(-3.4*m, 0, 0), // translation
		    logicConcreteWall,            // logic volume
		    "ConcreteWallLv1",            // name
		    logicWorld,                   // mother volume
		    false,                        // boolean operator
		    1,                            // copy number
		    checkOverlaps);               // check for overlaps

 new G4PVPlacement(0,                             // rotation
		   G4ThreeVector(3.4*m, 0, 0),   // translation
		    logicConcreteWall,            // logic volume
		    "ConcreteWallLv1",            // name
		    logicWorld,                   // mother volume
		    false,                        // boolean operator
		    2,                            // copy number
		    checkOverlaps);               // check for overlaps


  // Magnetit wall to cyclotron 
  G4Box* solidMagnetitWall = new G4Box("MagnetitWall", magnetitWallx, magnetitWally, magnetitWallz);
  G4LogicalVolume* logicMagnetitWall = new G4LogicalVolume(solidMagnetitWall,  // solid
							 MagnetitConcrete,    // material to be defined
							 "MagnetitWall",     // name
							 0,                  // field manager
							 0,                  // sensitive detektor
							 0);                 // user limits
  new G4PVPlacement(0,                                // rotation
		    G4ThreeVector(0, 0, -2.9*m),        // translation
		    logicMagnetitWall,                  // logic volume
		    "MagnetitWall",                     // name
		    logicWorld,                         // mother volume
                    false,                              // boolean operator
		    0,                                  // copy number
		    checkOverlaps);                     // check for overlaps

  
  // Magnetit shielding
  G4Box* solidMagnetitShielding = new G4Box("MagnetitShielding",magnetitShieldingx , magnetitShieldingy, magnetitShieldingz);

  G4LogicalVolume* logicMagnetitShielding = new G4LogicalVolume(solidMagnetitShielding,  // solid
							 MagnetitConcrete,               // material to be defined
							 "MagnetitShieldingLV",           // name
							 0,                             // field manager
							 0,                             // sensitive detektor
							 0);                            // user limits
 
  new G4PVPlacement(0,                                  // rotation
		    G4ThreeVector(-0.35*m, 0, 1.55*m),  // translation
		    logicMagnetitShielding,             // logic volume
		    "MagnetitShielding1",               // name
		    logicWorld,                         // mother volume
                    false,                              // boolean operator
		    1,                                  // copy number
		    checkOverlaps);                     // check for overlaps

  new G4PVPlacement(0,                                  // rotation
		    G4ThreeVector(-1.95*m, 0, 2.9*m),    // translation
		    logicMagnetitShielding,             // logic volume
		    "MagnetitShielding2",               // name
		    logicWorld,                         // mother volume
                    false,                              // boolean operator
		    2,                                  // copy number
		    checkOverlaps);                     // check for overlaps

  // Water shielding 1  
  G4Box* solidWaterShielding1 = new G4Box("WaterShielding1", 0.15*m, magnetitShieldingy, 0.70*m);
  G4LogicalVolume* logicWaterShielding1 = new G4LogicalVolume(solidWaterShielding1, //solid
							      water,                // material
							      "WaterShieldingLv1",  // name
							      0,                    // field manager
							      0,                    // sensitive detector
                                                              0);                   // user limits
  
  position = G4ThreeVector(-0.70*m, 0, -1.24*m);
  G4RotationMatrix rotm1 = G4RotationMatrix();
  rotm1.rotateY(-20*deg);
  transform = G4Transform3D(rotm1,position);
  new G4PVPlacement(transform,              // rotatio + translation
		    logicWaterShielding1,   // logical volume
		    "WaterShielding1",      // name
		    logicWorld,             // mother volyme
		    false,                  // boolean operator
		    0,                      // copy number
		    checkOverlaps);         // check for overlaps
  
// Water shielding 2  
  G4Box* solidWaterShielding2 = new G4Box("WaterShielding1", 0.15*m, magnetitShieldingy, 0.53*m);
  G4LogicalVolume* logicWaterShielding2 = new G4LogicalVolume(solidWaterShielding2, //solid
							      water,                // material
							      "WaterShieldingLv2",  // name
							      0,                    // field manager
							      0,                    // sensitive detector
                                                              0);                   // user limits
  G4RotationMatrix rotm2 = G4RotationMatrix();
  position = G4ThreeVector(-1.25*m, 0, -0.25*m);
  rotm2.rotateY(-76*deg);
  transform = G4Transform3D(rotm2,position);
  new G4PVPlacement(transform,              // rotatio + translation
		    logicWaterShielding2,   // logical volume
		    "WaterShielding2",      // name
		    logicWorld,             // mother volyme
		    false,                  // boolean operator
		    0,                      // copy number
		    checkOverlaps);         // check for overlaps

  //Beamline
  G4Tubs* solidBeamline = new G4Tubs("Beamline", innerRadius, beamlineRadius, beamlineLength, startAngle, spanningAngle);
  G4LogicalVolume* logicBeamline = new G4LogicalVolume(solidBeamline,     // solid
						       G4vacuum,          // material
						       "Beamline",        // name
						       0,                 // field manager
						       0,                 // sensitive detector
						       0);                // user limits

  // -1.3 , 0, -1.25
  G4RotationMatrix rotm3 = G4RotationMatrix();
  position = G4ThreeVector(-1.3*m, -0.70*m, -1.35*m);
  rotm3.rotateY(-30*deg);
  transform = G4Transform3D(rotm3,position);
  new G4PVPlacement(transform,                                   // rotation and translation
		    // G4ThreeVector(-1.3*m, 0, -1.65*m),           // translation
		    logicBeamline,                               // logic volume
		    "Beamline",                                  // name
		    logicWorld,                                  // mother volume
		    false,                                       // boolean operation
		    0,                                           // copy number
		    checkOverlaps);                              // check for overlaps
  
  // Target chamber a silver cylinder
  G4Tubs* solidTargetChamber = new G4Tubs("TargetChamber", innerRadius, outerRadius, chamberLength, startAngle, spanningAngle); 
  G4LogicalVolume* logicTargetChamber = new G4LogicalVolume(solidTargetChamber, // solid
							    Ag,                 // material
							    "TargetChamber",    // name
							    0,                  // field manager
							    0,                  // sensitive detector
							    0);                 // user limits
 
  new G4PVPlacement(0,                                           // rotation
		    G4ThreeVector(0, 0, 0),                      // translation
		    logicTargetChamber,                          // logical volume
		    "TargetChamber",                             // name
		    logicBeamline,                               // mother volume
		    false,                                       // no boolean operator
		    0,                                           // copy number
		    checkOverlaps);                              // check for overlaps
                                                           
// Target chamber seal, nickel
  G4Tubs* solidTargetSeal = new G4Tubs("TargetSeal", innerRadius, outerRadius, nickelfoilLength, startAngle, spanningAngle); 
  G4LogicalVolume* logicTargetSeal = new G4LogicalVolume(solidTargetSeal,    // solid
							 Ni,                 // material
							 "TargetSeal",       // name
							 0,                  // field manager
							 0,                  // sensitive detector
							 0);                 // user limits
 
  new G4PVPlacement(0,                                                   // rotation
		    G4ThreeVector(0, 0, -1.715e-3*m),                    // translation
		    logicTargetSeal,                                     // logical volume
		    "TargetSeal",                                        // name
		    logicBeamline,                                       // mother volume
		    false,                                               // no boolean operator
		    0,                                                   // copy number
		    checkOverlaps);                                      // check for overlaps
                                                           
// Water target H2 18-O
  G4Tubs* solidWaterTarget = new G4Tubs("WaterTarget", innerRadius, targetRadius, targetLength, startAngle, spanningAngle); 
  G4LogicalVolume* logicWaterTarget = new G4LogicalVolume(solidWaterTarget,   // solid
							  O18Water,           // material DEFINE!!
							  "WaterTarget",      // name
							  0,                  // field manager
							  0,                  // sensitive detector
							  0);                 // user limits
 
  new G4PVPlacement(0,                                                  // rotation
		    G4ThreeVector(0, 0, -0.25e-3*m),                    // translation
		    logicWaterTarget,                                   // logical volume
		    "WaterTarget",                                      // name
		    logicTargetChamber,                                 // mother volume
		    false,                                              // no boolean operator
		    0,                                                  // copy number
		    checkOverlaps);                                     // check for overlaps

  // Scoring volume water target
  fScoringVolume = logicWaterTarget;

  return physicalWorld;
}
