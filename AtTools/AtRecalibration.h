/**********************************************************************
 *   MassInformation class reads information from the amdc table      *
 *   for all the particles                                            *
 *   KinematicsCalculation obtains the kinematic lines and corrects   *
 *   the energy deviation                                             *
 *   Author: A. Munoz-Ramos alicia.munoz@rai.usc.es                   *
 *   Last Reviwed: 18/06/24                                           *
 *                                                                    *
 **********************************************************************/

// HEAVY == RECOIL
// LIGHT == SCATTERED

#ifndef ATRECALIBRATION_H
#define ATRECALIBRATION_H

#include <TF1.h>
#include <TSystem.h>

#include "TGraph.h"
#include "TH2.h"

#include <string>
#include <vector>

namespace AtTools {

class MassInformation {
public:
   MassInformation(int ABeam, std::string beamName, int ATarget, std::string targetName);
   MassInformation(int ABeam, int ZBeam, int ATarget, int ZTarget);

   MassInformation(int AElement, int ZElement);
   MassInformation(int AElement, std::string elementName);

   void SetLight(int ALight, std::string lightName);
   void SetLight(int ALight, int ZLight);
   void SetHeavy(int AHeavy, std::string heavyName);
   void SetHeavy(int AHeavy, int ZHeavy);

   double GetBeamMass() { return fMassBeam; }; // in MeV
   double GetTargetMass() { return fMassTarget; };
   double GetLightMass() { return fMassLight; };
   double GetHeavyMass() { return fMassHeavy; };
   double GetMass() { return fMassElement; }

   std::pair<int, int> GetBeamAZ() { return {fABeam, fZBeam}; } // first=A; second=Z
   std::pair<int, int> GetTargetAZ() { return {fATarget, fZTarget}; }
   std::pair<int, int> GetLightAZ() { return {fALight, fZLight}; }
   std::pair<int, int> GetHeavyAZ() { return {fAHeavy, fZHeavy}; }
   std::pair<int, int> GetAZ() { return {fAElement, fZElement}; }

   std::string GetBeamName() { return fBeamName; };
   std::string GetTargetName() { return fTargetName; };
   std::string GetLightName() { return fLightName; };
   std::string GetHeavyName() { return fHeavyName; };
   std::string GetName() { return fElementName; };

   double uma2Mev(double value) { return value * 931.494061; }
   double Mev2uma(double value) { return value / 931.494061; }

   void ReadMassTable();
   void Print();

private:
   int fABeam = -1;
   int fZBeam = -1;
   std::string fBeamName = "beam";
   double fMassBeam = -1; // MeV ->Already converted from the table

   int fATarget = -1;
   int fZTarget = -1;
   std::string fTargetName = "target";
   double fMassTarget = -1;

   int fALight = -1;
   int fZLight = -1;
   std::string fLightName = "light";
   double fMassLight = -1;

   int fAHeavy = -1;
   int fZHeavy = -1;
   std::string fHeavyName = "heavy";
   double fMassHeavy = -1;

   int fAElement = -1;
   int fZElement = -1;
   std::string fElementName = "element";
   double fMassElement = -1;

   std::string fMassFile = std::string(gSystem->Getenv("VMCWORKDIR")) + "/AtTools/massTable.txt";
   std::vector<int> fNeutronNumber, fProtonNumber, fANumber;
   std::vector<std::string> fElementShortName;
   std::vector<double> fMassUma;
};

class KinematicsCalculation {
public:
   KinematicsCalculation(MassInformation reactionInfo);
   ~KinematicsCalculation();

   // Set Excitation Energy for each particle
   void SetBeamEnergy(double value)
   {
      fBeamEnergyA = value;
      fBeamEnergyMeV = fBeamEnergyA * fReactionInfo.GetBeamAZ().first;
   } // beamEnergy units= MeV/A

   void SetExEnergyBeam(double value)
   {
      fExBeam = value;
      fwMassBeam = fReactionInfo.GetBeamMass() + fExBeam;
   };
   void SetExEnergyTarget(double value)
   {
      fExTarget = value;
      fwMassTarget = fReactionInfo.GetTargetMass() + fExTarget;
   };
   void SetExEnergyLight(double value)
   {
      fExLight = value;
      fwMassLight = fReactionInfo.GetLightMass() + fExLight;
   };
   void SetExEnergyHeavy(double value)
   {
      fExHeavy = value;
      fwMassHeavy = fReactionInfo.GetHeavyMass() + fExHeavy;
   };

   // Get Theoretical Kinematic Values
   std::vector<double> GetHeavyEnergy() { return fEnergyHeavy; };
   std::vector<double> GetLightEnergy() { return fEnergyLight; };
   std::vector<double> GetHeavyLabAngle() { return fAngleLABHeavy; };
   std::vector<double> GetLightLabAngle() { return fAngleLABLight; };
   std::vector<double> GetCMAngle() { return fThetaCM; };

   // Theoretical Kinematic Lines
   TGraph *lightGraph{};
   TGraph *heavyGraph{};
   TGraph *relationGraph{};
   TGraph *lightLabCMGraph{};

   // Plot of kinematic lines
   void Plot();
   void CalculateKinematics();

   // Real data from recalibration
   double CorrectEnergy(double beamCorrectedEnergy, double thetaLab, double originalEnergy);
   void FitEnergyCorrection(TH2F *hOriginalCorrectedEnergy, bool draw = false);
   std::pair<double, double> GetFitParameters() { return {fFuncPol1->GetParameter(0), fFuncPol1->GetParameter(1)}; };
   double Recalibration(double originalEnergy);
   std::pair<double, double> ExAndThetaCM(double energyLab, double thetaLAB, double TBeamCorrected);

private:
   double fBeamEnergyA{0.}, fBeamEnergyMeV{0.0};
   double fExTarget{0.}, fExBeam{0.}, fExLight{0.}, fExHeavy{0.};
   double fwMassBeam{0.}, fwMassTarget{0.}, fwMassLight{0.}, fwMassHeavy{0.};

   MassInformation fReactionInfo;

   std::vector<double> fEnergyLight, fEnergyHeavy, fThetaCM, fAngleLABLight, fAngleLABHeavy;

   TF1 *fFuncPol1 = new TF1("pol1", "pol1");
   bool fNoSolution = false;
};
} // namespace AtTools

#endif // ATRECALIBRATION_H
