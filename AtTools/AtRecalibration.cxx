#include "AtRecalibration.h"

#include <FairLogger.h>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TROOT.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <string>
#include <vector>

/*
 * Constructor for A and element names
 */

AtTools::MassInformation::MassInformation(int ABeam, std::string beamName, int ATarget, std::string targetName)
   : fABeam(ABeam), fBeamName(std::move(beamName)), fATarget(ATarget), fTargetName(std::move(targetName))
{
   if ((0 < fBeamName.size() && fBeamName.size() <= 2) && (0 < fTargetName.size() && fTargetName.size() <= 2)) {
      ReadMassTable();

      bool beamFound{false}, targetFound{false};

      for (int i = 0; i < fNeutronNumber.size(); i++) {
         if (fANumber[i] == fABeam && fElementShortName[i] == fBeamName) {
            fZBeam = fProtonNumber[i];
            fMassBeam = uma2Mev(fMassUma[i]);
            beamFound = true;
         }
         if (fANumber[i] == fATarget && fElementShortName[i] == fTargetName) {
            fZTarget = fProtonNumber[i];
            fMassTarget = uma2Mev(fMassUma[i]);
            targetFound = true;
         }
         if (targetFound && beamFound) {
            break;
         }
      }
      if (!beamFound) {
         LOG(fatal) << "Beam " << fBeamName << " wit: A=" << fABeam << " not found.";
      } else if (!targetFound) {
         LOG(fatal) << "Target " << fTargetName << " with A=" << fATarget << " not found.";
      }
   } else {
      LOG(fatal) << "Element names must be 1 or 2 characters";
   }
};

/*
 * Constructor for A and Z
 */
AtTools::MassInformation::MassInformation(int ABeam, int ZBeam, int ATarget, int ZTarget)
   : fABeam(ABeam), fZBeam(ZBeam), fATarget(ATarget), fZTarget(ZTarget)
{
   ReadMassTable();

   bool beamFound{false}, targetFound{false};
   for (int i = 0; i < fNeutronNumber.size(); i++) {
      if (fANumber[i] == fABeam && fProtonNumber[i] == fZBeam) {
         fBeamName = fElementShortName[i];
         fMassBeam = uma2Mev(fMassUma[i]);
         beamFound = true;
      }
      if (fANumber[i] == fATarget && fProtonNumber[i] == fZTarget) {
         fTargetName = fElementShortName[i];
         fMassTarget = uma2Mev(fMassUma[i]);
         targetFound = true;
      }
   }
   if (!targetFound) {
      LOG(fatal) << "Target " << fTargetName << " with A=" << fATarget << " not found.";
   } else if (!beamFound) {
      LOG(fatal) << "Beam " << fBeamName << " with A=" << fABeam << " not found.";
   }
};

/*
 * Constructor for an element A and name
 */
AtTools::MassInformation::MassInformation(int AElement, std::string elementName)
   : fAElement(AElement), fElementName(std::move(elementName))
{
   if ((0 < fElementName.size() && fElementName.size() <= 2)) {
      ReadMassTable();
      for (int i = 0; i < fNeutronNumber.size(); i++) {
         if (fANumber[i] == fAElement && fElementShortName[i] == fElementName) {
            fZElement = fProtonNumber[i];
            fMassElement = uma2Mev(fMassUma[i]);
            break;
         }
      }
   } else {
      LOG(fatal) << "Element name must be 1 or 2 characters";
   }
};

/*
 * Constructor for an element A and Z
 */
AtTools::MassInformation::MassInformation(int AElement, int ZElement) : fAElement(AElement), fZElement(ZElement)
{
   ReadMassTable();
   for (int i = 0; i < fNeutronNumber.size(); i++) {
      if (fANumber[i] == fAElement && fProtonNumber[i] == fZElement) {
         fElementName = fElementShortName[i];
         fMassElement = uma2Mev(fMassUma[i]);
         break;
      }
   }
};

/*
 * Reads amdc mass table information:
 * n Z A name mass[uma]
 * */
void AtTools::MassInformation::ReadMassTable()
{
   int neutronNumber, protonNumber, massNumber;
   std::string elementShortName;
   double massUma;

   std::ifstream streamer(fMassFile);
   while (streamer >> neutronNumber >> protonNumber >> massNumber >> elementShortName >> massUma) {
      fNeutronNumber.push_back(neutronNumber);
      fProtonNumber.push_back(protonNumber);
      fANumber.push_back(massNumber);
      fElementShortName.push_back(elementShortName);
      fMassUma.push_back(massUma);
   }
   streamer.close();
};

/*
 * Set light particle with A and name
 */
void AtTools::MassInformation::SetLight(int ALight, std::string lightName)
{
   fALight = ALight;
   fLightName = lightName;
   for (int i = 0; i < fNeutronNumber.size(); i++) {
      if (fANumber[i] == fALight && fElementShortName[i] == fLightName) {
         fZLight = fProtonNumber[i];
         fMassLight = uma2Mev(fMassUma[i]);
         break;
      }
   }

   fAHeavy = fABeam + fATarget - fALight;
   fZHeavy = fZBeam + fZTarget - fZLight;
   for (int i = 0; i < fNeutronNumber.size(); i++) {
      if (fANumber[i] == fAHeavy && fProtonNumber[i] == fZHeavy) {
         fHeavyName = fElementShortName[i];
         fMassHeavy = uma2Mev(fMassUma[i]);
         break;
      }
   }
};

/*
 * Set light particle with A and Z
 */
void AtTools::MassInformation::SetLight(int ALight, int ZLight)
{
   fALight = ALight;
   fZLight = ZLight;

   for (int i = 0; i < fNeutronNumber.size(); i++) {
      if (fANumber[i] == fALight && fProtonNumber[i] == fZLight) {
         fLightName = fElementShortName[i];
         fMassLight = uma2Mev(fMassUma[i]);
         break;
      }
   }

   fAHeavy = fABeam + fATarget - fALight;
   fZHeavy = fZBeam + fZTarget - fZLight;
   for (int i = 0; i < fNeutronNumber.size(); i++) {
      if (fANumber[i] == fAHeavy && fProtonNumber[i] == fZHeavy) {
         fHeavyName = fElementShortName[i];
         fMassHeavy = uma2Mev(fMassUma[i]);
         break;
      }
   }
};

/*
 * Set heavy particle with A and name
 */
void AtTools::MassInformation::SetHeavy(int AHeavy, std::string heavyName)
{
   fAHeavy = AHeavy;
   fHeavyName = heavyName;

   for (int i = 0; i < fNeutronNumber.size(); i++) {
      if (fANumber[i] == fAHeavy && fElementShortName[i] == fHeavyName) {
         fZHeavy = fProtonNumber[i];
         fMassHeavy = uma2Mev(fMassUma[i]);
         break;
      }
   }

   fALight = fABeam + fATarget - fAHeavy;
   fZLight = fZBeam + fZTarget - fZHeavy;
   for (int i = 0; i < fNeutronNumber.size(); i++) {
      if (fANumber[i] == fALight && fProtonNumber[i] == fZLight) {
         fLightName = fElementShortName[i];
         fMassLight = uma2Mev(fMassUma[i]);
         break;
      }
   }
}

/*
 * Set heavy with A and Z
 */
void AtTools::MassInformation::SetHeavy(int AHeavy, int ZHeavy)
{
   fAHeavy = AHeavy;
   fZHeavy = ZHeavy;

   for (int i = 0; i < fNeutronNumber.size(); i++) {
      if (fANumber[i] == fAHeavy && fProtonNumber[i] == fZHeavy) {
         fHeavyName = fElementShortName[i];
         fMassHeavy = uma2Mev(fMassUma[i]);
         break;
      }
   }

   fALight = fABeam + fATarget - fAHeavy;
   fZLight = fZBeam + fZTarget - fZHeavy;
   for (int i = 0; i < fNeutronNumber.size(); i++) {
      if (fANumber[i] == fALight && fProtonNumber[i] == fZLight) {
         fLightName = fElementShortName[i];
         fMassLight = uma2Mev(fMassUma[i]);
         break;
      }
   }
}

/*
 * Prints the reaction information
 * */
void AtTools::MassInformation::Print()
{
   std::cout << std::left;
   std::cout << "*--------- Reaction Information ---------*"
             << "\n\n";
   std::cout << std::setw(11) << "" << fABeam << fBeamName << " + " << fATarget << fTargetName << " -> " << fAHeavy
             << fHeavyName << " + " << fALight << fLightName << "\n\n";

   std::cout << std::setw(9) << "" << std::setw(4) << "Name" << std::setw(5) << "" << std::setw(3) << "Z"
             << std::setw(3) << "" << std::setw(3) << "A" << std::setw(2) << ""
             << "Mass[MeV]"
             << "\n";
   std::cout << std::setw(9) << "" << std::setw(4) << "----" << std::setw(4) << "" << std::setw(3) << "---"
             << std::setw(3) << "" << std::setw(3) << "---" << std::setw(3) << ""
             << "---------"
             << "\n";
   std::cout << std::setw(10) << "Beam: " << std::setw(5) << fBeamName << std::setw(2) << "" << std::setw(3) << fZBeam
             << std::setw(3) << "" << std::setw(3) << fABeam << std::setw(4) << "" << fMassBeam << "\n";
   std::cout << std::setw(10) << "Target: " << std::setw(5) << fTargetName << std::setw(2) << "" << std::setw(3)
             << fZTarget << std::setw(3) << "" << std::setw(3) << fATarget << std::setw(4) << "" << fMassTarget << "\n";
   std::cout << std::setw(10) << "Light: " << std::setw(5) << fLightName << std::setw(2) << "" << std::setw(3)
             << fZLight << std::setw(3) << "" << std::setw(3) << fALight << std::setw(4) << "" << fMassLight << "\n";
   std::cout << std::setw(10) << "Heavy: " << std::setw(5) << fHeavyName << std::setw(2) << "" << std::setw(3)
             << fZHeavy << std::setw(3) << "" << std::setw(3) << fAHeavy << std::setw(4) << "" << fMassHeavy << "\n";
   std::cout << "*----------------------------------------*"
             << "\n";
}

/******************************************************************/

/**
 * Setting the reaction information with the mass information given by MassInformation class
 * @param reactionInfo class containing all the components information
 * @param beamEnergy beam Energy in Mev/A
 **/
AtTools::KinematicsCalculation::KinematicsCalculation(MassInformation reactionInfo)
   : fReactionInfo(std::move(reactionInfo))
{

   fwMassBeam = fReactionInfo.GetBeamMass();
   fwMassTarget = fReactionInfo.GetTargetMass();
   fwMassLight = fReactionInfo.GetLightMass();
   fwMassHeavy = fReactionInfo.GetHeavyMass();
}

AtTools::KinematicsCalculation::~KinematicsCalculation()
{
   lightGraph->Delete();
   heavyGraph->Delete();
   relationGraph->Delete();
   lightLabCMGraph->Delete();
}

/**
 * Calculation of the Relativistic Kinematics in cm from 0.1 to 179.9 degrees TODO :: revisar angulos para evitar nan y
 * valores raros
 * TODO:: Revisar toda la cinematica
 * */
void AtTools::KinematicsCalculation::CalculateKinematics()
{
   // Reseting vectors
   fThetaCM.clear();
   fEnergyLight.clear();
   fEnergyHeavy.clear();
   fAngleLABLight.clear();
   fAngleLABHeavy.clear();

   // Initial Values
   double EBeamLAB = fBeamEnergyMeV + fwMassBeam;
   double beamMomentum =
      TMath::Sqrt(pow(EBeamLAB, 2) - pow(fwMassBeam, 2)); // p=sqrt(E² -m²)=sqrt(beamEnergy²+2*beamEnergy*beamMass)

   double ETotal = EBeamLAB + fwMassTarget;
   double ETotalCM = TMath::Sqrt(pow(ETotal, 2) - pow(beamMomentum, 2));

   // Check reaction is possible
   double TCM = ETotalCM - fwMassLight - fwMassHeavy;
   if (TCM < 0.0) {
      LOG(warning) << "Kinetic energy value is: " << TCM;
      fNoSolution = true;
      return;
   }

   // Out particles
   double ELightCM = 0.5 * (ETotalCM + (pow(fwMassLight, 2) - pow(fwMassHeavy, 2)) / ETotalCM);
   double lightMomentumCM = TMath::Sqrt(pow(ELightCM, 2) - pow(fwMassLight, 2));
   double EHeavyCM = 0.5 * (ETotalCM + (pow(fwMassHeavy, 2) - pow(fwMassLight, 2)) / ETotalCM);
   double heavyMomentumCM = TMath::Sqrt(pow(EHeavyCM, 2) - pow(fwMassHeavy, 2));

   // Beta and gamma
   double gammaC = ETotal / ETotalCM;
   double beta = beamMomentum / ETotal;

   // Lorentz transformations to lab
   for (double angle = 0.1; angle <= 179.9; angle = angle + 0.1) { // angle must be in radians for cos and sin
      fThetaCM.push_back(180 - angle);
      double angleRad = angle * TMath::DegToRad();

      // Light particle
      double ELight = gammaC * (ELightCM + beta * lightMomentumCM * cos(angleRad));
      double lightMomentumLabX = gammaC * (lightMomentumCM * cos(angleRad) + beta * ELightCM);
      double lightMomentumLabY = lightMomentumCM * sin(angleRad);
      double lightMomentumLab = TMath::Sqrt(pow(lightMomentumLabX, 2) + pow(lightMomentumLabY, 2));
      double cosThetaLight = (ELight * ETotal - ELightCM * ETotalCM) / (beamMomentum * lightMomentumLab);
      if (cosThetaLight > 1.0) {
         fAngleLABLight.push_back(0.0);
      } else {
         double thetaLight = acos(cosThetaLight);
         fAngleLABLight.push_back(thetaLight * TMath::RadToDeg());
      }
      fEnergyLight.push_back((ELight - fwMassLight));

      // Heavy particle
      double tg_thetaHeavy =
         heavyMomentumCM * sin(angleRad) /
         (gammaC * (heavyMomentumCM * cos(angleRad) +
                    beta * TMath::Sqrt(heavyMomentumCM * heavyMomentumCM + fwMassHeavy * fwMassHeavy)));
      fAngleLABHeavy.push_back(atan(tg_thetaHeavy) * TMath::RadToDeg());
      double EHeavy = ETotal - ELight;
      fEnergyHeavy.push_back((EHeavy - fwMassHeavy));
   }

   // filling tgraphs
   int nPoints = fThetaCM.size();
   lightGraph = new TGraph(nPoints, &fAngleLABLight[0], &fEnergyLight[0]);
   heavyGraph = new TGraph(nPoints, &fAngleLABHeavy[0], &fEnergyHeavy[0]);
   relationGraph = new TGraph(nPoints, &fThetaCM[0], &fAngleLABLight[0]);
   lightLabCMGraph = new TGraph(nPoints, &fAngleLABHeavy[0], &fAngleLABLight[0]);
}

/**
 * Calculates the corrected energy for a given Lab angle and the initial energy. Correction event by event.
 * @param beamCorrectedEnergy
 * @param thetaLab
 * @param originalEnergy
 **/
double AtTools::KinematicsCalculation::CorrectEnergy(double beamCorrectedEnergy, double thetaLab, double originalEnergy)
{
   // Setting Reaction masses
   double wMassBeam = fReactionInfo.GetBeamMass() + fExBeam;
   double wMassTarget = fReactionInfo.GetTargetMass() + fExTarget;
   double wMassLight = fReactionInfo.GetLightMass() + fExLight;
   double wMassHeavy = fReactionInfo.GetHeavyMass() + fExHeavy;

   // Initial Values
   double TBeam = beamCorrectedEnergy * fReactionInfo.GetBeamAZ().first;

   double EBeamLAB = TBeam + wMassBeam;
   double beamMomentumCM = TMath::Sqrt(2 * wMassBeam * TBeam + pow(TBeam, 2));
   double momentumLab = beamMomentumCM * cos(thetaLab * TMath::DegToRad());

   // Members for the second order energy equation
   double C1 = pow(wMassBeam, 2) + pow(wMassTarget, 2) + pow(wMassLight, 2) + 2 * EBeamLAB * wMassTarget;
   double X1 = pow(wMassHeavy, 2) - C1 + 2 * (EBeamLAB + wMassTarget) * wMassLight;

   double a = pow((EBeamLAB + wMassTarget) / (momentumLab), 2) - 1;
   double b = (((EBeamLAB + wMassTarget) * X1) / pow(momentumLab, 2)) - 2 * wMassLight;
   double c = pow(X1 / (2 * momentumLab), 2);

   // Solving equation taking into account both solutions
   double correctedEnergyPositive = (-b + TMath::Sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);
   double correctedEnergyNegative = (-b - TMath::Sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);

   // Selecting the corrected energy
   double diffPositive = TMath::Sqrt(pow(correctedEnergyPositive - originalEnergy, 2));
   double diffNegative = TMath::Sqrt(pow(correctedEnergyNegative - originalEnergy, 2));
   if (diffPositive < diffNegative) {
      return correctedEnergyPositive;
   } else {
      return correctedEnergyNegative;
   }
}

/**
 *
 */
void AtTools::KinematicsCalculation::FitEnergyCorrection(TH2F *hOriginalCorrectedEnergy, bool draw)
{
   fFuncPol1->SetRange(
      hOriginalCorrectedEnergy->GetXaxis()->GetBinLowEdge(hOriginalCorrectedEnergy->FindFirstBinAbove(0, 1)),
      hOriginalCorrectedEnergy->GetXaxis()->GetBinUpEdge(hOriginalCorrectedEnergy->FindLastBinAbove(0, 1)));
   hOriginalCorrectedEnergy->Fit(fFuncPol1, "RNQF");

   if (draw) {
      TCanvas *cEnergyRecalibration = new TCanvas("COriginalvsCorrectedE", "COriginalvsCorrectedE", 1000, 1000);
      cEnergyRecalibration->cd();
      hOriginalCorrectedEnergy->Draw("colz");
      fFuncPol1->Draw("SAME");
      fFuncPol1->SetLineWidth(1);
      fFuncPol1->SetLineColor(kRed);
      hOriginalCorrectedEnergy->SetTitle(" E corrections; E (MeV); E_{cor} (MeV)");
   }
}

/**
 * Recalibrates energy values event by event.
 */
double AtTools::KinematicsCalculation::Recalibration(double originalEnergy)
{
   return fFuncPol1->GetParameter(0) + fFuncPol1->GetParameter(1) * originalEnergy;
}

/**
 *Obtains Ex and thetaCM for a given ELab and thetaLAB
 */
std::pair<double, double>
AtTools::KinematicsCalculation::ExAndThetaCM(double energyLab, double thetaLAB, double TBeamCorrected)
{
   double Ex, thetaCM;
   double massBeam = fReactionInfo.GetBeamMass();
   double TBeam = TBeamCorrected + massBeam;
   double EBeam = massBeam + TBeamCorrected;

   double momentumBeam = TMath::Sqrt(2 * massBeam * TBeamCorrected + pow(TBeamCorrected, 2));
   double momentumLight = TMath::Sqrt(2 * fReactionInfo.GetLightMass() * energyLab + pow(energyLab, 2));
   double C1 = pow(massBeam, 2) + pow(fReactionInfo.GetTargetMass(), 2) + pow(fReactionInfo.GetLightMass(), 2) +
               2 * EBeam * fReactionInfo.GetTargetMass();
   Ex = TMath::Sqrt(C1 - 2 * (EBeam + fReactionInfo.GetTargetMass()) * (fReactionInfo.GetLightMass() + energyLab) +
                    2 * momentumBeam * momentumLight * cos(thetaLAB * TMath::DegToRad())) -
        fReactionInfo.GetHeavyMass();

   double massHeavy = fReactionInfo.GetHeavyMass() + Ex;
   double massLight = fReactionInfo.GetLightMass();
   double ELight = energyLab + massLight;
   double ETotal = EBeam + fReactionInfo.GetTargetMass();
   double ETotalCM = TMath::Sqrt(2 * EBeam * fReactionInfo.GetTargetMass() + pow(fReactionInfo.GetTargetMass(), 2) +
                                 pow(massHeavy, 2));
   double ELightCM =
      0.5 * (ETotalCM + (pow(fReactionInfo.GetLightMass(), 2) - pow(fReactionInfo.GetHeavyMass(), 2)) / ETotalCM);
   double momentumLightCM = TMath::Sqrt(pow(ELightCM, 2) - pow(fReactionInfo.GetLightMass(), 2));

   double gammaC = ETotal / ETotalCM;
   double beta = momentumBeam / ETotal;

   double cosThetaCM = ((ELight / gammaC) - ELightCM) / (beta * momentumLightCM);
   thetaCM = 180. - acos(cosThetaCM) * TMath::RadToDeg();
   return {Ex, thetaCM};
}

void AtTools::KinematicsCalculation::Plot()
{
   // TODO :si no se corre la cinematica el plot saldra vacio, comprobar si se ha corrido para los plots
   if (fThetaCM.empty()) {
      AtTools::KinematicsCalculation::CalculateKinematics();
   }

   lightGraph->GetXaxis()->SetTitle("#theta_{lab} (deg)");
   lightGraph->GetYaxis()->SetTitle("Energy (MeV)");
   heavyGraph->GetXaxis()->SetTitle("#theta_{lab} (deg)");
   heavyGraph->GetYaxis()->SetTitle("Energy (MeV)");
   relationGraph->GetXaxis()->SetTitle("#theta_{CM} (deg)");
   relationGraph->GetYaxis()->SetTitle("#theta_{Lab} (deg)");
   lightLabCMGraph->GetXaxis()->SetTitle("#theta_{Lab} Scatter (deg)");
   lightLabCMGraph->GetYaxis()->SetTitle("#theta_{Lab} Recoil (deg)");

   TCanvas *cKinematics = new TCanvas("KinematicLines", "Kinematic Lines", 1000, 1000);
   cKinematics->Divide(2, 2);

   cKinematics->Draw();
   cKinematics->cd(1);
   lightGraph->DrawClone("AC");

   cKinematics->cd(2);
   heavyGraph->DrawClone("AC");

   cKinematics->cd(3);
   relationGraph->DrawClone("AC");

   cKinematics->cd(4);
   lightLabCMGraph->DrawClone("AC");
}
