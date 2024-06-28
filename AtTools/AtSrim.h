/**************************************************************************
 *   Read and manipulate the data of SRIM generated tables                *
 *   from version SRIM-2013.00                                            *
 *   Author: A. Munoz-Ramos alicia.munoz@rai.usc.es                       *
 *   Last Reviwed: 18/6/24                                                *
 *                                                                        *
 **************************************************************************/

#ifndef ATSRIM_H
#define ATSRIM_H

#include <regex>
#include <string>
#include <vector>

namespace AtTools {

class SRIM {
public:
   SRIM(std::string fileName, std::string decimalSeparator = ".");

   // Getters for basic raw data
   std::vector<double> GetEnergyTable() { return fEnergy; };
   std::vector<double> GetRangeTable() { return fRange; }
   std::vector<double> GetStoppingPowerTable() { return fStoppingPower; }
   std::vector<double> GetLongitudinalStragglingTable() { return fLongitudinalStraggling; };
   std::vector<double> GetLateralStragglingTable() { return fLateralStraggling; }

   double GetRange(double energyValue);
   double GetEnergy(double rangeValue);
   double GetEnergyLeft(double reactionPoint, double initialEnergy);
   double GetInitialEnergy(double reactionPoint, double energyLeft);

   bool IsTargetGas() { return fTargetIsGas; };
   double GetDensity() { return fTargetDensity; };

   void Draw();
   void Print();
   /*
    double GetEnergyLoss;

    GetEfromEleftandTrackLength

    GetElossTable->PunchTroughPoint
    GetEfrom dE
 */
private:
   std::string fFileName;
   std::string fDecimalSeparator = ".";

   // Ion information
   std::string fIonName = "Ion";
   int fIonZ = -1;
   int fIonA = -1;

   // Target information
   bool fTargetIsGas = false;
   double fTargetDensity = -1;
   std::vector<std::string> fTargetAtom;
   std::vector<int> fTargetZ;
   std::vector<double> fTargetAtomicPercent;
   std::vector<double> fTargetMassPercent;

   // Stopping power information
   std::string fUnitsStoppingPower = "";
   std::map<std::string, double> fConversionUnitsStoppingPower;

   // Data information
   std::vector<double> fEnergy, fStoppingPower, fRange, fLongitudinalStraggling, fLateralStraggling;

   // Processing file
   void ReadFile();
   bool fIsHeader = true;
   void ReadHeaderInformation(std::string &line);
   void ReadTailInformation(std::string &line);
   double ConvertToDouble(std::string &value, const std::string &unitValue);
   void IonProcessor(std::smatch matchingElements);
   void DensityProcessor(std::smatch matchingElements);
   void CompositionProcessor(std::smatch matchingelements);
};
} // namespace AtTools

#endif // ATSTRIM_H
