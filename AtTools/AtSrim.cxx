#include "AtSrim.h"

#include <FairLogger.h>

#include <TGraph.h>
#include <TSpline.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <regex>
#include <string>

/**
 * Constructor to read a SRIM generated table with no modification needed
 * @param fileName string with the name of the SRIM file to read
 * @param decimalSeparator string of the decimal Separator
 **/
AtTools::SRIM::SRIM(std::string fileName, std::string decimalSeparator)
   : fDecimalSeparator(std::move(decimalSeparator)), fFileName(std::move(fileName))
{
   ReadFile();
};

/**
 * Converts string format data to double format and changes all units to MeV/mm.
 * Supported units to convert from:
 * Energy: eV, keV, MeV, GeV
 * Lenght: A, um, mm, cm, m, km
 * @param value
 * @param unitValue
 **/
double AtTools::SRIM::ConvertToDouble(std::string &value, const std::string &unitValue)
{
   double data;
   std::string newValue = value;
   std::map<std::string, double> conversion{{"eV", 1E-6}, {"keV", 1E-3}, {"MeV", 1}, {"GeV", 1E3},
                                            {"A", 1E-7},  {"um", 1E-3},  {"mm", 1},  {"cm", 1E1},
                                            {"m", 1E3},   {"km", 1E6},   {"None", 1}};
   // Check if decimal is a point or a comma and replace it
   if (fDecimalSeparator != ".") {
      newValue = std::regex_replace(value, std::regex(fDecimalSeparator), ".");
   }

   // Convert string to float
   try {
      data = std::stod(newValue);
   } catch (std::exception &e) {
      LOG(fatal) << "SRIM::ConvertToDouble() could not convert " << value << " to double";
   }
   try {
      data *= conversion.at(unitValue);
   } catch (std::exception &e) {
      LOG(fatal) << "SRIM::ConvertToDouble() could not convert the follow unit: " << unitValue;
   }
   return data;
};

/**
 * Extracts the information of the ion from the header
 * @param matchingElements matching string found in the header
 **/
void AtTools::SRIM::IonProcessor(std::smatch matchingElements)
{
   if (matchingElements.size() != 4)
      LOG(WARNING) << "Ion could not be read properly!";

   fIonName = matchingElements[1].str();
   fIonZ = stoi(matchingElements[2].str());
   fIonA = stoi(matchingElements[3].str());
}

/**
 * Extracts the target density from the header
 * @param matchingElements matching string found in the header
 **/
void AtTools::SRIM::DensityProcessor(std::smatch matchingElements)
{
   if (matchingElements.size() != 2)
      LOG(WARNING) << "Density could not be read properly!";
   std::ssub_match densityMatch = matchingElements[1];
   std::string densityStr = densityMatch.str();
   fTargetDensity = ConvertToDouble(densityStr, "None");
}

/**
 * Extracts the target composition from the header
 * @param matchingElements matching string found in the header
 **/
void AtTools::SRIM::CompositionProcessor(std::smatch matchingElements)
{
   if (matchingElements.size() != 5)
      LOG(WARNING) << "Composition could not be read properly!";
   fTargetAtom.push_back(matchingElements[1]);

   std::ssub_match ZMatch = matchingElements[2];
   std::string ZStr = ZMatch.str();
   fTargetZ.push_back(std::stoi(ZStr));

   std::ssub_match atomicMatch = matchingElements[3];
   std::string atomicStr = atomicMatch.str();
   fTargetAtomicPercent.push_back(ConvertToDouble(atomicStr, "None"));

   std::ssub_match massMatch = matchingElements[4];
   std::string massStr = massMatch.str();
   fTargetMassPercent.push_back(ConvertToDouble(massStr, "None"));
}

/**
 * Reads SRIM header, obtaining the following information:
 * @param line string of each line of the srim file
 **/

void AtTools::SRIM::ReadHeaderInformation(std::string &line)
{
   std::regex ionRegex(R"(^Ion(?=\s*=)\s*=\s*(\w+)\s*\[(\d+)\].*=\s*(\d+))");
   std::regex densityRegex(R"(^Target (?=Density)\w*\s*=\s*(\d+.{1}\d+E[+|-]\d+))");
   std::regex gasTarget(R"(^Target\s*(?=is))");
   std::regex stoppingUnits(R"(^Stopping Units\s*=\s*(.*))");
   std::regex compositionTarget(R"((^\w{1,2})\s*(\d+)\s*(\d+.{1}\d+)\s*(\d+.{1}\d+))");

   std::smatch matchingElements = {};

   // Removing leading and trailing spaces in the lines
   std::string newLine = std::regex_replace(std::move(line), std::regex(R"(^\s+|\s+$)"), "");

   // Finding Ion
   if (std::regex_search(newLine, matchingElements, std::regex(R"(^-{6,})"))) {
      fIsHeader = false;
   } else if (std::regex_search(newLine, matchingElements, ionRegex)) {
      IonProcessor(matchingElements);
   } else if (std::regex_search(newLine, matchingElements, densityRegex)) {
      DensityProcessor(matchingElements);
   } else if (std::regex_search(newLine, matchingElements, gasTarget)) {
      fTargetIsGas = true;
   } else if (std::regex_search(newLine, matchingElements, stoppingUnits)) {
      fUnitsStoppingPower = matchingElements[1].str();
   } else if (std::regex_search(newLine, matchingElements, compositionTarget)) {
      CompositionProcessor(matchingElements);
   }
};

/**
 * Reads SRIM tail that contains different stopping power units for conversion
 * @param line string of each line of the srim file
 **/
void AtTools::SRIM::ReadTailInformation(std::string &line)
{
   std::regex stoppingUnitsRegex(R"(^(\d+[,|.]\d+E[+|-]\d+)\s*(.*))");

   std::smatch matchingElements = {};

   std::string newLine = std::regex_replace(std::move(line), std::regex(R"(^\s+|\s+$)"), "");
   if (std::regex_search(newLine, matchingElements, stoppingUnitsRegex)) {
      std::ssub_match unit = matchingElements[1];
      std::string unitStr = unit.str();
      fConversionUnitsStoppingPower.insert({matchingElements[2].str(), ConvertToDouble(unitStr, "None")});
   }
}

/**
 * Reads a SRIM file with header
 **/
void AtTools::SRIM::ReadFile()
{
   std::ifstream streamer(fFileName);
   if (!streamer)
      LOG(fatal) << "Failed to open SRIM file " << fFileName;

   // Read lines
   std::string line{};
   LOG(info) << "Reading SRIM file: " << fFileName;
   LOG(info) << "Default separator is dot(.), make sure to specify other separator if it is the case. Faulty "
                "conversion will happen otherwise.";

   while (std::getline(streamer, line)) {

      // Reading header
      if (fIsHeader) {
         ReadHeaderInformation(line);
         continue;
      }

      // Reading data
      std::string energy, unitsEnergy, electronicStopoingPower, nuclearStoppingPower, range, unitsRange,
         longitudinalStraggling, unitsLongitudinalStraggling, lateralStraggling, unitsLateralStraggling;
      std::istringstream lineStreamer{line};

      if (lineStreamer >> energy >> unitsEnergy >> electronicStopoingPower >> nuclearStoppingPower >> range >>
          unitsRange >> longitudinalStraggling >> unitsLongitudinalStraggling >> lateralStraggling >>
          unitsLateralStraggling) {

         // Converting data to units MeV/mm
         fEnergy.push_back(ConvertToDouble(energy, unitsEnergy));
         fStoppingPower.push_back(ConvertToDouble(electronicStopoingPower, "None") +
                                  ConvertToDouble(nuclearStoppingPower, "None"));
         fRange.push_back(ConvertToDouble(range, unitsRange));
         fLongitudinalStraggling.push_back(ConvertToDouble(longitudinalStraggling, unitsLongitudinalStraggling));
         fLateralStraggling.push_back(ConvertToDouble(lateralStraggling, unitsLateralStraggling));

         continue;
      }

      // Reading unit conversion values
      // TO-DO: CONVERTIR Stopping power siempre a Mev/mm
      ReadTailInformation(line);
   }
   streamer.close();
   LOG(info) << "File closed";
   // Finish reading lines
};

/**
 * TSpline3 calculates the correspondent range value for the given energy value in MeV
 * @param energyValue energy in MeV for the evaluation
 **/
double AtTools::SRIM::GetRange(double energyValue)
{
   if (energyValue < fEnergy.front() || energyValue > fEnergy.back()) {
      LOG(warning) << "Value out of energy range. Calculated value by TSpline3::Eval(x) is not an interpolation.";
   }
   std::unique_ptr<TSpline3> energy2RangeSpline =
      std::make_unique<TSpline3>("energy2Range", &fEnergy[0], &fRange[0], fEnergy.size(), "b2e2", 0., 0.);
   return energy2RangeSpline->Eval(energyValue);
};

/**
 * TSpline3 calculates the correspondent energy value for the given range value
 * @param rangeValue range in mm for the evaluation
 * */
double AtTools::SRIM::GetEnergy(double rangeValue)
{
   if (rangeValue < fRange.front() || rangeValue > fRange.back()) {
      std::cout << fRange[-1] << "\n";
      LOG(warning) << "Value out of range range. Calculated value by TSpline3::Eval(x) is not an interpolation.";
   }
   std::unique_ptr<TSpline3> range2EnergySpline =
      std::make_unique<TSpline3>("range2Energy", &fRange[0], &fEnergy[0], fRange.size(), "b2e2", 0., 0.);
   return range2EnergySpline->Eval(rangeValue);
};

/**
 * Energy left for a given initial energy and distance
 * @param reactionPoint distance to calculate the energy left
 * @param initialEnergy initial energy in MeV
 * */
double AtTools::SRIM::GetEnergyLeft(double reactionPoint, double initialEnergy)
{
   double leftDistance, initialDistance;
   double leftEnergy;

   initialDistance = GetRange(initialEnergy);
   leftDistance = initialDistance - reactionPoint;
   leftEnergy = GetEnergy(leftDistance);
   return leftEnergy / fIonA;
}
/**
 * Reconstructs the energy value after a energy loss in a given point
 * @param reactionPoint
 * @param energyLeft
 */
double AtTools::SRIM::GetInitialEnergy(double reactionPoint, double energyLeft)
{
   double rangeLeft, rangeTotal;
   double initialEnergy;

   rangeLeft = GetRange(energyLeft);
   rangeTotal = rangeLeft + reactionPoint;

   initialEnergy = GetEnergy(rangeTotal);
   return initialEnergy / fIonA;
}

void AtTools::SRIM::Draw()
{
   std::vector<double> result;
   for (int i = 0; i < fEnergy.size(); ++i) {
      // result.push_back(std::abs(fEnergy[i]-fStoppingPower))
      std::cout << fStoppingPower[i] << "\n";
   }
   auto gr = new TGraph(fEnergy.size(), &fRange[0], &fStoppingPower[0]);
   gr->Draw("AL");
};
