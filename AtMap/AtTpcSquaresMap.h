/*********************************************************************
 *   Squares pad plane for the attpc Mapping Class                   *
 *   Author: A. Munoz-Ramos alicia.munoz@rai.usc.es                  *
 *   Last Reviwed: 09/04/2023                                        *
 *                                                                   *
 *********************************************************************/

#ifndef ATTPCSQUARESMAP_H
#define ATTPCSQUARESMAP_H

#include "AtMap.h"
#include <Math/Point2Dfwd.h>
#include <Rtypes.h>
#include <TString.h>

#include <map>
#include <vector>


#include <unordered_map>


class TBuffer;
class TClass;
class TFile;
class TMemberInspector;


class AtTpcSquaresMap : public AtMap {

public:
   AtTpcSquaresMap(float padSize = 2); // in mm
   ~AtTpcSquaresMap();

   void Dump() override;                                         
   void GeneratePadPlane() override;                             
   ROOT::Math::XYPoint CalcPadCenter(Int_t PadRef) override;     
   Int_t BinToPad(Int_t binval) override { return binval - 1; }; //check the return value

protected:
   float planeSize = 140.; // in mm
   int numPads = 0; //total number of pads
   float sizePad; //size of the pad in mm
   std::vector<int> squaresPerRow; //squares per row in one quadrant adapting the mid-point circle drawing algorithm
   ClassDefOverride(AtTpcSquaresMap, 1); 
};

#endif
