/*********************************************************************
 *   Triangles pad plane for the attpc Mapping Class                 *
 *   Author: A. Munoz-Ramos alicia.munoz@rai.usc.es                  *
 *   Last Reviwed: 09/01/2024                                        *
 *                                                                   *
 *********************************************************************/

#ifndef ATTPCTRIANGLESMAP_H
#define ATTPCTRIANGLESMAP_H

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


class AtTpcTrianglesMap : public AtMap {

public:
   AtTpcTrianglesMap(float padSize = 2.); 
   ~AtTpcTrianglesMap();

   void Dump() override;    
   void GeneratePadPlane() override;    
   ROOT::Math::XYPoint CalcPadCenter(Int_t PadRef) override;    
   Int_t BinToPad(Int_t binval) override { return binval - 1; }; //check the return value

protected:
   float planeSize = 140.; // in mm
   int numPads = 0; //total number of pads
   float sizePad; //size of the pad in mm (size is taken as the side of the triangle)
   float heightPad;
   std::vector<int> trianglesPerRow; //triangles per row in one sextant 
   ClassDefOverride(AtTpcTrianglesMap, 1); 
};

#endif

