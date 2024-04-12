#include "AtTpcTrianglesMap.h"

#include <FairLogger.h>

#include <Math/Point2D.h>
#include <TMath.h>
#include <Rtypes.h>
#include <TH2Poly.h>

#include <boost/multi_array/base.hpp>
#include <boost/multi_array/extent_gen.hpp>
#include <boost/multi_array/multi_array_ref.hpp>
#include <boost/multi_array/subarray.hpp>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

constexpr auto cRED = "\033[1;31m";
constexpr auto cYellow = "\033[1;33m";
constexpr auto cNORMAL = "\033[0m";

using XYPoint = ROOT::Math::XYPoint;


/*****************************************************************
 * Generates the number of triangles per row in one sextant inside
   the circle. 
 
 * The initial number of triangles in the tip is always one and 
   next rows are always incremented by 2:

                                 |n=R-1
   trianglesPerRow = t0*R + 2*SUM| (n)  ; R>0
                                 |n=0
*****************************************************************/
std::vector<int> RasterizeTriangle(int totalNPads){
    std::vector<int> values;

    int triangles = 1;

    while (totalNPads>0 && totalNPads>triangles){
        values.push_back(triangles);
        totalNPads -= triangles;
        triangles += 2;
    }

    return values;
}


/********************************************************
  * Constructor for the ATTPC Map shaped as an hexagon
    of 14 cm (radius) with pads shaped as triangles.


  * First, it calculates a equilateral triangle plane and
    the amount of small triangles that fit inside of it.
    Then it "rasterizes" the triangle to know how many small 
    triangles fit in each row.
 
  @param padSize: size of the pad in mm 
********************************************************/
AtTpcTrianglesMap::AtTpcTrianglesMap(float padSize) : AtMap() {

    sizePad = padSize;
    float padHeight = sizePad*TMath::Sqrt(3)/2.;
    heightPad = padHeight;

    float heightHexagonTriangle = planeSize*TMath::Sqrt(3)/2.;    
    int numPadsTriangles = (planeSize/sizePad)*(heightHexagonTriangle/heightPad); // Number of triangles in a equilateral triangle forming an hexagon inscribed in the padplane

    trianglesPerRow = RasterizeTriangle(numPadsTriangles); // Vector with the triangles in each row in a sextant 

    // Calculating how many pads are in the pad plane
    for (int i = 0; i<trianglesPerRow.size(); i++){
        numPads += trianglesPerRow[i];
    }

    numPads *=6;

    AtPadCoord.resize(boost::extents[numPads][3][2]); // pads [pindex] [point in the triangles] [x:0 or y:1]
    std::fill(AtPadCoord.data(), AtPadCoord.data() + AtPadCoord.num_elements(), 0);
    std::cout << " ATTPC Hexagon Map initialized with " << numPads << " triangles of size" << sizePad << std::endl;
    std::cout << " ATTPC Triangles Pad Coordinates container initialized " << std::endl;
    fNumberPads = numPads;
}


/*********************
 * Default destructor
*********************/
AtTpcTrianglesMap::~AtTpcTrianglesMap() = default;



/**************************************************************
 * Generates a txt file with the pad index and pad coordinates
**************************************************************/
void AtTpcTrianglesMap::Dump() {

    std::ofstream coordmap;
    coordmap.open("coordmap_Hexagon.txt");

    for (index i=0; i < 25; ++i)
    {
     coordmap << i << " ";
     coordmap << AtPadCoord[i][3][0] << " " << AtPadCoord[i][3][1];
     coordmap << std::endl;
    }

    coordmap.close();

}

/*********************************************************************************
 * Generates the pad plane. Fills the AtPadCoord[numPads][3][2] with the
   information of all the triangles as follows:

                                           | a is each corner of the triangle
   AtPadCoord[index of the pad][a][b] where|
                                           | b is 0 for x and 1 for y coordinates
 
*********************************************************************************/
void AtTpcTrianglesMap::GeneratePadPlane()
{

    if (fPadPlane) {
        LOG(error) << "Skipping generation of pad plane, it is already parsed!";
        return;
    }

    std::cout << "ATTPC Hexagon Map: Generating the map geometry of the ATTPC" << std::endl;

  int pad_num = 0;
    for (auto irow = 0; irow < trianglesPerRow.size(); ++irow) {
        for (auto ipad = 0; ipad < trianglesPerRow[irow]; ++ipad) {


         //1 sextant
         AtPadCoord[ipad + pad_num][0][0] = (-irow - 1 + ipad)*0.5*sizePad;
         AtPadCoord[ipad + pad_num][0][1] = (irow + ((ipad + 1)%2))*heightPad;
         AtPadCoord[ipad + pad_num][1][0] = (-irow + ipad)*0.5*sizePad;
         AtPadCoord[ipad + pad_num][1][1] = (irow + (ipad%2))*heightPad;
         AtPadCoord[ipad + pad_num][2][0] = (-irow + 1 + ipad)*0.5*sizePad;
         AtPadCoord[ipad + pad_num][2][1] = (irow + ((ipad + 1)%2))*heightPad;

         //2 sextant
         AtPadCoord[ipad + pad_num + numPads/6][0][0] = (-irow - 1 + ipad)*0.5*sizePad*TMath::Cos(TMath::Pi()/3) - (irow + ((ipad + 1)%2))*heightPad*TMath::Sin(TMath::Pi()/3);
         AtPadCoord[ipad + pad_num + numPads/6][0][1] = (-irow - 1 + ipad)*0.5*sizePad*TMath::Sin(TMath::Pi()/3) + (irow + ((ipad + 1)%2))*heightPad*TMath::Cos(TMath::Pi()/3);
         AtPadCoord[ipad + pad_num + numPads/6][1][0] = (-irow + ipad)*0.5*sizePad*TMath::Cos(TMath::Pi()/3) - (irow + (ipad%2))*heightPad*TMath::Sin(TMath::Pi()/3);
         AtPadCoord[ipad + pad_num + numPads/6][1][1] = (-irow + ipad)*0.5*sizePad*TMath::Sin(TMath::Pi()/3) + (irow + (ipad%2))*heightPad*TMath::Cos(TMath::Pi()/3);
         AtPadCoord[ipad + pad_num + numPads/6][2][0] = (-irow + 1 + ipad)*0.5*sizePad*TMath::Cos(TMath::Pi()/3) - (irow + ((ipad + 1)%2))*heightPad*TMath::Sin(TMath::Pi()/3);
         AtPadCoord[ipad + pad_num + numPads/6][2][1] = (-irow + 1 + ipad)*0.5*sizePad*TMath::Sin(TMath::Pi()/3) + (irow + ((ipad + 1)%2))*heightPad*TMath::Cos(TMath::Pi()/3);

         //3 sextant
         AtPadCoord[ipad + pad_num + 2*numPads/6][0][0] = (-irow - 1 + ipad)*0.5*sizePad*TMath::Cos(TMath::Pi()*2/3) - (irow + ((ipad + 1)%2))*heightPad*TMath::Sin(TMath::Pi()*2/3);
         AtPadCoord[ipad + pad_num + 2*numPads/6][0][1] = (-irow - 1 + ipad)*0.5*sizePad*TMath::Sin(TMath::Pi()*2/3) + (irow + ((ipad + 1)%2))*heightPad*TMath::Cos(TMath::Pi()*2/3);
         AtPadCoord[ipad + pad_num + 2*numPads/6][1][0] = (-irow + ipad)*0.5*sizePad*TMath::Cos(TMath::Pi()*2/3) - (irow + (ipad%2))*heightPad*TMath::Sin(TMath::Pi()*2/3);
         AtPadCoord[ipad + pad_num + 2*numPads/6][1][1] = (-irow + ipad)*0.5*sizePad*TMath::Sin(TMath::Pi()*2/3) + (irow + (ipad%2))*heightPad*TMath::Cos(TMath::Pi()*2/3);
         AtPadCoord[ipad + pad_num + 2*numPads/6][2][0] = (-irow + 1 + ipad)*0.5*sizePad*TMath::Cos(TMath::Pi()*2/3) - (irow + ((ipad + 1)%2))*heightPad*TMath::Sin(TMath::Pi()*2/3);
         AtPadCoord[ipad + pad_num + 2*numPads/6][2][1] = (-irow + 1 + ipad)*0.5*sizePad*TMath::Sin(TMath::Pi()*2/3) + (irow + ((ipad + 1)%2))*heightPad*TMath::Cos(TMath::Pi()*2/3);

         //4 sextant
         AtPadCoord[ipad + pad_num + 3*numPads/6][0][0] = (-irow - 1 + ipad)*0.5*sizePad;
         AtPadCoord[ipad + pad_num + 3*numPads/6][0][1] = -1*((irow + ((ipad + 1)%2))*heightPad);
         AtPadCoord[ipad + pad_num + 3*numPads/6][1][0] = (-irow + ipad)*0.5*sizePad;
         AtPadCoord[ipad + pad_num + 3*numPads/6][1][1] = -1*((irow + ipad%2)*heightPad);
         AtPadCoord[ipad + pad_num + 3*numPads/6][2][0] = (-irow + 1 + ipad)*0.5*sizePad;
         AtPadCoord[ipad + pad_num + 3*numPads/6][2][1] = -1*((irow + (ipad + 1)%2)*heightPad);

         //5 sextant
         AtPadCoord[ipad + pad_num + 4*numPads/6][0][0] = (-irow - 1 + ipad)*0.5*sizePad*TMath::Cos(-TMath::Pi()/3) - (irow + ((ipad + 1)%2))*heightPad*TMath::Sin(-TMath::Pi()/3);
         AtPadCoord[ipad + pad_num + 4*numPads/6][0][1] = (-irow - 1 + ipad)*0.5*sizePad*TMath::Sin(-TMath::Pi()/3) + (irow + ((ipad + 1)%2))*heightPad*TMath::Cos(-TMath::Pi()/3);
         AtPadCoord[ipad + pad_num + 4*numPads/6][1][0] = (-irow + ipad)*0.5*sizePad*TMath::Cos(-TMath::Pi()/3) - (irow + (ipad%2))*heightPad*TMath::Sin(-TMath::Pi()/3);
         AtPadCoord[ipad + pad_num + 4*numPads/6][1][1] = (-irow + ipad)*0.5*sizePad*TMath::Sin(-TMath::Pi()/3) + (irow + (ipad%2))*heightPad*TMath::Cos(-TMath::Pi()/3);
         AtPadCoord[ipad + pad_num + 4*numPads/6][2][0] = (-irow + 1 + ipad)*0.5*sizePad*TMath::Cos(-TMath::Pi()/3) - (irow + ((ipad + 1)%2))*heightPad*TMath::Sin(-TMath::Pi()/3);
         AtPadCoord[ipad + pad_num + 4*numPads/6][2][1] = (-irow + 1 + ipad)*0.5*sizePad*TMath::Sin(-TMath::Pi()/3) + (irow + ((ipad + 1)%2))*heightPad*TMath::Cos(-TMath::Pi()/3);

         //6 sextant
         AtPadCoord[ipad + pad_num + 5*numPads/6][0][0] = (-irow - 1 + ipad)*0.5*sizePad*TMath::Cos(-TMath::Pi()*2/3) - (irow + ((ipad + 1)%2))*heightPad*TMath::Sin(-TMath::Pi()*2/3);
         AtPadCoord[ipad + pad_num + 5*numPads/6][0][1] = (-irow - 1 + ipad)*0.5*sizePad*TMath::Sin(-TMath::Pi()*2/3) + (irow + ((ipad + 1)%2))*heightPad*TMath::Cos(-TMath::Pi()*2/3);
         AtPadCoord[ipad + pad_num + 5*numPads/6][1][0] = (-irow + ipad)*0.5*sizePad*TMath::Cos(-TMath::Pi()*2/3) - (irow + (ipad%2))*heightPad*TMath::Sin(-TMath::Pi()*2/3);
         AtPadCoord[ipad + pad_num + 5*numPads/6][1][1] = (-irow + ipad)*0.5*sizePad*TMath::Sin(-TMath::Pi()*2/3) + (irow + (ipad%2))*heightPad*TMath::Cos(-TMath::Pi()*2/3);
         AtPadCoord[ipad + pad_num + 5*numPads/6][2][0] = (-irow + 1 + ipad)*0.5*sizePad*TMath::Cos(-TMath::Pi()*2/3) - (irow + ((ipad + 1)%2))*heightPad*TMath::Sin(-TMath::Pi()*2/3);
         AtPadCoord[ipad + pad_num + 5*numPads/6][2][1] = (-irow + 1 + ipad)*0.5*sizePad*TMath::Sin(-TMath::Pi()*2/3) + (irow + ((ipad + 1)%2))*heightPad*TMath::Cos(-TMath::Pi()*2/3);

      }

      pad_num += trianglesPerRow[irow];
   }
  // Creating a TH2Poly to plot the pad plane
   kIsParsed = true;
   pad_num = pad_num*6;
   fPadPlane = new TH2Poly(); // NOLINT
   for (auto ipad = 0; ipad < pad_num; ++ipad) {
      Double_t px[] = {AtPadCoord[ipad][0][0], AtPadCoord[ipad][1][0], AtPadCoord[ipad][2][0], AtPadCoord[ipad][0][0]};
      Double_t py[] = {AtPadCoord[ipad][0][1], AtPadCoord[ipad][1][1], AtPadCoord[ipad][2][1], AtPadCoord[ipad][0][1]};
      fPadPlane->AddBin(4, px, py);
   }

   fPadPlane->SetName("Hexagone_Plane");
   fPadPlane->SetTitle("Hexagone_Plane");
   fPadPlane->ChangePartition(500, 500);

}


// Copied from GadgetMap
ROOT::Math::XYPoint  AtTpcTrianglesMap::CalcPadCenter(Int_t PadRef) {
   if (!kIsParsed) {
      LOG(error) << " AtGHexagonMap::CalcPadCenter Error : Pad plane has not been generated or parsed";
      return {-9999, -9999};
   }

   if (PadRef == -1) { // Boost multi_array crashes with a negative index
      LOG(debug) << " AtHexagonMap::CalcPadCenter Error : Pad not found";
      return {-9999, -9999};
   }

   Float_t x = (AtPadCoord[PadRef][0][0] + AtPadCoord[PadRef][1][0]) / 2.0;
   Float_t y = (AtPadCoord[PadRef][1][1] + AtPadCoord[PadRef][2][1]) / 2.0;
   return {x, y};
}

ClassImp(AtTpcTrianglesMap)

