#include "AtTpcSquaresMap.h"

#include <FairLogger.h>

#include <Math/Point2D.h>
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
 * Generates the number of squares per row in one quadrant inside
   the circle. (Based on the mid-point circle drawing algorithm)
 
 * Since the algorithm does not calculate the full quadrant and
   stops when y values are higher than x, a helper vector is used 
   to get all the information in one round.
*****************************************************************/
std::vector<int> RasterizeCircle(int radius){
    std::vector<int> values; // Main vector
    std::vector<int> helper; // Auxiliar vector
    
    int x = radius, y = 0;
    int counter = 0; // This is used to complete the cuadrant
    int P = 1 - radius;
    
    while (x>y){
        values.push_back(x); // Adds the value of squares (=x)
        y++;
        // Checks if it is outside or inside the circle
        if (P<=0){ 
            P = P + 2*y + 1;
            counter +=1; 
        }
        else{
            x--;
            P = P + 2*y - 2*x +1;
            counter += 1;
            helper.push_back(counter);
        }   
        if (x<y)
            break;
    }

    // Reverse the auxiliar vector to add it to all the values
    reverse(helper.begin(), helper.end());
    for (auto i = helper.begin(); i<helper.end(); i++) {
        values.push_back(*i);
    }
    return values;
}


/*******************************************************
  * Constructor for the ATTPC Map of a plane size 
    of 14 cm (radius) with pads shaped as squares. 
  * First calculates a square shaped plane and the 
    the amount of squares that fit inside. Then it
    "rasterizes" the circle using the amount of squares
    that fit in the first row for a quadrant.
 
  @param padSize: size of the pad in mm
*******************************************************/
AtTpcSquaresMap::AtTpcSquaresMap(float padSize) : AtMap() {

    sizePad = padSize;
    int numPadsSquare = (planeSize/sizePad)*(planeSize/sizePad); // Number of squares in a quadrant in a perfect square
    int padsInRow1 = planeSize/sizePad; // How many squares in the first row = longest row
    squaresPerRow = RasterizeCircle(padsInRow1); // Vector with the squares in each row in a quadrant 
 
    // Calculating how many pads are in the pad plane
    for (auto i = 0; i<squaresPerRow.size(); i++){
        numPads += squaresPerRow[i];
    }
    numPads = numPads*4;

    AtPadCoord.resize(boost::extents[numPads][4][2]); // pads [pindex] [point in the squares] [x:0 or y:1]
    std::fill(AtPadCoord.data(), AtPadCoord.data() + AtPadCoord.num_elements(), 0);
    std::cout << " ATTPC Squares Map initialized with " << numPads << " squares of size" << sizePad << std::endl;
    std::cout << " ATTPC Squardes Pad Coordinates container initialized " << std::endl;
    fNumberPads = numPads;
}


/*********************
 * Default destructor
*********************/
AtTpcSquaresMap::~AtTpcSquaresMap() = default;


/**************************************************************
 * Generates a txt file with the pad index and pad coordinates
**************************************************************/
void AtTpcSquaresMap::Dump() {

    std::ofstream coordmap;
    coordmap.open("coordmap_Squares.txt");
    
    for (index i=0; i < 25; ++i)
    {
     coordmap << i << " ";
     coordmap << AtPadCoord[i][3][0] << " " << AtPadCoord[i][3][1];
     coordmap << std::endl;
    }

    coordmap.close();

}



/*********************************************************************************
 * Generates the pad plane. Fills the AtPadCoord[numPads][4][2] with the
   information of all the squares as follows:

                                           | a is each corner of the square
   AtPadCoord[index of the pad][a][b] where|
                                           | b is 0 for x and 1 for y coordinates
 
*********************************************************************************/
void AtTpcSquaresMap::GeneratePadPlane()
{

    if (fPadPlane) {
        LOG(error) << "Skipping generation of pad plane, it is already parsed!";
        return;
    }

    std::cout << "ATTPC Squares Map: Generating the map geometry of the ATTPC" << std::endl;
    
    int pad_num = 0;
    for (auto irow = 0; irow < squaresPerRow.size(); ++irow) { 
        for (auto ipad = 0; ipad < squaresPerRow[irow]; ++ipad) {
         //First quadrant
         AtPadCoord[ipad + pad_num][0][0] = ipad * sizePad;
         AtPadCoord[ipad + pad_num][0][1] = -irow * sizePad;
         AtPadCoord[ipad + pad_num][1][0] = ipad * sizePad + sizePad;
         AtPadCoord[ipad + pad_num][1][1] = -irow * sizePad;
         AtPadCoord[ipad + pad_num][2][0] = ipad * sizePad + sizePad;
         AtPadCoord[ipad + pad_num][2][1] = -sizePad - irow * sizePad;
         AtPadCoord[ipad + pad_num][3][0] = ipad * sizePad;
         AtPadCoord[ipad + pad_num][3][1] = -sizePad - irow * sizePad;

         //Second quadrant
         AtPadCoord[ipad + pad_num + numPads/4][0][0] = -1*ipad * sizePad;
         AtPadCoord[ipad + pad_num + numPads/4][0][1] = -irow * sizePad;
         AtPadCoord[ipad + pad_num + numPads/4][1][0] = -1*(ipad * sizePad + sizePad);
         AtPadCoord[ipad + pad_num + numPads/4][1][1] = -irow * sizePad;
         AtPadCoord[ipad + pad_num + numPads/4][2][0] = -1*(ipad * sizePad + sizePad);
         AtPadCoord[ipad + pad_num + numPads/4][2][1] = -sizePad - irow * sizePad;
         AtPadCoord[ipad + pad_num + numPads/4][3][0] = -1*(ipad * sizePad);
         AtPadCoord[ipad + pad_num + numPads/4][3][1] = -sizePad - irow * sizePad;

         //Third quadrant
         AtPadCoord[ipad + pad_num + 2*numPads/4][0][0] = -1*(ipad * sizePad);
         AtPadCoord[ipad + pad_num + 2*numPads/4][0][1] = -1*(-irow * sizePad);
         AtPadCoord[ipad + pad_num + 2*numPads/4][1][0] = -1*(ipad * sizePad + sizePad);
         AtPadCoord[ipad + pad_num + 2*numPads/4][1][1] = -1*(-irow * sizePad);
         AtPadCoord[ipad + pad_num + 2*numPads/4][2][0] = -1*(ipad * sizePad + sizePad);
         AtPadCoord[ipad + pad_num + 2*numPads/4][2][1] = -1*(-sizePad - irow * sizePad);
         AtPadCoord[ipad + pad_num + 2*numPads/4][3][0] = -1*(ipad * sizePad);
         AtPadCoord[ipad + pad_num + 2*numPads/4][3][1] = -1*(-sizePad - irow * sizePad);

         //Fourth quadrant
         AtPadCoord[ipad + pad_num + 3*numPads/4][0][0] = ipad * sizePad;
         AtPadCoord[ipad + pad_num + 3*numPads/4][0][1] = -1*(-irow * sizePad);
         AtPadCoord[ipad + pad_num + 3*numPads/4][1][0] = ipad * sizePad + sizePad;
         AtPadCoord[ipad + pad_num + 3*numPads/4][1][1] = -1*(-irow * sizePad);
         AtPadCoord[ipad + pad_num + 3*numPads/4][2][0] = ipad * sizePad + sizePad;
         AtPadCoord[ipad + pad_num + 3*numPads/4][2][1] = -1*(-sizePad - irow * sizePad);
         AtPadCoord[ipad + pad_num + 3*numPads/4][3][0] = ipad * sizePad;
         AtPadCoord[ipad + pad_num + 3*numPads/4][3][1] = -1*(-sizePad - irow * sizePad);

      }

      pad_num += squaresPerRow[irow];
   }

   // Creating a TH2Poly to plot the pad plane
   kIsParsed = true;
   pad_num = pad_num*4;
   fPadPlane = new TH2Poly(); // NOLINT
   for (auto ipad = 0; ipad < pad_num; ++ipad) {
      Double_t px[] = {AtPadCoord[ipad][0][0], AtPadCoord[ipad][1][0], AtPadCoord[ipad][2][0], AtPadCoord[ipad][3][0],
                       AtPadCoord[ipad][0][0]};
      Double_t py[] = {AtPadCoord[ipad][0][1], AtPadCoord[ipad][1][1], AtPadCoord[ipad][2][1], AtPadCoord[ipad][3][1],
                       AtPadCoord[ipad][0][1]};
      fPadPlane->AddBin(5, px, py);
   }

   fPadPlane->SetName("Squares_Plane");
   fPadPlane->SetTitle("Squares_Plane");
   fPadPlane->ChangePartition(500, 500);

}


// Copied from GadgetMap
ROOT::Math::XYPoint  AtTpcSquaresMap::CalcPadCenter(Int_t PadRef) {
   if (!kIsParsed) {
      LOG(error) << " AtSquaresMap::CalcPadCenter Error : Pad plane has not been generated or parsed";
      return {-9999, -9999};
   }   

   if (PadRef == -1) { // Boost multi_array crashes with a negative index
      LOG(debug) << " AtSquaresMap::CalcPadCenter Error : Pad not found";
      return {-9999, -9999};
   }   

   Float_t x = (AtPadCoord[PadRef][0][0] + AtPadCoord[PadRef][1][0]) / 2.0;
   Float_t y = (AtPadCoord[PadRef][1][1] + AtPadCoord[PadRef][2][1]) / 2.0;
   return {x, y}; 
}

ClassImp(AtTpcSquaresMap)

