#ifndef ATHDF4WRITETASK_H
#define ATHDF4WRITETASK_H

#include <FairTask.h>

#include <Rtypes.h>  // for THashConsistencyHolder, Bool_t, ClassDefOverride
#include <TString.h> // for TString

#include <H5Cpp.h>

#include <memory> // for unique_ptr

class TBuffer;
class TClass;
class TClonesArray;
class TMemberInspector;

class AtHDF5WriteTask : public FairTask {

protected:
   TString fOutputFileName;
   TString fInputBranchName;
   bool fSpyral; // Parameter to create a hdf5 file base on Spyral clouds format

   int fMinEvent; // Attribute of the run indicating the min event Spyral will run from
   int fMaxEvent; // Attribute of the run indicating the max event Spyral will run to

   // Attributes of the event set to -1 by default
   float fICAmplitude;
   float fICIntegral;
   float fICCentroid;
   float fICMultiplicity;

   std::unique_ptr<H5::H5File> fFile{nullptr}; //!
   TClonesArray *fEventArray{nullptr};

   Bool_t fIsPersistence{false};
   /// If true events are indexed by ATTPCROOT event number. If false then use internal index [0-NumEventsInFile).
   Bool_t fUseEventNum{false};

   Int_t fEventNum{0};

public:
   AtHDF5WriteTask(TString fileName, TString branchName = "AtEventH", bool spyral = false);

   void SetEventNumbers(int minValue, int maxValue)
   {
      fMinEvent = minValue;
      fMaxEvent = maxValue;
   }
   void SetICAmplitude(float ICAmplitude) { fICAmplitude = ICAmplitude; };
   void SetICIntegral(float ICIntegral) { fICIntegral = ICIntegral; };
   void SetICCentroid(float ICCentroid) { fICCentroid = ICCentroid; };
   void SetICMultiplicity(float ICMultiplicity) { fICMultiplicity = ICMultiplicity; };

   void SetPersistence(bool val) { fIsPersistence = val; }
   void SetUseEventNum(bool val) { fUseEventNum = val; }
   virtual InitStatus Init() override;
   virtual void Exec(Option_t *opt) override;

   ClassDefOverride(AtHDF5WriteTask, 1);
};

#endif // #ifndef ATHDF4WRITETASK_H
