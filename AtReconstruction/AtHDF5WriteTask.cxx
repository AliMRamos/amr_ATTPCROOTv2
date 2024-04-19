#include "AtHDF5WriteTask.h"

#include "AtEvent.h"
#include "AtHit.h"

#include <FairLogger.h>      // for LOG
#include <FairRootManager.h> // for FairRootManager

#include <Math/Point3D.h> // for PositionVector3D
#include <TClonesArray.h>
#include <TObject.h> // for TObject

#include <H5Cpp.h>

#include <array>   // for array
#include <utility> // for move

AtHDF5WriteTask::AtHDF5WriteTask(TString fileName, TString branchName, bool spyral)
   : fOutputFileName(std::move(fileName)), fInputBranchName(std::move(branchName)), fSpyral(spyral), fMinEvent(0),
     fMaxEvent(1000), fICAmplitude(-1.), fICIntegral(-1.), fICCentroid(-1.), fICMultiplicity(-1.)
{
}

InitStatus AtHDF5WriteTask::Init()
{

   FairRootManager *ioMan = FairRootManager::Instance();
   if (ioMan == nullptr) {
      LOG(ERROR) << "Cannot find RootManager!";
      return kERROR;
   }

   fEventArray = dynamic_cast<TClonesArray *>(ioMan->GetObject(fInputBranchName));

   fFile = std::make_unique<H5::H5File>(fOutputFileName, H5F_ACC_TRUNC);

   // ---- HDF5 file format to use Spyral code ----
   if (fSpyral) {
      // ---- Creating group ----
      auto eventGroup = std::make_unique<H5::Group>(fFile->createGroup("/cloud"));
      H5::Attribute groupAttributeMin = eventGroup->createAttribute("min_event", H5::PredType::NATIVE_INT, H5S_SCALAR);
      groupAttributeMin.write(H5::PredType::NATIVE_INT, &fMinEvent);

      // ---- Adding attribute to the group ----
      H5::Attribute groupAttributeMax = eventGroup->createAttribute("max_event", H5::PredType::NATIVE_INT, H5S_SCALAR);
      groupAttributeMax.write(H5::PredType::NATIVE_INT, &fMaxEvent);
   }
   return kSUCCESS;
}

void AtHDF5WriteTask::Exec(Option_t *opt)
{
   auto *event = dynamic_cast<AtEvent *>(fEventArray->At(0));
   if (!event->IsGood())
      return;

   Int_t nHits = event->GetNumHits();
   const auto &traceEv = event->GetMesh();

   // ---- HDF5 file format to use Spyral code
   if (fSpyral) {

      hsize_t hitdim[] = {(hsize_t)nHits, 8}; // dimension of the DataSet
      H5::DataSpace hitSpace(2, hitdim);      // (rank, dimension) of the DataSet
      float hitsInformation[nHits][8];

      // ---- Filling variable array ----
      for (Int_t iHit = 0; iHit < nHits; iHit++) {
         const AtHit &hit = event->GetHit(iHit);

         auto hitPos = hit.GetPosition();
         hitsInformation[iHit][0] = hitPos.X();
         hitsInformation[iHit][1] = hitPos.Y();
         hitsInformation[iHit][2] = hitPos.Z();
         hitsInformation[iHit][3] = hit.GetCharge();
         hitsInformation[iHit][4] = hit.GetTraceIntegral();
         hitsInformation[iHit][5] = 1.0;
         hitsInformation[iHit][6] = hit.GetTimeStamp();
         hitsInformation[iHit][7] = 1.0;
      }

      int eventNum = fUseEventNum ? event->GetEventID() : fEventNum;

      // ---- Creating dataset for each event ----
      H5::DataSet cloudSet =
         fFile->createDataSet(TString::Format("/cloud/cloud_%d", eventNum), H5::PredType::NATIVE_FLOAT, hitSpace);

      cloudSet.write(hitsInformation, H5::PredType::NATIVE_FLOAT);

      // ---- Adding attributtes to the datasets ----
      H5::Attribute amplitudeAttr = cloudSet.createAttribute("ic_amplitude", H5::PredType::NATIVE_FLOAT, H5S_SCALAR);
      amplitudeAttr.write(H5::PredType::NATIVE_FLOAT, &fICAmplitude);

      H5::Attribute integralAttr = cloudSet.createAttribute("ic_integral", H5::PredType::NATIVE_FLOAT, H5S_SCALAR);
      integralAttr.write(H5::PredType::NATIVE_FLOAT, &fICIntegral);

      H5::Attribute centroidAttr = cloudSet.createAttribute("ic_centroid", H5::PredType::NATIVE_FLOAT, H5S_SCALAR);
      centroidAttr.write(H5::PredType::NATIVE_FLOAT, &fICCentroid);

      H5::Attribute multiplicityAttr =
         cloudSet.createAttribute("ic_multiplicity", H5::PredType::NATIVE_FLOAT, H5S_SCALAR);
      multiplicityAttr.write(H5::PredType::NATIVE_FLOAT, &fICMultiplicity);

   } else {

      hsize_t hitdim[] = {(hsize_t)nHits};            // dimension of the hits
      hsize_t tracedim[] = {(hsize_t)traceEv.size()}; // dimension of the trace
      H5::DataSpace hitSpace(1, hitdim);              // rank and dim
      H5::DataSpace traceSpace(1, tracedim);          // rank and dim

      /**
         @TODO At some point we may have to think about how to generalize this so it can also write
         derived types of AtHit.
      **/
      AtHit_t hits[nHits];
      auto hdf5Type = AtHit().GetHDF5Type();

      for (Int_t iHit = 0; iHit < nHits; iHit++) {
         const AtHit &hit = event->GetHit(iHit);

         auto hitPos = hit.GetPosition();
         hits[iHit].x = hitPos.X();
         hits[iHit].y = hitPos.Y();
         hits[iHit].z = hitPos.Z();
         hits[iHit].t = hit.GetTimeStamp();
         hits[iHit].A = hit.GetCharge();
      }

      int eventNum = fUseEventNum ? event->GetEventID() : fEventNum;

      auto eventGroup = std::make_unique<H5::Group>(fFile->createGroup(TString::Format("/Event_[%d]", eventNum)));
      H5::DataSet hitset = fFile->createDataSet(TString::Format("/Event_[%d]/HitArray", eventNum), hdf5Type, hitSpace);
      hitset.write(hits, hdf5Type);

      H5::DataSet traceset =
         fFile->createDataSet(TString::Format("/Event_[%d]/Trace", eventNum), H5::PredType::NATIVE_FLOAT, traceSpace);
      traceset.write(traceEv.data(), H5::PredType::NATIVE_FLOAT);
   }

   ++fEventNum;
}

ClassImp(AtHDF5WriteTask);
