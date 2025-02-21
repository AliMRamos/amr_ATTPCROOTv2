#ifndef AtGAS_H
#define AtGAS_H

#include <Rtypes.h>
// ROOT classes
#include <TString.h>

class TBuffer;
class TClass;
class TMemberInspector;

class AtGas {
public:
   // Constructor and Destructor
   AtGas(TString);
   ~AtGas() = default;

   // Why do we speefiy a assignment operator with a non-standard signature?
   void operator=(const AtGas &GasToCopy);

   // Getter
   Double_t GetEIonize();
   Double_t GetDriftVelocity();
   Double_t GetCoefAttachment();
   Double_t GetCoefDiffusionLong();
   Double_t GetCoefDiffusionTrans();
   Int_t GetGain();
   UInt_t GetRandomCS();

private:
   TString fGasFileName;
   void InitializeParameters();

   Double_t fEIonize{};            //!< effective ionization energy [eV]
   Double_t fDriftVelocity{};      //!< drift velocity [cm/ns]
   Double_t fCoefAttachment{};     //!< attachment coefficient
   Double_t fCoefDiffusionLong{};  //!< longitudinal diffusion coefficient
   Double_t fCoefDiffusionTrans{}; //!< transversal diffusion coefficient
   Double_t fGain{};               //!< gain factor from wire plane

   ClassDef(AtGas, 1)
};

#endif
