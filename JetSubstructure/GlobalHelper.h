#ifndef GET_GLOBALHELPERS_H
#define GET_GLOBALHELPERS_H
#include <TROOT.h>
#include <TH1D.h>
#include "xAODEventInfo/EventInfo.h"
//CaloCluster include
#include  "xAODCaloEvent/CaloCluster.h" 
#include  "xAODCaloEvent/CaloClusterContainer.h"
#define private public
#include "xAODHIEvent/HIEventShapeAuxContainer.h"
#undef private
#include "xAODHIEvent/HIEventShapeContainer.h" 

using namespace std;

//@CODE_begin
int GetGlobalBin(Int_t centralityScheme, float FCal_Et, bool isMC=false);
int GetCentralityBin(Int_t centralityScheme, float FCal_Et, bool isMC=false);
void SetRejectionHistogram(TH1D* h);
int GetCentralityNBins(Int_t centralityScheme);
float GetEventPlane(const xAOD::CaloClusterContainer *hiclus);
float GetEventPlane(const xAOD::HIEventShapeContainer* calos);
void SetupBinning(Int_t scheme, string variable, Double_t array[1000], Int_t &num);
Float_t GetAveragePsi(Float_t psi1, Float_t psi2);
//@CODE_end

#endif
