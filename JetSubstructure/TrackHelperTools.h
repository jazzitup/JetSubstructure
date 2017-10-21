#ifndef __TRACKHELPERTOOLS_H__
#define __TRACKHELPERTOOLS_H__

#include <string>
#include <vector>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <iostream>
#include "BaseClass.h"

using namespace std;

namespace TrackHelperTools
{
  void SetTrackCuts(int& selection,int& NCuts, bool& Dod0Param,int& nSISHits_cut, int& nSIHits_cut, int& nSIholes_cut, int& nPixholes_cut, int& nSCTHits_cut, int& nPixHits_cut, int& nBLHits_cut, int& nIBLHits_cut, double& d0_cut, double& z0sintheta_cut,int &nTRTHits_cut, float &sig_cut);
  //int getType(int barcode, int pdg, int nparent, int nparent_pdg);
  int getTypeTruth(int barcode, int pdg, int status, float charge);
  bool isStableParticle(int barcode, int pdg, int status);
  int getTypeReco(int barcode, int pdg, int status, float charge, float mc_prob, float mc_prob_cut);
  EL::StatusCode SetCutLevel(InDet::InDetTrackSelectionTool *trackSelectionTool, string cutlevel);

}
#endif


