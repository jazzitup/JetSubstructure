#ifndef __JETHELPERTOOLS_H__
#define __JETHELPERTOOLS_H__

#include <string>
#include <vector>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <iostream>
#include <TH3.h>

using namespace std;

namespace JetHelperTools
{
  enum SortType{ SortByPt=0, SortByEt, SortByEta, SortByAbsEta};

  double DeltaPhi(double phi1, double phi2);
  double PhiInPI(double phi);
  double DeltaEta(double eta1, double eta2);
  double DeltaR(double phi1, double eta1, double phi2, double eta2);
  double DeltaPsi(double phi, double psi);

  //in terms of TLorentzVectors for convenience, internally use above versions
  double DeltaPhi(const TLorentzVector& v1, const TLorentzVector& v2);
  double DeltaEta(const TLorentzVector& v1, const TLorentzVector& v2);
  double DeltaR(const TLorentzVector& v1, const TLorentzVector& v2);
  int GetJetFlavour(float jet_phi,float jet_eta, vector<float> *mc_unstable_gen_pt,vector<float> *mc_unstable_gen_eta,vector<float> *mc_unstable_gen_phi,vector<int> *mc_unstable_pdg);
  float GetJetWeight(double pt, double eta, double phi, TH3D* weights);
	std::vector<bool> GetIsolation(std::vector<double>& jetpT, std::vector<double>& jetEta, std::vector<double>& jetPhi, int iso_R, double iso_pT);


}
#endif
