#include "JetSubstructure/JetSubstructure.h"
using namespace std;

double JetHelperTools::DeltaPhi(double phi1, double phi2)
{
	Float_t diff;
	diff = phi1-phi2;
	while (diff < -TMath::Pi() ) diff += 2*TMath::Pi();
	while (diff >  TMath::Pi() ) diff -= 2*TMath::Pi();

	return diff;
}

double JetHelperTools::DeltaEta(double eta1, double eta2)
{
	return (eta1 - eta2);
}


double JetHelperTools::DeltaR(double phi1, double eta1, double phi2, double eta2)
{
	double deta=DeltaEta(eta1,eta2);
	double dphi=DeltaPhi(phi1, phi2);
	return std::sqrt(deta*deta+dphi*dphi);
}


double JetHelperTools::DeltaPsi (double phi, double psi) {

	double diff;
	diff = fabs(phi - psi);
	while (diff > TMath::Pi()/2. ) diff = TMath::Pi() - diff;
	return fabs(diff);
}


double JetHelperTools::DeltaPhi(const TLorentzVector& v1, const TLorentzVector& v2)
{
	return DeltaPhi(v1.Phi(),v2.Phi());
}

double JetHelperTools::DeltaEta(const TLorentzVector& v1, const TLorentzVector& v2)
{
	return DeltaEta(v1.Eta(),v2.Eta());
}


double JetHelperTools::DeltaR(const TLorentzVector& v1, const TLorentzVector& v2)
{
	double deta=DeltaEta(v1,v2);
	double dphi=DeltaPhi(v1,v2);
	return std::sqrt(deta*deta+dphi*dphi);
}

int JetHelperTools::GetJetFlavour(float jet_phi,float jet_eta, vector<float> *mc_unstable_gen_pt,vector<float> *mc_unstable_gen_eta,vector<float> *mc_unstable_gen_phi,vector<int> *mc_unstable_pdg) {
	int flavour=-1;              // this moment will identify the flavor (just an initialization)
	Float_t ptMax = 0;
	Int_t thePdg = 0;
	for (unsigned int j=0; j<mc_unstable_gen_eta->size(); j++)
	{
		int pdg = fabs( mc_unstable_pdg->at(j) );
		if ((pdg!=1) && (pdg!=2) && (pdg!=3) && (pdg!=4) && (pdg!=5) && (pdg!=21)) continue;
		Float_t deltaR = sqrt( pow(JetHelperTools::DeltaPhi(jet_phi,mc_unstable_gen_phi->at(j)),2) + pow(jet_eta-mc_unstable_gen_eta->at(j),2));
		if (deltaR > 0.4) continue;
		if (mc_unstable_gen_pt->at(j) > ptMax )
		{thePdg = pdg;
			ptMax = mc_unstable_gen_pt->at(j);
		}
	}
	if      ((thePdg>0) && (thePdg<=3)) flavour=1;    // quark jet
	else if (thePdg==21)                flavour=0;    // gluon jet
	else if ((thePdg>3) && (thePdg<=5)) flavour=2;   // heavy flavor jet
	else                                flavour=-1;  // unrecognized

	return flavour; //in units of GeV
}


std::vector<bool> JetHelperTools::GetIsolation(std::vector<double>& jetpT, std::vector<double>& jetEta, std::vector<double>& jetPhi, int iso_R, double iso_pT)
{
	bool isIso = true;
	std::vector<bool> Isolated;
	double pT_cut = iso_pT;

	for (int i = 0; i<jetpT.size(); i++)
	{
		if (iso_pT < 0 ) pT_cut = jetpT.at(i); // -1 is for isolating pT requirement is 100% of jet pT

		for (int j = 0; j < jetpT.size(); j++)
		{
			if (j == i) continue;

			double R = JetHelperTools::DeltaR(jetPhi.at(i),jetEta.at(i),jetPhi.at(j),jetEta.at(j));
			if (R < iso_R && jetpT.at(j) > pT_cut) isIso = false;
		}

		Isolated.push_back(isIso);
	}

	return Isolated;
}
