#ifndef JetCorrector_H
#define JetCorrector_H

#include "TFile.h"
#include "TH2F.h"
#include "TF1.h"
#include "JetSubstructure/JetHelperTools.h"
#include <fstream>
#include <TSystem.h>
#include <TAxis.h>

using namespace std;
using namespace JetHelperTools;

class JetCorrector
{
private:
	// meaning of these quantities is explained in the constructor where the default values are assigned
	TFile * weight_file;

	TH3D *event_weight_histo;
	TH1D *cent_weight_histo;
	vector<double> range_lo;
	vector<double> range_hi;
	vector<int> centiles;

	TFile * _f_JER;
	TH1F* _h_JER[8][8];
	Int_t _nEta_JER;
	Int_t _nCent_JER;
	Int_t _nCent_reweighting;

	TFile * _f_reweighting;
	TAxis * jet_pt_binning;
	TF1 *jet_spectra_weight[8][8];
	TH1D *FF_weight[8][8][20];
	TH1D *CHPS_weight[8][8][20];
	int etabin;

public:
	int nJetYBins;
	float JERcut;
	std::string _weight_file;
	std::string _centrality_weight_file;
	float min_jet_pt;
	float max_jet_pt;

	JetCorrector()
	{

		nJetYBins = 8;			// Default number of y bins
		JERcut = 999; //Cut in sigma of JER
					  //Track-to-jet balance
		_nEta_JER = 8;			// number of eta bins
		_nCent_JER = 8;			// number of centrality bins
		TString xfn = gSystem->GetFromPipe("echo $ROOTCOREBIN");
		_f_JER = new TFile(xfn + "/../JetSubstructure/data/jet_perf_histos.root","read");
		//Shifted by unity in Laura's file
		for(int i=0;i<_nEta_JER;i++){
			for(int j=0;j<_nCent_JER;j++){
				_h_JER[i][j]=(TH1F*)_f_JER->Get(Form("h1_JER_pT_y%i_c%i",i,j+1));
			}
		}

		//Reweighting
		_nCent_reweighting = 8;			// number of centrality bins
		_f_reweighting = new TFile(xfn + "/../JetSubstructure/data/spectra_weights_PbPb.root","read");
		jet_pt_binning = (TAxis*) ((TH1F*)_f_reweighting->Get("h_reco_jet_spectrum_7_cent0_system_0_PbPb"))->GetXaxis();
		//No inclusive in eta
		//for(int etabin=0;etabin<_nEta_JER;etabin++){
		int etabin=7;
		for(int j=0;j<_nCent_reweighting;j++){
			//Spectra weights
			jet_spectra_weight[etabin][j]=(TF1*)_f_reweighting->Get(Form("jet_spectra_weight_%i_cent%i_PbPb",etabin,j));
			for(int k=1;k<=jet_pt_binning->GetNbins();k++){
				//FF weights
				FF_weight[etabin][j][k] = (TH1D*)_f_reweighting->Get(Form("ff_weight_%i_cent%i_system_0_PbPb_jet_pt%i",etabin,j,k));
				CHPS_weight[etabin][j][k] = (TH1D*)_f_reweighting->Get(Form("CHPS_weight_%i_cent%i_system_0_PbPb_jet_pt%i",etabin,j,k));
			}
		}
		//}

		//Event weights
		_weight_file="Powheg.reweight.root";
		_centrality_weight_file="MB_FCal_Normalization.txt";

		TString event_weight_file = xfn + "/../JetSubstructure/data/"+ _weight_file;
		TString fcal_weight_file = xfn +"/../JetSubstructure/data/"+ _centrality_weight_file;

		weight_file = new TFile(event_weight_file.Data());
		event_weight_histo = (TH3D*)weight_file->Get("h3_pT_y_phi_rw");
		cent_weight_histo = (TH1D*)weight_file->Get("h1_cent_rw");

		std::ifstream ifs (fcal_weight_file.Data(), std::ifstream::in);
		if(!ifs)
		{
			cout << endl << "Failed to open file " << fcal_weight_file.Data();
		}

		int tmp1; double tmp2, tmp3;
		while (ifs >> tmp1 >> tmp2 >> tmp3)
		{
			range_lo.push_back(tmp2);
			range_hi.push_back(tmp3);
			centiles.push_back(tmp1);
			if(!ifs.good())break;
		}


	}// end of constructor

	float GetFCalWeight(float FCalEt);
	float GetJetWeight(double pt, double eta, double phi);
	int GetJetYBin(float y);
	Int_t GetJetEtaBin(float eta);
	bool MCJetJERClean(float truth_jet_pt,float reco_jet_pt, float reco_jet_eta, int cent);
	float GetJER(float reco_jet_pt, float reco_jet_eta, int cent);
	float GetJetReweightingFactor(double pt, double eta, int cent);
	float GetFFReweightingFactor(double z, double jet_pt, double jet_eta, int cent);
	float GetCHPSReweightingFactor(double pt, double jet_pt, double jet_eta, int cent);
	~JetCorrector() {}
};


#endif
