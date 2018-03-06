#ifndef TrackCorrector_H
#define TrackCorrector_H

#include "TFile.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "JetSubstructure/JetHelperTools.h"
#include <TSystem.h>

using namespace std;
using namespace JetHelperTools;

class TrackCorrector
{
private:
	// meaning of these quantities is explained in the constructor where the default values are assigned
	TFile * _f_JER;
	TH1F* _h_JER[8][8];
	Int_t _nEta_JER;
	Int_t _nCent_JER;

	TFile * _f_TMR;
	TF1* _tf1_TMR[8];

	TFile * _f_eff;
	TF1 *_tf1_eff[6][10];
	Int_t _nEta_eff;
	Int_t _nEta_eff_jetcorr;

	TFile * _f_eff_jetpt;
	TFile * _f_eff_trketa;

	TGraphAsymmErrors *_g_eff_jept[6][10][20];
//	TGraphAsymmErrors *_g_eff_trketa[25][10];
	TH1 *_h_eff_trketa[25][10];
	TAxis *new_trk_eta_binning;

	TAxis *jetpt_bins;

	TFile * _f_eff_eta;
	TH2D *_h_eff_eta[10];
	TAxis *trkpt_bins;

	TAxis *jety_bins;

	TFile* f_sagitta;
	TH2 * h_sagitta;
	
	double _max_eff_jet_pt;

public:
	string _cutlevel;
	float drmax;
	int _nCent;
	float _max_eff_pt;
	bool _eff_jety;
	int ndRBins;
	double dRrange[50];
	
	double trkpTThreshold;

	//For mutliplicity cuts
	int nMultThresholds;
	float MultThresholds[10];
	
	TrackCorrector(string cutlevel, int nCentBins, bool eff_jety):
	_cutlevel(cutlevel),
	_nCent(nCentBins),
	_eff_jety(eff_jety)
	{
		drmax = 1.0;
		TString xfn = gSystem->GetFromPipe("echo $ROOTCOREBIN");

		//Efficiency correction
		_nEta_eff_jetcorr=5;
		_nEta_eff=5;
		_max_eff_pt=350; //Maximum pt from fit, than constant
		_f_eff = new TFile(xfn + "/../pPbFragmentation/data/mc_eff_fits_" + _cutlevel + ".root","read");

		printf("Loading weights and efficiencies...\n");
		cout << "Using: " << xfn + "/../pPbFragmentation/data/mc_eff_fits_" + _cutlevel + ".root" << endl;
		//Fits are only use in "track-eta" based correction and for UE
		for(int cent=0;cent<_nCent;cent++){
			for(int e=0;e<_nEta_eff;e++){
				_tf1_eff[e][cent]=(TF1*)_f_eff->Get(Form("fit_eff_eta%i_cent%i",e,cent));
			}
		}

		//jet pt correction to the efficiency
		string main_eff_file = "mc_eff_jet_pt_corr";
		if (_eff_jety) main_eff_file = "mc_eff_fits_pt_exclusive";
		_f_eff_jetpt = new TFile(xfn + "/../pPbFragmentation/data/" + main_eff_file + "_" + _cutlevel + ".root","read");
		_f_eff_trketa = new TFile(xfn + "/../pPbFragmentation/data/mc_eff_trketa_jetptinc_" + _cutlevel + ".root","read");

		cout << "Using: " << _f_eff_jetpt->GetName() << endl;

		jetpt_bins=(TAxis*)_f_eff_jetpt->Get("jet_pt_binning");
		_max_eff_jet_pt = jetpt_bins->GetBinUpEdge(jetpt_bins->GetNbins());
		jety_bins = 0;
		if (_eff_jety) { 
			jety_bins = (TAxis*) _f_eff_jetpt->Get("jet_y_binning"); 
			_nEta_eff_jetcorr = jety_bins->GetNbins();
		}
		for(int cent=0;cent<_nCent;cent++){
			for(int e=0;e<_nEta_eff_jetcorr;e++){
				for(int jetb=1;jetb<=jetpt_bins->GetNbins();jetb++){
					if (!_eff_jety) _g_eff_jept[e][cent][jetb]=(TGraphAsymmErrors*)_f_eff_jetpt->Get(Form("ratio_eff_pt_eta%i_cent%i_pt%i",e,cent,jetb));
					else _g_eff_jept[e][cent][jetb]=(TGraphAsymmErrors*)_f_eff_jetpt->Get(Form("graph_eff_eta%i_cent%i_pt%i",e,cent,jetb));
				}
			}
		}

		//trk pt, trk eta based correction
		cout << "Using: " << _f_eff_trketa->GetName() << endl;
		new_trk_eta_binning = (TAxis*)_f_eff_trketa->Get("new_trk_eta_binning");
		int n_trketa_bins = new_trk_eta_binning->GetNbins();
		for(int cent=0;cent<_nCent;cent++)
		{
			for(int e=0;e<n_trketa_bins;e++)
			{
//				_g_eff_trketa[e][cent]=(TGraphAsymmErrors*)_f_eff_trketa->Get(Form("graph_eff_eta%i_cent%i",e,cent));
				_h_eff_trketa[e][cent]=(TH1*)_f_eff_trketa->Get(Form("hist_eff_eta%i_cent%i",e,cent));
			}
		}
		//Residual eta correction to the efficiency
		_f_eff_eta = new TFile(xfn + "/../pPbFragmentation/data/eta_corr_" + _cutlevel + ".root","read");
		//trkpt_bins=(TAxis*)_f_eff_eta->Get("trk_pt_binning");
		for(int cent=0;cent<_nCent;cent++){
			_h_eff_eta[cent]=(TH2D*)_f_eff_eta->Get(Form("eta_corr_cent%i",cent));
		}

		//Track-to-jet balance
		//JER
		_f_JER = new TFile(xfn + "/../pPbFragmentation/data/jet_perf_histos.root","read");
		_nEta_JER = 8;			// number of eta bins
		_nCent_JER = 8;			// number of centrality bins
								//Shifted by unity in Laura's file
		for(int i=0;i<_nEta_JER;i++){
			for(int j=0;j<_nCent_JER;j++){
				_h_JER[i][j]=(TH1F*)_f_JER->Get(Form("h1_JER_pT_y%i_c%i",i,j+1));
			}
		}

		//momentum resolution
		_f_TMR = new TFile(xfn + "/../pPbFragmentation/data/trk_resolution.root","read");
		//Shifted by unity in Laura's file
		for(int i=0;i<_nEta_eff;i++){
			_tf1_TMR[i]=(TF1*)_f_TMR->Get(Form("pp_fit_eta_%i",i));
		}


		//alignement correction
		f_sagitta = new TFile(xfn + "/../pPbFragmentation/data/5TeVHI2015_sagittaBias_pTmethod.root");
		h_sagitta = (TH2*)f_sagitta->Get("h_deltaSagittaMap");

		//For mutliplicity cuts
		nMultThresholds=6;
		MultThresholds[0]=1;MultThresholds[1]=4;MultThresholds[2]=6;MultThresholds[3]=10;MultThresholds[4]=20;MultThresholds[5]=50;

	}// end of constructor


	Int_t GetJetEtaBin(float eta);
	Int_t GetTrackEtaBin(float eta);
	Int_t GetxBin(float x, TAxis* axis);
	Int_t GetJetYBin(float y);
	Int_t GetdRBin(float R);
	void InitdRBinRange();
	//Int_t TrackCorrector::GetCentBin(int cent);
	bool PassTracktoJetBalance(float trk_pt,float jet_pt,float trk_eta,float jet_eta,int cent);
	float get_effcorr(float pt, float eta, int centrality, float uncert, float jet_pt, float jet_eta, float R);
	float get_effcorr(float pt, float eta, int centrality, float uncert);
	void correctChTrackpT(float &pt, int charge);
	void correctChTrackpT(float &pt, float eta, float phi, int charge);
	~TrackCorrector() {}
};


#endif
