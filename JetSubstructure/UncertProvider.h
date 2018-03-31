#ifndef UncertProvider_H
#define UncertProvider_H

#include "TFile.h"
#include "TH2F.h"
#include "TF1.h"
#include "TRandom3.h"
#include "JetSubstructure/JetHelperTools.h"
#include "JetUncertainties/JetUncertaintiesTool.h"
#include <fstream>
#include <TSystem.h>
#include <TAxis.h>
#include <TLorentzVector.h>
#include "xAODJet/JetContainer.h"
#include "JetCorrector.h"
#include "TrackCorrector.h"
#include "HIJESUncertaintyProvider.h"


using namespace std;
using namespace JetHelperTools;

class UncertProvider
{
  private:
  
	int GetSysShift(int uncert);
	int GetJESSysComponent(int uncert);
	void UncerJERIntrinsic_Gaus(xAOD::Jet* recon, xAOD::Jet* truth, int cent);
	void UncerJESIntrinsic(xAOD::Jet* recon);
	int GetEtaUJERBin(float eta);
	void UncerHIJESIntrinsic(xAOD::Jet* recon);
	void GetTrackUncert();
	void UncerHIJESCentrality(xAOD::Jet* recon, float FCalEt);
	int GetFinCentrality(float FCalEt);
	
	
	float UncerEffMaterial(float pt, float eta);
	//float UncerFit(float pt,TH1 * uncert){
	float UncerDense(float dR);
	
	//HI specific JES and JER
	TFile *f_HI_JER;
	TFile *f_HI_JES;
	TH1 * h_sJER_eta[7];
	TH1 * h_JER_eta[7];
	TH1 * h_sJES_eta[7];
	TRandom3 r;
	JetCorrector JERhelper;
	//	TrackCorrector Trackhelper;
	JetUncertaintiesTool jesProv;
	HIJESUncertaintyProvider HIJESProvider;
	
	//Tracking
	TFile *f_eff_uncert_2015_material;
	TH2F *h_eff_uncert_2015_material[4];
	TFile* f_sagitta;
	TH2 * h_sagitta;
	
	//TFile * _f_eff;
	//TH1F *_th1_eff_unc[6][10];
	//Int_t _nEta_eff;
	
	//Fine centrality for JES uncertainty
	int nFineCentBins;
	float fine_cent_bins[60];

  public:	
    int uncert_index;
    int uncert_class;
    float mcprobcut;
    string _cutlevel;
    int _nCent;
    
    UncertProvider(int ui, float mcprob, string cutlevel, int nCentBins, bool eff_jety):
     uncert_index(ui),
     mcprobcut(mcprob),
     _cutlevel(cutlevel),
     _nCent(nCentBins), 
     jesProv("JESProvider"),
     HIJESProvider("HIJESUncert_data15_5TeV.root"),
     fine_cent_bins{0.289595,0.308686,0.328744,0.349697,0.371561,0.394518,0.418573,0.443549,0.46959,0.49675,0.525092,0.554569,0.585275,0.617108,0.65018,0.684377,0.719896,0.756791,0.795018,0.834538,0.87541,0.917795,0.961609,1.0068,1.05367,1.10211,1.15214,1.20373,1.25693,1.31197,1.36875,1.42719,1.48744,1.55005,1.61434,1.68058,1.74932,1.81997,1.89316,1.96859,2.04651,2.12711,2.21002,2.29572,2.38468,2.47658,2.57162,2.66999,2.77237,2.87864,2.98931,3.10407,3.22397,3.34945,3.48077,3.61844,3.7635,3.91763,4.08137,4.26258}

       //     Trackhelper(cutlevel.c_str(),nCentBins, eff_jety)

     {
     	cout << "Inicialization of UncertProvider" << endl;
     	if (uncert_index==0) cout << "UncertProvider: No systematic variation is used...." << endl;
     	else cout << "Systematic uncertainty: " << GetSysName(uncert_index) << endl;
     	
     	
     	//***ppJES***
		//jesProv = new JetUncertaintiesTool("JESProvider");
		jesProv.setProperty("JetDefinition","AntiKt4EMTopo");
		jesProv.setProperty("MCType","MC15");
		jesProv.setProperty("ConfigFile","JES_2015/ICHEP2016/JES2015_19NP.config");
		// Initialize the tool
		jesProv.initialize();
     	
     	//***HI jet specific uncertainties***
     	TString xfn = gSystem->GetFromPipe("echo $ROOTCOREBIN");
     	
     	HIJESProvider.UseJESTool(true);
     	HIJESProvider.UseGeV(false);
     	
     	//HI JES <-> crosscalibration
     	TFile *f_HI_JES = new TFile(xfn + "/../JetSubstructure/data/cc_sys_090816.root","read");
     	nFineCentBins=60;
     	
     	//JER
     	f_HI_JER  = new TFile(xfn + "/../JetSubstructure/data/hi_jer.root","read");
   		
   		for(int ybin=0;ybin<7;ybin++) {
		  	h_sJER_eta[ybin] =(TH1*)f_HI_JER->Get(Form("delta_sigma_hi_eta%i",ybin));
		  	h_JER_eta[ybin] =(TH1*)f_HI_JER->Get(Form("sigma_hi_eta%i",ybin));
		  	h_sJES_eta[ybin] =(TH1*)f_HI_JES->Get(Form("fsys_rel_%i",ybin));  
		}
     	//***Tracking uncertainties
     	//Efficiency <->fit not needed now
     	/*
		  _nEta_eff=5;	  
		  _f_eff = new TFile(xfn + "/../pPbFragmentation/data/mc_eff_fits_" + _cutlevel + ".root","read");
	        
		  for(int cent=0;cent<_nCent;cent++){	  
			  for(int e=0;e<_nEta_eff;e++){	
				_th1_eff_unc[e][cent]=(TH1F*)_f_eff->Get(Form("conf_int_eta%i_cent%i",e,cent));
			  }
		  }
     	*/
     	//Material
     	f_eff_uncert_2015_material = new TFile(xfn + "/../JetSubstructure/data/TrackingEfficiencyRecommendations_20.7rel.root","read");
   		h_eff_uncert_2015_material[0]	= (TH2F*)f_eff_uncert_2015_material->Get("OneMinusRatioEfficiencyVSEtaPt_AfterRebinning_Old_Nominal_MCVSOld_5%Extra_MC_TightPrimary");
  		h_eff_uncert_2015_material[1]	= (TH2F*)f_eff_uncert_2015_material->Get("OneMinusRatioEfficiencyVSEtaPt_AfterRebinning_Old_Nominal_MCVSNew_Nominal_MC_TightPrimary");
   		h_eff_uncert_2015_material[2]	= (TH2F*)f_eff_uncert_2015_material->Get("OneMinusRatioEfficiencyVSEtaPt_AfterRebinning_Old_Nominal_MCVSOld_50%ExtraPP0_MC_TightPrimary");
   		h_eff_uncert_2015_material[3]	= (TH2F*)f_eff_uncert_2015_material->Get("OneMinusRatioEfficiencyVSEtaPt_AfterRebinning_Old_Nominal_MCVSOld_FTF_BIC_MC_TightPrimary");
     	
     	//Momentum
     	//Preliminary
     	//f_sagitta = new TFile(xfn + "/../pPbFragmentation/data/correctionmaps_HighGran_IBLon_NoGRL_INDET_2015_datareproAll25ns_correctedEp.root");
     	//h_sagitta = (TH2*)f_sagitta->Get("LambdaCorrectionVsEtaPhi_reweightedToEP");
     	//Statistical ucnertainty
     	f_sagitta = new TFile(xfn + "/../pPbFragmentation/data/5TeVHI2015_sagittaBias_pTmethod_statUncertainty.root");
	h_sagitta = (TH2*)f_sagitta->Get("h_deltaSagittaMap_statErr");
     	
	//     	GetTrackUncert();
     	
     }// end of constructor
		
	
	string GetSysName(int uncert);
	void CorrectJet(xAOD::Jet * reco, xAOD::Jet * truth, int cent, float FCalEt);
	float CorrectTrackEff(float pt, float eta, float dR, int cent);
	void UncerTrackMomentum(float &pT, float eta, float phi, int charge);
	//float UncerEffFit(float pt, float eta, int centrality);
	float GetMCProb();
	//CorrectTrack(xAODJet * reco);
    ~UncertProvider() {
    }
};


#endif
