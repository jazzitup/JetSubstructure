#ifndef JetSubstructure_H
#define JetSubstructure_H
#include "BaseClass.h"
#include <EventLoop/Algorithm.h>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TMath.h>


#include "BaseClass.h"

#include <EventLoop/Algorithm.h>
#include "xAODEventInfo/EventInfo.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "TrigConfxAOD/xAODConfigTool.h"
#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "JetCalibTools/JetCalibrationTool.h"
#include <boost/regex.hpp>
#include "JetSelectorTools/JetCleaningTool.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/VertexContainer.h"


#define private public
#include "xAODHIEvent/HIEventShapeAuxContainer.h"
#undef private
#include "xAODHIEvent/HIEventShapeContainer.h"

#include "JetSubstructure/JetHelperTools.h"
#include "JetSubstructure/GlobalHelper.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertexContainer.h"
#include "JetSubstructure/JetCorrector.h"

#include "fastjet/JetDefinition.hh"
#include "JetSubstructure/TrackHelperTools.h"
#include "JetSubstructure/UEEstimator.h"
#include "TMath.h"
#include "ZdcAnalysis/ZdcAnalysisTool.h"
#include "HIEventUtils/HIPileupTool.h"
#include <iostream>

#include "JetRecTools/ConstituentSubtractorTool.h"
#include "JetSubstructure/UncertProvider.h"

//#include "fastjet/contrib/SoftDrop.hh"

class JetSubstructure : public BaseClass
{
	// put your configuration variables here as public variables.
	// that way they can be set directly from CINT and python.
public:

	int nCentbins;


	// variables that don't get filled at submission time should be
	// protected from being send from the submission node to the worker
	// node (done by the //!)
public:
	// Tree *myTree; //!
	// TH1 *myHist; //!


	bool event_isTriggered[10]; //!
	
	UEEstimator* uee; //!
	
	vector<TH3D*> h_resp_cent; //!
	JetCorrector* jetcorr; //!
	TH2D * h_respR2R4; //!

	TH1D* hPtGenRaw; //!
	TH1D* hPtGenWgt; //!


	TTree *treeOut; //!
	TH2F* eventDisRecTow; //!
	TH2F* eventDisRecTow1; //!
	TH2F* eventDisRecTow2; //!
	TH2F* eventDisGen; //!
	TH2F* eventDisGen1; //!
	TH2F* eventDisGen2; //!
	TH2F* eventDisTrk; //!
	TH2F* eventDisTrk1; //!
	TH2F* eventDisTrk2; //!
	TH2F* eventDisChg; //!
	TH2F* eventDisChg1; //!
	TH2F* eventDisChg2; //!



	TH2D* hGenNcam; //!
        TH2D* hGenNchCam; //!
	TH2D* hRecoNcam; //!
        TH2D* hRecoNchCam; //!

        TH2D* hGenSdStat; //!
        TH2D* hGenSdChStat; //!
        TH2D* hRecoSdStat; //!
        TH2D* hRecoSdChStat; //!

	
        //Uncert tool
        vector<UncertProvider*> vUncertprovider; //!

	vector<TH2D*> hTrkPtEta_preCS_cent; //!
	vector<TH2D*> hTrkPtEta_postCS_cent; //!
	vector<TH2D*> hTrkPtEta_genMatch_cent; //!


        vector<TH3D*> h_trkGen_pt_dphi_cent; //!
        vector<TH3D*> h_allGen_pt_dphi_cent; //!
        vector<TH3D*> h_trkGen_pt_drap_cent; //!
        vector<TH3D*> h_allGen_pt_drap_cent; //!
	
	vector<TH3D*> h_reco_jet_cent; //!
	vector<TH3D*> h_reco_jet_cent_matched; //!
	vector<TH3D*> h_reco_jet_cent_unmatched; //!
	vector<TH3D*> h_truth_jet_cent; //!
	vector<TH3D*> h_truth_jet_cent_matched; //!

	vector<TH3D*> h_gen_reclst_ratio_cent; //!
	vector<TH3D*> h_reco_reclst_ratio_cent; //!
	
	vector<TH2D*> h_trkPt_trkBkgPt_cent; //!
	vector<TH2D*> h_bkgSubt_prePt_postPt_cent; //!
	vector<TH2D*> h_dRSubt_trkPt_cent; //!

	vector<TH2D*> h_trkPt_trkBkgPt_jetCone_cent; //!
        vector<TH2D*> h_bkgSubt_prePt_postPt_jetCone_cent; //!
        vector<TH2D*> h_dRSubt_trkPt_jetCone_cent; //!


	//Prescale sets
	TFile * f_trigger_RunNumber_prescale; //!
	TH2F * h2_trigger_RunNumber_prescale; //!


	// this is a standard constructor
	JetSubstructure ();

	// these are the functions inherited from Algorithm
	virtual EL::StatusCode setupJob (EL::Job& job);
	virtual EL::StatusCode fileExecute ();
	virtual EL::StatusCode histInitialize ();
	virtual EL::StatusCode changeInput (bool firstFile);
	virtual EL::StatusCode initialize ();
	virtual EL::StatusCode execute ();
	virtual EL::StatusCode postExecute ();
	virtual EL::StatusCode finalize ();
	virtual EL::StatusCode histFinalize ();
	
	HI::HIPileupTool                      *m_hiPileup;    //!
	ZDC::ZdcAnalysisTool *m_zdcTools; //!

	// this is needed to distribute the algorithm to the workers
	ClassDef(JetSubstructure, 1);
};

#endif
