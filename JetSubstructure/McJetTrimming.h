#ifndef McJetTrimming_H
#define McJetTrimming_H
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

//#include "fastjet/contrib/SoftDrop.hh"

class McJetTrimming : public BaseClass
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

	UEEstimator* uee; //!
	

        TH1D* h_pt; //!
        vector<TH2D*> h_mass_trimMassKt; //!
        vector<TH2D*> h_mass_trimMassAk; //!
        vector<TH2D*> h_mass_trimRateKt; //!
        vector<TH2D*> h_mass_trimRateAk; //!
	TH2D* h_pt_trimRateKt; //!
	TH2D* h_pt_trimRateAk; //!
	// this is a standard constructor
	McJetTrimming ();

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
	ClassDef(McJetTrimming, 1);
};

#endif
