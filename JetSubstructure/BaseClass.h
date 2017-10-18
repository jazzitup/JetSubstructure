#ifndef BaseClass_h
#define BaseClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
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
#include "xAODCore/ShallowAuxContainer.h"
#include "xAODCore/ShallowCopy.h"

#define private public
#include "xAODHIEvent/HIEventShapeAuxContainer.h"
#undef private
#include "xAODHIEvent/HIEventShapeContainer.h"

#include "JetSubstructure/JetHelperTools.h"
#include "JetSubstructure/GlobalHelper.h"
#include <iostream>


using namespace std;

class BaseClass : public EL::Algorithm
{

public:

	//Configuration from main macro
	int _isMC; //
	string _dataset; //
	int _isMB; //
	int _isHerwig; //
	int _jet_radius; //
	string _reco_jet_collection; //
	string _test_reco_jet_collection; //
	string _truth_jet_collection; //
	string _GRL; //
	string _cut_level; //
	int _centrality_scheme; //
	float _dR_truth_matching; //
	float _etaJetCut; //
	float _pTjetCut; //
	float _truthpTjetCut; //
	std::string _outputName; //
	double _pt_iso; //
	bool _truth_iso; //
	bool _reco_iso; //
	bool _ReclusterCA;
	float _ReclusterRadius;
	bool _saveLog;
	bool _pTtrkCut;
	//trigger
	string _trigger_collection; //
	vector <string> trigger_chains; //TODO
	vector <float> trigger_thresholds; //TODO
	vector <vector<float>> jet_pt_trig; //TODO
	int _nTriggers; //
	bool _applyReweighting; //
	bool _hasHI; //
	
	
	bool event_isTriggered[10]; //!	//TODO
	bool trigger[10]; //!	//TODO
	int jet_isTriggered[10]; //! //TODO
	float trig_prescale[10]; //! //TODO

	TH1D *h_FCal_Et; //!
	TH1D *h_RejectionHisto; //!
	TH1D * h_centrality; //!
	TH3D *hET_ETsub; //!
	TH2D *h_triggercounter; //!

	//Evnets
	int m_eventCounter; //!

	//GRL
	GoodRunsListSelectionTool *m_grl; //!

	//Trigger tools member variables
	Trig::TrigDecisionTool *m_trigDecisionTool; //!
	TrigConf::xAODConfigTool *m_trigConfigTool; //!
	vector<const Trig::ChainGroup*> _chainGroup; //!
	vector<const Trig::ChainGroup*> _referenceChainGroup; //!
	int _first_trigger;

	void SetTrigger_chains();
	void SetTrigger_hist(TH2D* h);

	//Calibration tools
	JetCalibrationTool * m_jetCalibration; //!
	JetCalibrationTool * m_jetCalibration_val; //!

	// this is a standard constructor
	BaseClass ();
	BaseClass(const BaseClass& base);

	// this is needed to distribute the algorithm to the workers
	ClassDef(BaseClass, 1);
};

#endif
