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

#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"
#include "InDetTrackSystematicsTools/InDetTrackSmearingTool.h"

#include "AsgTools/AnaToolHandle.h"
#include "JetInterface/IJetSelector.h"


#include "JetSubstructure/JetHelperTools.h"
#include "JetSubstructure/GlobalHelper.h"
#include "JetSubstructure/UncertProvider.h"


//Pileup tool
using namespace std;


namespace InDet { class InDetTrackSmearingTool; }


class BaseClass : public EL::Algorithm
{

public:

	//Configuration from main macro
	int _isPP; //
	int _isMC; //
	bool _doJES; //
	string _dataset; //
	int _isMB; //
	int _towerBkgKill; //
	bool _doTrimming; //
	int _defTrimAlgo; //
	float _fCut; //
	float _rSub; //
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
	std::unique_ptr<InDet::InDetTrackSmearingTool> m_trkSmearingTool; //! 

	double _pt_iso; //
	bool _truth_iso; //
	bool _reco_iso; //
	float _JetRadiusAna; //
	bool _saveLog; //
	bool _saveNtuple; //
	bool _saveEvtDisplay; //
	float _pTtrkCutReco; //
	float _pTtrkCutTruth; //
	float _etaTrkCut; //
	float _ptCutPostCS; //
	string _trk_cut_level; //
	//trigger
	string _trigger_collection; //
	vector <string> trigger_chains; //TODO
	vector <float> trigger_thresholds; //TODO
	vector <vector<float>> jet_pt_trig; //TODO
	int _nTriggers; //
	bool _applyReweighting; //
	bool _hasHI; //
	
	int _defJetRecl; //
	float _csMaxR; //
	float _Rktjet_bkg; //
	float _ghost_area; //
	float _alphaSubtr; //
	float _ptCutJetConeExc; //

	float _beta; //
	float _z_cut; //
	
	float _chParticleMassMeV; //


	bool event_isTriggered[10]; //!	//TODO
	bool trigger[10]; //!	//TODO
	int jet_isTriggered[10]; //! //TODO
	float trig_prescale[10]; //! //TODO

	float ptSysHI[50]; //! 
	float ptSysPP[50]; //! 
	int intrinsicComponent[50]; //!
	float intSignificance[50]; //!

	float trkJetMass4; //! 
	float trkJetMass6; //! 
	float trkJetMass8; //! 
	float trkJetMass10; //! 
	float trkJetPt4; //! 
	float trkJetPt6; //! 
	float trkJetPt8; //! 
	float trkJetPt10; //! 

	float trkJetMassRcSub2; //!
	float trkJetMassRcSub4; //!
	float trkJetMassRcSubEV; //!
	float trkJetMassRcSubMV; //!



	TH1D *h_FCal_Et; //!
	TH1D *h_RejectionHisto; //!
	TH1D *h_centrality; //!
	TH3D *hET_ETsub; //!
	TH2D *h_triggercounter; //!

        TF1* f_d0_cut; //!

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
	JetCalibrationTool * m_jetCalibration_insitu; //!

	//Track Selection Tool
	InDet::InDetTrackSelectionTool * m_trackSelectorTool; //!
	
	// Jet cleaner 
	asg::AnaToolHandle<IJetSelector> m_jetCleaningToolHandle; //!

	// this is a standard constructor
	BaseClass ();
	BaseClass(const BaseClass& base);

	// this is needed to distribute the algorithm to the workers
	ClassDef(BaseClass, 1);
};

#endif

