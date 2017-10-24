#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "JetSubstructure/JetSubstructure.h"
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/Message.h"
#include "xAODJet/JetContainer.h"
#include "xAODTrigger/JetRoIContainer.h"
#include "xAODTrigger/JetRoIAuxContainer.h"
#include <TFile.h>
#include <TSystem.h>
#include <AsgTools/MessageCheck.h>


// this is needed to distribute the algorithm to the workers
ClassImp(BaseClass)

BaseClass :: BaseClass ()
{
	// Here you put any code for the base initialization of variables,
	// e.g. initialize all pointers to 0.  Note that you should only put
	// the most basic initialization here, since this method will be
	// called on both the submission and the worker node.  Most of your
	// initialization code will go into histInitialize() and
	// initialize().
}

BaseClass :: BaseClass (const BaseClass& base)  {

	_isMC = base._isMC;
	_dataset = base._dataset;
	_isMB = base._isMB;
	_isHerwig = base._isHerwig;
	_jet_radius = base._jet_radius;
	_reco_jet_collection = base._reco_jet_collection;
	_test_reco_jet_collection = base._test_reco_jet_collection;
	_truth_jet_collection=base._truth_jet_collection;
	_GRL=base._GRL;
	_centrality_scheme=base._centrality_scheme;
	_dR_truth_matching = base._dR_truth_matching;
	_outputName=base._outputName;
	_etaJetCut=base._etaJetCut;
	_pTjetCut=base._pTjetCut;
	_truthpTjetCut=base._truthpTjetCut;
	_pt_iso=base._pt_iso;
	_applyReweighting=base._applyReweighting;
	_hasHI=base._hasHI;
	_truth_iso=base._truth_iso;
	_reco_iso=base._reco_iso;
	_ReclusterCA=base._ReclusterCA;
	_ReclusterRadius=base._ReclusterRadius;
        _trk_cut_level=base._trk_cut_level;
}













