#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <JetSubstructure/JetSubstructure.h>

#include <JetCalibTools/JetCalibrationTool.h>
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/Message.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetTypes.h"
#include "xAODHIEvent/HIEventShapeAuxContainer.h"
#include "xAODHIEvent/HIEventShapeContainer.h"

#include <AsgTools/MessageCheck.h>
#include <TSystem.h>
#include <TMath.h>
#include <TStyle.h>
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include <JetSubstructure/JetSubstructure.h>
#include <TFile.h>

#include "xAODTrigger/JetRoIContainer.h"
#include "xAODTrigger/JetRoIAuxContainer.h"
#include "xAODForward/ZdcModuleContainer.h"

using namespace std;
using namespace JetHelperTools;


// this is needed to distribute the algorithm to the workers
ClassImp(JetSubstructure)



JetSubstructure :: JetSubstructure ()
{
	// Here you put any code for the base initialization of variables,
	// e.g. initialize all pointers to 0.  Note that you should only put
	// the most basic initialization here, since this method will be
	// called on both the submission and the worker node.  Most of your
	// initialization code will go into histInitialize() and
	// initialize().
}



EL::StatusCode JetSubstructure :: setupJob (EL::Job& job)
{
	// Here you put code that sets up the job on the submission object
	// so that it is ready to work with your algorithm, e.g. you can
	// request the D3PDReader service or add output files.  Any code you
	// put here could instead also go into the submission script.  The
	// sole advantage of putting it here is that it gets automatically
	// activated/deactivated when you add/remove the algorithm from your
	// job, which may or may not be of value to you.
	xAOD::TReturnCode::enableFailure();
	job.useXAOD ();
	ANA_CHECK_SET_TYPE (EL::StatusCode); // set type of return code you are expecting (add to top of each function once)
	ANA_CHECK(xAOD::Init());

	nCentbins = GetCentralityNBins(_centrality_scheme);
	cout << "Number of centrality bins: " << nCentbins <<endl; //all bins + 1 inclusive

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode JetSubstructure :: histInitialize ()
{
	// Here you do everything that needs to be done at the very
	// beginning on each worker node, e.g. create histograms and output
	// trees.  This method gets called before any input files are
	// connected.

	cout << " Setting  histograms" << endl;


	int etaBinsN, phiBinsN, respBinsN, jetPtBinsN;
	double etaBins[1000], phiBins[1000], respBins[1000], jetPtBins[1000];

	SetupBinning(0, "eta-fine", etaBins, etaBinsN);
	SetupBinning(0, "phi-jet", phiBins, phiBinsN);
	SetupBinning(0, "resp", respBins, respBinsN);
	SetupBinning(0, "pt-jet-PbPb", jetPtBins, jetPtBinsN);

	//Basic histograms
	h_FCal_Et = new TH1D("h_FCal_Et",";FCal E_{T};N",100,0,5);
	h_FCal_Et->Sumw2();

	h_RejectionHisto = new TH1D("RejectionHisto","RejectionHisto",8,0,8);
	SetRejectionHistogram(h_RejectionHisto);

	h_centrality = new TH1D("Centrality","Centrality",10,0,10);
	h_centrality->Sumw2();

	hET_ETsub = new TH3D("hET_ETsub","hET_ETsub",200,0,200,100,-5,5,200,-50,50);
	hET_ETsub->Sumw2();
	
	hPtGenRaw = new TH1D("hPtGenRaw",";p_{T};Entries",5000,0,1000);
	hPtGenRaw->Sumw2();
	hPtGenWgt = (TH1D*)hPtGenRaw->Clone("hPtGenWgt");
	hPtGenWgt->Sumw2();
		

	//	h_triggercounter = new TH2D("h_triggercounter","h_triggercounter",_nTriggers,0,_nTriggers,2,-0.5,1.5);
	//	SetTrigger_hist(h_triggercounter);
	
	h_respR2R4 = new TH2D("h_respR2R4","h_respR2R4",jetPtBinsN, jetPtBins,respBinsN,respBins);
	wk()->addOutput(h_respR2R4);
	wk()->addOutput(hPtGenRaw);
	wk()->addOutput(hPtGenWgt);
	
	wk()->addOutput (h_FCal_Et);
	wk()->addOutput (h_RejectionHisto);
	wk()->addOutput (h_centrality);
	wk()->addOutput (hET_ETsub);
	//	wk()->addOutput (h_triggercounter);

	TH3D* temphist_3D = nullptr;
	TH2D* temphist_2D = nullptr;
	TH1D* temphist_1D = nullptr;
	

	for (int i=0;i<nCentbins;i++)
	{
		temphist_3D = new TH3D(Form("h_resp_cent%i",i),Form("h_resp_cent%i",i),jetPtBinsN, jetPtBins,respBinsN,respBins,etaBinsN,etaBins);
		h_resp_cent.push_back(temphist_3D);
		h_resp_cent.at(i)->Sumw2();


		temphist_3D = new TH3D(Form("h_truth_jet_cent%i",i),Form("h_truth_jet_cent%i",i),jetPtBinsN, jetPtBins,etaBinsN, etaBins, phiBinsN, phiBins);
		h_truth_jet_cent.push_back(temphist_3D);
		h_truth_jet_cent.at(i)->Sumw2();
		
		temphist_3D = new TH3D(Form("h_truth_jet_cent%i_matched",i),Form("h_truth_jet_cent%i_matched",i),jetPtBinsN, jetPtBins,etaBinsN, etaBins, phiBinsN, phiBins);
		h_truth_jet_cent_matched.push_back(temphist_3D);
		h_truth_jet_cent_matched.at(i)->Sumw2();
		
		temphist_3D = new TH3D(Form("h_reco_jet_cent%i_matched",i),Form("h_reco_jet_cent%i_matched",i),jetPtBinsN, jetPtBins,etaBinsN, etaBins, phiBinsN, phiBins);
		h_reco_jet_cent_matched.push_back(temphist_3D);
		h_reco_jet_cent_matched.at(i)->Sumw2();
		
		temphist_3D = new TH3D(Form("h_reco_jet_cent%i_unmatched",i),Form("h_reco_jet_cent%i_unmatched",i),jetPtBinsN, jetPtBins,etaBinsN, etaBins, phiBinsN, phiBins);
		h_reco_jet_cent_unmatched.push_back(temphist_3D);
		h_reco_jet_cent_unmatched.at(i)->Sumw2();
		
		temphist_3D = new TH3D(Form("h_reco_jet_cent%i",i),Form("h_reco_jet_cent%i",i),jetPtBinsN, jetPtBins,etaBinsN, etaBins, phiBinsN, phiBins);
		h_reco_jet_cent.push_back(temphist_3D);
		h_reco_jet_cent.at(i)->Sumw2();
		
		if (_isMC){
			wk()->addOutput (h_resp_cent.at(i));
			wk()->addOutput (h_reco_jet_cent.at(i));
			wk()->addOutput (h_truth_jet_cent.at(i));
			wk()->addOutput (h_truth_jet_cent_matched.at(i));
			wk()->addOutput (h_reco_jet_cent_matched.at(i));
			wk()->addOutput (h_reco_jet_cent_unmatched.at(i));
		}
		wk()->addOutput (h_reco_jet_cent.at(i));
	}

	cout << " Histograms ready" << endl;


	return EL::StatusCode::SUCCESS;
}



EL::StatusCode JetSubstructure :: fileExecute ()
{
	// Here you do everything that needs to be done exactly once for every
	// single file, e.g. collect a list of all lumi-blocks processed
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode JetSubstructure :: changeInput (bool firstFile)
{
	// Here you do everything you need to do when we change input files,
	// e.g. resetting branch addresses on trees.  If you are using
	// D3PDReader or a similar service this method is not needed.
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode JetSubstructure :: initialize ()
{
	// Here you do everything that you need to do after the first input
	// file has been connected and before the first event is processed,
	// e.g. create additional histograms based on which variables are
	// available in the input files.  You can also create all of your
	// histograms and trees in here, but be aware that this method
	// doesn't get called if no events are processed.  So any objects
	// you create here won't be available in the output if you have no
	// input events.
	ANA_CHECK_SET_TYPE (EL::StatusCode); // set type of return code you are expecting (add to top of each function once)
	xAOD::TEvent* event = wk()->xaodEvent();
	const xAOD::EventInfo* eventInfo = 0;
	m_eventCounter = 0;

	ANA_CHECK(event->retrieve( eventInfo, "EventInfo") );

	Info("initialize()", "Number of events = %lli", event->getEntries() ); // print long long int


	//Calibration tool
	const std::string name = "JetSubstructure"; //string describing the current thread, for logging
	TString jetAlgo = "AntiKt4HI"; //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo (see below)
	TString config = "JES_MC15c_HI_Nov2016.config"; //Path to global config used to initialize the tool (see below)
	TString calibSeq = "EtaJES"; //String describing the calibration sequence to apply (see below)

	//Call the constructor. The default constructor can also be used if the arguments are set with python configuration instead
	m_jetCalibration = new JetCalibrationTool (name);
	ANA_CHECK(m_jetCalibration->setProperty("JetCollection",jetAlgo.Data()));
	ANA_CHECK(m_jetCalibration->setProperty("ConfigFile",config.Data()));
	ANA_CHECK(m_jetCalibration->setProperty("CalibSequence",calibSeq.Data()));
	ANA_CHECK(m_jetCalibration->setProperty("IsData",!_isMC));
	ANA_CHECK(m_jetCalibration->initializeTool(name));
	
	//Calibration tool
	
	const std::string name_val = "JetSubstructure_val"; //string describing the current thread, for logging
	TString jetAlgo_val = "AntiKt4EMTopo"; //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo (see below)
	TString config_val = "JES_Full2012dataset_Preliminary_Jan13.config"; //Path to global config used to initialize the tool (see below)
	TString calibSeq_val = "AbsoluteEtaJES"; //String describing the calibration sequence to apply (see below)
	
	//Call the constructor. The default constructor can also be used if the arguments are set with python configuration instead
	
	m_jetCalibration_val = new JetCalibrationTool (name_val);
	ANA_CHECK(m_jetCalibration_val->setProperty("JetCollection",jetAlgo_val.Data()));
	ANA_CHECK(m_jetCalibration_val->setProperty("ConfigFile",config_val.Data()));
	ANA_CHECK(m_jetCalibration_val->setProperty("CalibSequence",calibSeq_val.Data()));
	ANA_CHECK(m_jetCalibration_val->setProperty("IsData",!_isMC));
	ANA_CHECK(m_jetCalibration_val->initializeTool(name_val));
	
	jetcorr = new JetCorrector();
	
	//Pileup tool
	
	// ZDCAnalysisTool
	m_zdcTools = new ZDC::ZdcAnalysisTool("ZdcAnalysisTool");
	// HIPileupTool
	m_hiPileup = new HI::HIPileupTool("PileupTool");
	
	ANA_CHECK(m_hiPileup->initialize());
	ANA_CHECK(m_zdcTools->initializeTool());
	
	// GRL
	TString xfn = gSystem->GetFromPipe("echo $ROOTCOREBIN");
	TString xmlfile = xfn + "/../JetSubstructure/data/"+ _GRL;

	m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
	std::vector<std::string> vecStringGRL;
	vecStringGRL.push_back(xmlfile.Data());
	ANA_CHECK(m_grl->setProperty( "GoodRunsListVec", vecStringGRL));
	ANA_CHECK(m_grl->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
	ANA_CHECK(m_grl->initialize());
	
	// Initialize and configure trigger tools
	if (_isMC==0)
	{
		m_trigConfigTool = new TrigConf::xAODConfigTool("xAODConfigTool"); // gives us access to the meta-data
		m_trigConfigTool->msg().setLevel( MSG::ERROR );
		ANA_CHECK(m_trigConfigTool->initialize());
		ToolHandle< TrigConf::ITrigConfigTool > trigConfigHandle( m_trigConfigTool );

		m_trigDecisionTool = new Trig::TrigDecisionTool("TrigDecisionTool");
		m_trigDecisionTool->msg().setLevel( MSG::ERROR );
		ANA_CHECK(m_trigDecisionTool->setProperty( "ConfigTool", trigConfigHandle ));
		ANA_CHECK(m_trigDecisionTool->setProperty( "TrigDecisionKey", "xTrigDecision"));
		ANA_CHECK(m_trigDecisionTool->initialize() );

		cout << "Adding following " << _nTriggers << " triggers: ";
		for (int i=0;i<_nTriggers;i++){
			cout << trigger_chains.at(i) << ", ";
			_chainGroup.push_back(m_trigDecisionTool->getChainGroup(trigger_chains.at(i)));
		}
		
		TString xfn = gSystem->GetFromPipe("echo $ROOTCOREBIN");
		f_trigger_RunNumber_prescale = new TFile(xfn + "/../pPbFragmentation/data/TriggerPrescales.root","READ");
		h2_trigger_RunNumber_prescale = (TH2F*)f_trigger_RunNumber_prescale->Get("h2_Trig_RunNumber_prescale");
		
		cout << endl << "Initialize triggers finished" << endl;
	}

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode JetSubstructure :: execute ()
{
	// Here you do everything that needs to be done on every single
	// events, e.g. read input variables, apply cuts, and fill
	// histograms and trees.  This is where most of your actual analysis
	// code will go.
	ANA_CHECK_SET_TYPE (EL::StatusCode); // set type of return code you are expecting (add to top of each function once)
	xAOD::TEvent* event = wk()->xaodEvent();

	int statSize=1;
	if(m_eventCounter!=0)
	{
		double power=std::floor(log10(m_eventCounter));
		statSize=(int)std::pow(10.,power);
	}
	if(m_eventCounter%statSize==0) std::cout << "Event: " << m_eventCounter << std::endl;
	m_eventCounter++;

	//All events
	bool keep = true;
	h_RejectionHisto->Fill(0.5);


	const xAOD::EventInfo* eventInfo = 0;
	ANA_CHECK(event->retrieve( eventInfo, "EventInfo"));

	// check if the event is data overlay or MC
	bool isHIJING = false; //For centrality
	if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) )
	{
		isHIJING = true;
	}
	else
	{
		const xAOD::TruthParticleContainer * particles = 0;
		if( event->xAOD::TVirtualEvent::retrieve(particles, "TruthParticles", true) )
		{
			isHIJING = false;
		}
	}
	//Get centrality bin and centile. Centile used for MB weighting (from MB_FCal_Normalization.txt)
	double FCalEt = 0;
	int cent_bin = 0;
	double event_weight_fcal = 1;
	//Centrality
	const xAOD::HIEventShapeContainer* calos=0;
	ANA_CHECK(event->retrieve( calos, "CaloSums"));
	if (_centrality_scheme>1)
	{
		FCalEt=calos->at(5)->et()*1e-6;
		cent_bin = GetCentralityBin(_centrality_scheme, FCalEt, isHIJING);

		event_weight_fcal = jetcorr->GetFCalWeight(FCalEt);

		if (_isMC && isHIJING) event_weight_fcal = 1;
	}

	if (cent_bin < 0)
	{
		h_RejectionHisto->Fill(1.5);
		keep = false;
	}

	// GRL
	if(!_isMC)
	{
		if(!m_grl->passRunLB(*eventInfo))
		{
			h_RejectionHisto->Fill(2.5);
			keep = false;
		}
	}
	
	//Vertex requirement
	const xAOD::VertexContainer * vertices = 0;
	if ( !event->retrieve( vertices, "PrimaryVertices" ).isSuccess() )
	{
		Error("execute()", "Failed to retrieve VertexContainer container. Exiting." );
		return EL::StatusCode::FAILURE;
	}
	if(vertices->size()<2)
	{
		h_RejectionHisto->Fill(3.5);
		keep = false;
	}
	
	//DAQ errors
	if(!_isMC){
		if(   (eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::SCT)==xAOD::EventInfo::Error ) || (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) ) ){
			h_RejectionHisto->Fill(4.5);
			keep = false;
		}
	}

	//Pileup		
	bool m_is_pileup = false;
	if (!_isMC) {	
		const xAOD::ZdcModuleContainer* zdcMod = 0;
		ANA_CHECK(event->retrieve( zdcMod, "ZdcModules")); 
		//const xAOD::HIEventShapeContainer* hiev = 0;
		//EL_RETURN_CHECK("execute",event->retrieve( hiev, "HIEventShape"));
		// ZDC
		m_zdcTools->reprocessZdc();
		// is Pileup
		m_is_pileup = m_hiPileup->is_pileup( *calos, *zdcMod); // SAVE pileup Decision HERE 0 = NO pileup, 1 = pileup
	}
	else m_is_pileup = (FCalEt > 4.8); //Remove pileup in MC
			
	if (m_is_pileup){
		h_RejectionHisto->Fill(6.5);
		keep = false;
	}
	
	
	if (!keep) return EL::StatusCode::SUCCESS; // go to the next event
	h_RejectionHisto->Fill(7.5);
	
	// trigger
	if (_isMC==0)
	{
		int event_passed_trigger=0;

		for (int i=0;i<_nTriggers;i++){

			event_isTriggered[i] = false;

			event_isTriggered[i] =  _chainGroup.at(i)->isPassed();
			h_triggercounter->Fill(i, (Double_t) event_isTriggered[i]);
			if(event_isTriggered[i]) event_passed_trigger=1;
		}

		if(!event_passed_trigger) return EL::StatusCode::SUCCESS; // go to next event
		else h_RejectionHisto->Fill(8.5);
	}
	
	h_FCal_Et->Fill(FCalEt, event_weight_fcal); //filled here to get proper event weight
	h_centrality->Fill(cent_bin,event_weight_fcal);
	
	//**Get Truth jets ***
	double event_weight = 1;
	double max_pt = 1;
	
	vector <double> truth_jet_pt_vector;
	vector <double> truth_jet_eta_vector;
	vector <double> truth_jet_phi_vector;
	vector<bool> truth_jet_isolation_vector;

	if (_isMC && _ReclusterTruthJets){
	  if(m_eventCounter%statSize==0) 	    cout << " re cluster!" << endl;
	  //**************** reclustering truth ***********
	
		std::vector<fastjet::PseudoJet> inputConst;
		const xAOD::TruthParticleContainer * particles = 0;
		ANA_CHECK(event->retrieve( particles, "TruthParticles"));
		xAOD::TruthParticleContainer::const_iterator truth_itr = particles->begin();
		xAOD::TruthParticleContainer::const_iterator truth_end = particles->end();
		for( ; truth_itr!=truth_end; ++truth_itr){
			int ty=TrackHelperTools::getTypeTruth((*truth_itr)->barcode(),(*truth_itr)->pdgId(),(*truth_itr)->status());
			if(ty!=1 && ty!=5) continue;
			//Pt cut on minimum track pt for reclustering;
			if ((*truth_itr)->pt() < 6000.) continue;
			inputConst.push_back( (*truth_itr)->p4() );		
		}
	
		fastjet::Selector jselector = fastjet::SelectorAbsRapRange(0.0,6.);
		fastjet::JetAlgorithm algo = fastjet::antikt_algorithm;
		float jetR = 0.4;
		fastjet::JetDefinition jetDef(algo, jetR);
		fastjet::ClusterSequence cs(inputConst, jetDef);
		vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
		// print out some infos
	  	//cout << "Clustering with " << jetDef.description() << endl;

	  	for (unsigned i = 0; i < jets.size(); i++) {
			
			double phi = jets[i].phi();
			if (phi> TMath::Pi()) phi = phi - 2*TMath::Pi();
			double pt = jets[i].pt()*0.001;
			double eta = jets[i].pseudorapidity();
			
			//if (jets[i].pt()*0.001 > 20.) cout << "jet " << i << ": "<< jets[i].pt()*0.001 << " " << jets[i].pseudorapidity() << " " << phi << endl;
			
			if (pt < _truthpTjetCut) continue;
			
			if (pt > max_pt)  // && _centrality_scheme>1) //  Removed centrlaity condition  (YS)     //event weight from leading truth jet (only if PbPb) 
			{

				event_weight = jetcorr->GetJetWeight(pt, eta, phi);
				max_pt = pt;
				//				if (_isMC && isHIJING) event_weight = 1;
			}

			truth_jet_pt_vector.push_back(pt);
			truth_jet_eta_vector.push_back(eta);
			truth_jet_phi_vector.push_back(phi);
		}
		inputConst.clear();
	}  
    
	else if (_isMC){
	  if(m_eventCounter%statSize==0)	  cout << "We don't do re cluster truth jets!" << endl;
		//**************** Getting truth ****************
	
		const xAOD::JetContainer* truth_jets = 0;
		ANA_CHECK(event->retrieve(truth_jets, _truth_jet_collection.c_str() ));

		xAOD::JetContainer::const_iterator truth_jet_itr = truth_jets->begin();
		xAOD::JetContainer::const_iterator truth_jet_end = truth_jets->end();

		for( ; truth_jet_itr != truth_jet_end; ++truth_jet_itr )
		{
			xAOD::JetFourMom_t jet_truth_4mom = (*truth_jet_itr)->jetP4();

			double pt     = (jet_truth_4mom.pt() * 0.001 );
			double eta    = (jet_truth_4mom.eta());
			double phi    = (jet_truth_4mom.phi());
		
			//if (pt > 20.) cout << "jet : "<< pt << " " << eta << " " << phi << endl;

                        if (pt > max_pt)  // && _centrality_scheme>1) //  Removed centrlaity condition  (YS)     //event weight from leading truth jet (only if PbPb)
			  {
			    event_weight = jetcorr->GetJetWeight(pt, eta, phi);
			    max_pt = pt;
			    //			    if (_isMC && isHIJING) event_weight = 1;
			  }
			
			//Not filling truth map here because don't know event weight yet

			if (pt < _truthpTjetCut) continue;

			truth_jet_pt_vector.push_back(pt);
			truth_jet_eta_vector.push_back(eta);
			truth_jet_phi_vector.push_back(phi);
		}
		truth_jet_isolation_vector = GetIsolation(truth_jet_pt_vector,truth_jet_eta_vector,truth_jet_phi_vector, 1.0, _pt_iso); // -1 is for pT_iso == jet pT
	}	
	


	event_weight = event_weight*event_weight_fcal; //event weight is only set if MC. Otherwise default is 1.
	
	//cout << "event_weight " << event_weight << endl;
	
	//Reco jets
	vector <double> reco_jet_pt_vector;
	vector <double> reco_jet_eta_vector;
	vector <double> reco_jet_phi_vector;
	vector<bool> reco_jet_isolation_vector;
	vector<bool> reco_jet_TM_vector;


	xAOD::TStore *store = new xAOD::TStore; //For calibration
	const xAOD::JetContainer* reco_jets = 0;
	ANA_CHECK(event->retrieve( reco_jets, _reco_jet_collection.c_str() ));

	xAOD::JetContainer* updatedjets = new xAOD::JetContainer();
	xAOD::AuxContainerBase* updatedjetsAux = new xAOD::AuxContainerBase();
	updatedjets->setStore( updatedjetsAux );

	xAOD::JetContainer::const_iterator jet_itr = reco_jets->begin();
	xAOD::JetContainer::const_iterator jet_end = reco_jets->end();

	ANA_CHECK(store->record( updatedjets, "updatedjets" ));
	ANA_CHECK(store->record( updatedjetsAux, "updatedjetsAux" ));

	for( ; jet_itr != jet_end; ++jet_itr )
	{
		xAOD::Jet newjet;
		newjet.makePrivateStore( **jet_itr );
		
		xAOD::Jet newjet_val;
		newjet_val.makePrivateStore( **jet_itr );

		xAOD::JetFourMom_t jet_4mom = newjet.jetP4("JetSubtractedScaleMomentum");
		if(m_eventCounter%statSize==0)		cout << " indexCali = " << _indexCali << endl;
		if (_indexCali == 0)  
		  jet_4mom = newjet.jetP4("JetSubtractedScaleMomentum");
		if (_indexCali == 1)  
		  jet_4mom = newjet.jetP4();


		//		double uncalib_jet_phi2  = (jet_4mom_cal2.phi());
		
		const xAOD::JetFourMom_t jet_4mom_unsubtracted  = newjet.jetP4("JetUnsubtractedScaleMomentum");
                if (_indexCali == 0) {
		  newjet.setJetP4("JetPileupScaleMomentum",jet_4mom); //Setting PileupScale and ConstitScale because they are not in DFAntiKt4HI
		  newjet.setJetP4("JetConstitScaleMomentum",jet_4mom_unsubtracted);
		  newjet.setJetP4("JetEMScaleMomentum",jet_4mom);
		  newjet_val.setJetP4("JetPileupScaleMomentum",jet_4mom); //Setting PileupScale and ConstitScale because they are not in DFAntiKt4HI
		  newjet_val.setJetP4("JetConstitScaleMomentum",jet_4mom_unsubtracted);
		  newjet_val.setJetP4("JetEMScaleMomentum",jet_4mom);
		  ANA_CHECK(m_jetCalibration->applyCalibration( newjet ) );
		  ANA_CHECK(m_jetCalibration_val->applyCalibration( newjet_val ) );
		}
		
		const xAOD::JetFourMom_t jet_4mom_xcalib = newjet.jetP4();
		const xAOD::JetFourMom_t jet_4mom_xcalib_val = newjet_val.jetP4();

		double jet_pt  = (newjet.pt() * 0.001);
		double jet_pt_val  = (newjet_val.pt() * 0.001);
		double jet_eta = newjet.eta();
		double jet_phi = newjet.phi();
		
		h_respR2R4->Fill(jet_pt,(jet_pt-jet_pt_val)/jet_pt);
	
		if (jet_pt < _pTjetCut) continue;

		reco_jet_pt_vector.push_back(jet_pt);
		reco_jet_phi_vector.push_back(jet_phi);
		reco_jet_eta_vector.push_back(jet_eta);
		reco_jet_TM_vector.push_back(false);
		
		//		cout << "test 1" << endl;
		const xAOD::JetConstituentVector constituents = (*jet_itr)->getConstituents();
		//		cout << "n constituents" <<(*jet_itr)->numConstituents() << endl;
		//		cout << "test 2" << endl;
		double sumETConst=0;
		int iterator=0;
		/*
		for (xAOD::JetConstituentVector::iterator itr = constituents.begin(); itr != constituents.end(); ++itr){ 
			iterator++; 
			//cout << "et " << (*itr)->e() << endl;
			//cout << "test 3" << endl;
			const xAOD::CaloCluster* cl=static_cast<const xAOD::CaloCluster*>(itr->rawConstituent());
			//cout << "test 4" << endl; 
			sumETConst+= cl->rawE() * 0.001 / std::cosh( cl->rawEta() ); 
		}
		*/

	}

	reco_jet_isolation_vector = GetIsolation(reco_jet_pt_vector,reco_jet_eta_vector,reco_jet_phi_vector, 1.0, _pt_iso); // -1 is for pT_iso == jet pT


	store->clear();
	delete store;
	
	if (_isMC){
		for (int i = 0; i<truth_jet_pt_vector.size(); i++)
		{
			double truth_pt  = truth_jet_pt_vector.at(i);
			double truth_eta  = truth_jet_eta_vector.at(i);
			double truth_phi  = truth_jet_phi_vector.at(i);

			if (_truth_iso) { if (!truth_jet_isolation_vector.at(i)) continue;}

			//	cout << " event weight ="  << event_weight << endl; // (YS)
			h_truth_jet_cent.at(cent_bin)->Fill(truth_pt,truth_eta, truth_phi, event_weight);
			h_truth_jet_cent.at(nCentbins-1)->Fill(truth_pt,truth_eta, truth_phi, event_weight);
			
			hPtGenRaw->Fill( truth_pt);
			hPtGenWgt->Fill( truth_pt, event_weight);
			
			double d_R_min = 9999;
			double d_R_itr = 0;

			double matched_reco_pt = 0;
			double matched_reco_eta = 0;
			double matched_reco_phi = 0;
			double weight = 0;

			for (int j = 0; j<reco_jet_pt_vector.size(); j++)
			{
				double reco_pt  = reco_jet_pt_vector.at(j);
				double reco_eta = reco_jet_eta_vector.at(j);
				double reco_phi = reco_jet_phi_vector.at(j);

				if (_reco_iso) { if (!reco_jet_isolation_vector.at(j)) continue;}

				double d_R_itr = DeltaR(reco_phi, reco_eta, truth_phi, truth_eta);

				if (d_R_itr < d_R_min)
				{
					matched_reco_pt = reco_pt;
					matched_reco_eta = reco_eta;
					matched_reco_phi = reco_phi;
					d_R_min = d_R_itr;
					reco_jet_TM_vector.at(j)=true;
				}
			}

			if(d_R_min < _dR_truth_matching) {

			  //Fill plots
			  h_resp_cent.at(cent_bin)->Fill(truth_pt,(matched_reco_pt - truth_pt)/truth_pt,truth_eta, event_weight);
			  h_resp_cent.at(nCentbins-1)->Fill(truth_pt,(matched_reco_pt - truth_pt)/truth_pt,truth_eta, event_weight);
			  
			  h_truth_jet_cent_matched.at(cent_bin)->Fill(truth_pt,truth_eta, truth_phi, event_weight);
			  h_truth_jet_cent_matched.at(nCentbins-1)->Fill(truth_pt,truth_eta, truth_phi, event_weight);
			}	
		}
	}	
	
	for (int j = 0; j<reco_jet_pt_vector.size(); j++){
		double reco_pt  = reco_jet_pt_vector.at(j);
		double reco_eta = reco_jet_eta_vector.at(j);
		double reco_phi = reco_jet_phi_vector.at(j);
		if (_reco_iso) { if (!reco_jet_isolation_vector.at(j)) continue;}
		if (_isMC){	
			if (!reco_jet_TM_vector.at(j)) {
				h_reco_jet_cent_unmatched.at(cent_bin)->Fill(reco_pt,reco_eta, reco_phi, event_weight);
			}
			else{
				h_reco_jet_cent_matched.at(cent_bin)->Fill(reco_pt,reco_eta, reco_phi, event_weight);
			}
		}	
		h_reco_jet_cent.at(cent_bin)->Fill(reco_pt,reco_eta, reco_phi, event_weight);
	}

	// Clear vectors
	reco_jet_pt_vector.clear();
	reco_jet_eta_vector.clear();
	reco_jet_phi_vector.clear();
	reco_jet_isolation_vector.clear();
	reco_jet_TM_vector.clear();

	truth_jet_pt_vector.clear();
	truth_jet_eta_vector.clear();
	truth_jet_phi_vector.clear();
	truth_jet_isolation_vector.clear();
	
	

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode JetSubstructure :: postExecute ()
{
	// Here you do everything that needs to be done after the main event
	// processing.  This is typically very rare, particularly in user
	// code.  It is mainly used in implementing the NTupleSvc.
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode JetSubstructure :: finalize ()
{
	// This method is the mirror image of initialize(), meaning it gets
	// called after the last event has been processed on the worker node
	// and allows you to finish up any objects you created in
	// initialize() before they are written to disk.  This is actually
	// fairly rare, since this happens separately for each worker node.
	// Most of the time you want to do your post-processing on the
	// submission node after all your histogram outputs have been
	// merged.  This is different from histFinalize() in that it only
	// gets called on worker nodes that processed input events.
	ANA_CHECK_SET_TYPE (EL::StatusCode); // set type of return code you are expecting (add to top of each function once)

	xAOD::TEvent* event = wk()->xaodEvent();
	cout << "Total counts = " << m_eventCounter << endl;
	//cleaning cleaning :)
	if( m_jetCalibration ) delete m_jetCalibration; m_jetCalibration = 0;
	if( m_jetCalibration_val ) delete m_jetCalibration_val; m_jetCalibration_val = 0;


	return EL::StatusCode::SUCCESS;
}



EL::StatusCode JetSubstructure :: histFinalize ()
{
	// This method is the mirror image of histInitialize(), meaning it
	// gets called after the last event has been processed on the worker
	// node and allows you to finish up any objects you created in
	// histInitialize() before they are written to disk.  This is
	// actually fairly rare, since this happens separately for each
	// worker node.  Most of the time you want to do your
	// post-processing on the submission node after all your histogram
	// outputs have been merged.  This is different from finalize() in
	// that it gets called on all worker nodes regardless of whether
	// they processed input events.
	return EL::StatusCode::SUCCESS;
}
