#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

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

#include "JetRec/JetSoftDrop.h"

using namespace std;
using namespace JetHelperTools;

struct jetSubStr {
  Int_t cent;
  float recoPt, recoRawPt, recoEta, recoRcPt, recoSdPt, recoSdmass, recoNrc, nTow, nTowPos;
  float genPt,  genEta, genRcPt,  genSdPt, genSdmass, genNrc, matchDr;
  float weight;
};
jetSubStr myJetSub;
TString myJetSubText = "cent/I:pt/F:rawPt:eta:rcpt:sdpt:sdmass:nrc:ntow:ntowp:genpt:geneta:genrcpt:gensdpt:gensdmass:gennrc:dr:weight";


// this is needed to distribute the algorithm to the workers
ClassImp(JetSubstructure)


void resetSubstr (jetSubStr &jetsub)
{
  jetsub.cent = -1 ;
  jetsub.recoPt= -1;
  jetsub.recoRawPt= -1;
  jetsub.recoEta=-1;
  jetsub.recoRcPt=-1;
  jetsub.recoSdmass=-1;
  jetsub.recoSdPt=-1;
  jetsub.recoNrc=0;
  jetsub.nTow =0;
  jetsub.nTowPos =0;
  
  jetsub.genPt=-1;
  jetsub.genEta=-1;
  jetsub.genRcPt=-1;
  jetsub.genSdmass=-1;
  jetsub.genSdPt=-1;
  jetsub.matchDr=-1;
  jetsub.genNrc=0;
  jetsub.weight = 1;

}

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
	double etaBins[1000], phiBins[1000], respBins[1000], jetPtBins[1000], ratioBins[1000];
	int ratioBinN = 500;
	
	for ( int i= 0 ; i<=ratioBinN ; i++) 
	  ratioBins[i] = 2./ratioBinN * i; 

	//	SetupBinning(0, "eta-fine", etaBins, etaBinsN);
	//	SetupBinning(0, "phi-jet", phiBins, phiBinsN);
	//	SetupBinning(0, "resp", respBins, respBinsN);
	//	SetupBinning(0, "pt-jet-PbPb", jetPtBins, jetPtBinsN);

	//Basic histograms
	h_FCal_Et = new TH1D("h_FCal_Et",";FCal E_{T};N",100,0,5);
	h_FCal_Et->Sumw2();

	h_RejectionHisto = new TH1D("RejectionHisto","RejectionHisto",8,0,8);
	SetRejectionHistogram(h_RejectionHisto);

	h_centrality = new TH1D("Centrality","Centrality",10,0,10);
	h_centrality->Sumw2();

	hET_ETsub = new TH3D("hET_ETsub","hET_ETsub",200,0,200,100,-5,5,200,-50,50);
	hET_ETsub->Sumw2();
	
	
        treeOut = new TTree("tr","new tree");
	treeOut->Branch("jets",&myJetSub,myJetSubText.Data());
	
	//	h_triggercounter = new TH2D("h_triggercounter","h_triggercounter",_nTriggers,0,_nTriggers,2,-0.5,1.5);
	//	SetTrigger_hist(h_triggercounter);
	
	//	h_respR2R4 = new TH2D("h_respR2R4","h_respR2R4",jetPtBinsN, jetPtBins,respBinsN,respBins);
	//	wk()->addOutput(h_respR2R4);

	//	wk()->addOutput (h_FCal_Et);
	wk()->addOutput (h_RejectionHisto);
	wk()->addOutput (h_centrality);
	wk()->addOutput (hET_ETsub);
	wk()->addOutput (treeOut);

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
	// input events./ here
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
		
	cout << "" <<endl;
	//Calibration tool
	
	
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
	ANA_CHECK(m_grl->setProperty("PassThrough", false));	// if true (default) will ignore result of GRL and will just pass all events
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
  float beta = 0 ;
  float z_cut = 0.1;
  

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
	if (_centrality_scheme>1)	  {
	  FCalEt=calos->at(5)->et()*1e-6;
	  cent_bin = GetCentralityBin(_centrality_scheme, FCalEt, isHIJING);
	  
	  event_weight_fcal = jetcorr->GetFCalWeight(FCalEt);
	  
	  if (_isMC && isHIJING) event_weight_fcal = 1;
	}
	if (cent_bin < 0) {
	  h_RejectionHisto->Fill(1.5);
	  keep = false;
	}
	// GRL
	if(!_isMC) {
	  if(!m_grl->passRunLB(*eventInfo)) {
	    h_RejectionHisto->Fill(2.5);
	    keep = false;
	  }}
	
	//Vertex requirement
	const xAOD::VertexContainer * vertices = 0;
	if ( !event->retrieve( vertices, "PrimaryVertices" ).isSuccess() ) {
	  Error("execute()", "Failed to retrieve VertexContainer container. Exiting." );
	  return EL::StatusCode::FAILURE;
	}
	if(vertices->size()<2) {
	  h_RejectionHisto->Fill(3.5);
	  keep = false;
	}
	
	//DAQ errors
	if(!_isMC) {
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
	  m_zdcTools->reprocessZdc();	  // ZDC
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
	if (_isMC==0) {
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
	


	vector <double> vpt_reco;
	vector <double> vptRaw_reco;  
	vector <double> veta_reco;
	vector <double> vphi_reco;
	vector <double> vptRc_reco;
	vector <double> vSdmass_reco;
	vector <double> vSdpt_reco;
	vector <int> vNrc_reco;
	vector <int> vNtow_reco;
	vector <int> vNtowp_reco;

	vector <double> vpt_gen;
	vector <double> veta_gen;
	vector <double> vphi_gen;
	vector <double> vptRc_gen;
	vector <double> vSdmass_gen;
	vector <double> vSdpt_gen;
	vector <int> vNrc_gen;

	
	/////////////   Reco jets /////////////////////////////////////////
	
	xAOD::TStore *store = new xAOD::TStore; //For calibration
	const xAOD::JetContainer* reco_jets = 0;
	ANA_CHECK(event->retrieve( reco_jets, _reco_jet_collection.c_str() ));

	xAOD::JetContainer::const_iterator jet_itr = reco_jets->begin();
	xAOD::JetContainer::const_iterator jet_end = reco_jets->end();
	for( ; jet_itr != jet_end; ++jet_itr ) {
	  
	  xAOD::Jet theRecoJet;
	  theRecoJet.makePrivateStore( **jet_itr );

	  const xAOD::JetFourMom_t jet_4momUnCal = theRecoJet.jetP4("JetSubtractedScaleMomentum"); // uncalib
	  const xAOD::JetFourMom_t jet_4momCalib = theRecoJet.jetP4();   // CALIBRATED!!! 
	  
	  double jet_pt  = jet_4momCalib.pt() * 0.001 ;
	  double jet_eta = jet_4momCalib.eta();
	  double jet_phi = jet_4momCalib.phi();
	  double jet_ptRaw = jet_4momUnCal.pt() * 0.001;

	  if (jet_pt < _pTjetCut) continue;
	  
	  const xAOD::JetConstituentVector constituents_tmp = (*jet_itr)->getConstituents();
	  //	  cout <<" number of RECO constituent = " << (*jet_itr)->numConstituents() << endl;
	  //	  cout <<" size of constituents = " << constituents_tmp.size() << endl;
	  xAOD::JetConstituentVector::iterator itCnst = constituents_tmp.begin();
	  xAOD::JetConstituentVector::iterator itCnst_E = constituents_tmp.end();
	  vector<fastjet::PseudoJet>  nonZeroConsts;
	  cout <<" Jet pT = " << jet_pt << ", raw pT = "<<jet_ptRaw<<endl;
	  cout <<" Constituent's pT: " << endl;
	  for( ; itCnst != itCnst_E; ++itCnst ) {
	    float thePt = (*itCnst)->pt();
	    fastjet::PseudoJet thisConst = fastjet::PseudoJet( (*itCnst)->Px(), (*itCnst)->Py(), (*itCnst)->Pz(), (*itCnst)->E() );
	    nonZeroConsts.push_back(thisConst);
	    cout << (*itCnst)->Pt()*0.001 <<" ("<< (*itCnst)->E() <<"), " ;
	  }
	  cout << endl;
	  // recluster by A/C
	  
	  fastjet::JetDefinition jetDefRe(fastjet::cambridge_algorithm, _ReclusterRadius);
	  fastjet::ClusterSequence csRe(nonZeroConsts, jetDefRe);
	  vector<fastjet::PseudoJet> jetsRe = fastjet::sorted_by_pt(csRe.inclusive_jets()); // return a vector of jets sorted into decreasing energy

	  if ( jetsRe.size() > 0  ) { 
	    cout << "Number of AcJets = " << jetsRe.size() << endl;
	    for ( int ic=0; ic< jetsRe.size() ; ic++) { 
	      cout <<"   "<<ic<<"th jet (pt,eta,phi) = " <<  jetsRe[ic].pt() *0.001<<", "<< jetsRe[ic].eta()<<", " << jetsRe[ic].phi()<<")"<<endl;
	    }
	  }
	  cout <<endl << endl;

	  fastjet::contrib::SoftDrop sd(beta, z_cut);
	  //	  fastjet::contrib::SoftDrop sd(beta, z_cut);
	  
	  
	  int theNrc = 0;
	  double thePtrc = 0;
	  double thesdpt = 0 ;
	  double thesdm = 0;
	  if ( jetsRe.size() > 0 )   {
	    theNrc = jetsRe.size();
	    thePtrc = jetsRe[0].pt() * 0.001;
	    fastjet::PseudoJet sd_jet = sd(jetsRe[0]);
	    thesdpt = sd_jet.pt() * 0.001; 
	    thesdm =  sd_jet.m() * 0.001; 
	  }  
	  vpt_reco.push_back(jet_pt);
	  vptRaw_reco.push_back(jet_ptRaw);
	  veta_reco.push_back(jet_eta);
	  vphi_reco.push_back(jet_phi);
	  vptRc_reco.push_back(thePtrc);
	  vNrc_reco.push_back(theNrc);
	  vSdmass_reco.push_back(thesdm);
	  vSdpt_reco.push_back(thesdpt);
	  vNtow_reco.push_back(nonZeroConsts.size() );
	  vNtowp_reco.push_back(nonZeroConsts.size() );

	  jetsRe.clear();
	  nonZeroConsts.clear();
	}
	
	
	//**Get Truth jets ***
	double event_weight = 1;
	double max_pt = 0;
	bool useReAntiKt = true; 
	
	if (_isMC){
	  if ( !useReAntiKt )  {
	    cout << " Don't use useReAntiKt=0 option yet.." << endl;
	    /*    const xAOD::JetContainer* truth_jets = 0;
		  ANA_CHECK(event->retrieve(truth_jets, _truth_jet_collection.c_str() ));
		  
		  xAOD::JetContainer::const_iterator truth_jet_itr = truth_jets->begin();
		  xAOD::JetContainer::const_iterator truth_jet_end = truth_jets->end();
	    
		  for( ; truth_jet_itr != truth_jet_end; ++truth_jet_itr ) {
		  xAOD::JetFourMom_t jet_truth_4mom = (*truth_jet_itr)->jetP4();
		  double pt     = (jet_truth_4mom.pt() * 0.001 );
		  double eta    = (jet_truth_4mom.eta());
		  double phi    = (jet_truth_4mom.phi());
		  if (pt > max_pt)  { 
		  event_weight = jetcorr->GetJetWeight(pt, eta, phi);
		  max_pt = pt;
		  //			    if (_isMC && isHIJING) event_weight = 1;
		  }
	      if (pt < _truthpTjetCut) continue;
	    */ 
	  }
	  else {  // if useReAntiKt
	    ////////// On-the-fly jet finder ///////////////////////////////////////////////
	    std::vector<fastjet::PseudoJet> inputConst;
	    const xAOD::TruthParticleContainer * particles = 0;
	    ANA_CHECK(event->retrieve( particles, "TruthParticles"));
	    xAOD::TruthParticleContainer::const_iterator truth_itr = particles->begin();
	    xAOD::TruthParticleContainer::const_iterator truth_end = particles->end();
	    for( ; truth_itr!=truth_end; ++truth_itr){
	      int ty=TrackHelperTools::getTypeTruth((*truth_itr)->barcode(),(*truth_itr)->pdgId(),(*truth_itr)->status());
	      if(ty!=1 && ty!=5) continue;
	      inputConst.push_back( (*truth_itr)->p4() );
	    }
	    
	    fastjet::JetAlgorithm algo = fastjet::antikt_algorithm;
	    fastjet::JetDefinition jetDef(algo, _ReclusterRadius);
	    fastjet::ClusterSequence cs(inputConst, jetDef);
	    vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
	    // Now, new genJet is ready!
	    
	    //////////////// event weight /////////////////////////////////////
	    max_pt = 0;
	    for (unsigned i = 0; i < jets.size(); i++) {
	      double jet_pt = jets[i].pt()*0.001;
	      double jet_eta = jets[i].eta();
	      double jet_phi = jets[i].phi();
	      if (jet_pt > max_pt) {
		max_pt = jet_pt;
		event_weight = jetcorr->GetJetWeight( max_pt, jet_eta, jet_phi);
	      }
	    }
	    event_weight = event_weight*event_weight_fcal; //event weight is only set if MC.
	    //////////////////////////////////////////////////////////////////

	    for (unsigned i = 0; i < jets.size(); i++) {
	      double jet_pt = jets[i].pt()*0.001;
              if (jet_pt < _truthpTjetCut) continue;
	      
	      double jet_eta = jets[i].pseudorapidity();
	      double jet_phi = jets[i].phi();
	      if (jet_phi> TMath::Pi()) jet_phi = jet_phi - 2*TMath::Pi();

	      // recluster by A/C
	      vector<fastjet::PseudoJet > constiRe = jets[i].constituents();
	      fastjet::JetDefinition jetDefRe(fastjet::cambridge_algorithm, _ReclusterRadius);
	      fastjet::ClusterSequence csRe(constiRe, jetDefRe);
              vector<fastjet::PseudoJet> jetsRe = fastjet::sorted_by_pt(csRe.inclusive_jets()); // return a vector of jets sorted into decreasing energy

    	      fastjet::contrib::SoftDrop sd(beta, z_cut);
	      int theNrc = 0;
	      double thesdpt = 0 ;
	      double thesdm = 0;
	      double thePtrc = 0;
	      if ( jetsRe.size() > 0 )   {
		theNrc = jetsRe.size();
		thePtrc = jetsRe[0].pt() * 0.001;
		fastjet::PseudoJet sd_jet = sd(jetsRe[0]);
		thesdpt = sd_jet.pt() * 0.001; 
		thesdm =  sd_jet.m() * 0.001;
	      }  
	      vpt_gen.push_back(jet_pt);
	      veta_gen.push_back(jet_eta);
	      vphi_gen.push_back(jet_phi);
	      vptRc_gen.push_back(thePtrc);
	      vNrc_gen.push_back(theNrc);
	      vSdmass_gen.push_back(thesdm);
	      vSdpt_gen.push_back(thesdpt);
	      
	      inputConst.clear();
	      jetsRe.clear();
	    }
	  }
	}
       

	// back to reco loop 
	for ( int ri = 0 ; ri< vpt_reco.size() ; ri++) { 
	  
	  int matchId = -1;
	  double drMin = 0.2;    // dR cut
	  for (int gi = 0; gi<vpt_gen.size() ; gi++) { 
	    double dr_itr = DeltaR(vphi_reco[ri], veta_reco[ri], vphi_gen[gi], veta_gen[gi]);
	    //	    cout <<"dR_itr = " << dr_itr << endl;
	    if ( dr_itr < drMin ) { 
	      matchId = gi;  
	      drMin = dr_itr;
	      //	      cout <<"found! dRmin = " << drMin << endl;
	    }
	  }
	  
	  resetSubstr(myJetSub);
	  //	  cout << " myJetsub is reset" << endl;
	  //	  cout << " myJetSub.matchDr  = " << myJetSub.matchDr << endl;
	  myJetSub.cent = cent_bin; 
	  myJetSub.recoPt = vpt_reco[ri];
	  myJetSub.recoRawPt = vptRaw_reco[ri];
	  myJetSub.recoEta = veta_reco[ri];
	  myJetSub.recoRcPt = vptRc_reco[ri];
	  myJetSub.recoNrc = vNrc_reco[ri];
	  myJetSub.nTow = vNtow_reco[ri];
	  myJetSub.nTowPos = vNtowp_reco[ri];
	  myJetSub.weight = event_weight;
	  myJetSub.recoSdPt = vSdpt_reco[ri];
	  myJetSub.recoSdmass = vSdmass_reco[ri];

	  if (matchId != -1) {
	    myJetSub.matchDr = drMin;
	    myJetSub.genPt = vpt_gen[matchId];
	    myJetSub.genEta = veta_gen[matchId];
	    myJetSub.genRcPt = vptRc_gen[matchId];
	    myJetSub.genSdmass = vSdmass_gen[matchId];
	    myJetSub.genSdPt = vSdpt_gen[matchId];
	    myJetSub.genNrc   =  vNrc_gen[matchId];
	  }
	  
	  treeOut->Fill();
	}
	
	
	

	
	
	// Clear vectors
	store->clear();
	delete store;
	
	vpt_reco.clear();
	vptRaw_reco.clear();
	veta_reco.clear();
	vphi_reco.clear();
	vptRc_reco.clear();
	vSdmass_reco.clear();
	vNrc_reco.clear();
	vNtow_reco.clear();
	vNtowp_reco.clear();
	vSdpt_reco.clear();

	vpt_gen.clear();
	veta_gen.clear();
	vphi_gen.clear();
	vptRc_gen.clear();
	vSdmass_gen.clear();
	vSdpt_gen.clear();
	vNrc_gen.clear();	

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

