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


#define EL_RETURN_CHECK( CONTEXT, EXP )			  \
  do {                                                    \
    if( ! EXP.isSuccess() ) {				  \
      Error( CONTEXT,					  \
	     XAOD_MESSAGE( "Failed to execute: %s" ),	  \
	     #EXP );					  \
      return EL::StatusCode::FAILURE;			  \
    }							  \
  } while( false )


struct jetSubStr {
  Int_t cent;
  float weight;
  float recoPt, recoRawPt, recoEta, recoRcPt, recoSdPt, recoSdMass, recoSdZ, recoSdTheta;
  float reChSdPt, reChSdMass, reChSdZ, reChSdTheta;
  float matchDr, genPt,  genEta, genRcPt,  genSdPt, genSdMass, genSdz, genSdtheta;
  float genChSdPt, genChSdMass, genChSdZ, genChSdTheta; 
};
jetSubStr myJetSub;

TString evtText = "cent/I:weight/F:";
TString recoText = "pt:rawPt:eta:rcPt:sdPt:sdMass:zg:theta:";
TString reChText = "chSdPt:chSdMass:chZg:chTheta:";
TString genText = "dr:genPt:genEta:genRcPt:genSdPt:genSdMass:genSdZg:genSdtheta:";
TString genChText = "genChSdPt:genChSdMass:genChSdZg:genChSdTheta";

TString myJetSubText = evtText+recoText+reChText+genText+genChText;

// this is needed to distribute the algorithm to the workers
ClassImp(JetSubstructure)


void resetSubstr (jetSubStr &jetsub)
{
  jetsub.cent = -1 ;
  jetsub.recoPt= -1;
  jetsub.recoRawPt= -1;
  jetsub.recoEta=-1;
  jetsub.recoRcPt=-1;
  jetsub.recoSdPt=-1;
  jetsub.recoSdMass=-1;
  jetsub.recoSdZ = 100;
  jetsub.recoSdTheta = 100;

  jetsub.genPt=-1;
  jetsub.genEta=-1;
  jetsub.genRcPt=-1;
  jetsub.genSdPt=-1;
  jetsub.genSdMass=-1;
  jetsub.matchDr=-1;
  jetsub.weight = 1;
  jetsub.genSdz = 100;
  jetsub.genSdtheta = 100;

}

int getMaxPtIndex ( vector<fastjet::PseudoJet>& jets ) {
  double max_pt = 0;
  int maxIndex = -1;
  for (unsigned i = 0; i < jets.size(); i++) {
    double jet_pt = jets[i].pt()*0.001;
    double jet_eta = jets[i].eta();
    double jet_phi = PhiInPI ( jets[i].phi() );  
    if (jet_pt > max_pt) {
      max_pt = jet_pt;
      maxIndex = i ; 
    }
  }
  return maxIndex;
}

void showLadder (fastjet::PseudoJet akJet, fastjet::PseudoJet camJet, fastjet::PseudoJet sd_jet) { 
  cout << "*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*" << endl;
  cout << "[ pT of antikT -> Cambridge -> SoftDrop  =  " << akJet.pt()*0.001 << " -> " << camJet.pt()*0.001 <<" -> "<< sd_jet.pt()*0.001 <<"]" << endl;
  cout << " Number of C/A constituents: " << camJet.constituents().size()  << ",  SoftDrop constituents: "<< sd_jet.constituents().size() << endl;
  cout << " Full consticuents lists (pt, eta, phi)" <<endl;
  vector<fastjet::PseudoJet> camConst = camJet.constituents() ;
  for ( int ic = 0 ; ic< camConst.size() ; ic++){
    cout << "    ( " << camConst[ic].pt() *0.001 << ", " << camConst[ic].eta() << ", "<< camConst[ic].phi() << ") "<< endl;
  }
  cout << "*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*" << endl;
  camConst.clear();
}

void drawTreeHistory ( fastjet::ClusterSequence &clustSeq ) {
  // C/A tree of hierarchy
  vector<fastjet::PseudoJet> vJetsCA = clustSeq.jets();
  vector<fastjet::ClusterSequence::history_element> vHistCA = clustSeq.history();
  int maxJetInd = -1;
  float maxPt = -1;
  for ( int ij= 0 ; ij< vJetsCA.size() ; ij++) {
    fastjet::PseudoJet childJet;
    if ( vJetsCA[ij].pt() > maxPt ) {
      maxPt = vJetsCA[ij].pt();
      maxJetInd = ij;
    }
  }
  if ( maxJetInd == -1 )    {    
    cout << "STRANGE!!!!  maxJetInd is -1.  drawTreeHistory is broken" << endl;
    return; 
  }
  
  const int nSteps = 30;
  int indGen[nSteps][200];
  int nGen[nSteps];
  for ( int step = 0 ; step<nSteps ; step++)
    nGen[step] = 0;

  nGen[0] = 1;
  indGen[0][0] = maxJetInd;

  cout << endl << "~~~~Start of tree for jet pT, eta, phi: " << vJetsCA[maxJetInd].pt()*0.001 << ", "<< vJetsCA[maxJetInd].eta()<<", "<< vJetsCA[maxJetInd].phi()  << "~~~~" << endl;
  for ( int step = 0 ; step< nSteps ; step++) {
    cout << "  step "<<step << endl;
    for ( int ijet = 0 ; ijet< nGen[step] ; ijet++) {
      int theIndex = indGen[step][ijet] ;
      int par1 = vHistCA[ theIndex ].parent1;
      int par2 = vHistCA[ theIndex ].parent2;
      if ( par1 < 0 )  {
	cout <<" par1 < 0 !!"<<endl;
	continue;
      }
      float dij = vHistCA[ theIndex ].dij;
      //              cout << "par1, par2 indice = " << par1<<", "<<par2<<endl;

      cout << "     Child jet pt=" << vJetsCA[theIndex].pt() * 0.001 << " GeV,  dij=" << dij<<endl;
      cout << "     Parent 1: jet pt= " << vJetsCA[par1].pt() * 0.001<< " GeV " << endl;
      cout << "     Parent 2: jet pt= " << vJetsCA[par2].pt() * 0.001 << " GeV "<< endl;

      indGen[step+1][ nGen[step+1] ] = par1;
      indGen[step+1][ nGen[step+1] +1 ] = par2;
      nGen[step+1] = nGen[step+1]+ 2;
    }
  }

  vJetsCA.clear();
  vHistCA.clear();
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl ;
}


void logSubJets( fastjet::PseudoJet &parent1, fastjet::PseudoJet &parent2 ) { 
  vector<fastjet::PseudoJet> constSub1 =  parent1.constituents() ;
  vector<fastjet::PseudoJet> constSub2 =  parent2.constituents() ;
  cout << "DR of subjets: " << DeltaR( parent1.phi(), parent1.rapidity(), parent2.phi(), parent2.rapidity() ) << endl;
  cout << "   * Subjet0: nConst, pt, eta, phi = [" <<  parent1.constituents().size() << ",   " <<parent1.pt()*0.001 << ", "<< parent1.eta() << ", "<< parent1.phi() << "] "<< endl;
  cout << "       constituents lists (pt,eta,phi)" <<endl;
  for ( int ic = 0 ; ic< constSub1.size() ; ic++){
    cout << "      ( " << constSub1[ic].pt() *0.001 << ", " << constSub1[ic].eta() << ", "<< constSub1[ic].phi() << ") " << endl;
  }
  cout << "   * Subjet1: nConst, pt, eta, phi = [" <<  parent2.constituents().size() << ",   " <<parent2.pt()*0.001 << ", "<< parent2.eta() << ", "<< parent2.phi() << "] "<< endl;
  cout << "       constituents lists (pt,eta,phi)" <<endl;
  for ( int ic = 0 ; ic< constSub2.size() ; ic++){
    cout << "      ( " << constSub2[ic].pt() *0.001<< ", " << constSub2[ic].eta() << ", "<< constSub2[ic].phi() << ")" << endl;
  }
  constSub1.clear();
  constSub2.clear();
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
       
	//	jetTree = new TTree("tr2","jet branching tree");
	//	jetTree->Branch("bran","");

	
	hGenNcam = new TH2D("hGenNcam","; p_{T} GeV/c ; Number of C/A clusters",100,0,1000,3,-0.5,2.5);
	hRecoNcam = (TH2D*)hGenNcam->Clone("hRecoNcam");
	hGenSdStat = new TH2D("hGenSdStat","; p_{T} GeV/c ; 0:fail, 1;pass",100,0,1000,2,-0.5,1.5);
	hRecoSdStat = (TH2D*)hGenSdStat->Clone("hRecoSdStat");
	
	
	
	wk()->addOutput (h_RejectionHisto);
	wk()->addOutput (h_centrality);
	wk()->addOutput (hET_ETsub);
	wk()->addOutput (treeOut);
	wk()->addOutput (hGenNcam);
	wk()->addOutput (hRecoNcam);
	wk()->addOutput (hGenSdStat);
	wk()->addOutput (hRecoSdStat);
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

  if (m_eventCounter%statSize==0)  std::cout << "Event: " << m_eventCounter << std::endl;
  
  
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
    m_is_pileup = m_hiPileup->is_pileup( *calos, *zdcMod); // SAVE pileup Decision HERE 0 = NO pileup, 1 = pileup
  }
  else m_is_pileup = (FCalEt > 4.8); //Remove pileup in MC
  if (m_is_pileup){
    h_RejectionHisto->Fill(6.5);
    keep = false;
  }
  
  if (!keep) return EL::StatusCode::SUCCESS; // go to the next event
  h_RejectionHisto->Fill(7.5);
  
  h_FCal_Et->Fill(FCalEt, event_weight_fcal); //filled here to get proper event weight
  h_centrality->Fill(cent_bin,event_weight_fcal);
  
  cout << "Centrality of this event = " << cent_bin << endl;
  vector <double> vpt_reco;
  vector <double> vptRaw_reco;  
  vector <double> veta_reco;
  vector <double> vphi_reco;
  vector <double> vptRc_reco;
  vector <double> vSdmass_reco;
  vector <double> vSdpt_reco;
  vector <float> vSdz_reco;
  vector <float> vSdtheta_reco;
  
  vector <double> vpt_gen;
  vector <double> veta_gen;
  vector <double> vphi_gen;
  vector <double> vptRc_gen;
  vector <double> vSdmass_gen;
  vector <double> vSdpt_gen;
  vector <float> vSdz_gen;
  vector <float> vSdtheta_gen;
  vector <float> vGenChSdPt;
  vector <float> vGenChSdMass  ;
  vector <float> vGenChSdZ  ;
  vector <float> vGenChSdTheta  ;




  // algorithm definition 
  fastjet::JetDefinition jetDefCam(fastjet::cambridge_algorithm, _ReclusterRadius);
  fastjet::JetDefinition jetDefAk(fastjet::antikt_algorithm, _ReclusterRadius);
  float beta = 0 ;
  float z_cut = 0.1;
  fastjet::contrib::SoftDrop softdropper(beta, z_cut);
	
	
  ///////////// tracks ////////////////////////////////////////
  const xAOD::TrackParticleContainer* recoTracks = 0;
  EL_RETURN_CHECK("execute",event->retrieve( recoTracks, "InDetTrackParticles"));
  for (const auto& trk : *recoTracks) {
    float pt = trk->pt()/1000.;
    float eta = trk->eta();
    float phi = trk->phi();
    cout << " track pt, eta, phi = " << pt <<", "<<eta<<", "<<phi<<endl;
  }
  
  
  
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
    
    fastjet::PseudoJet unCaliFourVec = fastjet::PseudoJet ( jet_4momUnCal.px(), jet_4momUnCal.py(), jet_4momUnCal.pz(), jet_4momUnCal.energy() );
    double jet_pt  = jet_4momCalib.pt() * 0.001 ;
    double jet_eta = jet_4momCalib.eta();
    double jet_phi = PhiInPI ( jet_4momCalib.phi() );
    double jet_ptRaw = jet_4momUnCal.pt() * 0.001;
    if (jet_pt < _pTjetCut)          continue;
    if (fabs(jet_eta) > _etaJetCut)  continue;
    
    if ( _saveLog) cout << "*~*~*~*~*~*~ RECO ~*~*~*~*~*~*" << endl << "  Anti-kT  jet [pt, eta, phi] : " << jet_pt <<", "<<jet_eta<<", "<<jet_phi<<endl << " Raw pT: " << jet_ptRaw << " GeV" << endl;
    
    const xAOD::JetConstituentVector recoConsts = (*jet_itr)->getConstituents();
    xAOD::JetConstituentVector::iterator itCnst   = recoConsts.begin();
    xAOD::JetConstituentVector::iterator itCnst_E = recoConsts.end();
    vector<fastjet::PseudoJet>  nonZeroConsts;
    vector<fastjet::PseudoJet>  toBeSubtracted; // reverse the pT only and subtract
    vector<bool>  nFlag; // reverse the pT only and subtract
    
    double ghostE = 0.000001;
    for( ; itCnst != itCnst_E; ++itCnst ) {
      double theEta = (*itCnst)->Eta() ; 
      double thePhi = PhiInPI ( (*itCnst)->Phi() ) ;
      
      const fastjet::PseudoJet thisConst = fastjet::PseudoJet( (*itCnst)->Px(), (*itCnst)->Py(), (*itCnst)->Pz(), (*itCnst)->E() );
      
      if ( (*itCnst)->pt() > 0 ) { // normal tower 
	nonZeroConsts.push_back(thisConst);
	toBeSubtracted.push_back ( fastjet::PseudoJet ( 0.000001, 0.000001, 0.000001, 0.000002 )) ;
	nFlag.push_back (false);
      }
      else {  // negative tower 
	double gpx = ghostE/cosh(theEta) * cos(thePhi) ;
	double gpy = ghostE/cosh(theEta) * sin(thePhi) ;
	double gpz = ghostE*tanh(theEta) ;
	double posPt = - (*itCnst)->pt() ;
	double posNorm = posPt * cosh(theEta);    // p = pT * cosh(
		
	nonZeroConsts.push_back( fastjet::PseudoJet ( gpx,gpy,gpz,ghostE )  );
	toBeSubtracted.push_back ( fastjet::PseudoJet ( gpx*posNorm/ghostE, gpy*posNorm/ghostE, gpz*posNorm/ghostE, posNorm) );
	nFlag.push_back(true);
      }
    }

    // Cambridge reclustering 
    if ( _saveLog) cout << "Reco reclustering starts!" << endl;
    fastjet::ClusterSequence recoCamSeq(nonZeroConsts, jetDefCam);
    vector<fastjet::PseudoJet> cambridgeJet = recoCamSeq.inclusive_jets();
    vector<int> pIndex = recoCamSeq.particle_jet_indices( cambridgeJet);
    if ( _saveLog)  {
      drawTreeHistory(recoCamSeq) ;
    }

    if ( cambridgeJet.size() > 0 ) { 
      for ( int ij= 0 ; ij < cambridgeJet.size() ; ij++) {
	if ( _saveLog)	cout<< ij<<"th Cambridge jet pT, eta phi: "<<cambridgeJet[ij].pt() *0.001  << ", "<< cambridgeJet[ij].eta() << ", " << cambridgeJet[ij].phi()<< endl;
	double subPx = cambridgeJet[ij].px() ; 
	double subPy = cambridgeJet[ij].py() ; 
	double subPz = cambridgeJet[ij].pz() ; 
	double subE  = cambridgeJet[ij].E() ; 
	
	//	      cout << "Missing constituent pt,eta,phi : " << endl;
	int towerCount=0;
	int towerCountN=0;
	for ( int ip = 0 ; ip < pIndex.size() ; ip++) {
	  if   (pIndex[ip]==ij) {   
	    towerCount++;
	    if ( nFlag[ip] ) {
	      towerCountN++;
	      subPx = subPx - toBeSubtracted[ip].px() ;
	      subPy = subPy - toBeSubtracted[ip].py() ;
	      subPz = subPz - toBeSubtracted[ip].pz() ;
	      subE = subE - toBeSubtracted[ip].E() ;
	    }
	  }
	}
	cambridgeJet[ij].reset_momentum( subPx, subPy, subPz, subE ) ;
	if ( _saveLog)  {
	  cout <<"   after Negative correction: " << cambridgeJet[ij].pt() *0.001  << " GeV " << endl;
	  cout <<"   Number of negative towers / (total): " << towerCountN << "/ ("<<towerCount<<")"<<endl;
	}
      }
    }
    
    vector<fastjet::PseudoJet> corrRecCamJets = fastjet::sorted_by_pt(cambridgeJet); // return a vector of jets sorted into decreasing energy
    
    
    // Yongsun:  http://acode-browser2.usatlas.bnl.gov/lxr-rel21/source/atlas/Reconstruction/Jet/JetRec/Root/JetSoftDrop.cxx line 67.
    
    double thePtrc = 0;
    double thesdpt = 0 ;
    double thesdm = 0;
    float thesdtheta = 100;
    float thesdz = -1;
 
    if ( corrRecCamJets.size() > 0 )   {
      thePtrc = corrRecCamJets[0].pt() * 0.001;
      fastjet::PseudoJet sd_jet = softdropper(corrRecCamJets[0]);
      thesdpt = sd_jet.pt() * 0.001; 
      thesdm =  sd_jet.m() * 0.001; 
      
      cout << "RECO JET softdrop:" << endl;
      if ( _saveLog) 
	showLadder ( unCaliFourVec, corrRecCamJets[0], sd_jet );
      
      fastjet::PseudoJet parent1, parent2;
      if (  sd_jet.has_parents(parent1, parent2) ) {
	thesdtheta =  DeltaR( parent1.phi(), parent1.rapidity(), parent2.phi(), parent2.rapidity() ) ;
	thesdz  = std::min( parent1.pt(),  parent2.pt() ) / ( parent1.pt() +  parent2.pt() ) ;
	if ( _saveLog)   {
	  logSubJets( parent1, parent2 ) ;
	}
      }
      else  {
	cout << "This softdrop jet is not divided!" << endl; 
      }
      
    }


    vpt_reco.push_back(jet_pt);
    vptRaw_reco.push_back(jet_ptRaw);
    veta_reco.push_back(jet_eta);
    vphi_reco.push_back(jet_phi);
    vptRc_reco.push_back(thePtrc);
    vSdmass_reco.push_back(thesdm);
    vSdpt_reco.push_back(thesdpt);
    vSdz_reco.push_back(thesdz);
    vSdtheta_reco.push_back(thesdtheta);
    
    corrRecCamJets.clear();
    nonZeroConsts.clear();
    toBeSubtracted.clear(); 
  }
  // sfsg
  
  //**Get Truth jets ***
  double event_weight = 1;
  bool useReAntiKt = true; 
  

  // truth jet constituent test
  /*  
      cout << " Don't use useReAntiKt=0 option yet.." << endl;
      const xAOD::JetContainer* truth_jets = 0;
      ANA_CHECK(event->retrieve(truth_jets, _truth_jet_collection.c_str() ));
      
      xAOD::JetContainer::const_iterator truth_jet_itr = truth_jets->begin();
      xAOD::JetContainer::const_iterator truth_jet_end = truth_jets->end();
      
      for( ; truth_jet_itr != truth_jet_end; ++truth_jet_itr ) {
      xAOD::JetFourMom_t jet_truth_4mom = (*truth_jet_itr)->jetP4();
      double pt     = (jet_truth_4mom.pt() * 0.001 );
      double eta    = (jet_truth_4mom.eta());
      double phi    = (jet_truth_4mom.phi());
      
      xAOD::JetConstituentVector vec1 = (*truth_jet_itr)->getConstituents();
      cout << " n const = " << vec1.size() << endl; 
      cout << "pt0 = " << vec1[0]->pt() << endl;
      cout << "pt1 = " << vec1[1]->pt() << endl;
      }
      /////// end of test  
      
      */
  
    
  if (_isMC){
    if ( useReAntiKt )  {
      ////////// On-the-fly jet finder ///////////////////////////////////////////////
      cout << endl << "///////////////////////////////////" << endl;
      cout << " MC information " << endl;

      const xAOD::TruthParticleContainer * genParCont = 0;
      ANA_CHECK(event->retrieve( genParCont, "TruthParticles"));
      std::vector<fastjet::PseudoJet> truthParticles;
      std::vector<fastjet::PseudoJet> truthCharges;

      for( xAOD::TruthParticleContainer::const_iterator truth_itr = genParCont->begin() ; truth_itr!= genParCont->end() ; ++truth_itr) {
	int thebc = (*truth_itr)->barcode(); 
	if ( (thebc > 0) && ( thebc < 200000 ) )  {
	  truthParticles.push_back( (*truth_itr)->p4() );
	  if ( fabs((*truth_itr)->charge()) > 0 ) 
	    truthCharges.push_back( (*truth_itr)->p4() ) ;
	} 
      }
      
      
      fastjet::ClusterSequence cs(truthParticles, jetDefAk);
      vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
      // Now, new genJet is ready!
      
      //////////////// event weight /////////////////////////////////////
      int maxPtId = getMaxPtIndex ( jets ) ;
      event_weight = 0;  
      if (maxPtId > -1)  
	event_weight = jetcorr->GetJetWeight( jets[maxPtId].pt(), jets[maxPtId].eta(), jets[maxPtId].phi() );
      event_weight = event_weight*event_weight_fcal; //event weight is only set if MC.
      if ( _saveLog ) cout << " Max Truth pT = " << jets[maxPtId].pt()*0.001 << " GeV " << endl ;
      //////////////////////////////////////////////////////////////////
      
      for (unsigned i = 0; i < jets.size(); i++) {  // MC anti-kT jets
	double jet_pt = jets[i].pt()*0.001;
	double jet_eta = jets[i].pseudorapidity();  
	double jet_rap = jets[i].rapidity();
	double jet_phi = PhiInPI ( jets[i].phi() ) ;
	if (jet_pt < _truthpTjetCut) continue;
	if ( fabs(jet_eta) > _etaJetCut ) continue;
	
	// Truth recluster by C/A
	vector<fastjet::PseudoJet > akConsts = jets[i].constituents();
	fastjet::ClusterSequence reCam(akConsts, jetDefCam);
	vector<fastjet::PseudoJet> camJets = fastjet::sorted_by_pt(reCam.inclusive_jets()); // return a vector of jets sorted into decreasing energy
	
	if ( _saveLog)  {   	  
	  cout << "Truth Cambridge jet " << endl ;
	  drawTreeHistory(reCam) ;  
	}
	// Charged SoftDrop
	// First find the charged particle SD
	
	double thesdpt = 0 ;
	double thesdm = 0;
	double thePtrc = 0;
	float thesdtheta = 100;
	float thesdz = -1;
	
	if ( camJets.size() > 0 )   {
	  thePtrc = camJets[0].pt() * 0.001;
	  fastjet::PseudoJet sd_jet = softdropper(camJets[0]);
	  thesdm =  sd_jet.m() * 0.001;
	  thesdpt = sd_jet.pt() * 0.001;
	  if ( _saveLog) { 
	    cout << "GEN JET softdrop:" << endl;
	    showLadder ( jets[i], camJets[0], sd_jet );  
	  }
	  
	  fastjet::PseudoJet parent1, parent2;
	  if (  sd_jet.has_parents(parent1, parent2) ) {
	    thesdtheta =  DeltaR( parent1.phi(), parent1.rapidity(), parent2.phi(), parent2.rapidity() ) ;
	    thesdz  = std::min( parent1.pt(),  parent2.pt() ) / ( parent1.pt() +  parent2.pt() ) ;
	    
	    if ( _saveLog)   {
	      logSubJets( parent1, parent2 ) ;
	    }
	  }
	  else {
	    cout << "This softdrop jet is not divided!" << endl; 
	  }
	}
	
	vpt_gen.push_back(jet_pt);
	veta_gen.push_back(jet_eta);
	vphi_gen.push_back(jet_phi);
	vptRc_gen.push_back(thePtrc);
	vSdmass_gen.push_back(thesdm);
	vSdpt_gen.push_back(thesdpt);
	vSdz_gen.push_back(thesdz);
	vSdtheta_gen.push_back(thesdtheta);
	camJets.clear();
	
	
	
	// Charged Particle reclustering
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "Charged particle clustering starts" << endl;
        vector<fastjet::PseudoJet> chargeConsts;
        for ( int ic=0; ic< truthCharges.size() ; ic++) {
          if ( DeltaR ( jet_phi, jet_rap, truthCharges[ic].phi(), truthCharges[ic].rapidity() ) <  _ReclusterRadius ) {
            chargeConsts.push_back( truthCharges[ic]) ;
          }
        }
	fastjet::ClusterSequence reChCam(chargeConsts, jetDefCam);
        vector<fastjet::PseudoJet> camChJets = fastjet::sorted_by_pt(reChCam.inclusive_jets()); 
	if ( _saveLog)  {
	  cout << "Charged Cambridge jets" << endl;
	  drawTreeHistory(reChCam) ;
	}
	// Charged SoftDrop
	// First find the charged particle SD
	
	float t_genChSdPt =0 ;
	float t_genChSdMass=0;
	float t_genChSdZ=0;
	float t_genChSdTheta=0;
	if ( camChJets.size() > 0 )   {
	  fastjet::PseudoJet chSd_jet = softdropper(camChJets[0]);
	  t_genChSdMass = chSd_jet.m() * 0.001;
	  t_genChSdPt = chSd_jet.pt() * 0.001;
	  fastjet::PseudoJet parent1, parent2;
          if (  chSd_jet.has_parents(parent1, parent2) ) {
	    t_genChSdTheta =  DeltaR( parent1.phi(), parent1.rapidity(), parent2.phi(), parent2.rapidity() ) ;
	    t_genChSdZ     = std::min( parent1.pt(),  parent2.pt() ) / ( parent1.pt() +  parent2.pt() ) ;
            if ( _saveLog)   {
              cout << " charged subjet profile:"<<endl;
	      logSubJets( parent1, parent2 ) ;
            }
          }
          else {
            cout << "This Charged Softdrop jet is not divided!" << endl;
          }
        }
	
	vGenChSdPt.push_back(t_genChSdPt);
	vGenChSdMass.push_back(t_genChSdMass)  ;
	vGenChSdZ.push_back(t_genChSdZ)  ;
	vGenChSdTheta.push_back(t_genChSdTheta) ;
	
      }
      truthParticles.clear();
      truthCharges.clear();
    }
    
    else  { // If to use the AOD container
      cout << " Don't use useReAntiKt=0 option yet.." << endl; 
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
    myJetSub.weight = event_weight;
    myJetSub.recoSdPt = vSdpt_reco[ri];
    myJetSub.recoSdMass = vSdmass_reco[ri];
    myJetSub.recoSdTheta = vSdtheta_reco[ri];
    myJetSub.recoSdZ = vSdz_reco[ri];

    
    if (matchId != -1) {
      myJetSub.matchDr = drMin;
      myJetSub.genPt = vpt_gen[matchId];
      myJetSub.genEta = veta_gen[matchId];
      myJetSub.genRcPt = vptRc_gen[matchId];
      myJetSub.genSdMass = vSdmass_gen[matchId];
      myJetSub.genSdPt = vSdpt_gen[matchId];
      myJetSub.genSdtheta = vSdtheta_gen[matchId];
      myJetSub.genSdz = vSdz_gen[matchId];
      myJetSub.genChSdPt = vGenChSdPt[matchId];
      myJetSub.genChSdMass = vGenChSdMass[matchId];
      myJetSub.genChSdZ = vGenChSdZ[matchId];
      myJetSub.genChSdTheta = vGenChSdTheta[matchId];

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
  vSdpt_reco.clear();
  vSdtheta_reco.clear();
  vSdz_reco.clear();

  vpt_gen.clear();
  veta_gen.clear();
  vphi_gen.clear();
  vptRc_gen.clear();
  vSdmass_gen.clear();
  vSdpt_gen.clear();
  vSdtheta_gen.clear();
  vSdz_gen.clear();
  
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

