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



#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/BackgroundEstimatorBase.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"

#include "JetRec/JetSoftDrop.h"
#include "JetRecTools/ConstituentSubtractorTool.h"
#include "fastjet/tools/Filter.hh"
#include "JetRec/JetTrimmer.h"

using namespace std;
using namespace JetHelperTools;
using namespace TrackHelperTools;
using namespace fastjet;

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
  float weight, rhoCh;
  float recoPt, recoRawPt, recoEta, recoRap, recoPhi, recoRcPt, recoSdPt, recoSdMass, recoSdZ, recoSdTheta, recoSdPt1, recoSdRap1, recoSdPhi1, recoSdPt2, recoSdRap2, recoSdPhi2;
  float reChSdPt, reChSdMass, reChSdZ, reChSdTheta;
  float matchDr, genPt,  genEta, genRap, genPhi, genRcPt,  genSdPt, genSdMass, genSdz, genSdtheta, genSdPt1, genSdRap1, genSdPhi1, genSdPt2, genSdRap2, genSdPhi2;
  float genChSdPt, genChSdMass, genChSdZ, genChSdTheta; 
  float recoTrMass, recoChTrMass, genTrMass, genChTrMass, recoTrNsub, recoChTrNsub, genTrNsub, genChTrNsub; 
};
jetSubStr myJetSub;

TString evtText = "cent/I:weight/F:rhoCh";
TString recoText = "pt:rawPt:eta:y:phi:rcPt:sdPt:sdMass:zg:theta:spt1:sy1:sphi1:spt2:sy2:sphi2";
TString reChText = "chSdPt:chSdMass:chZg:chTheta";
TString genText = "dr:genPt:genEta:genRap:genPhi:genRcPt:genSdPt:genSdMass:genZg:genTheta:genSpt1:genSy1:genSphi1:genSpt2:genSy2:genSphi2";
TString genChText = "genChSdPt:genChSdMass:genChZg:genChTheta";
TString trimText  = "trMass:chTrMass:genTrMass:genChTrMass:trNsub:chTrNsub:genTrNsub:genChTrNsub";

TString myJetSubText = evtText+":"+recoText+":"+reChText+":"+genText+":"+genChText+":"+trimText;


// this is needed to distribute the algorithm to the workers
ClassImp(JetSubstructure)


void resetSubstr (jetSubStr &jetsub)
{


  jetsub.cent = -1 ;
  jetsub.recoPt= -1;
  jetsub.recoRawPt= -1;
  jetsub.recoEta=-100;
  jetsub.recoRap=-100;
  jetsub.recoPhi=-100;
  jetsub.recoRcPt=-1;
  jetsub.recoSdPt=-1;
  jetsub.recoSdMass=-1;
  jetsub.recoSdZ = -1;
  jetsub.recoSdTheta = -1;
  jetsub.recoSdPt1 = -1;
  jetsub.recoSdRap1 = -100;
  jetsub.recoSdPhi1 = -100;
  jetsub.recoSdPt2 = -1;
  jetsub.recoSdRap2 = -100;
  jetsub.recoSdPhi2 = -100;


  jetsub.genPt=-1;
  jetsub.genEta=-10;
  jetsub.genRap=-10;
  jetsub.genPhi=-10;
  jetsub.genRcPt=-1;
  jetsub.genSdPt=-1;
  jetsub.genSdMass=-1;
  jetsub.matchDr=-1;
  jetsub.weight = 1;
  jetsub.genSdz = 100;
  jetsub.genSdtheta = 100;
  jetsub.genSdPt1 = -1;
  jetsub.genSdRap1 = 100;
  jetsub.genSdPhi1 = 100;
  jetsub.genSdPt2 = -1;
  jetsub.genSdRap2 = -10;
  jetsub.genSdPhi2 = -10;

  jetsub.recoTrMass = -1;
  jetsub.recoChTrMass = -1;
  jetsub.genTrMass = -1;
  jetsub.genChTrMass = -1;
  jetsub.recoTrNsub = -1;
  jetsub.recoChTrNsub = -1;
  jetsub.genTrNsub = -1;
  jetsub.genChTrNsub = -1;
}


int getMaxPtIndex ( vector<fastjet::PseudoJet>& jets ) {
  double max_pt = 0;
  int maxIndex = -1;
  for (int i = 0; i < jets.size(); i++) {
    double jet_pt = jets[i].pt()*0.001;
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
  cout << " Full constituents lists (pt, eta, phi)" <<endl;
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

	if (_saveEvtDisplay) {
	  eventDisRecTow = new TH2F("eventDisRecTow","",20,-1,1,20,-1,1);
	  eventDisRecTow1 = (TH2F*)eventDisRecTow->Clone("eventDisRecTow1");
	  eventDisRecTow2 = (TH2F*)eventDisRecTow->Clone("eventDisRecTow2");
	  eventDisGen = new TH2F("eventDisGen","",20,-1,1,20,-1,1);
	  eventDisGen1 = (TH2F*)eventDisGen->Clone("eventDisGen1");
	  eventDisGen2 = (TH2F*)eventDisGen->Clone("eventDisGen2");
	  treeOut->Branch("rect","TH2F",&eventDisRecTow,256000,0);
	  treeOut->Branch("rect1","TH2F",&eventDisRecTow1,256000,0);
	  treeOut->Branch("rect2","TH2F",&eventDisRecTow2,256000,0);
	  treeOut->Branch("genp","TH2F",&eventDisGen,256000,0);
	  treeOut->Branch("genp1","TH2F",&eventDisGen1,256000,0);
	  treeOut->Branch("genp2","TH2F",&eventDisGen2,256000,0);
	}
	

	//	jetTree = new TTree("tr2","jet branching tree");
	//	jetTree->Branch("bran","");
	
	TH2D* temphist_2d;
	TH3D* temphist_3d;

	const int nJetPtBinForEff = 9;
	double jetPtBinForEff[nJetPtBinForEff+1] = {0,100,150,200,250,300,350,400,600,1000};
	const int nTrkPtBinForEff = 15;
	double trkPtBinForEff[nTrkPtBinForEff+1] = {0,0.5,0.7,0.8,0.9,1,1.5,2,4,6,8,14,20,40,100,200};
	const int nDphiDetaBinForEff = 20;
	double dPhiBinForEff[nDphiDetaBinForEff+1];
	double dEtaBinForEff[nDphiDetaBinForEff+1];
	for ( int i=0 ; i<= nDphiDetaBinForEff ; i++) {
	  dPhiBinForEff[i] = -1 + i*0.1;
	  dEtaBinForEff[i] = -1 + i*0.1;
	}
	
	for (int i=0;i<nCentbins;i++)  {
	  temphist_3d = new TH3D(Form("h_trkGen_dphi_cent%i",i),";pt;dphi",
				 nJetPtBinForEff, jetPtBinForEff, nTrkPtBinForEff, trkPtBinForEff, nDphiDetaBinForEff, dPhiBinForEff);
	  h_trkGen_pt_dphi_cent.push_back( temphist_3d);
	  h_trkGen_pt_dphi_cent.at(i)->Sumw2();

	  temphist_3d = new TH3D(Form("h_allGen_dphi_cent%i",i),";pt;dphi",
				 nJetPtBinForEff, jetPtBinForEff, nTrkPtBinForEff, trkPtBinForEff, nDphiDetaBinForEff, dPhiBinForEff);
	  h_allGen_pt_dphi_cent.push_back( temphist_3d);
	  h_allGen_pt_dphi_cent.at(i)->Sumw2();

	  temphist_3d = new TH3D(Form("h_trkGen_drap_cent%i",i),";pt;deta",
				 nJetPtBinForEff, jetPtBinForEff, nTrkPtBinForEff, trkPtBinForEff, nDphiDetaBinForEff, dEtaBinForEff);
	  h_trkGen_pt_drap_cent.push_back( temphist_3d);
	  h_trkGen_pt_drap_cent.at(i)->Sumw2();

	  temphist_3d = new TH3D(Form("h_allGen_drap_cent%i",i),";pt;deta",
				 nJetPtBinForEff, jetPtBinForEff, nTrkPtBinForEff, trkPtBinForEff, nDphiDetaBinForEff, dEtaBinForEff);
	  h_allGen_pt_drap_cent.push_back( temphist_3d);
	  h_allGen_pt_drap_cent.at(i)->Sumw2();


	  temphist_2d = new TH2D(Form("hTrkPtEta_preCS_cent%i",i),";pT;eta",200,0,200,20,-3,3);
	  hTrkPtEta_preCS_cent.push_back(temphist_2d);
	  hTrkPtEta_preCS_cent.at(i)->Sumw2();

	  temphist_2d = new TH2D(Form("hTrkPtEta_postCS_cent%i",i),";pT;eta",200,0,200,20,-3,3);
	  hTrkPtEta_postCS_cent.push_back(temphist_2d);
	  hTrkPtEta_postCS_cent.at(i)->Sumw2();
	}
	
	hGenNcam = new TH2D("genCam","; p_{T} GeV/c ; Number of C/A clusters",100,0,1000,6,-0.5,5.5);
	hGenNchCam = (TH2D*)hGenNcam->Clone("genChCam");
	hRecoNcam = (TH2D*)hGenNcam->Clone("recoCam");
	hRecoNchCam = (TH2D*)hGenNcam->Clone("recoChCam");

	hGenSdStat = new TH2D("genSd","; p_{T} GeV/c ; [SoftDrop]  0=fail,  1=pass",100,0,1000,2,-0.5,1.5);
	hGenSdChStat = (TH2D*)hGenSdStat->Clone("genChSd");
	hRecoSdStat = (TH2D*)hGenSdStat->Clone("recoSd");
	hRecoSdChStat = (TH2D*)hGenSdStat->Clone("recoChSd");
	
	
	wk()->addOutput (h_RejectionHisto);
	wk()->addOutput (h_centrality);
	wk()->addOutput (hET_ETsub);
	wk()->addOutput (treeOut);

	wk()->addOutput (hGenNcam);
	wk()->addOutput (hGenNchCam);
	wk()->addOutput (hGenSdStat);
	wk()->addOutput (hGenSdChStat);
	wk()->addOutput (hRecoNcam);
	wk()->addOutput (hRecoNchCam);
	wk()->addOutput (hRecoSdStat);
	wk()->addOutput (hRecoSdChStat);
	
        for (int i=0;i<nCentbins;i++)  {
	  wk()->addOutput (h_trkGen_pt_dphi_cent.at(i));
	  wk()->addOutput (h_allGen_pt_dphi_cent.at(i));
	  wk()->addOutput (h_trkGen_pt_drap_cent.at(i));
	  wk()->addOutput (h_allGen_pt_drap_cent.at(i));
	  wk()->addOutput (hTrkPtEta_preCS_cent.at(i));
	  wk()->addOutput (hTrkPtEta_postCS_cent.at(i));
	}
	
	
        cout << "=======================================" << endl;
	cout << "Parametrization of d0 cut" << endl;
	f_d0_cut = new TF1("f1", "[0]*exp([1]*x)+[2]*exp([3]*x)", 0.4, 500);
        f_d0_cut->SetParameters(0.472367, -0.149934, 0.193095, 0.000337765);

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
		
	
	///////////// track selection ///////////////////////////////////
	m_trackSelectorTool = new InDet::InDetTrackSelectionTool("InDetTrackSelectorTool");
	SetCutLevel(m_trackSelectorTool, _trk_cut_level.c_str());
	EL_RETURN_CHECK("initialize()",m_trackSelectorTool->initialize());
	

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
  bool useReAntiKt = true; 


  


  
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
  //  double event_weight_fcal = 1;
  //Centrality
  const xAOD::HIEventShapeContainer* calos=0;
  ANA_CHECK(event->retrieve( calos, "CaloSums"));
  if (_centrality_scheme>1)	  {
    FCalEt=calos->at(5)->et()*1e-6;
    cent_bin = GetCentralityBin(_centrality_scheme, FCalEt, isHIJING);
    
    //    event_weight_fcal = jetcorr->GetFCalWeight(FCalEt);
    //    if (_isMC && isHIJING) event_weight_fcal = 1;
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
  

  //Vertex : 
  const xAOD::VertexContainer * vertices = 0;
  if ( !event->retrieve( vertices, "PrimaryVertices" ).isSuccess() ){
    Error("execute()", "Failed to retrieve VertexContainer container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  if(vertices->size()<2) {
    h_RejectionHisto->Fill(3.5);
    return EL::StatusCode::SUCCESS;
  }
  // Find primary vertex
  xAOD::VertexContainer::const_iterator vtx_itr = vertices->begin();
  xAOD::VertexContainer::const_iterator vtx_end = vertices->end();
  const xAOD::Vertex* primaryVertex = 0;
  for(;vtx_itr!=vtx_end;++vtx_itr)	  {
    if((*vtx_itr)->vertexType()==xAOD::VxType::PriVtx) {
      primaryVertex = (*vtx_itr);
      break;
    }}
  
  if(primaryVertex)	  {
    if (fabs(primaryVertex->z())>150.)
      {
	return EL::StatusCode::SUCCESS;
      }}
  else {
    h_RejectionHisto->Fill(5.5);
    return EL::StatusCode::SUCCESS;
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
  
  //  h_FCal_Et->Fill(FCalEt, event_weight_fcal); //filled here to get proper event weight
  h_centrality->Fill(cent_bin,1); //  weight is set 1 event_weight_fcal);
  
  if ( _saveLog) cout << "Centrality of this event = " << cent_bin << endl;

  TH2F* t_recTow[20];
  TH2F* t_recTow1[20];
  TH2F* t_recTow2[20];
  TH2F* t_genTow[20];
  TH2F* t_genTow1[20];
  TH2F* t_genTow2[20];


  vector <double> vpt_reco;
  vector <double> vptRaw_reco;  
  vector <double> veta_reco;
  vector <double> vrap_reco;
  vector <double> vphi_reco;
  vector <double> vptRc_reco;
  vector <double> vSdmass_reco;
  vector <double> vSdpt_reco;
  vector <float> vSdz_reco;
  vector <float> vSdtheta_reco;
  vector <float> vSdpt1_reco;
  vector <float> vSdrap1_reco;
  vector <float> vSdphi1_reco;
  vector <float> vSdpt2_reco;
  vector <float> vSdrap2_reco;
  vector <float> vSdphi2_reco;
  
  vector <float> vRecoChSdPt; 
  vector <float> vRecoChSdMass; 
  vector <float> vRecoChSdZ; 
  vector <float> vRecoChSdTheta; 
  //trimming
  vector <double> vTrMass_reco;
  vector <double> vChTrMass_reco;
  vector <double> vTrNsub_reco;
  
  vector <double> vpt_gen;
  vector <double> veta_gen;
  vector <double> vrap_gen;
  vector <double> vphi_gen;
  vector <double> vptRc_gen;
  vector <double> vSdmass_gen;
  vector <double> vSdpt_gen;
  vector <float> vSdz_gen;
  vector <float> vSdtheta_gen;
  vector <float> vSdpt1_gen;
  vector <float> vSdrap1_gen;
  vector <float> vSdphi1_gen;
  vector <float> vSdpt2_gen;
  vector <float> vSdrap2_gen;
  vector <float> vSdphi2_gen;

  vector <float> vChSdPt_gen;
  vector <float> vChSdMass_gen  ;
  vector <float> vChSdZ_gen  ;
  vector <float> vChSdTheta_gen  ;
  //trimming
  vector <double> vTrMass_gen;
  vector <float> vTrNsub_gen;
  vector <double> vChTrMass_gen;



  // algorithm definition 
  fastjet::JetDefinition jetDefCam(fastjet::cambridge_algorithm, _ReclusterRadius);
  fastjet::JetDefinition jetDefAk(fastjet::antikt_algorithm, _ReclusterRadius);
  fastjet::JetDefinition jetDefAk04(fastjet::antikt_algorithm, 0.4);
  float beta = 0 ;
  float z_cut = 0.1;
  fastjet::contrib::SoftDrop softdropper(beta, z_cut);
  fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::kt_algorithm, _rSub), 
                          fastjet::SelectorPtFractionMin(_fCut));
  

  if ( _saveLog) { 
    if ( (_bkgKill !=0) && (_doTrimming) )
      cout <<endl<<endl<<" bkgKill is not zero while doTrimming is on!!!!!" << endl<<endl<<endl ;
    if ( _doTrimming ) cout << "Trimming is on" << endl;
    else cout << "Trimming is off" << endl;
    cout << "Trimmer: " << trimmer.description()  << endl;
  }

  // reference : http://acode-browser2.usatlas.bnl.gov/lxr-rel21/source/atlas/Reconstruction/Jet/JetRec/Root/JetTrimmer.cxx
  
  ////////////// MC truth particles /////////////
  
  std::vector<fastjet::PseudoJet> truthParticles;
  std::vector<fastjet::PseudoJet> truthChargesAna;
  if (_isMC){
    const xAOD::TruthParticleContainer * genParCont = 0;
    ANA_CHECK(event->retrieve( genParCont, "TruthParticles"));
    
    for( xAOD::TruthParticleContainer::const_iterator truth_itr = genParCont->begin() ; truth_itr!= genParCont->end() ; ++truth_itr) {
      //      cout <<"  (*truth_itr)->p4().Pt() = " <<  (*truth_itr)->p4().Pt() << endl; // MeV is confirmed
      double ptTrk = (*truth_itr)->p4().Pt() * 0.001 ; 
      int thebc = (*truth_itr)->barcode();
      if ( (thebc > 0) && ( thebc < 200000 ) )  {
	
	truthParticles.push_back( (*truth_itr)->p4() );   // To be used for anti-kT jet reconstrucutre 
	
	if ( (fabs((*truth_itr)->charge()) > 0 ) && ( ptTrk > _pTtrkCut ) && ( fabs((*truth_itr)->p4().Eta()) < _etaTrkCut )  ) 
	  truthChargesAna.push_back( (*truth_itr)->p4() ) ; // For softdrop
	
      }
    }
  }
  


  double event_weight = 1;
  if (_isMC){
    if ( useReAntiKt )  {
      ////////// On-the-fly jet finder ///////////////////////////////////////////////
      if (_saveLog) {
	cout << "   * Reweighting factor using R=0.4 jet" << endl;
      }
      fastjet::ClusterSequence csR04(truthParticles, jetDefAk04);
      vector<fastjet::PseudoJet> jets04 = fastjet::sorted_by_pt(csR04.inclusive_jets());
      int maxPtId = getMaxPtIndex ( jets04 ) ;
      event_weight = 0;
      cout <<"     * STRANGE! there is no R=0.4 jet" << endl;
      if (maxPtId > -1)   {
	if (_saveLog)  {
	  cout << "      * Max pT (and y,phi) for R=0.4 Jets = " << jets04[maxPtId].pt() * 0.001 << " GeV, " << jets04[maxPtId].rapidity() << ", "<<jets04[maxPtId].phi() << endl;
	}
	event_weight = jetcorr->GetJetWeight( jets04[maxPtId].pt() * 0.001, jets04[maxPtId].rapidity(), PhiInPI(jets04[maxPtId].phi()) ); // pt unit needs to be GeV    ToBeFixed  : this value is strange  // phi must be +/- 3.14,  
	//	cout << "A:  " << (jets04[maxPtId].phi()) << ",  B: " << PhiInPI(jets04[maxPtId].phi()) << endl;
      }
      jets04.clear();
    }
  }
  
  
  
  
  

  
  ///////////// tracks ////////////////////////////////////////
  std::vector<fastjet::PseudoJet> selectedTrks;
  std::vector<fastjet::PseudoJet> selGenMatchTrks;
  
  const xAOD::TrackParticleContainer* recoTracks = 0;
  EL_RETURN_CHECK("execute",event->retrieve( recoTracks, "InDetTrackParticles"));
  for (const auto& trk : *recoTracks) {
    float pt = trk->pt()/1000.;
    float eta = trk->eta();
    float phi = trk->phi();
    if (_saveLog)   cout << "   * Track pt, eta, phi = " << pt <<", "<<eta<<", "<<phi<<endl;
    
    if ( pt < _pTtrkCut ) 
      continue;
    
    if ( fabs(eta) > _etaTrkCut ) 
      continue;
    
    
    
    if(!m_trackSelectorTool->accept(*trk, *vtx_itr )) 
      continue;
    if (_saveLog)    cout << "     passed selection cut "<< endl;
    
    double d0 = trk->d0();
    double d0_cut = f_d0_cut->Eval(pt);  // Make sure the unit is GeV!!!!
    if(fabs(d0) > d0_cut) continue;
    if (_saveLog)    cout << "     passed the d0 cut " << endl;
    
    double pionMass = 139.570 ; // in MeV
    double pionx = trk->p4().Px(); // Oct 24th, Yongsun confiremd p4() is in MeV unit
    double piony = trk->p4().Py();
    double pionz = trk->p4().Pz();
    double pione = sqrt (pionx*pionx + piony*piony + pionz*pionz + pionMass*pionMass) ;
    //    cout << " pt, px,py,pz = " << trk->pt() <<"," <<pionx<<", "<<piony<<", "<<pionz <<endl; 
    
    // truth particle matching 
    bool isTruthMatched=false;
    fastjet::PseudoJet matchPar = fastjet::PseudoJet (0,0,0,0);

    if ( _isMC) { 
      //Truth matching
      ElementLink< xAOD::TruthParticleContainer > truthLink = trk->auxdata<ElementLink< xAOD::TruthParticleContainer > >("truthParticleLink");
      float mcprob=-1.;
      if(truthLink.isValid()) {
	mcprob = trk->auxdata<float>("truthMatchProbability");
	if ( ( mcprob > 0.3 ) && ( (*truthLink)->barcode() > 0 ) &&  ( (*truthLink)->barcode() < 200000 ) ) {  // 0.3 is hard-coded at the momnet.  ToBeFixed   
	  isTruthMatched = true;
	  matchPar = fastjet::PseudoJet( (*truthLink)->p4() );
	}
      }
    }
    
    selectedTrks.push_back( fastjet::PseudoJet ( pionx, piony, pionz, pione ) );
    hTrkPtEta_preCS_cent.at(cent_bin)->Fill( pt, eta, event_weight); 
    if (isTruthMatched) selGenMatchTrks.push_back( matchPar );
    
  }
  if (_saveLog) cout << " number of reconstructed tracks: " << selectedTrks.size() << endl;
  
  
  // Background for tracks 
  // See https://github.com/cms-externals/fastjet-contrib/blob/master/ConstituentSubtractor/example.cc  
  // algo for background.  
  double ghost_area=0.01;  
  double R_bkg = 0.4; 
  fastjet::JetDefinition jet_def_for_rho(fastjet::kt_algorithm, R_bkg);
  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts,fastjet::GhostedAreaSpec(_etaTrkCut,1)); // the area definiton is used only for the jet backgroud estimator. It is not important for the ConstituentSubtractor when subtracting the whole event - this is not true when subtracting the individual jets

  fastjet::Selector rho_range_trk =  fastjet::SelectorAbsEtaMax(_etaTrkCut - R_bkg);
  fastjet::JetMedianBackgroundEstimator bge_rho_trk(rho_range_trk, jet_def_for_rho, area_def);
  fastjet::BackgroundJetScalarPtDensity *scalarPtDensity = new fastjet::BackgroundJetScalarPtDensity();
  bge_rho_trk.set_jet_density_class(scalarPtDensity); // this changes the computation of pt of patches from vector sum to scalar sum. The scalar sum seems more reasonable.
  bge_rho_trk.set_particles(selectedTrks);

  fastjet::contrib::ConstituentSubtractor subtractor_trk; 
  subtractor_trk.set_max_standardDeltaR(0.25); 
  subtractor_trk.set_alpha(0);
  subtractor_trk.set_ghost_area(ghost_area);  // same 
  subtractor_trk.set_background_estimator(&bge_rho_trk);
  subtractor_trk.set_common_bge_for_rho_and_rhom(true); // for massless input particles it\
  
  if ( _saveLog)   cout << endl << subtractor_trk.description() << endl ; 
  
  vector<PseudoJet> corrected_selectedTrks = subtractor_trk.subtract_event(selectedTrks, _etaTrkCut);
  // Fill the histogram:
  for ( int ii =0 ; ii < corrected_selectedTrks.size() ; ii++ ) {
    hTrkPtEta_postCS_cent.at(cent_bin)->Fill( corrected_selectedTrks[ii].pt()*0.001,  corrected_selectedTrks[ii].eta(), event_weight);
  }
  
  // background subtraction 
  if (_saveLog) {
    cout << endl <<" number of reconstructed tracks: " << selectedTrks.size() << endl;
    cout <<"  pT,eta,phi" << endl;
    vector<fastjet::PseudoJet> sorted = fastjet::sorted_by_pt(selectedTrks);
    for ( int ii =0 ; ii <sorted.size() ; ii++ ) {
      cout <<"    "<<sorted[ii].pt()*0.001 <<",   "<< sorted[ii].eta() << ",   "<<sorted[ii].phi() << endl;
    }
    cout << endl <<" After Subtraction number of tracks: " << corrected_selectedTrks.size() << endl;
    cout <<"  pT,eta,phi" << endl;
    vector<fastjet::PseudoJet> sorted2 = fastjet::sorted_by_pt(corrected_selectedTrks);
    for ( int ii =0 ; ii <sorted2.size() ; ii++ ) {
      cout <<"    "<<sorted2[ii].pt()*0.001 <<",   "<< sorted2[ii].eta() << ",   "<<sorted2[ii].phi() << endl;
    }
    cout << endl;
    sorted.clear();
    sorted2.clear();
  }



  /////////////   Reco jets /////////////////////////////////////////
  
  xAOD::TStore *store = new xAOD::TStore; //For calibration
  const xAOD::JetContainer* reco_jets = 0;
  ANA_CHECK(event->retrieve( reco_jets, _reco_jet_collection.c_str() ));
  
  xAOD::JetContainer::const_iterator jet_itr = reco_jets->begin();
  xAOD::JetContainer::const_iterator jet_end = reco_jets->end();
  int nRecoJetCounter=0;
  for( ; jet_itr != jet_end; ++jet_itr ) {
    
    xAOD::Jet theRecoJet;
    theRecoJet.makePrivateStore( **jet_itr );
    
    const xAOD::JetFourMom_t jet_4momUnCal = theRecoJet.jetP4("JetSubtractedScaleMomentum"); // uncalib
    const xAOD::JetFourMom_t jet_4momCalib = theRecoJet.jetP4();   // CALIBRATED!!! 
    
    fastjet::PseudoJet CalibFourVec = fastjet::PseudoJet ( jet_4momCalib.px(), jet_4momCalib.py(), jet_4momCalib.pz(), jet_4momCalib.energy() );
    fastjet::PseudoJet unCaliFourVec = fastjet::PseudoJet ( jet_4momUnCal.px(), jet_4momUnCal.py(), jet_4momUnCal.pz(), jet_4momUnCal.energy() );
    double jet_pt  = jet_4momCalib.pt() * 0.001 ;
    double jet_eta = jet_4momCalib.eta();
    double jet_rap = CalibFourVec.rapidity();
    double jet_phi = PhiInPI ( jet_4momCalib.phi() );
    double jet_ptRaw = jet_4momUnCal.pt() * 0.001;
    if (_saveLog) { 
      cout <<" jet_pt = " << jet_pt << endl;
      cout <<" jet_eta = " << jet_eta << endl;
    }
    if (jet_pt < _pTjetCut)          continue;
    if (fabs(jet_eta) > _etaJetCut)  continue;
    
    // WARNING! THERE MUST BE NO CONTINUE COMMAND FROM NOW ON IN THIS LOOP!!!!! 
    //    nRecoJetCounter++;    will be written at the end of this loop.
    
    if (_saveEvtDisplay) {
      eventDisRecTow->Reset();
      eventDisRecTow1->Reset();
      eventDisRecTow2->Reset();
      t_recTow[nRecoJetCounter] = (TH2F*)eventDisRecTow->Clone(Form("recTow_i%d",nRecoJetCounter));
      t_recTow1[nRecoJetCounter] = (TH2F*)eventDisRecTow->Clone(Form("recTow1_i%d",nRecoJetCounter));
      t_recTow2[nRecoJetCounter] = (TH2F*)eventDisRecTow->Clone(Form("recTow2_i%d",nRecoJetCounter));
      t_recTow[nRecoJetCounter]->Reset();
      t_recTow1[nRecoJetCounter]->Reset();
      t_recTow2[nRecoJetCounter]->Reset();
    }

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
      
      
      
      
      if ( _bkgKill == -1 ) { 
	if ( (*itCnst)->pt() > 0 ) { // normal tower 
	  nonZeroConsts.push_back(thisConst);
	  toBeSubtracted.push_back ( fastjet::PseudoJet ( 0.000001, 0.000001, 0.000001, 0.000002 )) ; // place holder
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
      else if ( _bkgKill == 0 ) {  // Just ignore negative towers
	if ( (*itCnst)->pt() > 0 ) {
	  nonZeroConsts.push_back(thisConst);
	  toBeSubtracted.push_back ( fastjet::PseudoJet ( 0.000001, 0.000001, 0.000001, 0.000002 )) ; // place holder
          nFlag.push_back (false);
	}
      }
      else if ( _bkgKill == 1 ) { // soft kill
	cout << "No SoftKilling module deployed yet" << endl; 
      }
    
      // Fill the eventdisplay histogram first! 
      if (_saveEvtDisplay) 
	t_recTow[nRecoJetCounter]->Fill(theEta - jet_rap, DeltaPhi(thePhi, jet_phi), (*itCnst)->pt() *0.001 ) ;
      
    }

    // Cambridge reclustering 
    if ( _saveLog) cout << "Reco re-clustering starts!" << endl;
    fastjet::ClusterSequence recoCamSeq(nonZeroConsts, jetDefCam);
    vector<fastjet::PseudoJet> cambridgeJet = recoCamSeq.inclusive_jets();
    vector<int> pIndex = recoCamSeq.particle_jet_indices( cambridgeJet);
    if ( _saveLog)  {
      drawTreeHistory(recoCamSeq) ;
    }
    
    hRecoNcam->Fill(jet_pt,cambridgeJet.size());
    
    // restore negative towers 
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
  
    double thePtrc = -1;
    double thesdpt = -1 ;
    double thesdm = -1;
    float thesdtheta =-1 ;
    float thesdz = -1;
    float thesdpt1 = -1;
    float thesdrap1 = -100;
    float thesdphi1 = -100;
    float thesdpt2 = -1;
    float thesdrap2 = -100;
    float thesdphi2 = -100;
    
    float t_recoChSdPt = -1 ;
    float t_recoChSdMass= -1;
    float t_recoChSdZ= -1 ;
    float t_recoChSdTheta= -1;
    float t_recoTrMass = -1;
    float t_recoTrNsub = -1;
    float t_recoChTrMass = -1;
    if ( corrRecCamJets.size() > 0 )   {
      fastjet::PseudoJet originalJet;
      // First trim this!! 
      if ( _doTrimming) {
	originalJet = trimmer(corrRecCamJets[0]);
	if (_saveLog) {
	  cout << " cambridge jet pt = " << originalJet.pt()*0.001 << endl;
	  cout << "after trimming: " << originalJet.pt()*0.0001 <<endl;
	  cout <<" nSubJets = " << originalJet.pieces().size() << endl;
	}
	t_recoTrMass = originalJet.m() * 0.001; 
	t_recoTrNsub   = originalJet.pieces().size() ;
      }
      else { 
	originalJet = corrRecCamJets[0]; 
      }
      
      
      thePtrc = originalJet.pt() * 0.001;
      fastjet::PseudoJet sd_jet = softdropper(originalJet);
      thesdpt = sd_jet.pt() * 0.001; 
      thesdm =  sd_jet.m() * 0.001; 
      
      if (_saveLog) {
	cout << "RECO JET softdrop:" << endl;
	showLadder ( unCaliFourVec, originalJet, sd_jet );
      }
      fastjet::PseudoJet parent1, parent2;
      if (  sd_jet.has_parents(parent1, parent2) ) {
	hRecoSdStat->Fill(jet_pt, 1);
	thesdtheta =  DeltaR( parent1.phi(), parent1.rapidity(), parent2.phi(), parent2.rapidity() ) ;
	thesdz  = std::min( parent1.pt(),  parent2.pt() ) / ( parent1.pt() +  parent2.pt() ) ;

	thesdpt1 = parent1.pt()*0.001;
	thesdphi1 = parent1.phi();
	thesdrap1 = PhiInPI(parent1.rapidity());
	thesdpt2 = parent2.pt()*0.001;
	thesdphi2 = parent2.phi();
	thesdrap2 = PhiInPI(parent2.rapidity());


	if ( _saveLog)   {
	  logSubJets( parent1, parent2 ) ;
	}
	//Fill the histogram
	if (_saveEvtDisplay) { 
	  vector<fastjet::PseudoJet> sub1 =  parent1.constituents() ;
	  vector<fastjet::PseudoJet> sub2 =  parent2.constituents() ;
	  //	DeltaR( parent1.phi(), parent1.rapidity(), parent2.phi(), parent2.rapidity()
	  for ( int ic = 0 ; ic< sub1.size() ; ic++)
	    t_recTow1[nRecoJetCounter]->Fill( sub1[ic].eta() - jet_rap, DeltaPhi( sub1[ic].phi(), jet_phi), sub1[ic].pt() * 0.001 );
	  for ( int ic = 0 ; ic< sub2.size() ; ic++)
	    t_recTow2[nRecoJetCounter]->Fill( sub2[ic].eta() - jet_rap, DeltaPhi( sub2[ic].phi(), jet_phi), sub2[ic].pt() * 0.001 );
	  sub1.clear();
	  sub2.clear();
	}
      }
      else  {
	hRecoSdStat->Fill(jet_pt, 0);
	if (_saveLog)	cout << "This softdrop jet is not divided!" << endl; 
      }
      
    }

    // Charged Particle reclustering
    if (_saveLog)    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << "(RECO) Charged track clustering starts" << endl;
    vector<fastjet::PseudoJet> trkConsts;
    for ( int ic=0; ic< corrected_selectedTrks.size() ; ic++) {
      if ( DeltaR ( jet_phi, jet_rap, corrected_selectedTrks[ic].phi(), corrected_selectedTrks[ic].rapidity() ) <  _ReclusterRadius ) {
	if ( corrected_selectedTrks[ic].pt()*0.001 < _pTtrkCutForSD )  
	  continue;
	
	trkConsts.push_back( corrected_selectedTrks[ic]) ;	
      }
    }
  
    // Tracking efficiency calculation 
    if ( _isMC) {  
      for ( int ic=0; ic< selGenMatchTrks.size() ; ic++) {
	if ( DeltaR ( jet_phi, jet_rap, selGenMatchTrks[ic].phi(), selGenMatchTrks[ic].rapidity() ) <  _ReclusterRadius ) {
	  h_trkGen_pt_dphi_cent.at(cent_bin)->Fill( jet_pt, selGenMatchTrks[ic].pt()*0.001, DeltaPhi(selGenMatchTrks[ic].phi(), jet_phi) ) ;
	  h_trkGen_pt_drap_cent.at(cent_bin)->Fill( jet_pt, selGenMatchTrks[ic].pt()*0.001, selGenMatchTrks[ic].rapidity() - jet_rap );
	}
      }
      for ( int ic=0; ic< truthChargesAna.size() ; ic++) {
	if ( DeltaR ( jet_phi, jet_rap, truthChargesAna[ic].phi(), truthChargesAna[ic].rapidity() ) <  _ReclusterRadius ) {
	  h_allGen_pt_dphi_cent.at(cent_bin)->Fill( jet_pt, truthChargesAna[ic].pt()*0.001, DeltaPhi(truthChargesAna[ic].phi(), jet_phi) ) ;
	  h_allGen_pt_drap_cent.at(cent_bin)->Fill( jet_pt, truthChargesAna[ic].pt()*0.001, truthChargesAna[ic].rapidity() - jet_rap );
	}
      }
    }
    
    
    fastjet::ClusterSequence reChCam(trkConsts, jetDefCam);
    vector<fastjet::PseudoJet> camChJets = fastjet::sorted_by_pt(reChCam.inclusive_jets());
    if ( _saveLog)  {
      cout << "(RECO) Charged Cambridge jets" << endl;
      drawTreeHistory(reChCam) ;
    }
    
    hRecoNchCam->Fill(jet_pt,camChJets.size() ) ;
    if ( camChJets.size() > 0 )   {
      fastjet::PseudoJet chSd_jet = softdropper(camChJets[0]);
      t_recoChSdMass = chSd_jet.m() * 0.001;
      t_recoChSdPt = chSd_jet.pt() * 0.001;
      fastjet::PseudoJet parent1, parent2;
      if (  chSd_jet.has_parents(parent1, parent2) ) {
	hRecoSdChStat->Fill(jet_pt, 1);
	t_recoChSdTheta =  DeltaR( parent1.phi(), parent1.rapidity(), parent2.phi(), parent2.rapidity() ) ;
	t_recoChSdZ     = std::min( parent1.pt(),  parent2.pt() ) / ( parent1.pt() +  parent2.pt() ) ;
	if ( _saveLog)  logSubJets( parent1, parent2 ) ;
      }
      else {
	hRecoSdChStat->Fill(jet_pt, 0);
      }
    }
    
    
    vpt_reco.push_back(jet_pt);
    vptRaw_reco.push_back(jet_ptRaw);
    veta_reco.push_back(jet_eta);
    vrap_reco.push_back(jet_rap);
    vphi_reco.push_back(jet_phi);
    vptRc_reco.push_back(thePtrc);
    
    vSdmass_reco.push_back(thesdm);
    vSdpt_reco.push_back(thesdpt);
    vSdz_reco.push_back(thesdz);
    vSdtheta_reco.push_back(thesdtheta);
    vSdpt1_reco.push_back(thesdpt1);
    vSdrap1_reco.push_back(thesdrap1);
    vSdphi1_reco.push_back(thesdphi1);
    vSdpt2_reco.push_back(thesdpt2);
    vSdrap2_reco.push_back(thesdrap2);
    vSdphi2_reco.push_back(thesdphi2);
  
    vRecoChSdPt.push_back(t_recoChSdPt);
    vRecoChSdMass.push_back(t_recoChSdMass)  ;
    vRecoChSdZ.push_back(t_recoChSdZ)  ;
    vRecoChSdTheta.push_back(t_recoChSdTheta) ;
    
    vTrMass_reco.push_back(t_recoTrMass);
    vTrNsub_reco.push_back(t_recoTrNsub);
    vChTrMass_reco.push_back(t_recoChTrMass);
    
    
    corrRecCamJets.clear();
    nonZeroConsts.clear();
    toBeSubtracted.clear(); 
    
    nRecoJetCounter++;   // This must cone to the end of the loop!
  }
  
  if (_saveLog) { 
    if (  (nRecoJetCounter == vpt_reco.size() ) )
      cout << "Count numbers are consistent!!" << endl ;
    else  cout << "Counts are inconsistent!!"<<endl;
  }
  
  


  
  
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
      if (_saveLog) {
	cout << endl << "///////////////////////////////////" << endl;
	cout << " MC information " << endl;
      }
      
      fastjet::ClusterSequence cs(truthParticles, jetDefAk);
      vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
      // Now, new genJet is ready!
      
      int nGenJetCounter =0;
      for (unsigned i = 0; i < jets.size(); i++) {  // MC anti-kT jets
	double jet_pt = jets[i].pt()*0.001;
	double jet_eta = jets[i].pseudorapidity();  
	double jet_rap = jets[i].rapidity();
	double jet_phi = PhiInPI ( jets[i].phi() ) ;
	if (jet_pt < _truthpTjetCut) continue;
	if ( fabs(jet_eta) > _etaJetCut+0.2 ) continue;
      
	// IMPORTNAT!  There must be no more continue command in this loop! 
	if (_saveEvtDisplay) {	
	  eventDisGen->Reset();
	  eventDisGen1->Reset();
	  eventDisGen2->Reset();
	  t_genTow[nGenJetCounter] = (TH2F*)eventDisGen->Clone(Form("genp_i%d",nRecoJetCounter));
	  t_genTow1[nGenJetCounter] = (TH2F*)eventDisGen->Clone(Form("genp1_i%d",nRecoJetCounter));
	  t_genTow2[nGenJetCounter] = (TH2F*)eventDisGen->Clone(Form("genp2_i%d",nRecoJetCounter));
	  t_genTow[nGenJetCounter]->Reset();
	  t_genTow1[nGenJetCounter]->Reset();
	  t_genTow2[nGenJetCounter]->Reset();
	}


	// Truth recluster by C/A
	vector<fastjet::PseudoJet > akConsts = jets[i].constituents();
	if (_saveEvtDisplay) {
	  for ( int ic = 0 ; ic< akConsts.size() ; ic++)  {
	    double theRap  = akConsts[ic].rapidity();
	    double thePhi  = PhiInPI( akConsts[ic].phi() );
	    double thePt   = akConsts[ic].pt() * 0.001;
	    t_genTow[nGenJetCounter]->Fill(theRap - jet_rap, DeltaPhi(thePhi, jet_phi), thePt ) ;
	  }
	}
	fastjet::ClusterSequence reCam(akConsts, jetDefCam);
	vector<fastjet::PseudoJet> camJets = fastjet::sorted_by_pt(reCam.inclusive_jets()); // return a vector of jets sorted into decreasing energy
	
	
	
	if ( _saveLog)  {   	  
	  cout << "Truth Cambridge jet " << endl ;
	  drawTreeHistory(reCam) ;  
	}
	// Charged SoftDrop
	// First find the charged particle SD
	
	double thesdpt = -1 ;
	double thesdm = -1;
	double thePtrc = -1 ;
	float thesdtheta = -1 ; 
	float thesdz = -1;
	float thesdpt1 = -1 ; 
	float thesdrap1 = -10 ; 
	float thesdphi1 = -10 ; 
	float thesdpt2 = -1 ; 
	float thesdrap2 = -10 ; 
	float thesdphi2 = -10 ; 
        float t_genTrMass = -1;
        float t_genTrNsub = -1;

	hGenNcam->Fill( jet_pt, camJets.size() ) ;
	
	if ( camJets.size() > 0 )   {  // If cambridge jet is made .. 
	  fastjet::PseudoJet originalJet;
	  if ( _doTrimming) {
	    originalJet = trimmer(camJets[0]);
	    if (_saveLog) {
	      cout << "Gen cambridge jet pt = " << originalJet.pt()*0.001 << endl;
	      cout << "Gen after trimming: " << originalJet.pt()*0.0001 <<endl;
	      cout <<" Gen nSubJets = " << originalJet.pieces().size() << endl;
	    }
	    t_genTrMass = originalJet.m() * 0.001;
	    t_genTrNsub =  originalJet.pieces().size();
	  }
	  else {
	    originalJet = camJets[0];
	  }

	  
	  
	  thePtrc = originalJet.pt() * 0.001;
	  fastjet::PseudoJet sd_jet = softdropper(originalJet);
	  thesdm =  sd_jet.m() * 0.001;
	  thesdpt = sd_jet.pt() * 0.001;
	  if ( _saveLog) { 
	    cout << "GEN JET softdrop:" << endl;
	    showLadder ( jets[i], originalJet, sd_jet );  
	  }
	  
	  fastjet::PseudoJet parent1, parent2;
	  if (  sd_jet.has_parents(parent1, parent2) ) { // If softdrop worked
	    hGenSdStat->Fill( jet_pt, 1 ) ;
	    thesdtheta =  DeltaR( parent1.phi(), parent1.rapidity(), parent2.phi(), parent2.rapidity() ) ;
	    thesdz  = std::min( parent1.pt(),  parent2.pt() ) / ( parent1.pt() +  parent2.pt() ) ;
	    thesdpt1 =  parent1.pt() * 0.001; 
	    thesdrap1 =  parent1.rapidity(); 
	    thesdphi1 =  PhiInPI(parent1.phi()) ;
	    thesdpt2 =  parent2.pt() * 0.001; 
	    thesdrap2 =  parent2.rapidity(); 
	    thesdphi2 =  PhiInPI(parent2.phi()) ;

	    if ( _saveLog)  logSubJets( parent1, parent2 ) ;
	    
	    if (_saveEvtDisplay) {
	      vector<fastjet::PseudoJet> sub1 =  parent1.constituents() ;
	      vector<fastjet::PseudoJet> sub2 =  parent2.constituents() ;
	      //    DeltaR( parent1.phi(), parent1.rapidity(), parent2.phi(), parent2.rapidity()
	      for ( int ic = 0 ; ic< sub1.size() ; ic++)
		t_genTow1[nGenJetCounter]->Fill( sub1[ic].rapidity() - jet_rap, DeltaPhi( sub1[ic].phi(), jet_phi), sub1[ic].pt() * 0.001 );
	      for ( int ic = 0 ; ic< sub2.size() ; ic++)
		t_genTow2[nGenJetCounter]->Fill( sub2[ic].rapidity() - jet_rap, DeltaPhi( sub2[ic].phi(), jet_phi), sub2[ic].pt() * 0.001 );
	      sub1.clear();
	      sub2.clear();
	    }

	    
	  }
	  else  {
	    hGenSdStat->Fill( jet_pt, 0 ) ;	    
	  }
	}
      
	vpt_gen.push_back(jet_pt);
	veta_gen.push_back(jet_eta);
	vrap_gen.push_back(jet_rap);
	vphi_gen.push_back(jet_phi);
	vptRc_gen.push_back(thePtrc);
	vSdmass_gen.push_back(thesdm);
	vSdpt_gen.push_back(thesdpt);
	vSdz_gen.push_back(thesdz);
	vSdtheta_gen.push_back(thesdtheta);
	vSdpt1_gen.push_back(thesdpt1);
	vSdrap1_gen.push_back(thesdrap1);
	vSdphi1_gen.push_back(thesdphi1);
	vSdpt2_gen.push_back(thesdpt2);
	vSdrap2_gen.push_back(thesdrap2);
	vSdphi2_gen.push_back(thesdphi2);

	vTrMass_gen.push_back(t_genTrMass);
	vTrNsub_gen.push_back(t_genTrNsub);

	camJets.clear();
	
	// Charged Particle reclustering
	if (_saveLog)	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	if (_saveLog)	cout << "(Truth) Charged particle clustering starts" << endl;
        vector<fastjet::PseudoJet> chargeConsts;
        for ( int ic=0; ic< truthChargesAna.size() ; ic++) {
          if ( DeltaR ( jet_phi, jet_rap, truthChargesAna[ic].phi(), truthChargesAna[ic].rapidity() ) <  _ReclusterRadius ) {
            if ( truthChargesAna[ic].pt() * 0.001 < _pTtrkCutForSD ) 
	      continue;

	    chargeConsts.push_back( truthChargesAna[ic]) ;
	  }
        }
	fastjet::ClusterSequence reChCam(chargeConsts, jetDefCam);
        vector<fastjet::PseudoJet> camChJets = fastjet::sorted_by_pt(reChCam.inclusive_jets()); 
	if ( _saveLog)  {
	  cout << "(Truth) Charged Cambridge jets" << endl;
	  drawTreeHistory(reChCam) ;
	}
	// Charged SoftDrop
	// First find the charged particle SD
	
	float t_genChSdPt = -1 ;
	float t_genChSdMass= -1;
	float t_genChSdZ= -1;
	float t_genChSdTheta=-1;
	float t_genChTrMass = -1;

	hGenNchCam->Fill( jet_pt, camChJets.size() );
	if ( camChJets.size() > 0 )   {
	  fastjet::PseudoJet chSd_jet = softdropper(camChJets[0]);
	  t_genChSdMass = chSd_jet.m() * 0.001;
	  t_genChSdPt = chSd_jet.pt() * 0.001;
	  fastjet::PseudoJet parent1, parent2;
          if (  chSd_jet.has_parents(parent1, parent2) ) {
	    hGenSdChStat->Fill(jet_pt, 1);
	    t_genChSdTheta =  DeltaR( parent1.phi(), parent1.rapidity(), parent2.phi(), parent2.rapidity() ) ;
	    t_genChSdZ     = std::min( parent1.pt(),  parent2.pt() ) / ( parent1.pt() +  parent2.pt() ) ;
	    if (_saveLog)	    cout << " Charged Softdrop jet splits into:" << endl;
            if ( _saveLog)       logSubJets( parent1, parent2 ) ;
	  }
	  else {
	    hGenSdChStat->Fill(jet_pt, 0);
	    if (_saveLog)	    cout << "This Charged Softdrop jet is not divided!" << endl;
          }
        }
	
	vChSdPt_gen.push_back(t_genChSdPt);
	vChSdMass_gen.push_back(t_genChSdMass)  ;
	vChSdZ_gen.push_back(t_genChSdZ)  ;
	vChSdTheta_gen.push_back(t_genChSdTheta) ;
	
	vChTrMass_gen.push_back(t_genChTrMass);
	
	
	nGenJetCounter++; // THIS MUST BE AT THE END OF THE jets LOOP! 
      }
      truthParticles.clear();
      truthChargesAna.clear();
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
    if (_saveEvtDisplay) { 
      eventDisRecTow->Reset();
      eventDisRecTow1->Reset();
      eventDisRecTow2->Reset();
      eventDisGen->Reset();
      eventDisGen1->Reset();
      eventDisGen2->Reset();
    }
    //	  cout << " myJetsub is reset" << endl;
    //	  cout << " myJetSub.matchDr  = " << myJetSub.matchDr << endl;
    myJetSub.cent = cent_bin; 
    myJetSub.rhoCh = bge_rho_trk.rho() * 0.001;
    myJetSub.recoPt = vpt_reco[ri];
    myJetSub.recoRawPt = vptRaw_reco[ri];
    myJetSub.recoEta = veta_reco[ri];
    myJetSub.recoRap = vrap_reco[ri];
    myJetSub.recoPhi = vphi_reco[ri];
    myJetSub.recoRcPt = vptRc_reco[ri];
    myJetSub.weight = event_weight;
    myJetSub.recoSdPt = vSdpt_reco[ri];
    myJetSub.recoSdMass = vSdmass_reco[ri];
    myJetSub.recoSdTheta = vSdtheta_reco[ri];
    myJetSub.recoSdZ = vSdz_reco[ri];
    myJetSub.recoSdPt1 = vSdpt1_reco[ri];
    myJetSub.recoSdRap1 = vSdrap1_reco[ri];
    myJetSub.recoSdPhi1 = vSdphi1_reco[ri];
    myJetSub.recoSdPt2 = vSdpt2_reco[ri];
    myJetSub.recoSdRap2 = vSdrap2_reco[ri];
    myJetSub.recoSdPhi2 = vSdphi2_reco[ri];

    myJetSub.reChSdPt = vRecoChSdPt[ri];
    myJetSub.reChSdMass = vRecoChSdMass[ri];
    myJetSub.reChSdZ = vRecoChSdZ[ri];
    myJetSub.reChSdTheta = vRecoChSdTheta[ri];

    myJetSub.recoTrMass = vTrMass_reco[ri];
    myJetSub.recoTrNsub = vTrNsub_reco[ri];

    myJetSub.recoChTrMass = vChTrMass_reco[ri];

    if (_saveEvtDisplay) { 
      eventDisRecTow->Add(t_recTow[ri]);
      eventDisRecTow1->Add(t_recTow1[ri]);
      eventDisRecTow2->Add(t_recTow2[ri]);
      delete t_recTow[ri];
      delete t_recTow1[ri];
      delete t_recTow2[ri];
    }
    if (matchId != -1) {
      myJetSub.matchDr = drMin;
      myJetSub.genPt = vpt_gen[matchId];
      myJetSub.genEta = veta_gen[matchId];
      myJetSub.genRap = vrap_gen[matchId];
      myJetSub.genPhi = vphi_gen[matchId];
      myJetSub.genRcPt = vptRc_gen[matchId];
      myJetSub.genSdMass = vSdmass_gen[matchId];
      myJetSub.genSdPt = vSdpt_gen[matchId];
      myJetSub.genSdtheta = vSdtheta_gen[matchId];
      myJetSub.genSdz = vSdz_gen[matchId];
      myJetSub.genSdPt1 = vSdpt1_gen[matchId];
      myJetSub.genSdPhi1 = vSdphi1_gen[matchId];
      myJetSub.genSdRap1 = vSdrap1_gen[matchId];
      myJetSub.genSdPt2 = vSdpt2_gen[matchId];
      myJetSub.genSdPhi2 = vSdphi2_gen[matchId];
      myJetSub.genSdRap2 = vSdrap2_gen[matchId];

      myJetSub.genChSdPt = vChSdPt_gen[matchId];
      myJetSub.genChSdMass = vChSdMass_gen[matchId];
      myJetSub.genChSdZ = vChSdZ_gen[matchId];
      myJetSub.genChSdTheta = vChSdTheta_gen[matchId];

      myJetSub.genTrMass = vTrMass_gen[matchId];
      myJetSub.genTrNsub = vTrNsub_gen[matchId];
      myJetSub.genChTrMass = vChTrMass_gen[matchId];


      if (_saveEvtDisplay) {
	eventDisGen->Add(t_genTow[matchId]);
	eventDisGen1->Add(t_genTow1[matchId]);
	eventDisGen2->Add(t_genTow2[matchId]);
	delete t_genTow[matchId];
	delete t_genTow1[matchId];
	delete t_genTow2[matchId];
      }

    }
    
    treeOut->Fill();
  }
  
  
	
  
  
  
  // Clear vectors
  store->clear();
  delete store;

  corrected_selectedTrks.clear();
  selectedTrks.clear();
  vpt_reco.clear();
  vptRaw_reco.clear();
  veta_reco.clear();
  vphi_reco.clear();
  vptRc_reco.clear();
  vSdmass_reco.clear();
  vSdpt_reco.clear();
  vSdtheta_reco.clear();
  vSdz_reco.clear();
  
  vRecoChSdPt.clear();
  vRecoChSdMass.clear();
  vRecoChSdZ.clear();
  vRecoChSdTheta.clear();

  
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
	
	EL_RETURN_CHECK( "Finalize", m_trackSelectorTool->finalize() );
	if( m_trackSelectorTool ) delete m_trackSelectorTool; m_trackSelectorTool=0;


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

