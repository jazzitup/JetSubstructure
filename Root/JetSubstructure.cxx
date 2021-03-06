
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
#include <JetSubstructure/UEEstimator.h>
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



// this is needed to distribute the algorithm to the workers
ClassImp(JetSubstructure)


struct jetSubStr {
  Int_t cent;
  float weight, rhoCh;
  float recoMass, recoPt, recoRawPt, recoEta, recoRap, recoPhi, recoRcPt, recoSdPt, recoSdMass, recoSdZ, recoSdTheta, recoSdPt1, recoSdRap1, recoSdPhi1, recoSdPt2, recoSdRap2, recoSdPhi2, recoSdArea1, recoSdArea2;
  float recoChPt, recoChMass, recoChPtRaw, recoChMassRaw, recoChMassGm, recoChSdPt, recoChSdMass, recoChSdZ, recoChSdTheta;
  float matchDr, genMass, genPt, genPt2,  genEta, genRap, genPhi, genRcPt,  genSdPt, genSdMass, genSdz, genSdtheta, genSdPt1, genSdRap1, genSdPhi1, genSdPt2, genSdRap2, genSdPhi2, genSdArea1, genSdArea2 ;
  float genNch, genChSdPt, genChSdMass, genChSdZ, genChSdTheta;

  
  float nTrkRaw, nTrkBkg, nTrkBkgNoWgt,  recoChPtRcSubt, recoChMassRcSubt, drTrkJetBkg, maxTrkPt ;   // random cone
  float fcalet;
};
jetSubStr myJetSub;

TString evtText = "cent/I:weight/F:rhoCh";
TString recoText = "mass:pt:rawPt:eta:y:phi:rcPt:sdPt:sdMass:zg:theta:spt1:sy1:sphi1:spt2:sy2:sphi2:sda1:sda2";
TString reChText = "chPtCSubt:chMassCSubt:chPtRaw:chMassRaw:chMassGm:chSdPt:chSdMass:chZg:chTheta";
TString genText = "dr:genMass:genPt:genPt2:genEta:genRap:genPhi:genRcPt:genSdPt:genSdMass:genZg:genTheta:genSpt1:genSy1:genSphi1:genSpt2:genSy2:genSphi2:genSda1:genSda2";
TString genChText = "genNch:genChSdPt:genChSdMass:genChZg:genChTheta";
TString rcText  = "nTrkRaw:nTrkBkg:nTrkBkgNoWgt:chPtRcSubt:chMassRcSubt:drTrkJetBkg:maxTrkPt";
TString extraText = "fcal";
TString myJetSubText = evtText+":"+recoText+":"+reChText+":"+genText+":"+genChText+":"+rcText+":"+extraText;


void resetSubstr (jetSubStr &jetsub)
{
  
  jetsub.cent = -1 ;
  jetsub.recoMass= -1;
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
  jetsub.recoSdArea1 = -1;
  jetsub.recoSdArea2 = -1;

  jetsub.genPt=-1;
  jetsub.genPt2=-1;
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
  jetsub.genSdArea1 = -1;
  jetsub.genSdArea2 = -1;
    
  jetsub.recoChPt = -1;
  jetsub.recoChMass = -10;
  jetsub.recoChPtRaw = -10;
  jetsub.recoChMassRaw = -10;
  jetsub.recoChPtRcSubt = -10;
  jetsub.recoChMassRcSubt = -10;
  jetsub.recoChMassGm = -10;
  
  jetsub.genMass = -1;
  jetsub.fcalet = -1;
  
  
}

struct EvtInfo  { 
  int run;
  int lumi;
  int event;
};

EvtInfo myEvt;


int getAkshatTypeTruth(int barcode, int pdg, int status, float charge)
{
  // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TrackingCPMoriond2016#Truth_definitions
  // 0 - fake (should not happened) , 1 - truth primary, 2 - truth secondary, 3 - truth primary out-of-phase-space, 4 - truth secondary neutral, 5 - truth primary strange baryons
  
  if(0<barcode && barcode<200000)
    {
      if(fabs(pdg)==3112 || fabs(pdg)==3222 || fabs(pdg)==3312 || fabs(pdg)==3334) return 5;
      else
	{
	  if(status==1 && fabs(charge)>0.5) return 1;
	  else return 3;
	}
    }
  
  if(200000<=barcode)
    {
      if(fabs(charge)>0.5) return 2;
      else return 4;
    }
  
  return 0; // to be safe
}



bool getYsTruthType(int barcode, int pdg, int status) 
{
  // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TrackingCPMoriond2016#Truth_definitions
  // 0 - fake (should not happened) , 1 - truth primary, 2 - truth secondary, 3 - truth primary out-of-phase-space, 4 - truth secondary neutral, 5 - truth primary strange baryons
  
  if (0<barcode && barcode<200000) { 
    if (fabs(pdg)==3112 || fabs(pdg)==3222 || fabs(pdg)==3312 || fabs(pdg)==3334 )
      return true; 
    else
      {
	if (status==1)  
	  return true; 
      }
  }
  else 
    return false;
}




void resetEvt (EvtInfo &theEvt)  {
  theEvt.run = -1;
  theEvt.lumi = -1;
  theEvt.event = -1;
}
TString evtString = "run/I:lumi:event";

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


float GetEventPlaneUsingEventShape(const xAOD::HIEventShapeContainer* calos){
  float psi_2;
  {
    float FCal_Et=(calos->at(5)->et()*0.001*0.001);
    double qx_2=(calos->at(5)->etCos().at(1));
    double qy_2=(calos->at(5)->etSin().at(1));
    float N_psi_2 = std::atan2(qy_2,qx_2);
    psi_2 = N_psi_2/2.0;
    //break;
  }
  //}
  return psi_2;
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
	if (_doJES) {
	  for ( int si = 0 ; si<=6; si++) { 
	    ptSysHI[si] = 0;
	    treeOut->Branch(Form("ptSysHI%d",si), &ptSysHI[si]); 
	  }
	  for ( int si = 0 ; si<=21; si++) {  // 20 - 41
	    ptSysPP[si] = 0;
	    treeOut->Branch(Form("ptSysPP%d",si), &ptSysPP[si]); 
	  }
	}
	// for jet mass systematics 
	treeOut->Branch("trkJetMass4",&trkJetMass4);
	treeOut->Branch("trkJetMass6",&trkJetMass6);
	treeOut->Branch("trkJetMass8",&trkJetMass8);
	treeOut->Branch("trkJetMass10",&trkJetMass10);
	treeOut->Branch("trkJetPt4",&trkJetPt4);
	treeOut->Branch("trkJetPt6",&trkJetPt6);
	treeOut->Branch("trkJetPt8",&trkJetPt8);
	treeOut->Branch("trkJetPt10",&trkJetPt10);

	treeOut->Branch("trkJetMassRcSub2",&trkJetMassRcSub2);
	treeOut->Branch("trkJetMassRcSub4",&trkJetMassRcSub4);
	treeOut->Branch("trkJetMassRcSubEV",&trkJetMassRcSubEV);
	treeOut->Branch("trkJetMassRcSubMV",&trkJetMassRcSubMV);
	
	
	treeOut->Branch("evt",&myEvt,evtString.Data());
	
	
	
	if (_saveEvtDisplay) {
	  eventDisRecTow = new TH2F("eventDisRecTow","",20,-1,1,20,-1,1);
	  eventDisRecTow1 = (TH2F*)eventDisRecTow->Clone("eventDisRecTow1");
	  eventDisRecTow2 = (TH2F*)eventDisRecTow->Clone("eventDisRecTow2");
	  eventDisGen = new TH2F("eventDisGen","",20,-1,1,20,-1,1);
	  eventDisGen1 = (TH2F*)eventDisGen->Clone("eventDisGen1");
	  eventDisGen2 = (TH2F*)eventDisGen->Clone("eventDisGen2");
	  eventDisTrk = new TH2F("eventDisTrk","",20,-1,1,20,-1,1);
          eventDisTrk1 = (TH2F*)eventDisTrk->Clone("eventDisTrk1");
          eventDisTrk2 = (TH2F*)eventDisTrk->Clone("eventDisTrk2");
	  eventDisChg = new TH2F("eventDisChg","",20,-1,1,20,-1,1);
          eventDisChg1 = (TH2F*)eventDisChg->Clone("eventDisChg1");
          eventDisChg2 = (TH2F*)eventDisChg->Clone("eventDisChg2");

	  treeOut->Branch("rect","TH2F",&eventDisRecTow,256000,0);
	  treeOut->Branch("rect1","TH2F",&eventDisRecTow1,256000,0);
	  treeOut->Branch("rect2","TH2F",&eventDisRecTow2,256000,0);
	  treeOut->Branch("genp","TH2F",&eventDisGen,256000,0);
	  treeOut->Branch("genp1","TH2F",&eventDisGen1,256000,0);
	  treeOut->Branch("genp2","TH2F",&eventDisGen2,256000,0);
	  treeOut->Branch("trk","TH2F",&eventDisTrk,256000,0);
	  treeOut->Branch("trk1","TH2F",&eventDisTrk1,256000,0);
	  treeOut->Branch("trk2","TH2F",&eventDisTrk2,256000,0);
	  treeOut->Branch("chg","TH2F",&eventDisChg,256000,0);
	  treeOut->Branch("chg1","TH2F",&eventDisChg1,256000,0);
	  treeOut->Branch("chg2","TH2F",&eventDisChg2,256000,0);
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

	  temphist_2d = new TH2D(Form("hTrkPtEta_preCS_cent%i",i),";pT;eta",400,0,100,20,-3,3);
	  hTrkPtEta_preCS_cent.push_back(temphist_2d);
	  hTrkPtEta_preCS_cent.at(i)->Sumw2();

	  temphist_2d = new TH2D(Form("hTrkPtEta_postCS_cent%i",i),";pT;eta",400,0,100,20,-3,3);
	  hTrkPtEta_postCS_cent.push_back(temphist_2d);
	  hTrkPtEta_postCS_cent.at(i)->Sumw2();

	  temphist_2d = new TH2D(Form("hTrkPtEta_genMatch_cent%i",i),";pT;eta",400,0,100,20,-3,3);
	  hTrkPtEta_genMatch_cent.push_back(temphist_2d);
	  hTrkPtEta_genMatch_cent.at(i)->Sumw2();

	  temphist_2d = new TH2D(Form("h_trkPt_trkBkgPt_cent%i",i),"",
                                 200,0,100, 220, -2,20);
	  h_trkPt_trkBkgPt_cent.push_back(temphist_2d);
	  h_trkPt_trkBkgPt_cent.at(i)->Sumw2();

	  temphist_2d = new TH2D(Form("h_trkPt_trkBkgPt_jetCone_cent%i",i),"",
                                 200,0,100, 220, -2,20);
	  h_trkPt_trkBkgPt_jetCone_cent.push_back(temphist_2d);
	  h_trkPt_trkBkgPt_jetCone_cent.at(i)->Sumw2();


	  temphist_2d = new TH2D(Form("h_bkgSubt_prePt_postPt_cent%i",i),"",
                                 1000, 0,100, 1000,0,100);
	  h_bkgSubt_prePt_postPt_cent.push_back(temphist_2d);
	  h_bkgSubt_prePt_postPt_cent.at(i)->Sumw2();
	  
	  temphist_2d = new TH2D(Form("h_bkgSubt_prePt_postPt_jetCone_cent%i",i),"",
                                 1000, 0,100, 1000,0,100);
	  h_bkgSubt_prePt_postPt_jetCone_cent.push_back(temphist_2d);
	  h_bkgSubt_prePt_postPt_jetCone_cent.at(i)->Sumw2();


	  temphist_2d = new TH2D(Form("h_dRSubt_trkPt_cent%i",i),"",
				 110,0,0.0000011,1000,0,100);
	  h_dRSubt_trkPt_cent.push_back(temphist_2d);
	  h_dRSubt_trkPt_cent.at(i)->Sumw2();
	  
	  temphist_2d = new TH2D(Form("h_dRSubt_trkPt_jetCone_cent%i",i),"",
				 110,0,0.0000011,1000,0,100);
	  h_dRSubt_trkPt_jetCone_cent.push_back(temphist_2d);
	  h_dRSubt_trkPt_jetCone_cent.at(i)->Sumw2();
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

	if ( _saveNtuple ) wk()->addOutput (treeOut);

	//	wk()->addOutput (hGenNcam);
	//	wk()->addOutput (hGenNchCam);
	//	wk()->addOutput (hGenSdStat);
	//	wk()->addOutput (hGenSdChStat);
	//	wk()->addOutput (hRecoNcam);
	//	wk()->addOutput (hRecoNchCam);
	//	wk()->addOutput (hRecoSdStat);
	//	wk()->addOutput (hRecoSdChStat);
	
	if ( 1 == 0 ) {
	  for (int i=0;i<nCentbins;i++)  {
	    wk()->addOutput (h_trkGen_pt_dphi_cent.at(i));
	    wk()->addOutput (h_allGen_pt_dphi_cent.at(i));
	    wk()->addOutput (h_trkGen_pt_drap_cent.at(i));
	    wk()->addOutput (h_allGen_pt_drap_cent.at(i));
	    wk()->addOutput (hTrkPtEta_preCS_cent.at(i));
	    wk()->addOutput (hTrkPtEta_postCS_cent.at(i));
	    wk()->addOutput (hTrkPtEta_genMatch_cent.at(i));
	    wk()->addOutput (h_trkPt_trkBkgPt_cent.at(i));
	    wk()->addOutput (h_bkgSubt_prePt_postPt_cent.at(i));
	    wk()->addOutput (h_dRSubt_trkPt_cent.at(i));
	    wk()->addOutput (h_trkPt_trkBkgPt_jetCone_cent.at(i));
	    wk()->addOutput (h_bkgSubt_prePt_postPt_jetCone_cent.at(i));
	    wk()->addOutput (h_dRSubt_trkPt_jetCone_cent.at(i));
	  }
	}
	
        cout << "=======================================" << endl;
	cout << "Parametrization of d0 cut" << endl;
	f_d0_cut = new TF1("f1", "[0]*exp([1]*x)+[2]*exp([3]*x)", 0.4, 500);
        f_d0_cut->SetParameters(0.472367, -0.149934, 0.193095, 0.000337765);

	// *~*~*~*~*~UEE estimatior *~*~*~*~*~*~*~*~*~*~*~*~*~* //
	uee = new UEEstimator();
	uee->ptBkgrThreshold = 10; // cones with track with _trkptBkgrThreshold GeV excluded
	uee->jetptBkgrThreshold = 90; // cones with proximity of jet with _jetptBkgrThreshold GeV excluded
	uee->m_maxjetdeltaR = 0.8; // definition of "proximity" to a real jet




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
	//	const std::string name = "JetSubstructure"; //string describing the current thread, for logging
	//	TString jetAlgo = "AntiKt4HI"; //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo (see below)
	//	TString config = "JES_MC15c_HI_Nov2016.config"; //Path to global config used to initialize the tool (see below)
	//	TString calibSeq = "EtaJES"; //String describing the calibration sequence to apply (see below)

	//Call the constructor. The default constructor can also be used if the arguments are set with python configuration instead
	//	m_jetCalibration = new JetCalibrationTool (name);
	//	ANA_CHECK(m_jetCalibration->setProperty("JetCollection",jetAlgo.Data()));
	//	ANA_CHECK(m_jetCalibration->setProperty("ConfigFile",config.Data()));
	//	ANA_CHECK(m_jetCalibration->setProperty("CalibSequence",calibSeq.Data()));
	//	ANA_CHECK(m_jetCalibration->setProperty("IsData",!_isMC));
	//	ANA_CHECK(m_jetCalibration->initializeTool(name));
		
	//in-situ calib
	TString jetAlgo_insitu = "AntiKt4EMTopo"; //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo (see below)
	TString config_insitu = "JES_MC15cRecommendation_May2016_xCalib.config"; //Path to global config used to initialize the tool (see below)
	const std::string name_insitu = "insitu"; //string describing the current thread, for logging
	TString calibSeq_insitu = "Insitu_DEV"; //String describing the calibration sequence to apply (see below)
	
	m_jetCalibration_insitu = new JetCalibrationTool(name_insitu, jetAlgo_insitu, config_insitu, calibSeq_insitu, true);
	EL_RETURN_CHECK("initialize()",m_jetCalibration_insitu->initializeTool(name_insitu));

	
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
	
	
	// Jet cleaner for pp 
	m_jetCleaningToolHandle.setTypeAndName("JetCleaningTool/JetCleaning"); 
	ANA_CHECK(m_jetCleaningToolHandle.setProperty("CutLevel","LooseBad"));
	ANA_CHECK(m_jetCleaningToolHandle.setProperty("DoUgly", false));
	ANA_CHECK(m_jetCleaningToolHandle.retrieve());
	
	
	
	// Tracking Unc. 
	TFile* f_eff_uncert_2015_material = new TFile(xfn + "/../JetSubstructure/data/TrackingEfficiencyRecommendations_20.7rel.root");
	h_eff_uncert_2015_material[0] = (TH2F*)f_eff_uncert_2015_material->Get("OneMinusRatioEfficiencyVSEtaPt_AfterRebinning_Old_Nominal_MCVSOld_5%Extra_MC_TightPrimary"); 
	h_eff_uncert_2015_material[1]    = (TH2F*)f_eff_uncert_2015_material->Get("OneMinusRatioEfficiencyVSEtaPt_AfterRebinning_Old_Nominal_MCVSNew_Nominal_MC_TightPrimary");
	h_eff_uncert_2015_material[2]    = (TH2F*)f_eff_uncert_2015_material->Get("OneMinusRatioEfficiencyVSEtaPt_AfterRebinning_Old_Nominal_MCVSOld_50%ExtraPP0_MC_TightPrimary");
	h_eff_uncert_2015_material[3]    = (TH2F*)f_eff_uncert_2015_material->Get("OneMinusRatioEfficiencyVSEtaPt_AfterRebinning_Old_Nominal_MCVSOld_FTF_BIC_MC_TightPrimary");

	TFile* f_fake_uncert = new TFile(xfn + "/../JetSubstructure/data/fake_uncert_variation.root");
	h_fake_uncert[0] = (TH2D*)f_fake_uncert->Get("fake_uncert_jety0p3");
	h_fake_uncert[1] = (TH2D*)f_fake_uncert->Get("fake_uncert_jety0p3_0p8");
	h_fake_uncert[2] = (TH2D*)f_fake_uncert->Get("fake_uncert_jety0p8_1p2");
	h_fake_uncert[3] = (TH2D*)f_fake_uncert->Get("fake_uncert_jety1p2_2p1");
	
	TFile* f_sagitta = new TFile(xfn + "/../JetSubstructure/data/5TeVHI2015_sagittaBias_pTmethod_statUncertainty.root");
	h_sagitta = (TH2D*)f_sagitta->Get("h_deltaSagittaMap_statErr");
	
	hRandomUnit = new TH1D("hRandomUnit","",1,0,1);
	hRandomUnit->Fill(0.5);
	
	
	// *~*~*~*~*~JES/JER uncertainty *~*~*~*~*~*~*~* //
	float _mcProbCut = 0.5;
	bool _eff_jety = false;
	
	if ( _doJES ) {
	  vUncertIndex.push_back(1);
	  vUncertIndex.push_back(6);
	  vUncertIndex.push_back(7);
	  vUncertIndex.push_back(8);
	  vUncertIndex.push_back(9);
	  vUncertIndex.push_back(16);
	  vUncertIndex.push_back(17); 

	  cout << "number of uncertainty factors = " << vUncertIndex.size() << endl;
	  for ( int ii = 0 ; ii< vUncertIndex.size() ; ii++) {
	    UncertProvider *tempUncet = new UncertProvider(vUncertIndex.at(ii),_mcProbCut,"_cut_level.c_str()", 30 , _eff_jety);
	    vUncProvHI.push_back(tempUncet);
	  }
	
	  // pp intrinsic 
	  intrinsicComponent[0] = 0;  
	  intrinsicComponent[1] = 0;
	  intrinsicComponent[2] = 1;  
	  intrinsicComponent[3] = 1;
	  intrinsicComponent[4] = 8;  
	  intrinsicComponent[5] = 8;
	  intrinsicComponent[6] = 9;  
	  intrinsicComponent[7] = 9;
	  intrinsicComponent[8] = 10;  
	  intrinsicComponent[9] = 10;
	  intrinsicComponent[10] = 11;  
	  intrinsicComponent[11] = 11;
	  intrinsicComponent[12] = 12;  
	  intrinsicComponent[13] = 12;
	  intrinsicComponent[14] = 13;
	  intrinsicComponent[15] = 13;
	  intrinsicComponent[16] = 14;
	  intrinsicComponent[17] = 14;
	  intrinsicComponent[18] = 15;
	  intrinsicComponent[19] = 15;
	  intrinsicComponent[20] = 16;
	  intrinsicComponent[21] = 16;
	  for ( int ii = 0 ; ii<=21 ; ii++)  {
	    if ( (ii%2) == 0 )  intSignificance[ii] = 1. ;
	    else intSignificance[ii] = -1. ;
	  }
	  
	  
	}
	//	vUncProvHI = new UncertProvider(20,_mcProbCut,"_cut_level.c_str()", 30 , _eff_jety);
	
	
	return EL::StatusCode::SUCCESS;
}


EL::StatusCode JetSubstructure :: execute ()
{
  bool useReAntiKt = false; // Set false!! 

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
  int cent_bin_scheme30 = 0;
  //  double event_weight_fcal = 1;
  //Centrality
  const xAOD::HIEventShapeContainer* calos=0;
  ANA_CHECK(event->retrieve( calos, "CaloSums"));
  FCalEt=calos->at(5)->et()*1e-6;
  if (_centrality_scheme>1)	  {
    cent_bin = GetCentralityBin(_centrality_scheme, FCalEt, isHIJING);
    cent_bin_scheme30 = GetCentralityBin(30, FCalEt, isHIJING);
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
    //    m_zdcTools->reprocessZdc();	  // ZDC
    //    m_is_pileup = m_hiPileup->is_pileup( *calos, *zdcMod); // SAVE pileup Decision HERE 0 = NO pileup, 1 = pileup
    //    cout << "here 4 " << endl;
  }
  //  else m_is_pileup = (FCalEt > 4.8); //Remove pileup in MC
  //  if (m_is_pileup){
  //  h_RejectionHisto->Fill(6.5);
  //    keep = false;
  //  }
  
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

  TH2F* t_trkTow[20];
  TH2F* t_trkTow1[20];
  TH2F* t_trkTow2[20];
  TH2F* t_chgTow[20];
  TH2F* t_chgTow1[20];
  TH2F* t_chgTow2[20];


  vector <double> vmass_reco;
  vector <double> vpt_reco;
  vector <double> vptRaw_reco;  
  vector <double> vetaRaw_reco;  
  vector <double> vphiRaw_reco;  
  vector <double> vNconst_reco;
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
  
  vector <float> vChSdPt_reco; 
  vector <float> vChSdMass_reco; 
  vector <float> vChSdZ_reco; 
  vector <float> vChSdTheta_reco; 

  vector <float> vChPt_reco; 
  vector <float> vChMass_reco; 
  vector <float> vChPtRaw_reco; 
  vector <float> vChMassRaw_reco; 
  vector <float> vChPtRcSubt_reco; 
  vector <float> vChMassRcSubt_reco; 
  vector <float> vChMassRcSubt2_reco; 
  vector <float> vChMassRcSubt4_reco; 
  vector <float> vChMassRcSubtEV_reco; 
  vector <float> vChMassRcSubtMV_reco; 


  vector <float> vChMassGm_reco; 

  vector <double> vMass_gen;
  vector <double> vpt_gen;
  vector <double> vpt2_gen;
  vector <double> veta_gen;
  vector <double> vrap_gen;
  vector <double> vphi_gen;
  vector <double> vptRc_gen;
  vector <double> vSdmass_gen;
  vector <double> vSdPt_gen;
  vector <float> vSdz_gen;
  vector <float> vSdTheta_gen;
  vector <float> vSdpt1_gen;
  vector <float> vSdrap1_gen;
  vector <float> vSdphi1_gen;
  vector <float> vSdpt2_gen;
  vector <float> vSdrap2_gen;
  vector <float> vSdphi2_gen;



  vector <float> vNch_gen;
  vector <float> vChSdPt_gen;
  vector <float> vChSdMass_gen  ;
  vector <float> vChSdZ_gen  ;
  vector <float> vChSdTheta_gen  ;

  vector <float> vNChRaw;
  vector <float> vNChBkg;
  vector <float> vNChBkgNoWgt;
  vector <float> vDrChJetBkg;
  vector <float> vMaxTrkPt;


  // algorithm definition 
  fastjet::JetDefinition jetDefReclus(fastjet::cambridge_algorithm, _JetRadiusAna);
  if ( _defJetRecl == 0)   jetDefReclus.set_jet_algorithm( fastjet::cambridge_algorithm ) ;
  else if ( _defJetRecl == 1)   jetDefReclus.set_jet_algorithm( fastjet::kt_algorithm ) ;

  fastjet::JetDefinition jetDefAk(fastjet::antikt_algorithm, _JetRadiusAna);
  fastjet::JetDefinition jetDefAk04(fastjet::antikt_algorithm, 0.4);  // Used for event reweighting factor
  fastjet::contrib::SoftDrop softdropper(_beta, _z_cut);
  if ( _saveLog) cout << softdropper.description() << endl;
  
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
     
      //      if ( getYsTruthType( thebc, (*truth_itr)->pdgId(), (*truth_itr)->status() ) == false )
      //	continue; 

      //      if (  (*truth_itr)->status() !=1 ) continue;


      //      if ( status !=1 ) continue; 
      if  ( !( (0<thebc)  && (thebc < 200000) ) ) continue;
      
      truthParticles.push_back( (*truth_itr)->p4() );   // To be used for anti-kT jet reconstruction 
      
      if ( (fabs((*truth_itr)->charge()) > 0 ) && ( ptTrk > _pTtrkCutTruth ) && ( fabs((*truth_itr)->p4().Eta()) < _etaTrkCut )  ) 
	truthChargesAna.push_back( (*truth_itr)->p4() ) ; // For softdrop
      
    }
  }
  
  
  
  double event_weight = 1;
  if (_isMC){
    ////////// On-the-fly jet finder ///////////////////////////////////////////////
    if (_saveLog) {
      cout << "   * Reweighting factor using R=0.4 jet" << endl;
    }
    const xAOD::JetContainer* truth_jets = 0;
    ANA_CHECK(event->retrieve(truth_jets, _truth_jet_collection.c_str() ));
    xAOD::JetContainer::const_iterator truth_jet_itr = truth_jets->begin();
    xAOD::JetContainer::const_iterator truth_jet_end = truth_jets->end();

    double maxGenPt = 0;
    
    for( ; truth_jet_itr != truth_jet_end; ++truth_jet_itr ) {
      xAOD::JetFourMom_t jet_truth_4mom = (*truth_jet_itr)->jetP4();
      
      double pt     = (jet_truth_4mom.pt() * 0.001 );
      double eta    = (jet_truth_4mom.eta());
      double phi    = (jet_truth_4mom.phi());
      if ( pt > maxGenPt )  {
	maxGenPt = pt; 
	event_weight = jetcorr->GetJetWeight(pt, eta, phi);
      }
    }

    if (_saveLog)  {
      cout << "      * Max pT (and y,phi) for R=0.4 Jets = " << maxGenPt << " GeV, " << endl;
    }
    
  }
  
  
  
  //  find the jet cone for exclusion area in background estimation procedure 
  xAOD::TStore *store = new xAOD::TStore; //For calibration
  
  //  xAOD::JetContainer* updatedjets = new xAOD::JetContainer();
  //  xAOD::AuxContainerBase* updatedjetsAux = new xAOD::AuxContainerBase();
  //  updatedjets->setStore( updatedjetsAux );   
  //  store->record(updatedjets,"updatedjets");
  //  store->record(updatedjetsAux,"updatedjetsAux");

  const xAOD::JetContainer* reco_jets_forExclusion = 0;
  ANA_CHECK(event->retrieve( reco_jets_forExclusion, "DFAntiKt4HIJets" ) );

  xAOD::JetContainer::const_iterator jetExclu_itr = reco_jets_forExclusion->begin();
  xAOD::JetContainer::const_iterator jetExclu_end = reco_jets_forExclusion->end();
  std::vector<float> etaJetExclu;
  std::vector<float> phiJetExclu;

  std::vector<float> _vJetPtForRC;
  std::vector<float> _vJetEtaForRC;
  std::vector<float> _vJetPhiForRC;

  for( ; jetExclu_itr != jetExclu_end; ++jetExclu_itr ) {
    xAOD::Jet theRecoJet;
    theRecoJet.makePrivateStore( **jetExclu_itr );
    
    const xAOD::JetFourMom_t jet_4momCalib = theRecoJet.jetP4();   // CALIBRATED!!!
    float jet_pt =0; 
    float jet_eta =0;
    float jet_phi =0 ;

    if ( _isMC) {   
      jet_pt  = jet_4momCalib.pt() * 0.001 ;
      jet_eta = jet_4momCalib.eta();
      jet_phi = jet_4momCalib.phi();

      const xAOD::JetFourMom_t jet_4mom_unsubtracted = theRecoJet.jetP4("JetUnsubtractedScaleMomentum");
      theRecoJet.setJetP4("JetConstitScaleMomentum",jet_4mom_unsubtracted); //Required
      //      const xAOD::JetFourMom_t jet_4mom_subtracted_uncalib = theRecoJet.jetP4("JetPileupScaleMomentum");
      if ( _saveLog) { 
	cout << "jet_pt = " << jet_pt*0.001 << endl;
	cout << "Unsubtrcted & uncalibrated pt : " << jet_4mom_unsubtracted.pt() * 0.001 << endl;
	//	cout << "Subtrcted & uncalibrated pt : " << jet_4mom_subtracted_uncalib.pt() * 0.001 << endl;
      }
      
    }
    else {  // if (!_isMC)
      const xAOD::JetFourMom_t jet_4mom_unsubtracted = theRecoJet.jetP4("JetUnsubtractedScaleMomentum");
      theRecoJet.setJetP4("JetConstitScaleMomentum",jet_4mom_unsubtracted); //Required

      const xAOD::JetFourMom_t jet_4mom_EMJES = theRecoJet.jetP4(); //This is default, up to EM+JES calibrated four momentum
      theRecoJet.setJetP4("JetGSCScaleMomentum", jet_4mom_EMJES);



      EL_RETURN_CHECK("execute()", m_jetCalibration_insitu->applyCalibration( theRecoJet ) );
      const xAOD::JetFourMom_t jet_4momCroCal = theRecoJet.jetP4(); // uncalib
           

      jet_pt  = jet_4momCroCal.pt() * 0.001 ;
      jet_eta = jet_4momCroCal.eta();
      jet_phi = jet_4momCroCal.phi();
    }
    
    _vJetPtForRC.push_back(jet_pt);  // GeV!
    _vJetEtaForRC.push_back(jet_eta);
    _vJetPhiForRC.push_back(jet_phi);

    if (jet_pt < _ptCutJetConeExc)
      continue;
    if (fabs(jet_eta) > _etaTrkCut - 0.4)    // because we are using DFAntiKt4HIJets
      continue;  

    etaJetExclu.push_back(jet_eta);
    phiJetExclu.push_back(jet_phi);
  }
  if (_saveLog)  cout << endl << " There are " << etaJetExclu.size() << " exclusion cones (pT > " << _ptCutJetConeExc << " GeV)" << endl;

  //  find the reco jet cones for analysis 

  const xAOD::JetContainer* reco_jets = 0;
  
  ANA_CHECK(event->retrieve( reco_jets, _reco_jet_collection.c_str() ));
  xAOD::JetContainer::const_iterator jetcone_itr = reco_jets->begin();
  xAOD::JetContainer::const_iterator jetcone_end = reco_jets->end();

  std::vector<float> etaJetCone150;
  std::vector<float> phiJetCone150;
  for( ; jetcone_itr != jetcone_end; ++jetcone_itr ) {
    xAOD::Jet theRecoJet;
    theRecoJet.makePrivateStore( **jetcone_itr );
    const xAOD::JetFourMom_t jet_4momCalib = theRecoJet.jetP4();   // CALIBRATED!!!
    float jet_pt  = jet_4momCalib.pt() * 0.001 ;
    float jet_eta = jet_4momCalib.eta();
    float jet_phi = jet_4momCalib.phi();
    if (jet_pt < 150)     
      continue; // hard coded
    if (fabs(jet_eta) > _etaJetCut )  
      continue;
    
    etaJetCone150.push_back(jet_eta);
    phiJetCone150.push_back(jet_phi);
  }
  
  
  
  
  

  
  ///////////// tracks ////////////////////////////////////////
  std::vector<fastjet::PseudoJet> selectedTrks;
  std::vector<fastjet::PseudoJet> selAndExcldTrks; // after jet cone exclusion
  std::vector<fastjet::PseudoJet> selGenMatchTrks;
  std::vector<fastjet::PseudoJet> selectedTrksMV; // momentum varied
  
  std::vector<float> _vTrkPtForRC;
  std::vector<float> _vTrkEtaForRC;
  std::vector<float> _vTrkPhiForRC;


  const xAOD::TrackParticleContainer* recoTracks = 0;
  EL_RETURN_CHECK("execute",event->retrieve( recoTracks, "InDetTrackParticles"));
  for (const auto& trk : *recoTracks) {
    float pt = trk->pt()/1000.;
    float eta = trk->eta();
    float phi = trk->phi();
    if (_saveLog)   cout << " ==== Track pt, eta, phi = " << pt <<", "<<eta<<", "<<phi<<endl;
    
    if ( pt < _pTtrkCutReco ) 
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
    
    _vTrkPtForRC.push_back(pt);  // GeV!!!!
    _vTrkEtaForRC.push_back(eta);
    _vTrkPhiForRC.push_back(phi);
    
    double pionMass = _chParticleMassMeV; // 139.570 ; // in MeV
    double pionx = trk->p4().Px(); // Oct 24th, Yongsun confiremd p4() is in MeV unit
    double piony = trk->p4().Py();
    double pionz = trk->p4().Pz();
    double pione = sqrt (pionx*pionx + piony*piony + pionz*pionz + pionMass*pionMass) ;
    
    // Momentum varied tracks by Saggita bias 
    float charge_tmp= trk->charge();
    if (charge_tmp>0.) charge_tmp=1;
    else if (charge_tmp<0.) charge_tmp =-1;
    
    float eta_tmp=eta;
    if(eta_tmp>2.499 )  eta_tmp= 2.499;
    if(eta_tmp<-2.499)  eta_tmp=-2.499;

    double orgPt = trk->p4().Pt() * 0.001; // GeV!!!!
    double saggitaRatio = 1. / ( 1. + charge_tmp * orgPt * (h_sagitta->GetBinContent(h_sagitta->FindBin(eta_tmp, phi))) * 0.001 ); // 0.001 has nothing to do with GeV->MeV conversion.
    //    cout << "charge = " << charge_tmp << endl;  cout << "orgPt = " << orgPt << " GeV" << endl; cout << "ratio = " << saggitaRatio << endl;
    double pionxMV = pionx*saggitaRatio;
    double pionyMV = piony*saggitaRatio;
    double pionzMV = pionz*saggitaRatio;
    double pioneMV = sqrt (pionxMV*pionxMV + pionyMV*pionyMV + pionzMV*pionzMV + pionMass*pionMass) ;
    
    
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
    selectedTrksMV.push_back( fastjet::PseudoJet ( pionxMV, pionyMV, pionzMV, pioneMV ) );
    
    if (isTruthMatched) selGenMatchTrks.push_back( matchPar );
    

    
    /////// Jet cone excluded  //////
    bool isExcluded = false;
    for ( int jj =0 ; jj < etaJetExclu.size() ; jj++ ) {
      float drTrkJet = DeltaR( phiJetExclu[jj], etaJetExclu[jj], phi, eta );
      if ( drTrkJet < 0.4 )   
	isExcluded = true; 
    }
    if ( isExcluded == false ) {
      selAndExcldTrks.push_back(fastjet::PseudoJet ( pionx, piony, pionz, pione )) ;
    }
    ////////////////////  


}
  if (_saveLog) cout << " number of reconstructed tracks: " << selectedTrks.size() << endl;
  if (_saveLog) cout << "         Out of exclusion cone : " << selAndExcldTrks.size() << endl;

    
  //Get reaction plane (maybe not needed for first iteration):
  uee->Psi = GetEventPlaneUsingEventShape(calos);
  uee->ExcludeConesByJetandTrack(_vTrkPtForRC,_vTrkEtaForRC,_vTrkPhiForRC,_vJetPtForRC,_vJetEtaForRC,_vJetPhiForRC);
  
  
  // Background for tracks 
  // See https://github.com/cms-externals/fastjet-contrib/blob/master/ConstituentSubtractor/example.cc  
  // algo for background.  
  //  double _ghost_area=0.01;  
  //  double _Rktjet_bkg = 0.4; 
  fastjet::JetDefinition jet_def_for_rho(fastjet::kt_algorithm, _Rktjet_bkg);
  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts,fastjet::GhostedAreaSpec(_etaTrkCut,1)); // the area definiton is used only for the jet backgroud estimator. It is not important for the ConstituentSubtractor when subtracting the whole event - this is not true when subtracting the individual jets

  fastjet::Selector rho_range_trk =  fastjet::SelectorAbsEtaMax(_etaTrkCut - _Rktjet_bkg);
  fastjet::JetMedianBackgroundEstimator bge_rho_trk(rho_range_trk, jet_def_for_rho, area_def);
  fastjet::BackgroundJetScalarPtDensity *scalarPtDensity = new fastjet::BackgroundJetScalarPtDensity();
  bge_rho_trk.set_jet_density_class(scalarPtDensity); // this changes the computation of pt of patches from vector sum to scalar sum. The scalar sum seems more reasonable.
  bge_rho_trk.set_particles(selAndExcldTrks);

  fastjet::contrib::ConstituentSubtractor subtractor_trk; 
  subtractor_trk.set_max_standardDeltaR(_csMaxR); 
  subtractor_trk.set_alpha(_alphaSubtr);
  subtractor_trk.set_ghost_area(_ghost_area);  // same 
  subtractor_trk.set_background_estimator(&bge_rho_trk);
  subtractor_trk.set_common_bge_for_rho_and_rhom(true); // for massless input particles it\
  
  if ( _saveLog)   cout << endl << subtractor_trk.description() << endl ; 
  
  vector<PseudoJet> corrected_selectedTrks = subtractor_trk.subtract_event(selectedTrks, _etaTrkCut);
  
  
  
  // Fill the histograms:
  for ( int ii =0 ; ii < selectedTrks.size() ; ii++ ) {
    hTrkPtEta_preCS_cent.at(cent_bin)->Fill( selectedTrks[ii].pt()*0.001, selectedTrks[ii].eta(), event_weight);
  }
  for ( int ii =0 ; ii < corrected_selectedTrks.size() ; ii++ ) {
    hTrkPtEta_postCS_cent.at(cent_bin)->Fill( corrected_selectedTrks[ii].pt()*0.001,  corrected_selectedTrks[ii].eta(), event_weight);
  }
  for ( int ii =0 ; ii < selGenMatchTrks.size() ; ii++ ) {
    hTrkPtEta_genMatch_cent.at(cent_bin)->Fill( selGenMatchTrks[ii].pt()*0.001,  selGenMatchTrks[ii].eta(), event_weight);
  }


  for ( int ii =0 ; ii < selectedTrks.size() ; ii++ ) {
    float ipt = selectedTrks[ii].pt()*0.001;
    float ieta = selectedTrks[ii].eta();
    float iphi = selectedTrks[ii].phi();

    float minR = 0.0001; 
    int jMatch = -1;
    for ( int jj =0 ; jj < corrected_selectedTrks.size() ; jj++ ) {
      float jpt = corrected_selectedTrks[jj].pt()*0.001;
      float jeta = corrected_selectedTrks[jj].eta();
      float jphi = corrected_selectedTrks[jj].phi();
      float drij = DeltaR( iphi, ieta, jphi, jeta) ;
      if ( drij < minR )   {
	minR = drij;
	jMatch = jj ; 
      }  
    }
    float jmpt = 0;    // if no tracks are found in distane of 0.01
    if ( jMatch > -1 ) { 
      jmpt = corrected_selectedTrks[jMatch].pt()*0.001;
    }
    h_bkgSubt_prePt_postPt_cent.at(cent_bin)->Fill( ipt, jmpt, event_weight);
    h_trkPt_trkBkgPt_cent.at(cent_bin)->Fill( ipt, ipt-jmpt, event_weight);
    h_dRSubt_trkPt_cent.at(cent_bin)->Fill( minR, ipt, event_weight);
    
    for ( int ijc = 0 ; ijc < etaJetCone150.size() ; ijc++)   {  // If it is in jet cone    
      float drTrkJet = DeltaR( phiJetCone150[ijc], etaJetCone150[ijc], iphi, ieta) ;
      if (  drTrkJet < _JetRadiusAna ) { 
	h_bkgSubt_prePt_postPt_jetCone_cent.at(cent_bin)->Fill( ipt, jmpt, event_weight);
	h_trkPt_trkBkgPt_jetCone_cent.at(cent_bin)->Fill( ipt, ipt-jmpt,  event_weight);
	h_dRSubt_trkPt_jetCone_cent.at(cent_bin)->Fill( minR, ipt, event_weight);
      }
    }
    
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
  

  ////
  

  /////////////   Main Loop:  Reco jets /////////////////////////////////////////
  xAOD::JetContainer::const_iterator jet_itr = reco_jets->begin();
  xAOD::JetContainer::const_iterator jet_end = reco_jets->end();
  int nRecoJetCounter=0;

  vector<xAOD::Jet> selectedRecoJets; 
  for( ; jet_itr != jet_end; ++jet_itr ) {
    //    cout << " jet_itr = << " << jet_itr << endl;
    xAOD::Jet theRecoJet;
    theRecoJet.makePrivateStore( **jet_itr );
    
    const xAOD::JetFourMom_t jet_4momUnCal = theRecoJet.jetP4("JetSubtractedScaleMomentum"); // uncalib
    xAOD::JetFourMom_t jet_4momCalib = theRecoJet.jetP4();   // CALIBRATED!!! 

    
    fastjet::PseudoJet unCaliFourVec = fastjet::PseudoJet ( jet_4momUnCal.px(), jet_4momUnCal.py(), jet_4momUnCal.pz(), jet_4momUnCal.energy() );
    if ( !_isMC) { 
      const xAOD::JetFourMom_t jet_4mom_unsubtracted = theRecoJet.jetP4("JetUnsubtractedScaleMomentum");
      theRecoJet.setJetP4("JetConstitScaleMomentum",jet_4mom_unsubtracted); //Required
      const xAOD::JetFourMom_t jet_4mom_EMJES = theRecoJet.jetP4(); //This is default, up to EM+JES calibrated four momentum
      theRecoJet.setJetP4("JetGSCScaleMomentum", jet_4mom_EMJES);
      EL_RETURN_CHECK("execute()", m_jetCalibration_insitu->applyCalibration( theRecoJet ) );
      jet_4momCalib = theRecoJet.jetP4(); // uncalib
    }

    fastjet::PseudoJet CalibFourVec = fastjet::PseudoJet ( jet_4momCalib.px(), jet_4momCalib.py(), jet_4momCalib.pz(), jet_4momCalib.energy() );
  
    double jet_pt  = jet_4momCalib.pt() * 0.001 ;
    double jet_eta = jet_4momCalib.eta();
    double jet_phi = PhiInPI ( jet_4momCalib.phi() );
    double jet_mass = CalibFourVec.m() * 0.001 ;
    double jet_rap  = CalibFourVec.rapidity();
    double jet_ptRaw = jet_4momUnCal.pt() * 0.001;
    double jet_etaRaw = jet_4momUnCal.eta() ;
    double jet_phiRaw = jet_4momUnCal.phi() ;
    double jet_massRaw = unCaliFourVec.m() * 0.001;
    
    
    if (_saveLog) { 
      cout << "In the analysis collection : " << endl;
      cout <<" pt = " << jet_pt << endl;
      cout <<" pt/e = " << jet_pt *1000/ jet_4momCalib.energy() << endl;
      cout <<" m/e = " << jet_mass *1000/ jet_4momCalib.energy() << endl;
      
      cout << " In Uncal jet : " << endl; 
      cout <<" pt   =  " << jet_ptRaw  << endl;
      cout <<" pt/e = "  << jet_4momUnCal.pt() / jet_4momUnCal.energy() << endl;
      cout <<" m/e  = "  << jet_4momUnCal.mass() / jet_4momUnCal.energy() << endl;
      
    }
    
    if (jet_pt < _pTjetCut)          continue;
    if (fabs(jet_eta) > _etaJetCut)  continue;
    if( (_isPP) && (!m_jetCleaningToolHandle->keep( **jet_itr )) )   {
      if ( _saveLog) cout << " this jet is cleaned " << endl;
      continue;
    }

    
    vector<double> jet_TrigPresc_vector ;
    vector<bool> jet_IsTrig_vector ; 

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

      t_trkTow[nRecoJetCounter] = (TH2F*)eventDisTrk->Clone(Form("trk_i%d",nRecoJetCounter));
      t_trkTow1[nRecoJetCounter] = (TH2F*)eventDisTrk->Clone(Form("trk1_i%d",nRecoJetCounter));
      t_trkTow2[nRecoJetCounter] = (TH2F*)eventDisTrk->Clone(Form("trk2_i%d",nRecoJetCounter));
      t_trkTow[nRecoJetCounter]->Reset();
      t_trkTow1[nRecoJetCounter]->Reset();
      t_trkTow2[nRecoJetCounter]->Reset();

    }
  
    if ( _saveLog) cout << "*~*~*~*~*~*~ RECO ~*~*~*~*~*~*" << endl << "  Anti-kT  jet [pt, eta, phi] : " << jet_pt <<", "<<jet_eta<<", "<<jet_phi<<endl << " Raw pT: " << jet_ptRaw << " GeV" << endl;
    
    int icoutB =1 ;
    const xAOD::JetConstituentVector recoConsts = (*jet_itr)->getConstituents();
    xAOD::JetConstituentVector::iterator itCnst   = recoConsts.begin();
    xAOD::JetConstituentVector::iterator itCnst_E = recoConsts.end();
    vector<fastjet::PseudoJet>  nonZeroConsts;
    vector<fastjet::PseudoJet>  toBeSubtracted; // reverse the pT only and subtract
    vector<bool>  nFlag; // reverse the pT only and subtract
    
    // Test for jet mass reproducability  //this test!!  cross-check

    /*      if ( _saveLog ) { 
      double totalPx =0;
      double totalPy =0;
      double totalPz =0;
      double totalE =0; 
      for( ; itCnst != itCnst_E; ++itCnst ) {
      totalE = totalE + (*itCnst)->pt() * cosh ((*itCnst)->eta() ) ;
      totalPx = totalPx +  (*itCnst)->pt() * cos ( (*itCnst)->phi() ) ;
      totalPy = totalPy + (*itCnst)->pt() * sin ( (*itCnst)->phi() ) ;
      totalPz = totalPz + (*itCnst)->pt() * sinh ((*itCnst)->eta() ) ;
      //	cout << "(*itCnst)->pt() = " << (*itCnst)->pt()  << endl;
      //	cout << "(*itCnst)->eta() = " << (*itCnst)->eta()  << endl;
      //	cout << "(*itCnst)->phi() = " << (*itCnst)->phi()  << endl << endl;
      
      }
      double totolM = sqrt ( totalE*totalE - totalPx*totalPx - totalPy*totalPy - totalPz*totalPz ) ;
      cout  << " analysis jet pt = " << jet_pt << endl;
      cout << " Constituent measured pt = " << sqrt( totalPx*totalPx+ totalPy*totalPy) << endl;
      cout  << " analysis jet m/pt = " << jet_mass/jet_pt << endl;
      cout  << " Constituent measured m/pt = " << totolM/ (sqrt( totalPx*totalPx+ totalPy*totalPy)) << endl; 
      cout << "ratio =" << totolM/ (sqrt( totalPx*totalPx+ totalPy*totalPy)) / ( jet_mass/jet_pt ) << endl;
    }
    itCnst   = recoConsts.begin();
    */

    // For unsubtracted towers
    if ( _useUnbtMass )  {
      double totalPx =0;
      double totalPy =0;
      double totalPz =0;
      double totalE =0; 
      for( ; itCnst != itCnst_E; ++itCnst ) {
	const xAOD::CaloCluster* cl=static_cast<const xAOD::CaloCluster*>(itCnst->rawConstituent());
	if ( _saveLog) { 
	  cout << "-> cl->rawE() = " <<  cl->rawE()  << endl;
	  cout << "-> cl->altE() = " <<  cl->altE()  << endl;
	  cout << "(*itCnst)->pt() * cosh ((*itCnst)->eta() )" << (*itCnst)->pt() * cosh ((*itCnst)->eta() ) << endl;
	}
	totalE = totalE + cl->altE()  ;
	totalPx = totalPx + cl->altE() /cosh ((*itCnst)->eta()) * cos ( (*itCnst)->phi() ) ;
	totalPy = totalPy + cl->altE() /cosh ((*itCnst)->eta()) * sin ( (*itCnst)->phi() ) ;
	totalPz = totalPz + cl->altE()/cosh ((*itCnst)->eta())  * sinh ((*itCnst)->eta() ) ;
      }
      double jet_massUnSubt = sqrt ( totalE*totalE - totalPx*totalPx - totalPy*totalPy - totalPz*totalPz ) * 0.001;
      jet_mass = jet_massUnSubt;  // highjacking
      //  cout  << " analysis jet mass = " << jet_mass << endl;
      //      cout  << " Constituent measured mass = " << jet_massUnSubt*0.001 << endl;
      //      cout << "ratio =" << jet_massUnSubt*0.001/jet_mass << endl;
    }
    

    double ghostE = 0.00001;
    itCnst   = recoConsts.begin(); // VERY IMPORANT! 
    for( ; itCnst != itCnst_E; ++itCnst ) {
      int icoutC =1 ;
      double theEta = (*itCnst)->Eta() ; 
      double thePhi = PhiInPI ( (*itCnst)->Phi() ) ;
      const fastjet::PseudoJet thisConst = fastjet::PseudoJet( (*itCnst)->Px(), (*itCnst)->Py(), (*itCnst)->Pz(), (*itCnst)->E() );
      
      if ( _towerBkgKill == -1 ) { 
	if ( (*itCnst)->pt() > 0 ) { // normal tower 
	  nonZeroConsts.push_back(thisConst);
	  toBeSubtracted.push_back ( fastjet::PseudoJet (0,0,0,ghostE));
      nFlag.push_back (false);
	  //	  if ( _saveLog) cout << "positive pt =" << (*itCnst)->pt()*0.001<< endl;
	  //	  if ( _saveLog) cout << "positive E =" << (*itCnst)->E()*0.001<< endl;
	  //	  if ( _saveLog) cout << "positive eta =" << (*itCnst)->eta()<< endl;
	  //	  if ( _saveLog) cout << "positive phi =" << (*itCnst)->phi()<< endl;
	}
	else {  // negative tower 
	  double gpx = ghostE/cosh(theEta) * cos(thePhi) ;
	  double gpy = ghostE/cosh(theEta) * sin(thePhi) ;
	  double gpz = ghostE*tanh(theEta) ;
	  
	  //  if ( _saveLog) cout << "negative pt =" << (*itCnst)->pt()*0.001<< endl;
	  //	  if ( _saveLog) cout << "negative E =" << (*itCnst)->E()*0.001<< endl;
	  //	  if ( _saveLog) cout << "negative eta =" << (*itCnst)->eta()<< endl;
	  //	  if ( _saveLog) cout << "negative phi =" << (*itCnst)->phi()<< endl;
	  double posPt = - (*itCnst)->pt() ;
	  double posPx = posPt * cos(thePhi);
          double posPy = posPt * sin(thePhi);
          double posPz = posPt * sinh(theEta);
          double posE =  posPt * cosh(theEta);
	  nonZeroConsts.push_back(fastjet::PseudoJet ( gpx,gpy,gpz,ghostE )  );
	  toBeSubtracted.push_back(fastjet::PseudoJet ( posPx, posPy, posPz, posE) );
	  nFlag.push_back(true);
	}
      }
      else if ( _towerBkgKill == 0 ) {  // Just ignore negative towers
	if ( (*itCnst)->pt() > 0 ) {
	  nonZeroConsts.push_back(thisConst);
	  toBeSubtracted.push_back ( fastjet::PseudoJet ( 0.000001, 0.000001, 0.000001, 0.000002 )) ; // place holder
          nFlag.push_back (false);
	}
      }
      else if ( _towerBkgKill == 1 ) { // soft kill
	cout << "No SoftKilling module deployed yet" << endl; 
      }
      
      // Fill the eventdisplay histogram first! 
      if (_saveEvtDisplay) 
	t_recTow[nRecoJetCounter]->Fill(theEta - jet_rap, DeltaPhi(thePhi, jet_phi), (*itCnst)->pt() *0.001 ) ;
      
    }
    double jet_nConst = nonZeroConsts.size();

    
 
    // Cambridge reclustering 
    if ( _saveLog) cout << "Reco re-clustering starts!" << endl;
    fastjet::ClusterSequence recoCamSeq(nonZeroConsts, jetDefReclus);
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

    float t_recoChPt      = -1 ;
    float t_recoChMass    = -100;
    float t_recoChPtRaw = -100;
    float t_recoChMassRaw = -100;
    float t_recoChPtRcSubt = -100;
    float t_recoChMassRcSubt = -100;
    float t_recoChMassRcSubt2 = -100;
    float t_recoChMassRcSubt4 = -100;
    float t_recoChMassRcSubtEV = -100;
    float t_recoChMassRcSubtMV = -100;


    float t_recoChMassGm  = -100;

    float t_recoNChRaw = -1 ;  
    float t_recoNChBkg = -1 ;  
    float t_recoNChBkgNoWgt = -1 ;  
    float t_recoDrChJetBkg = -1 ;  
    float t_recoMaxTrkPt = -1 ;  


    vector<fastjet::PseudoJet> corrRecCamJets = fastjet::sorted_by_pt(cambridgeJet); // return a vector of jets sorted into decreasing energy

    if ( corrRecCamJets.size() > 0 )   {
      fastjet::PseudoJet originalJet = corrRecCamJets[0];
      
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
    // tracks in jet cone after CS
    vector<fastjet::PseudoJet> trkConsts;
    for ( int ic=0; ic< corrected_selectedTrks.size() ; ic++) {
      if ( DeltaR ( jet_phi, jet_rap, corrected_selectedTrks[ic].phi(), corrected_selectedTrks[ic].rapidity() ) <  _JetRadiusAna ) {
	if ( corrected_selectedTrks[ic].pt()*0.001 < _ptCutPostCS )  
	  continue;
	trkConsts.push_back( corrected_selectedTrks[ic]) ;	
	if (_saveEvtDisplay)
	  t_trkTow[nRecoJetCounter]->Fill(corrected_selectedTrks[ic].rapidity() - jet_rap, DeltaPhi(corrected_selectedTrks[ic].phi(), jet_phi), corrected_selectedTrks[ic].pt()*0.001 );
      }
    }
    // tracks in jet cone before CS
    vector<fastjet::PseudoJet> trkConstsRaw;
    vector<fastjet::PseudoJet> trkConstsRaw2;
    vector<fastjet::PseudoJet> trkConstsRaw4;
    vector<fastjet::PseudoJet> trkConstsRawEV;
    vector<fastjet::PseudoJet> trkConstsRawMV;
    for ( int ic=0; ic< selectedTrks.size() ; ic++) {
      if ( DeltaR ( jet_phi, jet_rap, selectedTrks[ic].phi(), selectedTrks[ic].rapidity() ) <  _JetRadiusAna ) {
	trkConstsRaw.push_back( selectedTrks[ic]) ;	
	
	if ( selectedTrks[ic].pt() * 0.001 > 2)  trkConstsRaw2.push_back( selectedTrks[ic]) ;	
	if ( selectedTrks[ic].pt() * 0.001 > 4)  trkConstsRaw4.push_back( selectedTrks[ic]) ;	

	// Efficiency variation 
	double densePart = 0;
	if ( DeltaR ( jet_phi, jet_eta, selectedTrks[ic].phi(), selectedTrks[ic].eta()) < 0.1) 
	  densePart =0.004;
	double modifiedTrkPt = 19;
	if ( selectedTrks[ic].pt() * 0.001 < 19 )   
	  modifiedTrkPt = selectedTrks[ic].pt() * 0.001; 
	
	double effUnc2 = densePart*densePart;
	for (int i=0;i<4;i++){
	  int iptbin = h_eff_uncert_2015_material[i]->GetXaxis()->FindBin(modifiedTrkPt);
	  int ietabin = h_eff_uncert_2015_material[i]->GetYaxis()->FindBin(selectedTrks[ic].eta());
	  effUnc2 = effUnc2 + pow(h_eff_uncert_2015_material[i]->GetBinContent(iptbin,ietabin),2);
	}
	double effUnc1 = sqrt(effUnc2);
	
	// fake rate : 
	int theJetEtaBin = 0;  
	if ( fabs(jet_eta) < 0.3 )  theJetEtaBin = 0;
	else if ( fabs(jet_eta) < 0.8 )  theJetEtaBin = 1;
	else if ( fabs(jet_eta) < 1.2 )  theJetEtaBin = 2;
	else                             theJetEtaBin = 3;
	int binInFakeHist = h_fake_uncert[theJetEtaBin]->FindBin( selectedTrks[ic].pt() * 0.001,  jet_pt) ; // both in GeV unit
	double fakeUnc1 = h_fake_uncert[theJetEtaBin]->GetBinContent(binInFakeHist);
		
	double totUnc = fabs(effUnc1) + fabs(fakeUnc1);
	//	cout  <<  effUnc1 <<"+ "<<fakeUnc1<<" = "<< totUnc << endl;
	// tracking eff/fake rate variation 
	if ( hRandomUnit->GetRandom() > totUnc )     // random number between 0 - 1 
	  trkConstsRawEV.push_back( selectedTrks[ic]) ;	
	
	// Momentum variation  // 
	trkConstsRawMV.push_back( selectedTrksMV[ic]) ; // ys <= variation 
	
      }
    }
    // tracks matched to gen particle
    vector<fastjet::PseudoJet> trkConstsGm;
    for ( int ic=0; ic< selGenMatchTrks.size() ; ic++) {
      if ( DeltaR (jet_phi, jet_rap, selGenMatchTrks[ic].phi(), selGenMatchTrks[ic].rapidity() ) <  _JetRadiusAna ) {
	trkConstsGm.push_back( selGenMatchTrks[ic]) ;	
      }
    }
    
    // Raw charged mass 
    fastjet::PseudoJet trkSumRaw = fastjet::PseudoJet(0,0,0,0);
    fastjet::PseudoJet trkSumRaw2 = fastjet::PseudoJet(0,0,0,0);
    fastjet::PseudoJet trkSumRaw4 = fastjet::PseudoJet(0,0,0,0);
    fastjet::PseudoJet trkSumRawEV = fastjet::PseudoJet(0,0,0,0);
    fastjet::PseudoJet trkSumRawMV = fastjet::PseudoJet(0,0,0,0);

    for ( int ii=0; ii< trkConstsRaw.size() ; ii++)       {
      trkSumRaw = trkSumRaw + trkConstsRaw[ii] ;
      if ( trkConstsRaw[ii].pt() >  t_recoMaxTrkPt )
	t_recoMaxTrkPt = trkConstsRaw[ii].pt() ;
    }

    for ( int ii=0; ii< trkConstsRaw2.size() ; ii++)       {
      trkSumRaw2 = trkSumRaw2 + trkConstsRaw2[ii] ;
    }
    for ( int ii=0; ii< trkConstsRaw4.size() ; ii++)       {
      trkSumRaw4 = trkSumRaw4 + trkConstsRaw4[ii] ;
    }
    for ( int ii=0; ii< trkConstsRawEV.size() ; ii++)       {
      trkSumRawEV = trkSumRawEV + trkConstsRawEV[ii] ;
    }
    for ( int ii=0; ii< trkConstsRawMV.size() ; ii++)       {
      trkSumRawMV = trkSumRawMV + trkConstsRawMV[ii] ;
    }


    t_recoMaxTrkPt = t_recoMaxTrkPt * 0.001;
    float rawPhiTrkSum = trkSumRaw.phi();
    float rawEtaTrkSum = trkSumRaw.eta();
  
    // Gne Massched track mass
    fastjet::PseudoJet trkSumGm = fastjet::PseudoJet(0,0,0,0);
    for ( int ii=0; ii< trkConstsGm.size() ; ii++)       {
      trkSumGm = trkSumGm + trkConstsGm[ii] ;
    }
    t_recoChMassGm =  trkSumGm.m() *0.001 ;

    // tracks in random cone :
    if ( _saveLog) 	cout << "_vTrkPtForRC.size() = " << _vTrkPtForRC.size() << endl; 
    if ( _saveLog) 	cout << "selectedTrks.size() = " << selectedTrks.size() << endl; 
    if ( _vTrkPtForRC.size() != selectedTrks.size() )   
      cout << " Error!!  _vTrkPtForRC.size() and selectedTrks.size() are not same!!! " << endl;


    int   nTrkBkgCounts = 0;
    float sumPxBkg = 0;
    float sumPyBkg = 0;
    float sumPzBkg = 0;
    float sumEBkg = 0;
    float nTrkBkg = 0;

    float sumPxBkg2 = 0;
    float sumPyBkg2 = 0;
    float sumPzBkg2 = 0;
    float sumEBkg2 = 0;

    float sumPxBkg4 = 0;
    float sumPyBkg4 = 0;
    float sumPzBkg4 = 0;
    float sumEBkg4 = 0;

    float sumPxBkgEV = 0;
    float sumPyBkgEV = 0;
    float sumPzBkgEV = 0;
    float sumEBkgEV = 0;

    float sumPxBkgMV = 0;
    float sumPyBkgMV = 0;
    float sumPzBkgMV = 0;
    float sumEBkgMV = 0;

    for ( int ic=0 ; ic< selectedTrks.size() ; ic++) { 
      float iTrkPt = selectedTrks[ic].pt() * 0.001;
      float iTrkEta = selectedTrks[ic].eta();
      float iTrkPhi = selectedTrks[ic].phi();

      float iTrkPtMV = selectedTrksMV[ic].pt() * 0.001;
      float iTrkEtaMV = selectedTrksMV[ic].eta();
      float iTrkPhiMV = selectedTrksMV[ic].phi();
      if ( fabs(selectedTrksMV[ic].eta() - selectedTrks[ic].eta() ) > 0.0001 )  {
	cout << "!!!!!! iTrkPt =! iTrPtMV " << endl; 
	cout << " iTrkPt, iTrkPtMV = " << iTrkPt <<", " << iTrkPtMV << endl;
	cout << " iTrkEta, iTrkEtaMV = " << iTrkEta <<", " << iTrkEtaMV << endl;
      }

      uee->FindCone(iTrkPt, iTrkEta, iTrkPhi);

      float deltaRBkgr = uee->GetDeltaRToConeAxis();
      float deltaEtaBkgr = uee->GetDeltaEtaToConeAxis();
      float deltaPhiBkgr = uee->GetDeltaPhiToConeAxis();
      if (deltaRBkgr <= 0.4) { 
	//around calo jet 
	float w_ncones = uee->GetNConesWeight();

	float w_eta  = uee->CalculateEtaWeight(iTrkPt, iTrkEta, jet_eta, cent_bin_scheme30);
	float w_flow = uee->CalculateFlowWeight( iTrkPt, iTrkEta, iTrkPhi, jet_phi,  FCalEt ); //Flow correction
	float w_bkgr = w_eta * w_ncones * w_flow;
	
	// around track jet 
	float w_eta2  = uee->CalculateEtaWeight(iTrkPt, iTrkEta, rawEtaTrkSum, cent_bin_scheme30);
	float w_flow2 = uee->CalculateFlowWeight( iTrkPt, iTrkEta, iTrkPhi, rawPhiTrkSum,  FCalEt ); //Flow correction
        float w_bkgr2 = w_eta2 * w_ncones * w_flow2;
	
	if ( _saveLog)      cout << "deta, dephi in random cone  = " << deltaEtaBkgr <<",  "<<deltaPhiBkgr << endl;
	float newEta = jet_eta + deltaEtaBkgr;
	float newPhi = jet_phi + deltaPhiBkgr;
	float newTrkPx = iTrkPt * cos(newPhi);
	float newTrkPy = iTrkPt * sin(newPhi);
	float newTrkPz = iTrkPt * sinh(newEta);
	float newTrkP  = iTrkPt * cosh(newEta);
	float newTrkE  = sqrt( newTrkP*newTrkP + _chParticleMassMeV*_chParticleMassMeV*0.000001);

	float newTrkPxMV = iTrkPtMV * cos(newPhi);
	float newTrkPyMV = iTrkPtMV * sin(newPhi);
	float newTrkPzMV = iTrkPtMV * sinh(newEta);
	float newTrkPMV  = iTrkPtMV * cosh(newEta);
	float newTrkEMV  = sqrt( newTrkPMV*newTrkPMV + _chParticleMassMeV*_chParticleMassMeV*0.000001);



	float newEta2 = rawEtaTrkSum + deltaEtaBkgr;
	float newPhi2 = rawPhiTrkSum + deltaPhiBkgr;
	float newTrkPx2 = iTrkPt * cos(newPhi2);
	float newTrkPy2 = iTrkPt * sin(newPhi2);
	float newTrkPz2 = iTrkPt * sinh(newEta2);
	float newTrkP2  = iTrkPt * cosh(newEta2);
	float newTrkE2  = sqrt( newTrkP2*newTrkP2 + _chParticleMassMeV*_chParticleMassMeV*0.000001);
	
	sumPxBkg = sumPxBkg + newTrkPx*w_bkgr ;  
	sumPyBkg = sumPyBkg + newTrkPy*w_bkgr ;  
	sumPzBkg = sumPzBkg + newTrkPz*w_bkgr ;  
	sumEBkg = sumEBkg + newTrkE*w_bkgr ;  

	/// momentum variation 
	sumPxBkgMV = sumPxBkgMV + newTrkPxMV*w_bkgr ;  
	sumPyBkgMV = sumPyBkgMV + newTrkPyMV*w_bkgr ;  
	sumPzBkgMV = sumPzBkgMV + newTrkPzMV*w_bkgr ;  
	sumEBkgMV = sumEBkgMV + newTrkEMV*w_bkgr ;  
      


	if ( iTrkPt > 2) { 
	  sumPxBkg2 = sumPxBkg2 + newTrkPx*w_bkgr ;  
	  sumPyBkg2 = sumPyBkg2 + newTrkPy*w_bkgr ;  
	  sumPzBkg2 = sumPzBkg2 + newTrkPz*w_bkgr ;  
	  sumEBkg2 = sumEBkg2 + newTrkE*w_bkgr ;  
	}
	if ( iTrkPt > 4) { 
	  sumPxBkg4 = sumPxBkg4 + newTrkPx*w_bkgr ;  
	  sumPyBkg4 = sumPyBkg4 + newTrkPy*w_bkgr ;  
	  sumPzBkg4 = sumPzBkg4 + newTrkPz*w_bkgr ;  
	  sumEBkg4 = sumEBkg4 + newTrkE*w_bkgr ;  
	}

	// Efficiency variation
	double densePart = 0;
	if ( deltaRBkgr < 0.1) 
          densePart =0.004;
	double modifiedTrkPt = 19;
        if ( iTrkPt < 19 )
          modifiedTrkPt = iTrkPt;

        double effUnc2 = densePart*densePart;
        for (int i=0;i<4;i++){
          int iptbin = h_eff_uncert_2015_material[i]->GetXaxis()->FindBin(modifiedTrkPt);
          int ietabin = h_eff_uncert_2015_material[i]->GetYaxis()->FindBin(iTrkEta);
          effUnc2 = effUnc2 + pow(h_eff_uncert_2015_material[i]->GetBinContent(iptbin,ietabin),2);
        }
	double effUnc1 = sqrt(effUnc2);
	// fake rate :
        int theJetEtaBin = 0;
        if ( fabs(jet_eta) < 0.3 )  theJetEtaBin = 0;
        else if ( fabs(jet_eta) < 0.8 )  theJetEtaBin = 1;
        else if ( fabs(jet_eta) < 1.2 )  theJetEtaBin = 2;
        else                             theJetEtaBin = 3;
        int binInFakeHist = h_fake_uncert[theJetEtaBin]->FindBin( iTrkPt, jet_pt); // in GeV
	double fakeUnc1   = h_fake_uncert[theJetEtaBin]->GetBinContent(binInFakeHist);
	// total rate :
	double totUnc = fabs(effUnc1) + fabs(fakeUnc1);
	//	cout  <<  effUnc1 <<"+ "<<fakeUnc1<<" = "<< totUnc << endl;

	/// eff/fake rate variation 
        if ( hRandomUnit->GetRandom() > totUnc )   {   // random number between 0 - 1
	  sumPxBkgEV = sumPxBkgEV + newTrkPx*w_bkgr ;  
	  sumPyBkgEV = sumPyBkgEV + newTrkPy*w_bkgr ;  
	  sumPzBkgEV = sumPzBkgEV + newTrkPz*w_bkgr ;  
	  sumEBkgEV = sumEBkgEV + newTrkE*w_bkgr ;  
	}



	//	else {  cout << "rejected" << endl; }
	
	nTrkBkg = nTrkBkg +  w_bkgr;
	nTrkBkgCounts = nTrkBkgCounts + 1;
      }
    }
    if ( _saveLog)
      cout << "Sum of charged particle background px,py,pz,E = " << sumPxBkg<<", " << sumPyBkg<<", " << sumPzBkg << ",   " << sumEBkg << endl;
    
    
    // Calculate charged pT sum and charged mass 
    fastjet::PseudoJet trkSum = fastjet::PseudoJet(0,0,0,0);
    for ( int ii=0; ii< trkConsts.size() ; ii++) {
      trkSum = trkSum + trkConsts[ii] ;   //  CS'ed track sum
    }
    fastjet::PseudoJet trkSumRcBkg = fastjet::PseudoJet(  sumPxBkg*1000., sumPyBkg*1000., sumPzBkg*1000., sumEBkg*1000.); // in MeV
    fastjet::PseudoJet trkSumRcSubt = trkSumRaw - trkSumRcBkg; 

    fastjet::PseudoJet trkSumRcBkg2 = fastjet::PseudoJet(  sumPxBkg2*1000., sumPyBkg2*1000., sumPzBkg2*1000., sumEBkg2*1000.); // in MeV
    fastjet::PseudoJet trkSumRcSubt2 = trkSumRaw2 - trkSumRcBkg2; 

    fastjet::PseudoJet trkSumRcBkg4 = fastjet::PseudoJet(  sumPxBkg4*1000., sumPyBkg4*1000., sumPzBkg4*1000., sumEBkg4*1000.); // in MeV
    fastjet::PseudoJet trkSumRcSubt4 = trkSumRaw4 - trkSumRcBkg4;

    fastjet::PseudoJet trkSumRcBkgEV = fastjet::PseudoJet(  sumPxBkgEV*1000., sumPyBkgEV*1000., sumPzBkgEV*1000., sumEBkgEV*1000.); // in MeV
    fastjet::PseudoJet trkSumRcSubtEV = trkSumRawEV - trkSumRcBkgEV;

    fastjet::PseudoJet trkSumRcBkgMV = fastjet::PseudoJet(  sumPxBkgMV*1000., sumPyBkgMV*1000., sumPzBkgMV*1000., sumEBkgMV*1000.); // in MeV
    fastjet::PseudoJet trkSumRcSubtMV = trkSumRawMV - trkSumRcBkgMV;
 

    t_recoChPtRaw =  trkSumRaw.pt() *0.001 ;   // raw mass
    t_recoChMassRaw =  trkSumRaw.m() *0.001 ;   // raw mass
    t_recoChPt =  trkSum.pt() *0.001 ;  // CS subtracted pt sum
    t_recoChMass =  trkSum.m() *0.001 ; // CS subtracted mass

    t_recoChPtRcSubt =  trkSumRcSubt.pt() *0.001 ; 
    t_recoChMassRcSubt =  trkSumRcSubt.m() *0.001 ; 
    
    t_recoChMassRcSubt2  =  trkSumRcSubt2.m() *0.001 ; 
    t_recoChMassRcSubt4  =  trkSumRcSubt4.m() *0.001 ; 
    t_recoChMassRcSubtEV =  trkSumRcSubtEV.m() *0.001 ; 
    t_recoChMassRcSubtMV =  trkSumRcSubtMV.m() *0.001 ; 
    

    t_recoNChRaw = trkConstsRaw.size() ; 
    t_recoNChBkg = nTrkBkg ;
    t_recoNChBkgNoWgt = nTrkBkgCounts ;

    t_recoDrChJetBkg =  DeltaR(trkSumRaw.phi(), trkSumRaw.eta(), trkSumRcBkg.phi(), trkSumRcBkg.eta()) ; 
    
    if (_saveLog) { 
      cout << "cent = " << cent_bin << endl;
      
      cout << "N^ch in Raw = " << trkConstsRaw.size() << endl;
      cout << "N^ch(weighted) in Bkg = " << nTrkBkg << endl;
      cout << "N^ch(number  ) in Bkg = " << nTrkBkgCounts << endl;
      cout << "trkSumRcSubt pT / jet pT = " << trkSumRcSubt.pt() / jet_pt /1000 << endl;
      cout << "Gen Match M   = " << t_recoChMassGm << endl;
      cout << "ChMassRcSubt1 = " << t_recoChMassRcSubt << endl;
      cout <<" t_recoMaxTrkPt = " << t_recoMaxTrkPt << endl;
      cout << " trkSumRaw pt,pz,E: " << trkSumRaw.pt()*0.001 <<",  " <<trkSumRaw.pz()*0.001 <<",  " << trkSumRaw.e()*0.001 << endl;
      cout << " DR jet axis and trk jet axis: " <<  DeltaR(jet_phi, jet_eta, trkSumRaw.phi(), trkSumRaw.eta()) << endl;
      cout << " DR jet axis and bkg1       : " <<  DeltaR(jet_phi, jet_eta, trkSumRcBkg.phi(), trkSumRcBkg.eta()) << endl;
      cout << " DR track jet axis and bkg1      : " <<  DeltaR(trkSumRaw.phi(), trkSumRaw.eta(), trkSumRcBkg.phi(), trkSumRcBkg.eta()) << endl;
      cout << " trkSumRcBkg pt,pz,E: " << trkSumRcBkg.pt()*0.001 <<",  " <<trkSumRcBkg.pz()*0.001 <<", " << trkSumRcBkg.e()*0.001 << endl;
      cout << " trkSumRcSubt pt,pz,E: " << trkSumRcSubt.pt()*0.001 <<",  " <<trkSumRcSubt.pz()*0.001 <<", " << trkSumRcSubt.e()*0.001 << endl;
      
    }
    



    // Tracking efficiency calculation 
    if ( _isMC) {  
      for ( int ic=0; ic< selGenMatchTrks.size() ; ic++) {
	if ( DeltaR ( jet_phi, jet_rap, selGenMatchTrks[ic].phi(), selGenMatchTrks[ic].rapidity() ) <  _JetRadiusAna ) {
	  h_trkGen_pt_dphi_cent.at(cent_bin)->Fill( jet_pt, selGenMatchTrks[ic].pt()*0.001, DeltaPhi(selGenMatchTrks[ic].phi(), jet_phi) ) ;
	  h_trkGen_pt_drap_cent.at(cent_bin)->Fill( jet_pt, selGenMatchTrks[ic].pt()*0.001, selGenMatchTrks[ic].rapidity() - jet_rap );
	}
      }
      for ( int ic=0; ic< truthChargesAna.size() ; ic++) {
	if ( DeltaR ( jet_phi, jet_rap, truthChargesAna[ic].phi(), truthChargesAna[ic].rapidity() ) <  _JetRadiusAna ) {
	  h_allGen_pt_dphi_cent.at(cent_bin)->Fill( jet_pt, truthChargesAna[ic].pt()*0.001, DeltaPhi( truthChargesAna[ic].phi(), jet_phi) ) ;
	  h_allGen_pt_drap_cent.at(cent_bin)->Fill( jet_pt, truthChargesAna[ic].pt()*0.001, truthChargesAna[ic].rapidity() - jet_rap );
	}
      }
    }

    fastjet::ClusterSequence reChCam(trkConsts, jetDefReclus);
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
	//Fill the histogram
	if (_saveEvtDisplay) {
	  vector<fastjet::PseudoJet> sub1 =  parent1.constituents() ;
	  vector<fastjet::PseudoJet> sub2 =  parent2.constituents() ;
	  for ( int ic = 0 ; ic< sub1.size() ; ic++)
	    t_trkTow1[nRecoJetCounter]->Fill( sub1[ic].eta() - jet_rap, DeltaPhi( sub1[ic].phi(), jet_phi), sub1[ic].pt() * 0.001 );
	  for ( int ic = 0 ; ic< sub2.size() ; ic++)
	    t_trkTow2[nRecoJetCounter]->Fill( sub2[ic].eta() - jet_rap, DeltaPhi( sub2[ic].phi(), jet_phi), sub2[ic].pt() * 0.001 );
	  sub1.clear();
	  sub2.clear();
	}    
	if ( _saveLog && (cent_bin==0)) {   
	  cout << " centBin = " << cent_bin << endl;
	  cout << " dR        = " << t_recoChSdTheta << endl;
	  cout << " Track nConst1 = " << parent1.constituents().size() << endl;
	  cout << " Track nConst2 = " << parent2.constituents().size() << endl;
	}
	
      }
      else {
	  hRecoSdChStat->Fill(jet_pt, 0);
	}
      
    }
  
    selectedRecoJets.push_back(theRecoJet);
    vmass_reco.push_back(jet_mass);
    vpt_reco.push_back(jet_pt);
    vptRaw_reco.push_back(jet_ptRaw);
    vetaRaw_reco.push_back(jet_etaRaw);
    vphiRaw_reco.push_back(jet_phiRaw);
    vNconst_reco.push_back(jet_nConst);
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
  
    vChSdPt_reco.push_back(t_recoChSdPt);
    vChSdMass_reco.push_back(t_recoChSdMass)  ;
    vChSdZ_reco.push_back(t_recoChSdZ)  ;
    vChSdTheta_reco.push_back(t_recoChSdTheta) ;
    
    vChPt_reco.push_back(t_recoChPt);
    vChMass_reco.push_back(t_recoChMass)  ;
    vChPtRaw_reco.push_back(t_recoChPtRaw)  ;
    vChPtRaw_reco.push_back(t_recoChPtRaw)  ;
    vChMassRaw_reco.push_back(t_recoChMassRaw)  ;
    vChPtRcSubt_reco.push_back(t_recoChPtRcSubt)  ;
    vChMassRcSubt_reco.push_back(t_recoChMassRcSubt)  ;

    vChMassRcSubt2_reco.push_back(t_recoChMassRcSubt2)  ;
    vChMassRcSubt4_reco.push_back(t_recoChMassRcSubt4)  ;
    vChMassRcSubtEV_reco.push_back(t_recoChMassRcSubtEV)  ;
    vChMassRcSubtMV_reco.push_back(t_recoChMassRcSubtMV)  ;

    vChMassGm_reco.push_back(t_recoChMassGm)  ;

    vNChRaw.push_back(t_recoNChRaw);
    vNChBkg.push_back(t_recoNChBkg);
    vNChBkgNoWgt.push_back(t_recoNChBkgNoWgt);
    vDrChJetBkg.push_back(t_recoDrChJetBkg);
    vMaxTrkPt.push_back(t_recoMaxTrkPt);

    
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


  
  vector<xAOD::Jet> selectedGenJets;

  if (_isMC){
    if (_saveLog) {
      cout << endl << " /////// MC information " << endl;
    }
    
    const xAOD::JetContainer* truth_jets = 0;
    ANA_CHECK(event->retrieve(truth_jets, _truth_jet_collection.c_str() ));
    xAOD::JetContainer::const_iterator truth_jet_itr = truth_jets->begin();
    xAOD::JetContainer::const_iterator truth_jet_end = truth_jets->end();
    
    //  On the fly gen jet
    fastjet::ClusterSequence cs(truthParticles, jetDefAk);
    vector<fastjet::PseudoJet> OtfGenJets = fastjet::sorted_by_pt(cs.inclusive_jets());
    
    
    int nGenJetCounter =0;
    for( ; truth_jet_itr != truth_jet_end; ++truth_jet_itr ) {
      xAOD::JetFourMom_t jet_truth_4mom = (*truth_jet_itr)->jetP4();
      fastjet::PseudoJet jet_truth_pj = fastjet::PseudoJet ( jet_truth_4mom.px(), jet_truth_4mom.py(), jet_truth_4mom.pz(), jet_truth_4mom.energy() );
      xAOD::Jet theGenJet;
      theGenJet.makePrivateStore( **truth_jet_itr );
      
      double jet_mass = jet_truth_pj.m()*0.001;
      if (jet_mass < 0 ){ 
	cout << " jet_mass =  " << jet_mass << endl;
	cout << " jet_mass0 = " << jet_truth_4mom.M() * 0.001 << endl;
      }
      double jet_pt = jet_truth_pj.pt()*0.001;
      double jet_eta = jet_truth_pj.eta();
      double jet_rap = jet_truth_pj.rapidity();
      double jet_phi = PhiInPI ( jet_truth_pj.phi() ) ;
      if (jet_pt < _truthpTjetCut) continue;
      if ( fabs(jet_eta) > _etaJetCut+0.2 ) continue;
      
      // find the matched Oft genjet 
      double jet_pt2 = -10;
      fastjet::PseudoJet matchedOtf = fastjet::PseudoJet(0,0,0,0);;
      for ( int ij =0  ; ij< OtfGenJets.size() ; ij++) { 
	if ( OtfGenJets[ij].pt() * 0.001 < _truthpTjetCut ) continue;
	if ( fabs(OtfGenJets[ij].eta()) > _etaJetCut+0.2 ) continue;
	
	if ( DeltaR(   OtfGenJets[ij].phi(),   OtfGenJets[ij].eta(),  jet_phi, jet_eta ) < 0.1 )  {
	  matchedOtf =  OtfGenJets[ij]; 
	}
      }
      if (_saveLog)  cout << " jet_pt2 = " << jet_pt2 << endl;
      if (_saveLog)  cout << " jet_pt:jet_eta = " << jet_pt <<", " << jet_eta << endl;
      if (_saveLog)  cout << " otf_pt:otf_eta = " <<  matchedOtf.pt() * 0.001 <<", " <<  matchedOtf.eta() << endl;

      jet_pt2 = matchedOtf.pt() * 0.001 ;

      /* sanity checks 
      if ( matchedOtf != fastjet::PseudoJet(0,0,0,0) ) { 
	vector<fastjet::PseudoJet> sub1 =  matchedOtf.constituents() ;
	if (_saveLog)	cout << " sub pt : " ;
	for ( int ij =0  ; ij< sub1.size() ; ij++) {
	cout << sub1[ij].pt()*0.001<<", " ;
	  }
	cout << endl;
	}
	else if (_saveLog)	
      cout << " no otf matching" << endl;

      */ 

      // IMPORTNAT!  There must be no more continue command in this loop! 
      if (_saveEvtDisplay) {	
	eventDisGen->Reset();
	eventDisGen1->Reset();
	eventDisGen2->Reset();
	t_genTow[nGenJetCounter] = (TH2F*)eventDisGen->Clone(Form("genp_i%d",nGenJetCounter));
	t_genTow1[nGenJetCounter] = (TH2F*)eventDisGen->Clone(Form("genp1_i%d",nGenJetCounter));
	t_genTow2[nGenJetCounter] = (TH2F*)eventDisGen->Clone(Form("genp2_i%d",nGenJetCounter));
	t_genTow[nGenJetCounter]->Reset();
	t_genTow1[nGenJetCounter]->Reset();
	t_genTow2[nGenJetCounter]->Reset();
	
	eventDisChg->Reset();
	eventDisChg1->Reset();
	eventDisChg2->Reset();
	t_chgTow[nGenJetCounter] = (TH2F*)eventDisChg->Clone(Form("chg_i%d",nGenJetCounter));
	t_chgTow1[nGenJetCounter] = (TH2F*)eventDisChg->Clone(Form("chg1_i%d",nGenJetCounter));    
	t_chgTow2[nGenJetCounter] = (TH2F*)eventDisChg->Clone(Form("chg2_i%d",nGenJetCounter));
	t_chgTow[nGenJetCounter]->Reset();
	t_chgTow1[nGenJetCounter]->Reset();
	t_chgTow2[nGenJetCounter]->Reset();
      }
      vector<fastjet::PseudoJet > akConsts;
      //	vector<ElementLink < xAOD::IParticleContainer > > truthLinkVector =  (*truth_jet_itr)->constituentLinks ();
      //	vector<fastjet::PseudoJet > akConsts = (*truth_jet_itr)->constituents();
      if ( matchedOtf != fastjet::PseudoJet(0,0,0,0) )
	akConsts =  matchedOtf.constituents() ;
      // cout << "number of constituents = " <<  truthLinkVector.size() << endl;
      //	for (int i = 0; i < truthLinkVector.size(); i++) {
      //	  fastjet::PseudoJet akElement = fastjet::PseudoJet ( (*(truthLinkVector[i]))->p4() );
      //	  akConsts.push_back(akElement);
      //	}
      //	cout << "after truth_jet_itr->constituents()" << endl;
      
      // Truth recluster by C/A
      if (_saveEvtDisplay) {
	for ( int ic = 0 ; ic< akConsts.size() ; ic++)  {
	  double theRap  = akConsts[ic].rapidity();
	  double thePhi  = PhiInPI( akConsts[ic].phi() );
	  double thePt   = akConsts[ic].pt() * 0.001;
	  t_genTow[nGenJetCounter]->Fill(theRap - jet_rap, DeltaPhi(thePhi, jet_phi), thePt ) ;
	}
      }
      fastjet::ClusterSequence reCam(akConsts, jetDefReclus);
      vector<fastjet::PseudoJet> camJets = fastjet::sorted_by_pt(reCam.inclusive_jets()); // return a vector of jets sorted into decreasing energy
      hGenNcam->Fill( jet_pt, camJets.size() ) ;
      
      
      
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
      
      if ( camJets.size() > 0 )   {  // If cambridge jet is made .. 
	fastjet::PseudoJet originalJet = camJets[0];
	
	thePtrc = originalJet.pt() * 0.001;
	fastjet::PseudoJet sd_jet = softdropper(originalJet);
	thesdm =  sd_jet.m() * 0.001;
	thesdpt = sd_jet.pt() * 0.001;
	
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
    
      selectedGenJets.push_back(theGenJet);
      vMass_gen.push_back(jet_mass);
      vpt_gen.push_back(jet_pt);
      vpt2_gen.push_back(jet_pt2);
      veta_gen.push_back(jet_eta);
      vrap_gen.push_back(jet_rap);
      vphi_gen.push_back(jet_phi);
      vptRc_gen.push_back(thePtrc);
      vSdmass_gen.push_back(thesdm);
      vSdPt_gen.push_back(thesdpt);
      vSdz_gen.push_back(thesdz);
      vSdTheta_gen.push_back(thesdtheta);
      vSdpt1_gen.push_back(thesdpt1);
      vSdrap1_gen.push_back(thesdrap1);
      vSdphi1_gen.push_back(thesdphi1);
      vSdpt2_gen.push_back(thesdpt2);
      vSdrap2_gen.push_back(thesdrap2);
      vSdphi2_gen.push_back(thesdphi2);
      
	
      camJets.clear();
      
      // Charged Particle reclustering
      if (_saveLog)	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
      if (_saveLog)	cout << "(Truth) Charged particle clustering starts" << endl;
      
      vector<fastjet::PseudoJet> chargeConsts;
      float t_genNch = 0 ;
      for ( int ic=0; ic< truthChargesAna.size() ; ic++) {
	if ( DeltaR ( jet_phi, jet_rap, truthChargesAna[ic].phi(), truthChargesAna[ic].rapidity() ) <  _JetRadiusAna ) {
	  chargeConsts.push_back( truthChargesAna[ic]) ;
	    t_genNch = t_genNch + 1;
	}
      }
      
      if (_saveEvtDisplay) {
	for ( int ic = 0 ; ic< chargeConsts.size() ; ic++)  {
	  double theRap  = chargeConsts[ic].rapidity();
	  double thePhi  = PhiInPI( chargeConsts[ic].phi() );
	  double thePt   = chargeConsts[ic].pt() * 0.001;
	  t_chgTow[nGenJetCounter]->Fill(theRap - jet_rap, DeltaPhi(thePhi, jet_phi), thePt ) ;
	}
      }
      
      
      fastjet::ClusterSequence reChCam(chargeConsts, jetDefReclus);
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
	  
	  if (_saveEvtDisplay) {
	    vector<fastjet::PseudoJet> sub1 =  parent1.constituents() ;
	    vector<fastjet::PseudoJet> sub2 =  parent2.constituents() ;
	    for ( int ic = 0 ; ic< sub1.size() ; ic++)
	      t_chgTow1[nGenJetCounter]->Fill( sub1[ic].rapidity() - jet_rap, DeltaPhi( sub1[ic].phi(), jet_phi), sub1[ic].pt() * 0.001 );
	    for ( int ic = 0 ; ic< sub2.size() ; ic++)
                t_chgTow2[nGenJetCounter]->Fill( sub2[ic].rapidity() - jet_rap, DeltaPhi( sub2[ic].phi(), jet_phi), sub2[ic].pt() * 0.001 );
	    sub1.clear();
	    sub2.clear();
	  }
	  
	  
	  
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
      vNch_gen.push_back(t_genNch);
	
      nGenJetCounter++; // THIS MUST BE AT THE END OF THE jets LOOP! 
    }
    truthParticles.clear();
    truthChargesAna.clear();
    
  }
  
  
  // Initialize uncertainty tool
  if (_doJES) { 
    jesProv = new JetUncertaintiesTool("JESProvider");
    jesProv->setProperty("JetDefinition","AntiKt4EMTopo");
    jesProv->setProperty("MCType","MC15");
    jesProv->setProperty("ConfigFile","JES_2015/ICHEP2016/JES2015_19NP.config");
    jesProv->initialize();
  }

  // back to reco loop 
  for ( int ri = 0 ; ri< vpt_reco.size() ; ri++) { 
    
    int matchId = -1;
    double drMin = 0.2;    // dR cut
    if ( _isMC) { 
      for (int gi = 0; gi<vpt_gen.size() ; gi++) { 
	double dr_itr = DeltaR(vphi_reco[ri], veta_reco[ri], vphi_gen[gi], veta_gen[gi]);
	if ( dr_itr < drMin ) { 
	  matchId = gi;  
	  drMin = dr_itr;
	}
      }
    }
    
    if ( matchId != -1 ) { 
      
      if ( _doJES ) { 
	if ( _saveLog)  cout << "======= systematics ==== " << endl;
	
	xAOD::Jet* recoJetSysPP = new xAOD::Jet( selectedRecoJets[ri] );
	double iniPt = recoJetSysPP->pt();
	
	for ( int si = 0 ; si<=21; si++) {  // 20 - 41
	  double theUncertainty = intSignificance[si]  * (jesProv->getUncertainty( intrinsicComponent[si], (*recoJetSysPP)));
	  ptSysPP[si] = iniPt * (1. + theUncertainty) * 0.001;
	  
	}

	for ( int ii=0 ; ii<vUncProvHI.size() ;ii++) { 
	  //	  xAOD::Jet* thisJet = new xAOD::Jet();
	  xAOD::Jet* recoJetSys = new xAOD::Jet( selectedRecoJets[ri] );
	  xAOD::Jet* genJetSys = new xAOD::Jet( selectedGenJets[matchId] );
	  
	  if  (  (_isPP)  && ( (vUncertIndex.at(ii) == 16)||(vUncertIndex.at(ii) == 17) ) ) {
	    //	    This applies only for HI
	    ptSysHI[ii] = recoJetSys->pt() * 0.001 ;
	  }
	  else {
	    vUncProvHI.at(ii)->CorrectJet ( recoJetSys, genJetSys, cent_bin, FCalEt ) ;
	    ptSysHI[ii] = recoJetSys->pt() * 0.001;
	  }
	  delete recoJetSys; 
	  delete genJetSys; 
	}
	
	
	delete recoJetSysPP;
      }
    }
    
    // Fill the track mass and pt  //ys
    fastjet::PseudoJet theTrkSum4 = fastjet::PseudoJet (0,0,0,0);
    fastjet::PseudoJet theTrkSum6 = fastjet::PseudoJet (0,0,0,0);
    fastjet::PseudoJet theTrkSum8 = fastjet::PseudoJet (0,0,0,0);
    fastjet::PseudoJet theTrkSum10 = fastjet::PseudoJet (0,0,0,0);
    for ( int ii =0 ; ii < selectedTrks.size() ; ii++ ) {
      double ipt = selectedTrks[ii].pt()*0.001;
      double ieta = selectedTrks[ii].eta();
      double iphi = selectedTrks[ii].phi();
      double theDr = DeltaR (  vphi_reco[ri],  veta_reco[ri], iphi, ieta);
      if ( theDr < _JetRadiusAna ) { 
	if ( ipt > 4 ) 
	  theTrkSum4  = theTrkSum4 +  selectedTrks[ii];
	if ( ipt > 6 ) 
	  theTrkSum6  = theTrkSum6 +  selectedTrks[ii];
	if ( ipt > 8 ) 
	  theTrkSum8  = theTrkSum8 +  selectedTrks[ii];
	if ( ipt > 10 ) 
	  theTrkSum10  = theTrkSum10 +  selectedTrks[ii];
      }
    }
    trkJetMass4 = theTrkSum4.m() * 0.001;      trkJetPt4 = theTrkSum4.pt() * 0.001;
    trkJetMass6 = theTrkSum6.m() * 0.001;      trkJetPt6 = theTrkSum6.pt() * 0.001;
    trkJetMass8 = theTrkSum8.m() * 0.001;      trkJetPt8 = theTrkSum8.pt() * 0.001;
    trkJetMass10 = theTrkSum10.m() * 0.001;      trkJetPt10 = theTrkSum10.pt() * 0.001;
    
    
    resetSubstr(myJetSub);
    resetEvt(myEvt);
    if (_saveEvtDisplay) { 
      eventDisRecTow->Reset();
      eventDisRecTow1->Reset();
      eventDisRecTow2->Reset();
      eventDisGen->Reset();
      eventDisGen1->Reset();
      eventDisGen2->Reset();
      eventDisChg->Reset();
      eventDisChg1->Reset();
      eventDisChg2->Reset();
      eventDisTrk->Reset();
      eventDisTrk1->Reset();
      eventDisTrk2->Reset();
    }
    
      
    //	  cout << " myJetsub is reset" << endl;
    //	  cout << " myJetSub.matchDr  = " << myJetSub.matchDr << endl;
    myEvt.run = eventInfo->runNumber();
    myEvt.lumi = eventInfo->lumiBlock();
    myEvt.event =  eventInfo->eventNumber() ;
    
    myJetSub.cent = cent_bin; 
    myJetSub.fcalet = FCalEt;
    myJetSub.rhoCh = bge_rho_trk.rho() * 0.001;
    myJetSub.recoMass = vmass_reco[ri];
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

    myJetSub.recoChSdPt = vChSdPt_reco[ri];
    myJetSub.recoChSdMass = vChSdMass_reco[ri];
    myJetSub.recoChSdZ = vChSdZ_reco[ri];
    myJetSub.recoChSdTheta = vChSdTheta_reco[ri];
  
    myJetSub.recoChPt = vChPt_reco[ri];
    myJetSub.recoChMass = vChMass_reco[ri];
    myJetSub.recoChPtRaw = vChPtRaw_reco[ri];
    myJetSub.recoChMassRaw = vChMassRaw_reco[ri];
    myJetSub.recoChPtRcSubt = vChPtRcSubt_reco[ri];
    myJetSub.recoChMassRcSubt = vChMassRcSubt_reco[ri];
    myJetSub.recoChMassGm = vChMassGm_reco[ri];
    
    trkJetMassRcSub2 = vChMassRcSubt2_reco[ri];
    trkJetMassRcSub4 = vChMassRcSubt4_reco[ri];
    trkJetMassRcSubEV = vChMassRcSubtEV_reco[ri];
    trkJetMassRcSubMV = vChMassRcSubtMV_reco[ri];
    
    myJetSub.nTrkRaw = vNChRaw[ri]; 
    myJetSub.nTrkBkg = vNChBkg[ri]; 
    myJetSub.nTrkBkgNoWgt = vNChBkgNoWgt[ri];   

    myJetSub.drTrkJetBkg = vDrChJetBkg[ri];
    myJetSub.maxTrkPt = vMaxTrkPt[ri];
    
    

    if (_saveEvtDisplay) { 
      eventDisRecTow->Add(t_recTow[ri]);
      eventDisRecTow1->Add(t_recTow1[ri]);
      eventDisRecTow2->Add(t_recTow2[ri]);
      eventDisTrk->Add(t_trkTow[ri]);
      eventDisTrk1->Add(t_trkTow1[ri]);
      eventDisTrk2->Add(t_trkTow2[ri]);
      delete t_recTow[ri];
      delete t_recTow1[ri];
      delete t_recTow2[ri];
      delete t_trkTow[ri];
      delete t_trkTow1[ri];
      delete t_trkTow2[ri];


    }
    if (matchId != -1) {
      myJetSub.matchDr = drMin;
      myJetSub.genMass = vMass_gen[matchId];
      myJetSub.genPt = vpt_gen[matchId];
      myJetSub.genPt2 = vpt2_gen[matchId];
      myJetSub.genEta = veta_gen[matchId];
      myJetSub.genRap = vrap_gen[matchId];
      myJetSub.genPhi = vphi_gen[matchId];
      myJetSub.genRcPt = vptRc_gen[matchId];
      myJetSub.genSdMass = vSdmass_gen[matchId];
      myJetSub.genSdPt = vSdPt_gen[matchId];
      myJetSub.genSdtheta = vSdTheta_gen[matchId];
      myJetSub.genSdz = vSdz_gen[matchId];
      myJetSub.genSdPt1 = vSdpt1_gen[matchId];
      myJetSub.genSdPhi1 = vSdphi1_gen[matchId];
      myJetSub.genSdRap1 = vSdrap1_gen[matchId];
      myJetSub.genSdPt2 = vSdpt2_gen[matchId];
      myJetSub.genSdPhi2 = vSdphi2_gen[matchId];
      myJetSub.genSdRap2 = vSdrap2_gen[matchId];

      myJetSub.genNch    = vNch_gen[matchId];
      myJetSub.genChSdPt = vChSdPt_gen[matchId];
      myJetSub.genChSdMass = vChSdMass_gen[matchId];
      myJetSub.genChSdZ = vChSdZ_gen[matchId];
      myJetSub.genChSdTheta = vChSdTheta_gen[matchId];


      if (_saveEvtDisplay) {
	eventDisGen->Add(t_genTow[matchId]);
	eventDisGen1->Add(t_genTow1[matchId]);
	eventDisGen2->Add(t_genTow2[matchId]);
	eventDisChg->Add(t_chgTow[matchId]);
	eventDisChg1->Add(t_chgTow1[matchId]);
	eventDisChg2->Add(t_chgTow2[matchId]);
	delete t_genTow[matchId];
	delete t_genTow1[matchId];
	delete t_genTow2[matchId];
	delete t_chgTow[matchId];
	delete t_chgTow1[matchId];
	delete t_chgTow2[matchId];
      }
      
    }
    
    treeOut->Fill();
  
    
    if ( _saveLog) {
      cout << " eventInfo->eventNumber() == "<< eventInfo->eventNumber() << endl;
      cout << " eventInfo->runNumber() == " <<  eventInfo->runNumber() << endl;
      cout << " eventInfo->lumiBlock() " << eventInfo->lumiBlock() << endl;
    }
    
  }    
    
  if (_doJES) {   // Very important!
    delete jesProv;		        
  }  
  delete scalarPtDensity;
  // Clear vectors
  //  store->clear();
  delete store;

  corrected_selectedTrks.clear();
  selectedTrks.clear();
  selAndExcldTrks.clear();
  vmass_reco.clear();
  vpt_reco.clear();
  vptRaw_reco.clear();
  veta_reco.clear();
  vphi_reco.clear();
  vptRc_reco.clear();
  vSdmass_reco.clear();
  vSdpt_reco.clear();
  vSdtheta_reco.clear();
  vSdz_reco.clear();
  
  vChSdPt_reco.clear();
  vChSdMass_reco.clear();
  vChSdZ_reco.clear();
  vChSdTheta_reco.clear();

  vChPt_reco.clear();
  vChMass_reco.clear();
  vChPtRaw_reco.clear();
  vChMassRaw_reco.clear();
  vChPtRcSubt_reco.clear();
  vChMassRcSubt_reco.clear();
  vChMassGm_reco.clear();
  
  vMass_gen.clear();
  vpt_gen.clear();
  veta_gen.clear();
  vphi_gen.clear();
  vptRc_gen.clear();
  vSdmass_gen.clear();
  vSdPt_gen.clear();
  vSdTheta_gen.clear();
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
	//	if( m_jetCalibration ) delete m_jetCalibration; m_jetCalibration = 0;
	
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

