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
#include <JetSubstructure/McJetTrimming.h>
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
ClassImp(McJetTrimming)




McJetTrimming :: McJetTrimming ()
{
	// Here you put any code for the base initialization of variables,
	// e.g. initialize all pointers to 0.  Note that you should only put
	// the most basic initialization here, since this method will be
	// called on both the submission and the worker node.  Most of your
	// initialization code will go into histInitialize() and
	// initialize().
}



EL::StatusCode McJetTrimming :: setupJob (EL::Job& job)
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



EL::StatusCode McJetTrimming :: histInitialize ()
{
	// Here you do everything that needs to be done at the very
	// beginning on each worker node, e.g. create histograms and output
	// trees.  This method gets called before any input files are
	// connected.

	cout << " Setting  histograms" << endl;

	//Basic histograms
	const int nPtBins = 9;
	double ptBin[nPtBins+1] = {0,100,150,200,250,300,350,400,600,1000};
	h_pt = new TH1D("h_pt",";pT;",nPtBins,ptBin);
	
	h_pt_trimRateKt = new TH2D("h_pt_trRate_kt",";p_{T} (GeV);Trimmed ratio",400,0,400,400,0,2);
	h_pt_trimRateAk = new TH2D("h_pt_trRate_ak",";p_{T} (GeV);Trimmed ratio",400,0,400,400,0,2);

	TH2D* temphist_2d;
	for (int i=0;i<=nPtBins;i++) { 
	  temphist_2d = new TH2D(Form("h_m_mTr_kt_pt%d",i),";mass(GeV);Trimmed mass(GeV)",400,0,400,400,0,400);
	  h_mass_trimMassKt.push_back(temphist_2d);
	  h_mass_trimMassKt.at(i)->Sumw2();
	
	  temphist_2d = new TH2D(Form("h_m_mTr_ak_pt%d",i),";mass(GeV);Trimmed mass(GeV)",400,0,400,400,0,400);
	  h_mass_trimMassAk.push_back(temphist_2d);
	  h_mass_trimMassAk.at(i)->Sumw2();

	  temphist_2d = new TH2D(Form("h_m_trRate_kt_pt%d",i),";mass(GeV);Trimmed ratio",400,0,400,400,0,2);
	  h_mass_trimRateKt.push_back(temphist_2d);
	  h_mass_trimRateKt.at(i)->Sumw2();

	  temphist_2d = new TH2D(Form("h_m_trRate_ak_pt%d",i),";mass(GeV);Trimmed ratio",400,0,400,400,0,2);
	  h_mass_trimRateAk.push_back(temphist_2d);
	  h_mass_trimRateAk.at(i)->Sumw2();
	}	
	
	wk()->addOutput (h_pt);
	wk()->addOutput (h_pt_trimRateKt);
	wk()->addOutput (h_pt_trimRateAk);
	for (int i=0;i<nPtBins;i++) { 
	  wk()->addOutput (h_mass_trimMassKt.at(i));
	  wk()->addOutput (h_mass_trimMassAk.at(i));
	  wk()->addOutput (h_mass_trimRateKt.at(i));
	  wk()->addOutput (h_mass_trimRateAk.at(i));
}
	
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode McJetTrimming :: fileExecute ()
{
	// Here you do everything that needs to be done exactly once for every
	// single file, e.g. collect a list of all lumi-blocks processed
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode McJetTrimming :: changeInput (bool firstFile)
{
	// Here you do everything you need to do when we change input files,
	// e.g. resetting branch addresses on trees.  If you are using
	// D3PDReader or a similar service this method is not needed.
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode McJetTrimming :: initialize ()
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
	const std::string name = "McJetTrimming"; //string describing the current thread, for logging
	// Initialize and configure trigger tools

	return EL::StatusCode::SUCCESS;
}


EL::StatusCode McJetTrimming :: execute ()
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
  
  // algorithm definition 
  fastjet::JetDefinition jetDefReclus(fastjet::cambridge_algorithm, _JetRadiusAna);
  if ( _defJetRecl == 0)        jetDefReclus.set_jet_algorithm( fastjet::cambridge_algorithm ) ;
  else if ( _defJetRecl == 1)   jetDefReclus.set_jet_algorithm( fastjet::kt_algorithm ) ;
  
  fastjet::JetDefinition jetDefAk(fastjet::antikt_algorithm, _JetRadiusAna);

  fastjet::contrib::SoftDrop softdropper(_beta, _z_cut);
  if ( _saveLog) cout << softdropper.description() << endl;
  
  std::vector<fastjet::PseudoJet> truthParticles;
  const xAOD::TruthParticleContainer * genParCont = 0;
  ANA_CHECK(event->retrieve( genParCont, "TruthParticles"));
  
  for( xAOD::TruthParticleContainer::const_iterator truth_itr = genParCont->begin() ; truth_itr!= genParCont->end() ; ++truth_itr) {
    //      cout <<"  (*truth_itr)->p4().Pt() = " <<  (*truth_itr)->p4().Pt() << endl; // MeV is confirmed
    double ptTrk = (*truth_itr)->p4().Pt() * 0.001 ; 
    int thebc = (*truth_itr)->barcode();
    if ( (thebc > 0) && ( thebc < 200000 ) )  {
      truthParticles.push_back( (*truth_itr)->p4() );   // To be used for anti-kT jet reconstrucutre 
    }
  }
  
  ////////// On-the-fly jet finder ///////////////////////////////////////////////
  fastjet::ClusterSequence cs(truthParticles, jetDefAk);
  vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
  
  int nGenJetCounter =0;
  for (unsigned i = 0; i < jets.size(); i++) {  // MC anti-kT jets
    double jet_mass = jets[i].m()*0.001;
    double jet_pt = jets[i].pt()*0.001;
    double jet_eta = jets[i].pseudorapidity();  
    double jet_rap = jets[i].rapidity();
    double jet_phi = PhiInPI ( jets[i].phi() ) ;
    if (jet_pt < 50 ) continue;
    if ( fabs(jet_eta) > _etaJetCut ) continue;

    int jetPtBin = h_pt->FindBin(jet_pt);
    
    vector<fastjet::PseudoJet > akConsts = jets[i].constituents();
    
    // Truth trimmer goes here:
    float t_genTrNsub  = 0;
    float t_genTrTheta = 0;
    float t_genTrDels  = 100;
    float t_genTrMass  = -1;
    float t_genSdDels  = 100;
    
    float radiusSubJet = _rSub;
    float thef = _fCut;
    
    fastjet::JetDefinition jetDefTrimAk = fastjet::JetDefinition(fastjet::antikt_algorithm, radiusSubJet);
    fastjet::JetDefinition jetDefTrimKt = fastjet::JetDefinition(fastjet::kt_algorithm, radiusSubJet);

    fastjet::ClusterSequence trimSeqAk(akConsts, jetDefTrimAk);
    vector<fastjet::PseudoJet> trimmedJetsAk = trimSeqAk.inclusive_jets();
    fastjet::ClusterSequence trimSeqKt(akConsts, jetDefTrimKt);
    vector<fastjet::PseudoJet> trimmedJetsKt = trimSeqKt.inclusive_jets();

    fastjet::PseudoJet sumKt = fastjet::PseudoJet(0,0,0,0);
    fastjet::PseudoJet sumAk = fastjet::PseudoJet(0,0,0,0);
    
    for ( int ij= 0 ; ij < trimmedJetsAk.size() ; ij++) {
      if ( ( trimmedJetsAk[ij].pt() * 0.001 / jet_pt ) > thef ) {
	sumAk = sumAk + trimmedJetsAk[ij] ;
      }
    }
    for ( int ij= 0 ; ij < trimmedJetsKt.size() ; ij++) {
      if ( ( trimmedJetsKt[ij].pt() * 0.001 / jet_pt ) > thef ) {
	sumKt = sumKt + trimmedJetsKt[ij] ;
      }
    }
    
    /*    vector<fastjet::PseudoJet> sortedTrimmedJets = fastjet::sorted_by_pt(trimmedJets);
    if ( t_genTrNsub > 1) {
    t_genTrDels    =  ( sortedTrimmedJets[0].pt() - sortedTrimmedJets[1].pt() )*0.001 /jet_pt ;
    t_genTrTheta =  DeltaR( sortedTrimmedJets[0].phi(), sortedTrimmedJets[0].eta(), sortedTrimmedJets[1].phi(), sortedTrimmedJets[1].eta() ) ;
    }*/
    h_pt_trimRateKt->Fill( jet_pt, sumKt.pt()*0.001/jet_pt);
    h_pt_trimRateAk->Fill( jet_pt, sumAk.pt()*0.001/jet_pt);
    h_mass_trimMassKt.at(jetPtBin)->Fill(jet_mass, sumKt.m() * 0.001); 
    h_mass_trimMassAk.at(jetPtBin)->Fill(jet_mass, sumAk.m() * 0.001); 
    h_mass_trimRateKt.at(jetPtBin)->Fill(jet_mass, sumKt.m() * 0.001/jet_mass);
    h_mass_trimRateAk.at(jetPtBin)->Fill(jet_mass, sumAk.m() * 0.001/jet_mass);

  }

  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode McJetTrimming :: postExecute ()
{
	// Here you do everything that needs to be done after the main event
	// processing.  This is typically very rare, particularly in user
	// code.  It is mainly used in implementing the NTupleSvc.
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode McJetTrimming :: finalize ()
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

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode McJetTrimming :: histFinalize ()
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

