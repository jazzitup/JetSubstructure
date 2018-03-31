#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoop/CondorDriver.h"
#include "EventLoopGrid/PrunDriver.h"
#include "EventLoopGrid/GridDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoop/OutputStream.h>
#include <TSystem.h>

#include "SampleHandler/ScanDir.h"

#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <exception>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include "JetSubstructure/BaseClass.h"
#include "JetSubstructure/JetSubstructure.h"
#include "JetSubstructure/McJetTrimming.h"

int main(int argc, char *argv[])
{

	//configuration
	int isPP;
	int isMC = 1;
	bool doJES = false;
	string dataset = ""; // 0 - pPb, 1 - 2p76TeV pp, 2 - 7TeV pp , 3 - 5p02TeV pp, 4 - 5p02TeV PbPb
	int isMB=0;
	int isHerwig=0;
	int num_evt;
	int jet_radius = 4;
	int DxAODMode;
	string reco_jet_collection;
	string test_reco_jet_collection;
	string truth_jet_collection;
	string grl;
	std::string input_directory;
	bool isGridJob;
	bool isCondor;
	std::string submitDir;
	std::string InDS;
	std::string OutDS;
	float dR_truth_matching;
	int centrality_scheme;
	int nFilesPerJob;
	float pTjetCut, truthpTjetCut, etaJetCut, ptCutJetConeExc;	

	float etaTrkCut;
	float pTtrkCutReco, pTtrkCutTruth; 
	float ptCutPostCS;
	std::string grid_configuration="";
	string weight_file;
	string centrality_weight;
	std::string output_file_name="ntuple";
	bool applyReweighting;
	float JetRadiusAna;
	bool saveLog;
	bool saveNtuple;
	bool saveEvtDisplay;
	string trk_cut_level;
	
	float chParticleMassMeV;

	float ghost_area;
	float Rktjet_bkg;
	float alphaSubtr; 

	int towerBkgKill; // -1: none, 0= 0 GeV, 1 = SoftKill

	bool doTrimming; // 
	int defTrimAlgo; // definition for jet algo used in trimming
	float fCut; //  Functioning only when doTrimming = true;
	float rSub; //  Functioning only when doTrimming = true;

	float csMaxR;
	
	float beta; 
	float z_cut;
	
	int defJetRecl;
	
	//Boost configuration
	//1) command line only: options can only be given on command line, not in config file
	boost::program_options::options_description cmd_only_options("command line only options");
	std::string config_file_name;

	cmd_only_options.add_options() //note unusual syntax when adding options!
	("help,h","produce help message")
	("config,c",boost::program_options::value<std::string>(&config_file_name),"name of configuration file");

	//2) main options: most likely to be set by user, can be specified both via command line or config
	//explaination is included in help message
	boost::program_options::options_description main_options("main options");

	main_options.add_options()
	("output_file",boost::program_options::value<std::string>(&output_file_name)->default_value("ntuple"),"name of output root file")
	("grl,g",boost::program_options::value<std::string>(&grl)->default_value("GRL_pPb_5p02TeV_2013.xml"),"grl file name")
	("trk_cut_level,C",boost::program_options::value<std::string>(&trk_cut_level)->default_value("ppTight"),"Trk cut level")
	("reco_jet_collection",boost::program_options::value<std::string>(&reco_jet_collection)->default_value("antikt4HIItrEM"),"Jet collection")
	("test_reco_jet_collection",boost::program_options::value<std::string>(&test_reco_jet_collection)->default_value("none"),"Test Jet collection")
	("truth_jet_collection",boost::program_options::value<std::string>(&truth_jet_collection)->default_value("antikt4Truth"),"Truth jet collection")
	("dataset",boost::program_options::value<string>(&dataset)->default_value("PbPb_5p02"),"Type of input data")
	("isMB",boost::program_options::value<int>(&isMB)->default_value(0),"MB or HP")
	("towerBkgKill",boost::program_options::value<int>(&towerBkgKill)->default_value(0),"Background subtraction method")
	("doTrimming",boost::program_options::value<bool>(&doTrimming)->default_value(false),"Trimming?")
	("defTrimAlgo",boost::program_options::value<int>(&defTrimAlgo)->default_value(2),"Trimming?")
	("fCut",boost::program_options::value<float>(&fCut)->default_value(0.09),"fCut")
	("rSub",boost::program_options::value<float>(&rSub)->default_value(0.2),"Reclustering ratius")
	("csMaxR",boost::program_options::value<float>(&csMaxR)->default_value(0.25),"CS max R")
        ("beta",boost::program_options::value<float>(&beta)->default_value(0.),"beta for SoftDrop")
        ("z_cut",boost::program_options::value<float>(&z_cut)->default_value(0.1),"z_cut for SoftDrop")
	("isHerwig",boost::program_options::value<int>(&isHerwig)->default_value(0),"Pythia or Herwig")
        ("isPP",boost::program_options::value<int>(&isPP)->default_value(0),"pp (1) or PbPb(0)")
	("isMC",boost::program_options::value<int>(&isMC)->default_value(1),"MC or data mode")
	("doJES",boost::program_options::value<bool>(&doJES)->default_value(false),"Do sys")
	("num_evt,n",boost::program_options::value<int>(&num_evt)->default_value(-1),"number of events, -1 runs all events")
	("saveLog",boost::program_options::value<bool>(&saveLog)->default_value(false),"Save the log?")
        ("saveNtuple",boost::program_options::value<bool>(&saveNtuple)->default_value(true),"Save Ntuple?")
	("saveEvtDisplay",boost::program_options::value<bool>(&saveEvtDisplay)->default_value(false),"Save the log?")
	("isGridJob",boost::program_options::value<bool>(&isGridJob)->default_value(0),"is it grid job?")
	("isCondor",boost::program_options::value<bool>(&isCondor)->default_value(0),"is it running on condor?")
	("JetRadiusAna",boost::program_options::value<float>(&JetRadiusAna)->default_value(1.0),"Jet Radius for analysis")
	("input_directory",boost::program_options::value<std::string>(&input_directory)->default_value("/afs/cern.ch/work/m/mrybar/xAOD/"),"name of input directory containing all files")
	("submit_directory",boost::program_options::value<std::string>(&submitDir)->default_value("submitDir"),"name of output directory")
	("InDS,i",boost::program_options::value<std::string>(&InDS)->default_value(""),"InDS for grid job")
	("OutDS,o",boost::program_options::value<std::string>(&OutDS)->default_value(""),"OutDS for grid job")
	("dR_truth_matching",boost::program_options::value<float>(&dR_truth_matching)->default_value(0.2),"dR truth matching parameter")
	("nFilesPerJob",boost::program_options::value<int>(&nFilesPerJob)->default_value(1),"Number of files per grid job")
	("jet_radius",boost::program_options::value<int>(&jet_radius)->default_value(4),"Jet radius")
	("centrality_scheme,s",boost::program_options::value<int>(&centrality_scheme)->default_value(1),"Centrality scheme")
	("jet_eta_cut",boost::program_options::value<float>(&etaJetCut)->default_value(1.2),"Jet eta cut")
	("jet_pT_cut",boost::program_options::value<float>(&pTjetCut)->default_value(10.),"Jet pT cut")
	("truth_jet_pT_cut",boost::program_options::value<float>(&truthpTjetCut)->default_value(10.),"Truth jet pT cut")
	("applyReweighting",boost::program_options::value<bool>(&applyReweighting)->default_value(0),"apply reweighting to match shape between data and MC?")
	("grid_configuration",boost::program_options::value<std::string>(&grid_configuration)->default_value(""),"Settings for grid configuration")
	("truth_track_pT_cut",boost::program_options::value<float>(&pTtrkCutTruth)->default_value(1),"Truth pT cut")
	("reco_track_pT_cut",boost::program_options::value<float>(&pTtrkCutReco)->default_value(1),"Reco pT cut")
	("reco_track_pT_cut_postCS",boost::program_options::value<float>(&ptCutPostCS)->default_value(1),"Track pT cut for SoftDrop")
	("track_eta_cut",boost::program_options::value<float>(&etaTrkCut)->default_value(2.4),"Track eta cut")
	("nFilesPerJob",boost::program_options::value<int>(&nFilesPerJob)->default_value(1),"Number of files per grid job")
	("ghost_area",boost::program_options::value<float>(&ghost_area)->default_value(0.005),"Ghost area")
	("Rktjet_bkg",boost::program_options::value<float>(&Rktjet_bkg)->default_value(0.4),"Rktjet_bkg")
	("alphaSubtr",boost::program_options::value<float>(&alphaSubtr)->default_value(1),"alphaSubtr")
	("ptCutJetConeExc",boost::program_options::value<float>(&ptCutJetConeExc)->default_value(80.),"Jet pT cut for exclusion area")
	("defJetRecl",boost::program_options::value<int>(&defJetRecl)->default_value(0),"defJetRecl")
        ("chParticleMassMeV",boost::program_options::value<float>(&chParticleMassMeV)->default_value(0),"Charged Particle Mass in MeV")

	  ;

//	if (!jet_performance_mode) doForward=false;
	//combine options types for parsing
	//all options may be specified on command line
	boost::program_options::options_description cmdline_options;
	cmdline_options.add(cmd_only_options).add(main_options);

	//all options except command line only may be specified in config file
	boost::program_options::options_description config_options;
	config_options.add(main_options);

	boost::program_options::variables_map vm;

	//first parse command line
	try
	{
		boost::program_options::store(boost::program_options::parse_command_line(argc, argv, cmdline_options), vm);
		boost::program_options::notify(vm);
	}
	catch(std::exception& e)
	{
		std::cerr << "Bad command line argument" << std::endl;
		std::cerr << e.what() << std::endl;
		return 1;
	}

	//if config was specified, also parse config file
	if(vm.count("config"))
	{
		ifstream config_stream(config_file_name.c_str());
		try
		{
			boost::program_options::store(boost::program_options::parse_config_file(config_stream,cmdline_options), vm);
			boost::program_options::notify(vm);
		}
		catch(std::exception& e)
		{
			std::cerr << "Bad config file argument" << std::endl;
			std::cerr << e.what() << std::endl;
			return 1;
		}
	}

	//Need to be in MC mode for performance
//	if (performance_mode==1) data_switch = 1;

	cout << endl << endl << "*********Configuration**********" << endl;
//	if (isMC==0) cout << "Run Number:" << run <<endl;
	if (dataset=="pPb_5p02") {cout << "Using p+Pb 5.02 TeV setup" << endl;}
	if (dataset=="pp_2p76") {cout << "Using pp 2.76 TeV setup" << endl;}
	if (dataset=="pp_7p00") {cout << "Using pp 7.00 TeV setup" << endl;}
	if (dataset=="pp_5p02") {cout << "Using pp 5.02 TeV setup" << endl;}
	if (dataset=="PbPb_5p02") {cout << "Using PbPb 5.02 TeV setup" << endl;}


	if (doJES) {cout << endl << "SYSTEMATICS MODE" << endl << endl << endl;}
	if (isHerwig) {cout << "Using Herwig" << endl;}
	if (isMB) {cout << "Using MB data" << endl;} else { cout << "Using HP data" << endl;}
	if (isMC) cout << "Running in MC mode" << endl; else { cout << "Running in data mode" << endl;}

	if (isPP) cout << "This is pp, so jet cleaing is running" << endl; 
	else { cout << "This is PbPb" << endl;}

	cout << "Input directory:          }" << input_directory << endl;
	cout << "Output directory: " << submitDir << endl;

	cout << "dR truth matching parameter: " << dR_truth_matching << endl;
	cout << "Jet collection: " << reco_jet_collection << endl;
	if (strcmp (test_reco_jet_collection.c_str(),"none") != 0) cout << "Using test jet collection: " << test_reco_jet_collection << endl;
	if (isMC==1) cout << "Truth jet collection: " << truth_jet_collection << endl;
	cout << "grl: " << grl << endl;
	cout << "Centrality scheme: "<< centrality_scheme << endl;
	cout << "jet pt cut: "<< pTjetCut << endl;
	cout << "jet |eta| cut: "<< etaJetCut << endl;
	cout <<" Nominal jet radius :  " << JetRadiusAna << endl;
        cout << "Track Selection: " << trk_cut_level << endl;
	cout << "============ Calo Tower Background  ============= " << endl;
	cout << "Tower fuctuation Kill(// -1: none, 0= 0 GeV, 1 = SoftKill): "<< towerBkgKill << endl;
	cout << "Trimming (1:yes,  0:no)               : " << doTrimming << endl;
	cout << "Trimming Algo (0=Cambridge,1=kT,2=ak) : " << defTrimAlgo << endl;
	cout << "    doTrimming*fCut: "   << fCut * doTrimming << endl;
	cout << "    doTrimming*rSub: "   << rSub * doTrimming << endl;
	cout << endl;
	cout <<"============== Constituent Subtraction =========== " <<endl;
	cout <<" Radius of kt jet: " <<  Rktjet_bkg << endl;
	cout <<" Ghost Area:       " <<  ghost_area << endl;
	cout <<" csMaxR: " << csMaxR << endl;
	cout <<" alpha (1=kT like): "<< alphaSubtr << endl;
	cout <<" Jet pT cut for Exclusion area alpha: " << ptCutJetConeExc << endl;
	cout <<"=================== Soft Drop ==================== " <<endl;
	cout <<" Reclustering algo (0=Cambridge, 1=kT): " << defJetRecl << endl;
	cout <<"     beta         : " << beta << endl;    
	cout <<"     z_cut        : " << z_cut<< endl;    
	cout << endl;
	cout <<"======== Track Selection  ====== " <<endl;
	cout << "eta cut: " << etaTrkCut << endl;
	cout << "Truth track pT cut          : " << pTtrkCutTruth << " GeV" <<endl;
	cout << "Reco  track pT cut          : " << pTtrkCutReco << " GeV"<< endl;
	cout << "Reco  track pT cut after CS : " << ptCutPostCS << " GeV" << endl;
	cout << "*_*_*_*_ Log & output _*_*_*_*_*_*_*_*_*_*_*_*" <<endl;
	cout << "Save log?                   : " << saveLog << endl; 
	cout << "Save event display?         : " << saveEvtDisplay << endl; 
	cout << "Save ntuple?                : " << saveNtuple << endl; 
	cout << "Charged Particle Mass (0MeV): " << chParticleMassMeV  << " MeV" << endl;
	if (strcmp (grid_configuration.c_str(),"") != 0)  cout << "Additional grid configuration: " << grid_configuration.c_str() << endl;

	cout << "********************************" << endl << endl << endl;

	// Set up the job for xAOD access:
	xAOD::Init().ignore();

	// Construct the samples to run on:
	SH::SampleHandler sh;

	// Get input file (! be careful about path -- not to include last folder !)
	if (!isGridJob){
		SH::ScanDir().filePattern("*").scan(sh, input_directory);
	}
	else {
		SH::scanDQ2 (sh, InDS.c_str());
		sh.setMetaString( "nc_grid_filter", "*AOD*");
	}
	sh.setMetaString( "nc_tree", "CollectionTree" );

	// Print what we found:
	sh.print();

	// Create an EventLoop job:
	cout << "Creating EventLoop job" << endl;
	EL::Job job;

	//Set outputFile
	EL::OutputStream output(output_file_name.c_str());
	job.outputAdd (output);
	//	EL::NTupleSvc *ntuple = new EL::NTupleSvc(output_file_name.c_str());
	//	job.algsAdd (ntuple);

	job.sampleHandler( sh );

	cout << "Seting maximum events to " << num_evt << endl;
	job.options()->setDouble (EL::Job::optMaxEvents, num_evt);

	// To automatically delete submitDir
	job.options()->setDouble(EL::Job::optRemoveSubmitDir, 1);


	// Add our analysis to the job:
	cout << "Add our analysis to the job" << endl;
	BaseClass* alg;
	alg = new JetSubstructure();

	//TODO write a copy const of base class

	//Set parameters
	alg->_isPP = isPP;
	alg->_isMC = isMC;
	alg->_doJES = doJES;
	alg->_dataset = dataset;
	alg->_GRL = grl;
	alg->_reco_jet_collection=reco_jet_collection.c_str();
	alg->_test_reco_jet_collection=test_reco_jet_collection.c_str();
	alg->_truth_jet_collection=truth_jet_collection.c_str();
	alg->_centrality_scheme = centrality_scheme;
	alg->_jet_radius = jet_radius;
	alg->_isMB = isMB;
	alg->_towerBkgKill = towerBkgKill;
	alg->_doTrimming = doTrimming;
	alg->_defTrimAlgo = defTrimAlgo;
	alg->_fCut = fCut;
	alg->_rSub = rSub;
	alg->_csMaxR = csMaxR;
	alg->_defJetRecl = defJetRecl;
	alg->_beta = beta;
	alg->_z_cut = z_cut;
	alg->_isHerwig = isHerwig;
	alg->_dR_truth_matching = dR_truth_matching;
	alg->_etaJetCut=etaJetCut;
	alg->_pTjetCut=pTjetCut;
	alg->_truthpTjetCut=truthpTjetCut;
	alg->_applyReweighting=applyReweighting;
	alg->_JetRadiusAna=JetRadiusAna;
	alg->_pTtrkCutReco = pTtrkCutReco;
	alg->_pTtrkCutTruth = pTtrkCutTruth;
	alg->_ptCutPostCS = ptCutPostCS;
	alg->_etaTrkCut = etaTrkCut;
        alg->_trk_cut_level = trk_cut_level;
        alg->_ghost_area = ghost_area;
        alg->_Rktjet_bkg = Rktjet_bkg;
        alg->_alphaSubtr = alphaSubtr;
	alg->_ptCutJetConeExc= ptCutJetConeExc;
	alg->_saveLog = saveLog;
	alg->_saveNtuple = saveNtuple;
	alg->_saveEvtDisplay = saveEvtDisplay;
	alg->_chParticleMassMeV = chParticleMassMeV;


	//Initialzie trigger
	//	alg->SetTrigger_chains();
	job.algsAdd( alg );
	alg->_outputName = output_file_name.c_str(); // give the name of the output to our algorithm
												 // Run the job using the local/direct driver:
	cout << "Run the job" << endl;
	//Split level protection
	job.options()->setString(EL::Job::optXaodAccessMode, EL::Job::optXaodAccessMode_athena);

	if(!isGridJob){
		EL::DirectDriver driver;
		driver.submit( job, submitDir );
	}
	else if(!isGridJob && isCondor){
		EL::CondorDriver driver;
		job.options()->setString(EL::Job::optCondorConf, "getenv = true\naccounting_group = group_atlas.uillu");
		driver.submit( job, submitDir );
	}else {
		EL::PrunDriver driver;
		driver.options()->setString("nc_outputSampleName",OutDS.c_str());
		if (nFilesPerJob > 0) driver.options()->setDouble("nc_nFilesPerJob", nFilesPerJob);
		driver.options()->setString("nc_cmtConfig", "x86_64-slc6-gcc49-opt");
		driver.options()->setDouble(EL::Job::optGridMergeOutput, 1); //run merging jobs for all samples before downloading (recommended)
		if (strcmp (grid_configuration.c_str(),"") != 0) job.options()->setString (EL::Job::optSubmitFlags, grid_configuration.c_str()); //allow task duplication
		driver.submitOnly( job, submitDir );

	}
	cout << "We are done!" << endl;

	std::cout << "done looping" << std::endl;
}


//int main( int argc, char* argv[] ) {
//
//	// Take the submit directory from the input if provided:
//	std::string submitDir = "submitDir";
//	if( argc > 1 ) submitDir = argv[ 1 ];
//
//	// Set up the job for xAOD access:
//	xAOD::Init().ignore();
//
//	// Construct the samples to run on:
//	SH::SampleHandler sh;
//
//	// use SampleHandler to scan all of the subdirectories of a directory for particular MC single file:
//	const char* inputFilePath = gSystem->ExpandPathName ("$ALRB_TutorialData/p2622/");
//	SH::ScanDir().filePattern("DAOD_SUSY1.08377960._000012.pool.root.1").scan(sh,inputFilePath);
//
//
//	// Set the name of the input TTree. It's always "CollectionTree"
//	// for xAOD files.
//	sh.setMetaString( "nc_tree", "CollectionTree" );
//
//	// Print what we found:
//	sh.print();
//
//	// Create an EventLoop job:
//	EL::Job job;
//	job.sampleHandler( sh );
//	job.options()->setDouble (EL::Job::optMaxEvents, 500);
//
//	// Add our analysis to the job:
//	MyxAODAnalysis* alg = new MyxAODAnalysis();
//	job.algsAdd( alg );
//
//	// Run the job using the local/direct driver:
//	EL::DirectDriver driver;
//	driver.submit( job, submitDir );
//
//	return 0;
//}
