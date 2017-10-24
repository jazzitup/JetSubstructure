#include "JetSubstructure/TrackHelperTools.h"
#include "TF1.h"
#include "xAODRootAccess/tools/Message.h"

#define EL_RETURN_CHECK( CONTEXT, EXP )			\
  do {							\
    if( ! EXP.isSuccess() ) {				\
      Error( CONTEXT,					\
	     XAOD_MESSAGE( "Failed to execute: %s" ),	\
	     #EXP );					\
      return EL::StatusCode::FAILURE;			\
    }							\
  } while( false )

using namespace std;

int TrackHelperTools::getTypeTruth(int barcode, int pdg, int status, float charge)
{
	// https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TrackingCPMoriond2016#Truth_definitions
	// 0 - fake (should not happened) , 1 - truth primary, 2 - truth secondary, 3 - truth primary out-of-phase-space, 4 - truth secondary neutral, 5 - truth primary strange baryons
	
	if(0<barcode && barcode<200000)
	{
	  if(fabs(pdg)==3112 || fabs(pdg)==3222 || fabs(pdg)==3312 || fabs(pdg)==3334) return 5;
	  else { 
	    //	if(status==1 && fabs(charge)>0.5) return 1;
	        if(status==1) return 1;
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

/*
int TrackHelperTools::getType(int barcode, int pdg, int nparent, int nparent_pdg)
{
	Warning("TrackHelperTools", "getType deprecated, use getTypeReco or getTypeTruth instead" );
	return getTypeTruth(barcode, pdg, 1, 1.0);
}*/

int TrackHelperTools::getTypeReco(int barcode, int pdg, int status, float charge, float mc_prob, float mc_prob_cut)
{
	// https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TrackingCPMoriond2016#Truth_definitions
	// 0 - fake, 1 - primary, 2 - secondary, 3 - primary out-of-phase-space, 4 - secondary neutral, 5 - truth primary strange baryons
	
	if(barcode<=0 || mc_prob<mc_prob_cut) return 0;
	else return getTypeTruth(barcode, pdg, status, charge);
	
	
	return 0; // to be safe
}

void TrackHelperTools::SetTrackCuts(int& selection, int& NCuts, bool& Dod0Param, int& nSISHits_cut, int& nSIHits_cut, int& nSIHoles_cut, int& nPixHoles_cut, int& nSCTHits_cut, int& nPixHits_cut, int& nBLHits_cut, int& nIBLHits_cut, double& d0_cut, double& z0sintheta_cut, int &nTRTHits_cut, float &sig_cut){
	
	//Default (Run2) = loose primary
	
	Dod0Param=true;
	NCuts = 11;
		
	nSCTHits_cut=0;
   	nPixHits_cut=0;
   	nIBLHits_cut=0;
   	nBLHits_cut=0;
   	
   	d0_cut=1.5;
   	z0sintheta_cut=1.5;
   	
   	nSISHits_cut=1;
   	
   	nSIHits_cut=10;
   	nSIHoles_cut=2;
   	nPixHoles_cut=1;
   	nTRTHits_cut = 0;
   	sig_cut = 999.;
   	
	//Loose
	if (selection==1){
		d0_cut=2.;
		z0sintheta_cut=2.;
		nSISHits_cut=999;
	   	nSIHits_cut=7;
	}
	//Tight		
	if (selection==2){
		d0_cut=1.;
		z0sintheta_cut=1.;
		nSISHits_cut=999;
		nSIHits_cut=9;
	   	nPixHoles_cut=0;
	   	NCuts = 13;
	}
	
	//Cut flow
	if (selection==3){
		d0_cut=2.;
		z0sintheta_cut=2.;
	}
	
	if (selection==4){
		d0_cut=1.;
		z0sintheta_cut=1.;
	}
	
	//Default without SharedHits cut
	if (selection==5){
		nSISHits_cut=999.;
	}
	
	//Tight	without SharedHits cut	+ TRT cuts
	if (selection==6){
		d0_cut=1.;
		z0sintheta_cut=1.;
		nSISHits_cut=999;
		nSIHits_cut=9;
	   	nPixHoles_cut=0;
	   	nTRTHits_cut =8;
	   	NCuts = 13;
	}
	
	//Tight	without SharedHits cut	+ sign cuts 
	if (selection==7){ //6
		d0_cut=1.;
		z0sintheta_cut=1.;
		nSISHits_cut=999;
		nSIHits_cut=9;
	   	nPixHoles_cut=0;
	   	sig_cut = 3.;
	   	NCuts = 13;
	}
	//Tight	without SharedHits cut	+ sign cuts + TRT cuts
	if (selection==8){ //7 bug
		d0_cut=1.;
		z0sintheta_cut=1.;
		nSISHits_cut=999;
		nSIHits_cut=9;
	   	nPixHoles_cut=0;
	   	sig_cut = 3.;
	   	nTRTHits_cut =8;
	   	NCuts = 13;
	}
	
	//ala Run1 cust
	//Default
	if (selection==20){
		d0_cut=1.5;
		z0sintheta_cut=1.5;
		nSISHits_cut=999;
	   	nSIHits_cut=0;
	   	nSIHoles_cut=999;
	   	nPixHoles_cut=999;
	   	
	   	nSCTHits_cut=6;
   	    nPixHits_cut=1;
   	    nIBLHits_cut=1;
	}
			
    cout << "Usign following tracking cuts with the selection: " << selection << " SCT Hits: "<< nSCTHits_cut << " Pixel Hits: "<< nPixHits_cut << " BL Hits: "<< nIBLHits_cut << " d0: "<< d0_cut << " zo sin: "<< z0sintheta_cut << endl;
}

EL::StatusCode TrackHelperTools::SetCutLevel(InDet::InDetTrackSelectionTool *trackSelectionTool, string cutlevel)
{

	if (!cutlevel.std::string::compare("HITight"))
	{
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("CutLevel","HITight"));
		return EL::StatusCode::SUCCESS;
	}
	
	else if (!cutlevel.std::string::compare("FJR"))
	{
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("CutLevel","HILoose"));
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("minPt",4000.));
		return EL::StatusCode::SUCCESS;
	}

	else if (!cutlevel.std::string::compare("ppTight"))
	{
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("CutLevel","TightPrimary"));
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("maxZ0SinTheta",1.0));
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("minPt",1.));
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("maxNSiSharedModules",100));
		return EL::StatusCode::SUCCESS;
	}

	else if (!cutlevel.std::string::compare("ppTight_manual"))
	{
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("CutLevel","NoCut"));
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("maxAbsEta",2.5));
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("minNSiHits", 9));
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("maxNSiHoles", 2));
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("maxNPixelHoles", 0));
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("minEtaForStrictNSiHitsCut", 1.65));
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("minNSiHitsAboveEtaCutoff", 11));
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("useMinBiasInnermostLayersCut",1)); // s timhle asi projdou tracky, ktere nemaji hity v IBL ani v BL a ktere nemaji expected hity v BL ani v IBL; podle me takove
		return EL::StatusCode::SUCCESS;
	}

	else if (!cutlevel.std::string::compare("ppTight_tight"))
	{
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("CutLevel","TightPrimary"));
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("maxZ0SinTheta",1.0));
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("minPt",1.));
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("maxD0overSigmaD0",3.0));
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("maxZ0SinThetaoverSigmaZ0SinTheta",3.0));
		return EL::StatusCode::SUCCESS;
	}

	else if (!cutlevel.std::string::compare("NoCuts"))
	{
		EL_RETURN_CHECK("initialize()",trackSelectionTool->setProperty("CutLevel","NoCut"));
		return EL::StatusCode::SUCCESS;
	}


	else cout << "****TRACK SELECTOR TOOL NOT INITIALIZED****" << endl;

	return EL::StatusCode::SUCCESS;

}



