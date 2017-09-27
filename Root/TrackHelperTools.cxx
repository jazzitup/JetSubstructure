#include "JetSubstructure/TrackHelperTools.h"
#include "TF1.h"

using namespace std;


int TrackHelperTools::getTypeTruth(int barcode, int pdg, int status)
{
	// https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TrackingCPMoriond2016#Truth_definitions
	// 0 - fake (should not happened) , 1 - truth primary, 2 - truth secondary, 3 - truth primary out-of-phase-space, 4 - truth secondary neutral, 5 - truth primary strange baryons
	
	//if (pdg==13 || pdg==-13) return 6;
	
	if(0<barcode && barcode<200000)
	{
		if(fabs(pdg)==3112 || fabs(pdg)==3222 || fabs(pdg)==3312 || fabs(pdg)==3334) return 5;
		else
		{
			if(status==1) return 1;
			else return 3;
		}
	}
	
	if(200000<=barcode)
	{
		//if(fabs(charge)>0.5) return 2;
		//else return 4;
		return 2;
	}
	
	return 0; // to be safe
}

