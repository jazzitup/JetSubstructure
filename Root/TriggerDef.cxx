#define TriggerDef_cxx
#include "JetSubstructure/JetSubstructure.h"

using namespace std;

void BaseClass::SetTrigger_chains(){

	if (_dataset.compare("pp_5p02") == 0)
	{
		if(_isMB)
		{
			_nTriggers=1;
			trigger_chains.push_back("HLT_mb_sptrk");
			
			trigger_thresholds.push_back(0);
			
			jet_pt_trig.resize(_nTriggers);
			(jet_pt_trig[0]).push_back(0); (jet_pt_trig[0]).push_back(10000.);
		}
		else
		{
			_nTriggers=6;
			_trigger_collection="a4tcemsubjesFS";
			trigger_chains.push_back("HLT_j10_320eta490");
			trigger_chains.push_back("HLT_j15_320eta490");
			trigger_chains.push_back("HLT_j25_320eta490");
			trigger_chains.push_back("HLT_j35_320eta490");
			trigger_chains.push_back("HLT_j45_320eta490");
			trigger_chains.push_back("HLT_j55_320eta490");
		
			trigger_thresholds.push_back(10);
			trigger_thresholds.push_back(15);
			trigger_thresholds.push_back(25);
			trigger_thresholds.push_back(35);
			trigger_thresholds.push_back(45);
			trigger_thresholds.push_back(55);
		
		
			jet_pt_trig.resize(_nTriggers);
			(jet_pt_trig[0]).push_back(20); (jet_pt_trig[0]).push_back(25);
			(jet_pt_trig[1]).push_back(25); (jet_pt_trig[1]).push_back(35);
			(jet_pt_trig[2]).push_back(35); (jet_pt_trig[2]).push_back(45);
			(jet_pt_trig[3]).push_back(45); (jet_pt_trig[3]).push_back(55);
			(jet_pt_trig[4]).push_back(55); (jet_pt_trig[4]).push_back(65);
			(jet_pt_trig[5]).push_back(65); (jet_pt_trig[5]).push_back(10000.);	
		}
	}
	//PbPb
	if (_dataset.compare("PbPb_5p02") == 0)
	{
		//Read prescales sets
		//Order of Y bins starting from bin #3: "HLT_j40_ion_L1TE20", "HLT_j50_ion_L1TE20", "HLT_j60_ion_L1TE50", "HLT_j75_ion_L1TE50", "HLT_j100_ion_L1TE50","HLT_j20", "HLT_j30_L1TE5", "HLT_j40_L1TE10",  "HLT_j50_L1J12","HLT_j60_L1J15","HLT_j75_L1J20","HLT_j85","HLT_mb_sptrk_ion_L1ZDC_A_C_VTE50","HLT_noalg_mb_L1TE50","HLT_mb_sptrk"]
		cout << "Setting triggers....";
		//First trigger for PbPb FF
		_first_trigger = 1; // <=> j40
		
		
		if(_isMB)
		{
			_nTriggers=2;
			_trigger_collection="a4ionemsubjesFS";
			trigger_chains.push_back("HLT_noalg_mb_L1TE50");
			trigger_chains.push_back("HLT_mb_sptrk_ion_L1ZDC_A_C_VTE50");
					
			trigger_thresholds.push_back(0);
			trigger_thresholds.push_back(0);
			
			jet_pt_trig.resize(_nTriggers);
			(jet_pt_trig[0]).push_back(0); (jet_pt_trig[0]).push_back(10000.);
			(jet_pt_trig[1]).push_back(0); (jet_pt_trig[1]).push_back(10000.);
		}
		else
		{
			_nTriggers=6;
			_trigger_collection="a4ionemsubjesFS";
			trigger_chains.push_back("HLT_j30_ion_L1TE20");//0
			trigger_chains.push_back("HLT_j40_ion_L1TE20");//1
			trigger_chains.push_back("HLT_j50_ion_L1TE20");//2
			trigger_chains.push_back("HLT_j60_ion_L1TE50");//3
			trigger_chains.push_back("HLT_j75_ion_L1TE50");//4
			trigger_chains.push_back("HLT_j100_ion_L1TE50");//5
			//trigger_chains.push_back("HLT_j100_ion_L1J10");
			//trigger_chains.push_back("HLT_j150_ion_L1TE50");
					
			trigger_thresholds.push_back(30);
			trigger_thresholds.push_back(40);
			trigger_thresholds.push_back(50);
			trigger_thresholds.push_back(60);
			trigger_thresholds.push_back(75);
			trigger_thresholds.push_back(100);
			//trigger_thresholds.push_back(100);
			//trigger_thresholds.push_back(150);
			
			
			jet_pt_trig.resize(_nTriggers);
			(jet_pt_trig[0]).push_back(45);   (jet_pt_trig[0]).push_back(58); //TODO: need to be reevaluated
			(jet_pt_trig[1]).push_back(58); (jet_pt_trig[1]).push_back(70); 
			(jet_pt_trig[2]).push_back(70); (jet_pt_trig[2]).push_back(76); 
			(jet_pt_trig[3]).push_back(76); (jet_pt_trig[3]).push_back(90); 
			(jet_pt_trig[4]).push_back(90); (jet_pt_trig[4]).push_back(116.);
			(jet_pt_trig[5]).push_back(116.); (jet_pt_trig[5]).push_back(10000.);
		}
	}
	
}

void BaseClass::SetTrigger_hist(TH2D* h){
	for(int i=1; i<=h->GetNbinsX();i++){
		//h->GetXaxis()->SetBinLabel(i,Form("j%.0f",trigger_thresholds.at(i-1)));
		h->GetXaxis()->SetBinLabel(i,Form("%s",trigger_chains.at(i-1).c_str()) );
	}
}
