#include "JetSubstructure/JetCorrector.h"


int JetCorrector::GetJetYBin(float y){
	int yBin=0;
	if (y < 1.6) yBin =6;
	if (y < 1.2) yBin =5;
	if (y < 0.8) yBin =4;
	if (y < 0.3) yBin =3;
	if (y < -0.3) yBin =2;
	if (y < -0.8) yBin =1;
	if (y < -1.2) yBin =0;
	return yBin;
}

Int_t JetCorrector::GetJetEtaBin(float eta) {
    	if (fabs(eta)<0.4) return 1;
    	if (fabs(eta)<0.8) return 2;
    	if (fabs(eta)<1.2) return 3;
    	if (fabs(eta)<1.6) return 4;
    	if (fabs(eta)<2.0) return 5;
    	if (fabs(eta)<2.4) return 6;
    	if (fabs(eta)<2.8) return 7;  
    	else return 0;
    }

float JetCorrector::GetJetWeight(double pt, double eta, double phi)
{
	int xb=event_weight_histo->GetXaxis()->FindBin(pt);
	int yb=event_weight_histo->GetYaxis()->FindBin(eta);
	int zb=event_weight_histo->GetZaxis()->FindBin(phi);
	float jet_weight=event_weight_histo->GetBinContent(xb,yb,zb);

	return jet_weight;

}


float JetCorrector::GetFCalWeight(float FCalEt) {
	int centile = -1;
	float event_weight_fcal=1;
	for (int i=0; i<centiles.size(); i++)
	{
		if ( FCalEt < range_hi.at(i) && FCalEt >= range_lo.at(i) ) centile = centiles.at(i);
	}
	
	event_weight_fcal = cent_weight_histo->GetBinContent(cent_weight_histo->FindBin(centile));
	return event_weight_fcal;			
}

bool JetCorrector::MCJetJERClean(float truth_jet_pt,float reco_jet_pt, float reco_jet_eta, int cent){
	bool pass = true;
	float JER =  GetJER(reco_jet_pt, reco_jet_eta, cent);
	if (fabs (( truth_jet_pt - reco_jet_pt) /  truth_jet_pt ) > JERcut * JER ) pass = false;
	return pass;
	
}

float JetCorrector::GetJER(float reco_jet_pt, float reco_jet_eta, int cent){
	float JER =  _h_JER[GetJetEtaBin(reco_jet_eta)][cent]->GetBinContent(_h_JER[GetJetEtaBin(reco_jet_eta)][cent]->FindBin(reco_jet_pt));
	return JER;	
}

float JetCorrector::GetJetReweightingFactor(double pt, double eta, int cent){
	//int jet_eta_bin = GetJetEtaBin(eta);
	//cout << " jet w: " << jet_spectra_weight[7][cent]->Eval(pt) << " pt: " << pt << " cent: " << cent << endl;
	return jet_spectra_weight[7][cent]->Eval(pt);
}

float JetCorrector::GetFFReweightingFactor(double z, double jet_pt, double jet_eta, int cent){
	//int jet_eta_bin = GetJetEtaBin(eta);
	if (jet_pt < min_jet_pt) jet_pt = min_jet_pt; //Minimum reco pT
	if (jet_pt > max_jet_pt) jet_pt = max_jet_pt; //Maximum reco pT  
	int jet_pt_bin = jet_pt_binning->FindBin(jet_pt);
	int z_bin = FF_weight[7][cent][jet_pt_bin]->FindBin(z);
	//if ( cent ==5 ) cout << " FF w: " << FF_weight[7][cent][jet_pt_bin]->GetBinContent(z_bin) << " z " << z <<" jet pt: " << jet_pt << " cent: " << cent << endl;
	return FF_weight[7][cent][jet_pt_bin]->GetBinContent(z_bin);
}

float JetCorrector::GetCHPSReweightingFactor(double pt, double jet_pt, double jet_eta, int cent){
	//int jet_eta_bin = GetJetEtaBin(eta);
	if (jet_pt < min_jet_pt) jet_pt = min_jet_pt; //Minimum reco pT
	if (jet_pt > max_jet_pt) jet_pt = max_jet_pt; //Maximum reco pT
	int jet_pt_bin = jet_pt_binning->FindBin(jet_pt);
	int pt_bin = CHPS_weight[7][cent][jet_pt_bin]->FindBin(pt);
	//cout << " CHPS w: " << CHPS_weight[7][cent][jet_pt_bin]->GetBinContent(pt_bin) << " jet pt: " << jet_pt << " cent: " << cent << endl;
	return CHPS_weight[7][cent][jet_pt_bin]->GetBinContent(pt_bin);
}
