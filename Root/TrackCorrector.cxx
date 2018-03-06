#include "JetSubstructure/TrackCorrector.h"
#include <boost/assign.hpp>

void TrackCorrector::InitdRBinRange()
{
	//	ndRBins = 20;	// Default number of dR bins
	//	float s = 0.0025; //0.08
	//	float e = drmax;

	cout << "Setting up dR binning: ";

	//	dRrange[0] = 0.;
	//	for (int i = 1; i < ndRBins+1 ; i++)
	//	{
	//		dRrange[i] = s*pow((e/s),(1./ndRBins)*i);
	//		dRrange[i] = dRrange[0]+(i*(e - s)/ndRBins);
	//		dRrange[i] = sqrt(i)*s;
	//		dRrange[i] = i*i*s;
	//	}

	vector<double> tmp_drrange;
	boost::assign::push_back( tmp_drrange )
	( 0 )( 0.05 )( 0.1 )( 0.15 )( 0.2 )( 0.25 )( 0.3 )( 0.4 )( 0.5 )( 0.6 )( 0.7 )( 0.8 )( 1.0 )( 1.2 );

	for (int i = 0; i < tmp_drrange.size(); i++)
	{
		dRrange[i] = tmp_drrange.at(i);
	}
	ndRBins = tmp_drrange.size()-1;	// Default number of dR bins

	for (int i = 0; i < ndRBins+1; i++) cout << Form("%f...",dRrange[i]);
	cout << Form("--- %i bins",ndRBins) << endl;
	cout << endl;

}

Int_t TrackCorrector::GetdRBin(float dR)
{
	int bin = -1;
	if (dR >= dRrange[ndRBins]) bin = ndRBins;
	else
	{
		for (int i = 0; i < ndRBins ; i++)
		{
			if (dR >= dRrange[i] && dR < dRrange[i+1]) bin = i;
		}
	}

	return bin;
}

Int_t TrackCorrector::GetJetEtaBin(float eta) {
    	if (fabs(eta)<0.4) return 1;
    	if (fabs(eta)<0.8) return 2;
    	if (fabs(eta)<1.2) return 3;
    	if (fabs(eta)<1.6) return 4;
    	if (fabs(eta)<2.0) return 5;
    	if (fabs(eta)<2.4) return 6;
    	if (fabs(eta)<2.8) return 7;  
    	else return 0;
    }

Int_t TrackCorrector::GetTrackEtaBin(float eta){
	int EtaBin=0;
	if (eta < 2.5) EtaBin =4;
	if (eta < 2.) EtaBin =3;
	if (eta < 1.) EtaBin =2;
	if (eta < -1.) EtaBin =1;
	if (eta < -2.) EtaBin =0;
	
	return EtaBin;
}

Int_t TrackCorrector::GetxBin(float x, TAxis* axis){

	int xBin = -1;
	for (int i = 0; i < axis->GetNbins(); i++)
	{
		if (x >= axis->GetBinLowEdge(i+1) &&
			x < axis->GetBinUpEdge(i+1)) xBin = i;
	}
	if (xBin == -1)
	{
		cout << "Warning: did not find bin for variable at " << x << ". Using last bin" << endl;
		xBin = axis->GetNbins() - 1;
	}

	return xBin;
}

Int_t TrackCorrector::GetJetYBin(float y){
	int yBin=0;
	if (fabs(y) < 2.1) yBin =3;
	if (fabs(y) < 1.2) yBin =2;
	if (fabs(y) < 0.8) yBin =1;
	if (fabs(y) < 0.3) yBin =0;
	return yBin;
}

//Int_t TrackCorrector::GetCentBin(int cent) {
//    	return cent+1; //Shifted by unity in Laura's file 
//    }

bool TrackCorrector::PassTracktoJetBalance(float trk_pt,float jet_pt,float trk_eta,float jet_eta,int cent)
//@brief: apply 3 sigma JER track-to-jet balance cut 
{
	bool pass = true;
	float JER =  _h_JER[GetJetEtaBin(jet_eta)][cent]->GetBinContent(_h_JER[GetJetEtaBin(jet_eta)][cent]->FindBin(jet_pt));
	float TMR =  _tf1_TMR[GetTrackEtaBin(trk_eta)]->Eval(trk_pt);
	//float res = sqrt(JER*JER + TMR*TMR);
	if (trk_pt > (jet_pt + 3 * sqrt(pow(JER * jet_pt,2) + pow(TMR * trk_pt,2))) ) pass = false;
	
	//ultimate cut for overflow bins
	if (trk_pt > 2* jet_pt) pass = false;
	
	return pass;
}

float TrackCorrector::get_effcorr(float pt, float eta, int centrality, float uncert, float jet_pt, float jet_y, float R){

  float weight;
  
  //eta weights
  int etabin =  GetTrackEtaBin(eta);
  int ybin =0;// =  GetJetYBin(jet_y);  
  if (pt>_max_eff_pt) pt = _max_eff_pt; 
  if ( _eff_jety && jet_pt > _max_eff_jet_pt ) jet_pt = _max_eff_jet_pt;
  
  //jet weight, applied only within jets and for !eff_jety scenario
  float eff_jet_corr = 1.;
  float eff_pt_for_jetc = 1.; //minimum value for jet-eff correction
  if (pt>eff_pt_for_jetc) eff_pt_for_jetc = pt;
  if (R<=0.4 && !_eff_jety) eff_jet_corr = _g_eff_jept[etabin][centrality][jetpt_bins->FindBin(jet_pt)]->Eval(eff_pt_for_jetc); //correct only track within jet (not in UE)
  //Protection for empty bins
  if (eff_jet_corr<0.1) eff_jet_corr = 1.;
  
  //eta weight, now disabled  
  float eff_eta_corr = 1;
  //eff_eta_corr = _h_eff_eta[centrality]->GetBinContent(_h_eff_eta[centrality]->FindBin(eta,pt));
  
  //main efficiency
  float eff=1.;
  if (!_eff_jety || R>0.4) eff = _tf1_eff[etabin][centrality]->Eval(pt);
  else {
  	ybin = jety_bins->FindBin(fabs(jet_y))-1; //axis starts from 1, naming from 0
  	//cout << "_max_eff_jet_pt " << _max_eff_jet_pt << " ybin " << ybin << " jetpt " << jet_pt << " jetpt_bins->FindBin(jet_pt) " << jetpt_bins->FindBin(jet_pt) << " centrality " << centrality << " pt " << pt << endl;
  	eff = _g_eff_jept[ybin][centrality][jetpt_bins->FindBin(jet_pt)]->Eval(pt);
  }    
  //pt weight
  weight = 1./((eff + uncert)*eff_jet_corr*eff_eta_corr);
  
  bool test = false;
  //bool test = true;
  if ((weight < 0.1 || weight>10.) || test)
  {
	  if (R<0.4) cout << "Probably efficiency problem, weight: " << weight << " pt  " << pt << " eta " << eta << " eta bin " << etabin <<  " centrality "<< centrality << " eff_jet_corr " << eff_jet_corr << " jet pt " << jet_pt << " jet y " << jet_y << " y bin " << ybin << " dR " << R <<endl;
  	///if (fabs(jet_y)>1.2)	cout << "Probably efficiency problem, weight: " << weight << " pt  " << pt << " eta " << eta << " eta bin " << etabin <<  " centrality "<< centrality << " eff_jet_corr " << eff_jet_corr << " jet pt " << jet_pt << " jet y " << jet_y << " y bin " << ybin << " dR " << R <<endl;
  	//if (fabs(jet_y)>2. && R<0.4 && pt >1.) cout << "Probably efficiency problem, weight: " << weight << " pt  " << pt << " eta " << eta << " eta bin " << etabin <<  " centrality "<< centrality << " eff_jet_corr " << eff_jet_corr << " jet pt " << jet_pt << " jet y " << jet_y << " y bin " << ybin << " dR " << R << " eff_eta_corr " << eff_eta_corr << " eff_jet_corr " << eff_jet_corr << " uncert " << uncert <<endl;
  	if (!test) cout << "Warning: Setting eff weight to 1" << endl;
  	if (!test) weight=1; 	
  	if (!test && pt > jet_pt)
	{
		eff = _g_eff_jept[ybin][centrality][jetpt_bins->FindBin(jet_pt)]->Eval(jet_pt);
		weight=1/(eff+uncert);
		cout << Form("Warning: edge effect with trk pt > jet pt. Setting eff weight to 1/efficiency evaluated at jetpt. Eff corr = %f", weight) << endl;
	}
  }
  return weight; 
}


float TrackCorrector::get_effcorr(float pt, float eta, int centrality, float uncert){

	float weight = -1;

	//eta weights
	int etabin =  GetxBin(eta, new_trk_eta_binning);
	if (pt > _max_eff_pt) pt = _max_eff_pt;

	//main efficiency
	float eff=1.;
//	cout << "Probably efficiency problem, weight: " << weight << " eff " << eff << " pt  " << pt << " eta " << eta << " eta bin " << etabin <<  " centrality "<< centrality << endl;

//	eff = _g_eff_trketa[etabin][centrality]->Eval(pt);
	eff = _h_eff_trketa[etabin][centrality]->GetBinContent(_h_eff_trketa[etabin][centrality]->FindBin(pt));

	//pt weight
	weight = 1./((eff + uncert));

	bool test = false;
//	bool test = true;
	if ((weight < 0.1 || weight>10.) || test)
	{

		cout << "Probably efficiency problem, weight: " << weight << " eff " << eff << " pt  " << pt << " eta " << eta << " eta bin " << etabin <<  " centrality "<< centrality << endl;

		if (!test) cout << "Warning: Setting eff weight to 1" << endl;
		if (!test) weight=1;
	}
	return weight;
}


//correction based on muon momentum
void TrackCorrector::correctChTrackpT(float &pt, int charge){
	float minpt=pt;
	if(minpt<25) minpt=25;
	float s=0.0171947*TMath::Log(minpt)-0.0353475;
	pt*=(1.0-charge*s);
}

void TrackCorrector::correctChTrackpT(float& pt, float eta, float phi, int charge)
//Aply uncertainty to track momentum
{
	float charge_tmp=charge;
	if (charge_tmp>0.) charge_tmp=1;
	else if (charge_tmp<0.) charge_tmp =-1;
	
	float eta_tmp=eta;
	if(eta_tmp>2.499 )  eta_tmp= 2.499;
	if(eta_tmp<-2.499)  eta_tmp=-2.499;
	pt/=(1+charge*pt*(h_sagitta->GetBinContent(h_sagitta->FindBin(eta_tmp, phi)))*1e-3 );
}

