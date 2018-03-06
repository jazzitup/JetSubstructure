#include "JetSubstructure/UncertProvider.h"


//Mapping of systematic
void UncertProvider::GetTrackUncert(){         
   if (uncert_index==1) uncert_class=3; //JER uncert
   else if (uncert_index>5 && uncert_index<10) uncert_class=2; //HI JES  
   else if (uncert_index>19) uncert_class=1; //intrincis JES
   else if (uncert_index>9 && uncert_index < 15) uncert_class=4; //tracking efficiency
   else if (uncert_index ==15) uncert_class=5; //Trk resolution
   else if (uncert_index ==16 || uncert_index ==17) uncert_class=6; //Trk resolution
   //else if (uncert_index ==18) uncert_class=7; //Trk charge scale for preliminary only
   else uncert_class = 0;
   cout << "Uncertainty class... " << uncert_class << endl; 
} 


void UncertProvider::CorrectJet(xAOD::Jet * reco, xAOD::Jet * truth = 0, int cent = 0, float FCalEt=0){
	switch (uncert_class){
		case 1:
		{	 
			 UncerJESIntrinsic(reco);
			 break;
		}	 	 	 	  	 
		case 2:
		{
			UncerHIJESIntrinsic(reco);
			break;
		}
		case 3:
		{
		 	 UncerJERIntrinsic_Gaus(reco, truth, cent);
		 	 break;
		}
		case 6:
		{
		 	 UncerHIJESCentrality(reco,FCalEt);
		 	 break;
		} 	 	 	 	  	 
	}
	return;
}


int UncertProvider::GetEtaUJERBin(float eta){
	int yBin=0;
	if (fabs(eta) < 0.3) yBin =0;
	if (fabs(eta) > 0.3 && fabs(eta) < 0.8) yBin =1;
	if (fabs(eta) > 0.8 && fabs(eta) < 1.2) yBin =2;
	if (fabs(eta) > 1.2 && fabs(eta) < 2.1) yBin =3;
	if (fabs(eta) > 2.1 && fabs(eta) < 2.8) yBin =4;
	if (fabs(eta) > 2.8 && fabs(eta) < 3.6) yBin =5;
	if (fabs(eta) > 3.6) yBin =6;

	return yBin;
}

void UncertProvider::UncerJERIntrinsic_Gaus(xAOD::Jet* recon, xAOD::Jet* truth, int cent)
{
       Float_t jetPtRecon = recon->pt();
  	   Float_t jetPtTruth = truth->pt();
   	   Float_t jetEtaRecon = recon->eta();
   	   Float_t jetEtaTruth = truth->eta();
   	   Float_t jetPhiRecon = recon->phi();
   	   Float_t jetMRecon = recon->m();
            
      Float_t uncertainty = 0;
      uncertainty = h_sJER_eta[GetEtaUJERBin(jetEtaRecon)]->Interpolate(jetPtRecon/1000.);

      //Float_t JER = JERhelper.GetJER(jetPtRecon/1000.,jetEtaRecon,cent);  centrality dependent JER
      Float_t JER =h_JER_eta[GetEtaUJERBin(jetEtaRecon)]->Interpolate(jetPtRecon/1000.);  
      Float_t smearingFactorSyst = sqrt(pow(JER+uncertainty,2)-pow(JER,2));
	  	  
	  Float_t correction = r.Gaus(0., smearingFactorSyst); 
	  //cout << "jet pt" << jetPtRecon << " truth " <<  jetPtTruth << " uncert: " << uncertainty << " JER " << JER << " correction " << correction << endl;
      Float_t jetPtRecoNew =jetPtRecon + jetPtTruth * correction;
      if (jetPtRecoNew<10) {jetPtRecoNew=1; jetMRecon=0;}
      recon->setJetP4( xAOD::JetFourMom_t(jetPtRecoNew,jetEtaRecon,jetPhiRecon,jetMRecon) );
}


//JES uncert
void UncertProvider::UncerJESIntrinsic(xAOD::Jet* recon)
//@brief: applies JES shift based on the JES uncertainty provider tool
//@note: significance should be +/- 1
{
   
   Float_t significance=GetSysShift(uncert_index);
   Int_t component = GetJESSysComponent(uncert_index); 
   Float_t uncertainty=0;
   Float_t jetPt = recon->pt();
   Float_t jetEta = recon->eta();
   Float_t jetPhi = recon->phi();
   Float_t jetM = recon->m();
   
   //cout << "jet pt " << jetPt << " jet eta " << jetEta << " jet phi " << jetPhi << " jet m" << jetM << endl;  
   uncertainty = 1+ significance * (jesProv.getUncertainty(component,(*recon)));
   //cout << "Uncert 0: " <<  (jesProv.getUncertainty(0,(*recon))) << endl;
   //cout << "Uncert 1: " <<  (jesProv.getUncertainty(1,(*recon))) << endl ;
   //cout << "Uncert 18: " <<  (jesProv.getUncertainty(18,(*recon))) << endl ;
   //cout << "Uncert 19: " <<  (jesProv.getUncertainty(19,(*recon))) << endl ;
   recon->setJetP4( xAOD::JetFourMom_t(jetPt*uncertainty,jetEta,jetPhi,jetM) ) ;
}


void UncertProvider::UncerHIJESIntrinsic(xAOD::Jet* recon)
//@brief: applies HI JES shift based on the JES uncertainty provider tool
//@note: significance should be +/- 1
{
   Float_t significance=GetSysShift(uncert_index);
   Int_t component = GetJESSysComponent(uncert_index);
   Float_t jetPt = recon->pt();
   Float_t jetEta = recon->eta();
   Float_t jetPhi = recon->phi();
   Float_t jetM = recon->m();
   
   Float_t uncertainty=0;
   Float_t HIJESuncertainty=0;
   if (component==1) HIJESuncertainty = sqrt(pow(HIJESProvider.GetUncertaintyComponent("flav_composition",jetPt, jetEta),2)+pow(HIJESProvider.GetUncertaintyComponent("flav_response",jetPt, jetEta),2));
   if (component==2) HIJESuncertainty = h_sJES_eta[GetEtaUJERBin(jetEta)]->Interpolate(jetPt/1000.)-1;
   //cout << "jet pt" << jetPt << " jet eta " << jetEta << " uncert " << HIJESuncertainty << endl;
   uncertainty = 1 + significance * HIJESuncertainty; 
   recon->setJetP4( xAOD::JetFourMom_t(jetPt*uncertainty,jetEta,jetPhi,jetM) );
}

void UncertProvider::UncerHIJESCentrality(xAOD::Jet* recon, float FCalEt)
//@brief: applies HI JES shift based on the JES uncertainty provider tool
//@note: significance should be +/- 1
{
   Float_t significance=GetSysShift(uncert_index);
   Float_t jetPt = recon->pt();
   Float_t jetEta = recon->eta();
   Float_t jetPhi = recon->phi();
   Float_t jetM = recon->m();
   
   Float_t uncertainty=0;
   Float_t HIJESuncertainty;
   int centrality = GetFinCentrality(FCalEt);
     
   if (centrality>60) HIJESuncertainty = 0;
   else HIJESuncertainty = (60.-(float)centrality)/60.*0.005;
   
   //cout << " FCal " << FCalEt << " cent " << centrality << " uncert " << HIJESuncertainty << endl;

   uncertainty = 1 + significance * HIJESuncertainty; 
   recon->setJetP4( xAOD::JetFourMom_t(jetPt*uncertainty,jetEta,jetPhi,jetM) );
}

int UncertProvider::GetFinCentrality(float FCalEt){
	int bin=nFineCentBins;
	
	for (int i=0;i<nFineCentBins;i++){ 
		if (fine_cent_bins[i]<FCalEt) bin = nFineCentBins - (i + 1);
		else break;
	}
	return bin;
}


string UncertProvider::GetSysName(int uncert){   
   string UncertLabel;
   UncertLabel="Non";
   int lastNonJES=19;
   switch(uncert){
		break;
		case 1: UncertLabel="JER_Intrinsic";
	    break;
	    case 2: UncertLabel="Tracking_Sign";
	    break;
	    case 3: UncertLabel="Unfolding";
	    break;
	    case 4: UncertLabel="Tracking_mcprob0p4";
	    break;
	    case 5: UncertLabel="Tracking_mcprob0p5";
	    break;
	    case 6: UncertLabel="JES_HI_1_P";
		break;
		case 7: UncertLabel="JES_HI_1_N";
		break;
		case 8: UncertLabel="JES_HI_2_P";
		break;
		case 9: UncertLabel="JES_HI_2_N";
		break;
		case 10: UncertLabel="Material_P";
		break;
		case 11: UncertLabel="Material_N";
		break;
		case 12: UncertLabel="eff_P";
		break;
		case 13: UncertLabel="eff_N";
		break;
		case 14: UncertLabel="effjet_N";
		break;
		case 15: UncertLabel="trk_resolution";
		break;
		case 16: UncertLabel="JES_HIC_P";
		break;
		case 17: UncertLabel="JES_HIC_N";
		break;
		case 18: UncertLabel="trk_ChScale";
		break;


	}
	if (uncert>lastNonJES) {
		if (uncert%2==1) UncertLabel=Form("JES_Intrinsic_%i_N",(int)(uncert-lastNonJES)/2);
		else UncertLabel=Form("JES_Intrinsic_%i_P",(int)(uncert-lastNonJES+1)/2);
	}
	return UncertLabel;	
}

int UncertProvider::GetSysShift(int uncert){   
    int shift=1;
    if (uncert%2==1) shift=-1; //odds are negativni
	return shift;
}

//Return JES components
int UncertProvider::GetJESSysComponent(int uncert){   
   int Component=0;
   int lastNonJES=19;
   if (uncert<lastNonJES){
	   switch(uncert){
			case 6: Component=1;
			break;
			case 7: Component=1;
			break;
			case 8: Component=2;
			break;
			case 9: Component=2;
			break;
		}
	}
	//insitu components needed for HI with JES2015_19NP.config (no pileup an b-jets)
	//TODO do it in better way 
	if (uncert>lastNonJES){   
	   switch(uncert){
			case 20: Component=0;//1
			break;
			case 21: Component=0;//1
			break;
			case 22: Component=1;//2
			break;
			case 23: Component=1;//2
			break;
			case 24: Component=8;//9
			break;
			case 25: Component=8;//9
			break;
			case 26: Component=9;//10
			break;
			case 27: Component=9;//10
			break;
			case 28: Component=10;//11
			break;
			case 29: Component=10;//11
			break;
			case 30: Component=11;//12
			break;
			case 31: Component=11;//12
			break;
			case 32: Component=12;//13
			break;
			case 33: Component=12;//13
			break;
			case 34: Component=13;//14
			break;
			case 35: Component=13;//14
			break;
			case 36: Component=14;//15
			break;
			case 37: Component=14;//15
			break;
			case 38: Component=15;//16
			break;
			case 39: Component=15;//16
			break;
			case 40: Component=16;//17
			break;
			case 41: Component=16;//17
			break;
		}
    }
   return Component;
}
