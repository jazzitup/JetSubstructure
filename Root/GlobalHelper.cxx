#include "JetSubstructure/GlobalHelper.h"
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <TMath.h>
//#include <xAODCaloEvent/​CaloClusterContainer.h>
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// there are following options currently implemented :
//
//  0 = do not distinguish centrality => out->iNum = 1
//  2 = Pb+Pb centrality              => out->iNum = 7, 0-10%=b, 10-20%=c, 20-30%=m, 30-40%=n, 40-50%=o, 50-60%=p, 60-80%=j
//  20 = p+Pb centrality              => out->iNum = 7, 0-10%=b, 10-20%=c, 20-30%=m, 30-40%=n, 40-50%=o, 50-60%=p, 60-100%=j
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


float GetEventPlane(const xAOD::CaloClusterContainer *hiclus)
// @brief: returns a event plane angle
{
	float totalEt_P=0;
	float totalEt_N=0;
	float m_Qnx_P;
	float m_Qny_P;
	float m_Qnx_N;
	float m_Qny_N;
	for(auto cl: *hiclus){
		double eta = cl->eta0();
		double phi = cl->phi0();
		if(fabs(eta) > 3.2 and fabs(eta) < 4.8){
			double et = cl->altE()/cosh(eta);
			if(eta > 0.) totalEt_P +=et;
			if(eta < 0.) totalEt_N +=et;
			if(eta > 0.){
				m_Qnx_P +=(et*TMath::Cos(2*phi));
				m_Qny_P +=(et*TMath::Sin(2*phi));
			}
			if(eta < 0.){
				m_Qnx_N +=(et*TMath::Cos(2*phi));
				m_Qny_N +=(et*TMath::Sin(2*phi));
			}
		}
	}

	float psiEP_N = std::atan2(m_Qny_N/totalEt_N, m_Qnx_N/totalEt_N);
	float psiEP_P =  std::atan2(m_Qny_P/totalEt_P, m_Qnx_P/totalEt_P);
	return GetAveragePsi(psiEP_N, psiEP_P);
}

float GetEventPlane(const xAOD::HIEventShapeContainer* calos){
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

Float_t GetAveragePsi(Float_t psi1, Float_t psi2)
//@brief: calculate the average event plane
{
	Float_t phase = (fabs(psi1-psi2)<(TMath::Pi()/2.))? 0:(TMath::Pi()/2.);
	return ( (psi1+psi2)/2. + phase );
}

int GetCentralityBin(Int_t centralityScheme, float FCal_Et, bool isMC)
// @brief: returns a centrality bin [0-9] based on centralityScheme and MTGlobalEvent
{
	return GetGlobalBin(centralityScheme,FCal_Et, isMC);
}

void SetRejectionHistogram(TH1D* h)
{
	h->GetXaxis()->SetBinLabel(1,"All");
	h->GetXaxis()->SetBinLabel(2,"centrality");
	h->GetXaxis()->SetBinLabel(3,"GRL");
	h->GetXaxis()->SetBinLabel(4,"Vertex");
	h->GetXaxis()->SetBinLabel(5,"LAr Quality");
	h->GetXaxis()->SetBinLabel(6,"vx_z");
	h->GetXaxis()->SetBinLabel(7,"Accepted");
	h->GetXaxis()->SetBinLabel(8,"Passed trigger");
}

int GetGlobalBin(Int_t centralityScheme, float FCal_Et, bool isMC)
// @brief: returns a global bin [0-9] based on centralityScheme, MTGlobalEvent, and MTEvent (e.g. position of leading jet wrt RP)
{
	Float_t centrality = FCal_Et;

	if (centralityScheme==1)
	{
		return 0;
	}
	else if (centralityScheme==30) // Pb+Pb 2015
	{
		// nominal 85%, full Fcal
		if ( 2.98931 	<=centrality && centrality< 6.00  ) return 0;		// 0-10%
		if ( 2.04651	<=centrality && centrality< 2.98931  ) return 1;	// 10-20%
		if ( 1.36875	<=centrality && centrality< 2.04651  ) return 2;	// 20-30%
		if ( 0.87541	<=centrality && centrality< 1.36875  ) return 3;	// 30-40%
		if ( 0.525092	<=centrality && centrality< 0.87541  ) return 4;	// 40-50%
		if ( 0.289595	<=centrality && centrality< 0.525092 ) return 5;	// 50-60%
		if ( 0.063719	<=centrality && centrality< 0.289595 ) return 6;	// 60-80%

		// Hijing doesn't have the same FCal distribution as data, let's keep everything
		if(isMC && centrality>= 6.0) return 0;
		if(isMC && 0.0<=centrality && centrality<0.063719) return 6;

		// you cannot do this -- if you want 60-70% or 70-80% you need to setup a different centrality scheme
		//if ( 0.144140 	<=centrality && centrality< 0.289595 ) return 7;		// 60-70%
		//if ( 0.063719 	<=centrality && centrality< 0.144140 ) return 8; 		// 70-80%

		return -1;
	}
	else if (centralityScheme==31) // Pb+Pb 2015, merged 40-60%
	{
		// nominal 85%, full Fcal
		if ( 2.98931 	<=centrality && centrality< 6.00  ) return 0;		// 0-10%
		if ( 2.04651	<=centrality && centrality< 2.98931  ) return 1;	// 10-20%
		if ( 1.36875	<=centrality && centrality< 2.04651  ) return 2;	// 20-30%
		if ( 0.87541	<=centrality && centrality< 1.36875  ) return 3;	// 30-40%
		if ( 0.525092	<=centrality && centrality< 0.87541  ) return 4;	// 40-50%
		if ( 0.289595	<=centrality && centrality< 0.525092 ) return 4;	// 50-60%
		if ( 0.063719	<=centrality && centrality< 0.289595 ) return 5;	// 60-80%

		// Hijing doesn't have the same FCal distribution as data, let's keep everything
		if(isMC && centrality>= 6.0) return 0;
		if(isMC && 0.0<=centrality && centrality<0.063719) return 5;

		// you cannot do this -- if you want 60-70% or 70-80% you need to setup a different centrality scheme
		//if ( 0.144140 	<=centrality && centrality< 0.289595 ) return 7;		// 60-70%
		//if ( 0.063719 	<=centrality && centrality< 0.144140 ) return 8; 		// 70-80%

		return -1;
	}
	else if (centralityScheme==32) // Pb+Pb 2015, merged bins for tests
	{
		// nominal 85%, full Fcal
		if ( 2.98931 	<=centrality && centrality< 6.00  ) return 0;		// 0-10%
		if ( 2.04651	<=centrality && centrality< 2.98931  ) return 0;	// 10-20%
		if ( 1.36875	<=centrality && centrality< 2.04651  ) return 0;	// 20-30%
		if ( 0.87541	<=centrality && centrality< 1.36875  ) return 0;	// 30-40%
		if ( 0.525092	<=centrality && centrality< 0.87541  ) return 1;	// 40-50%
		if ( 0.289595	<=centrality && centrality< 0.525092 ) return 1;	// 50-60%
		if ( 0.063719	<=centrality && centrality< 0.289595 ) return 1;	// 60-80%

		// Hijing doesn't have the same FCal distribution as data, let's keep everything
		if(isMC && centrality>= 6.0) return 0;
		if(isMC && 0.0<=centrality && centrality<0.063719) return 1;

		// you cannot do this -- if you want 60-70% or 70-80% you need to setup a different centrality scheme
		//if ( 0.144140 	<=centrality && centrality< 0.289595 ) return 7;		// 60-70%
		//if ( 0.063719 	<=centrality && centrality< 0.144140 ) return 8; 		// 70-80%

		return -1;
	}
	else if (centralityScheme==20)	// p+Pb centrality
	{
		centrality = FCal_Et;

		if ((10e9 >= centrality) && (53.74 <= centrality)) return 0; //0-10%
		if ( 40.04 <= centrality) return 1; //10-20%
		if ( 31.07 <= centrality) return 2; //20-30%
		if ( 24.10 <= centrality) return 3; //30-40%
		if ( 13.41 <= centrality) return 4; //40-60%
		if ((13.41 > centrality) && (5.585 <= centrality)) return 5;  //60-90%

		return -1;
	}
	else {
		centrality = 0;
		return -2;
	}
}

int GetCentralityNBins(Int_t centralityScheme)
{
	//Number + 1 to  include the inclusive bin
	if (centralityScheme==1) return 2;
	if (centralityScheme==2) return 7;
	if (centralityScheme==20) return 7;
	if (centralityScheme==30) return 8;
	if (centralityScheme==31) return 7;
	if (centralityScheme==32) return 3;

	else return 1;
}


void SetupBinning(Int_t scheme, string variable, Double_t array[1000], Int_t &num)
//@brief: to setup binning for histograms
{
	Float_t a, c;
	Int_t k;

	if ((scheme==0) && (variable=="pt-jet-rebin"))
	{num = 14;
		//def
		array[0] = 20; array[1] = 32;array[2] = 45;array[3] = 60;array[4] = 80;array[5] = 110;array[6] = 160;array[7] = 210;array[8] = 260;array[9] = 310;array[10] = 400;array[11] = 500; array[12] = 600; array[13] = 800; array[14] = 1000;
		//array[0] = 32; array[1] = 40;array[2] = 45;array[3] = 60;array[4] = 80;array[5] = 110;array[6] = 160;array[7] = 210;array[8] = 260;array[9] = 310;array[10] = 400;array[11] = 500; array[12] = 600; array[13] = 800; array[14] = 1000;
		printf("... pt-jet-binning : ");
		//num = 100;
		Float_t value=0;
		for (int i=0; i<=num; i++)
		{
			//array[i] = value;
			//value = value + 2;
			printf("%.0f, ", array[i]);
		}
		cout << " ... " << num << endl;
	}

	if ((scheme==0) && (variable=="pt-jet-PbPb"))
	{num = 16;
		//def
		//double bins[13]={63.096, 80.0, 100.000, 125.892,  158.488,  199.525,  251.186,  316.224,  398.101,  501.178,  630.944, 794.308, 999.970};
		double bins[17]={25.119, 31.623, 40.0, 50.119, 63.096, 79.433, 100.000, 125.892,  158.488,  199.525,  251.186,  316.224,  398.101,  501.178,  630.944, 794.308, 999.970};
		//double bins[12]={80.0, 100.000, 125.892,  158.488,  199.525,  251.186,  316.224,  398.101,  501.178,  630.944, 794.308, 999.970};
		printf("... pt-jet-binning : ");
		//num = 100;
		Float_t value=0;
		for (int i=0; i<=num; i++)
		{
			array[i]=bins[i];
			//value = value + 2;
			printf("%.0f, ", array[i]);
		}
		cout << " ... " << num << endl;
	}

	if ((scheme==0) && (variable=="pt-jet"))
	{
		num = 14;
		array[0] = 10; array[1] = 32;array[2] = 45;array[3] = 60;array[4] = 80;array[5] = 110;array[6] = 160;array[7] = 210;array[8] = 260;array[9] = 310;array[10] = 400;array[11] = 500;array[12] = 600;array[13] = 800;array[14] = 1000;
		printf("\n... pt-jet-binning : ");
		for (int i=0; i<=num; i++)
		{
			printf("%.0f, ", array[i]);
		}
		cout << " ... \n" << num << endl;
	}

	if ((scheme==0) && (variable=="eta-jet"))
	{
		printf("\n... eta-jet-binning : four simple bins ");
		array[0] = 0;
		array[1] = 0.3;
		array[2] = 0.8;
		array[3] = 1.2;
		array[4] = 2.1;
		num = 4;
	}

	if ((scheme==0) && (variable=="phi-jet"))
	{
		printf("\n... phi-jet binning ");
		num = 64;
		Float_t value=-TMath::Pi();
		for (int i=0; i<=num; i++)
		{
			array[i] = value;
			printf("%.4f, ", array[i]);
			value = value + 0.1;
		}
	}


	if ((scheme==0) && (variable=="dpT_fine"))
	{
		printf("\n... delta pT bins");
		num = 21;
		//Float_t value=0;
		array[0]=-1.;array[1]=-0.1;array[2]=-0.05;array[3]=-0.03;array[4]=-0.02;array[5]=-0.01;array[6]=-0.005;array[7]=-0.001;array[8]=-0.0005;array[9]=-0.0002;array[10]=-0.0001;array[11]=0.0;array[12]=0.0001;array[13]=0.0002;array[14]=0.0005;array[15]=0.001;array[16]=0.005;array[17]=0.01;array[18]=0.02;array[19]=0.03;array[20]=0.05;array[21]=0.1;array[21]=1;
		printf("\n... delta pt-binning : ");
		for (int i=0; i<=num; i++)
		{
			printf("%.0f, ", array[i]);
		}
		cout << " ... \n" << num << endl;
	}
	if ((scheme==0) && (variable=="dpT"))
	{
		printf("\n... delta pT bins ");
		num = 14;
		//Float_t value=0;
		array[0]=-100;array[1]=-10.;array[2]=-1.;array[3]=-0.1;array[4]=-0.01;array[5]=-0.005;array[6]=-0.001;array[7]=0.0;array[8]=0.001;array[9]=0.005;array[10]=0.01;array[11]=0.1;array[12]=1;array[13]=10;array[14]=100;
		printf("\n... delta pt-binning : ");
		for (int i=0; i<=num; i++)
		{
			printf("%.0f, ", array[i]);
		}
		cout << " ... \n" << num << endl;
	}

	if ((scheme==0) && (variable=="resp"))
	{
		printf("\n... jes-binning : ");
		num = 200;
		array[0] = -1.005;
		for (int i=1; i<=num; i++)
		{
			array[i] = array[0]+0.01*i;
			printf("%4.4f, ", array[i]);
		}
		cout << " ... " << num << endl;

	}

	if ((scheme==0) && (variable=="eta-fine"))
	{
		printf("\n... fine eta binning : ");
		num = 100;
		array[0] = -5;
		for (int i=1; i<=num; i++)
		{
			array[i] = array[0]+i*0.1;
			printf("%4.4f, ", array[i]);
		}
		cout << " ... " << num << endl;
		
	}
	
	
}


