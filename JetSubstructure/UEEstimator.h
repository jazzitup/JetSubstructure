#ifndef UEESTIMATOR_H
#define UEESTIMATOR_H

#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "JetSubstructure/JetHelperTools.h"
#include "TAxis.h"

using namespace std;
using namespace JetHelperTools;

class UEEstimator
{
  private:
    // meaning of these quantities is explained in the constructor where the default values are assigned
    Bool_t _bkgrCones[50];
    Int_t _nEta;
    Int_t _nPhi;
    Int_t _maxNCones; 
    Int_t _w_ncones;
    Float_t _etaOfCone;
    Float_t _phiOfCone;
    Float_t _deltaRToConeAxis;
    Float_t _maxConePt;
    Int_t _maxConeIndex;
    TFile * _f_weights;
    TH1F* _h_eta_w[80][7];
    Float_t _bkgrCones_hpT[50];
    TF1* _f1_trkEta[10];
    TAxis * ptaxis;
    
    TH3D * _h_v2_EP;    
    TFile * _f_flowv2;

  public:	
    // meaning of these quantities is explained in the constructor where the default values are assigned
    Float_t ptBkgrThreshold;
    Float_t jetptBkgrThreshold;
    Float_t Psi;
    Float_t m_maxjetdeltaR;

    UEEstimator() 
     {
      _nEta = 5;			// number of cones in eta direction 
      _nPhi = 8;			// number of cones in phi direction
      _maxNCones=_nEta*_nPhi;		// ID will be covered by _nEta*_nPhi=40 cones
      _w_ncones = _maxNCones;		// number of active cones (to be specified in ExcludeCones)
      for (int i=0; i<_maxNCones; i++) _bkgrCones[i] = 1;	// all cones can be used by default
      for (int i=0; i<_maxNCones; i++) _bkgrCones_hpT[i] = 0;	// maximal track pT
      ptBkgrThreshold = 10;		// by default tracks with pt>4*GeV can be coming from a hard scattering event
      jetptBkgrThreshold = 90;		// remove random cones that overlap with a jet
      m_maxjetdeltaR = 0.8;
      
      //parametrization from pPbCentrality 
      
      _f_weights = new TFile("$ROOTCOREBIN/../JetSubstructure/data/UE_eta_weight_PbPb_2015.root","read");
      ptaxis = (TAxis*) _f_weights->Get("UE_track_pt"); 
      for(int i=ptaxis->FindBin(1.);i<=ptaxis->FindBin(15.);i++){ //First bin at 1 GeV last at 15 GeV
		for(int j=0;j<7;j++){
			_h_eta_w[i][j]=(TH1F*)_f_weights->Get(Form("trk_eta_pt%i_cent%i",i,j)); //TODO update eta disitrbution
		}				
	  }
	  
	  //Parametrization from 2.76 TeV PbPb //TODO update eta disitrbution
	  /*
	  for (int i=0; i<10; i++)
        {_f1_trkEta[i] = new TF1( Form("f1_trkEta%i", i), "[0]*exp(-pow(x,2)/[1]) + [2]*exp(-pow(x,2)/[3])", -2.5, 2.5);
        }
   
      _f1_trkEta[0]->SetParameters( 234825.477990, 7.419620, 197604.101166, 571962746.539300);
      _f1_trkEta[1]->SetParameters( 150372.652701, 6.084889, 153719.937835, 1088114749.919683);
      _f1_trkEta[2]->SetParameters( 208082.808258, 13.145839, 7.746952, -0.806354);
      _f1_trkEta[3]->SetParameters( 73744.400684, 4.978877, 70348.910934, 1689597696.631848);
      _f1_trkEta[4]->SetParameters( 96950.006459, 9.767454, 1.325937, -0.713972);
      _f1_trkEta[5]->SetParameters( 47689.995455, 4.888658, 18614.531204, 414177156.622444);
      _f1_trkEta[6]->SetParameters( 90762.767213, 5.251247, 8.180272, -0.989233);
	   */
	   
	   _f_flowv2 = new TFile("$ROOTCOREBIN/../JetSubstructure/data/flowV2.root","read");
	   _h_v2_EP = (TH3D*)_f_flowv2->Get("V2_EP");
	           
      _etaOfCone = 0;	// eta position of a cone found for a given particle
      _phiOfCone = 0;	// phi position of a cone found for a given particle
      _deltaRToConeAxis = 0;	// distance between a position of a cone found for a given particle and that particle
      _maxConePt = 0; //maximum track pt asociated with the cone
      _maxConeIndex = -1;

     }// end of constructor

    // helping functions
    Float_t GetDeltaRToConeAxis() {return _deltaRToConeAxis;}
    Float_t GetetaOfConeAxis() {return _etaOfCone;}
    Float_t GetphiOfConeAxis() {return _phiOfCone;}
    Float_t GetMaxConepT() {return _maxConePt;}
    Int_t GetMaxConeIndex() {return _maxConeIndex;}     
    Float_t GetNConesWeight() {return (_w_ncones>0)? 1./_w_ncones : 0;}
    Float_t GetDeltaPsi(float phi, float psi);

    // functions to calculate weights -- explained in .cxx file
    void ExcludeConesByTrk(vector<float> &trk_pt,vector<float> &trk_eta,vector<float> &trk_phi);
    void ExcludeConesByJet(vector<float> &jet_pt,vector<float> &jet_eta,vector<float> &jet_phi);
    void ExcludeConesByJetandTrack(vector<float> &trk_pt,vector<float> &trk_eta,vector<float> &trk_phi,vector<float> &jet_pt,vector<float> &jet_eta,vector<float> &jet_phi);   
    void FindCone(float trk_pt,float trk_eta,float trk_phi);
    Float_t CalculateEtaWeight(float trk_pT, float trk_eta, float jet_eta, Int_t icent);
    Float_t CalculateFlowWeight(float trk_pt,float trk_eta,float trk_phi, float nearJetPhi, float FCalEt);
    
    Int_t GetTrackpTBin(float pt) {
    	Int_t bin=-1;
    		if (pt>0.2) bin=0; if (pt>1.) bin=1; if (pt>2.) bin=2; if (pt>3.) bin=3; if (pt>4.) bin=4; if (pt>5.) bin=5; 
    	return bin;
    }
    ~UEEstimator() {}
};


#endif
