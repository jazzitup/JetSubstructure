#include "JetSubstructure/UEEstimator.h"


void UEEstimator::ExcludeConesByTrk(vector<float> &trk_pt,vector<float> &trk_eta,vector<float> &trk_phi)
//@brief: this is to indentify which cones have to be excluded due to possibly containing a jet 
{
   for (int i=0; i<_maxNCones; i++) _bkgrCones[i] = 1;	// any cone can be used unless we prove otherwise
   for (int i=0; i<_maxNCones; i++) _bkgrCones_hpT[i] = 0.;	// reseting maximal pT

   for (unsigned int j1=0; j1<trk_pt.size(); j1++)
     {
      //if (trk.ApplyAllCuts(j1)==0) continue; //?????
      //if (trk_pt.at(j1)/1000. < ptBkgrThreshold ) continue;

      for (int iEta=0; iEta<_nEta; iEta++)
        {for (int iPhi=0; iPhi<_nPhi; iPhi++)
           {Float_t thePhi = (iPhi*2./_nPhi -1 + 1./_nPhi)*TMath::Pi();
            Float_t theEta = (iEta*2./_nEta -1 + 1./_nEta)*2.0;
            Float_t deltaR = DeltaR( trk_phi.at(j1), trk_eta.at(j1), thePhi, theEta );
            Int_t idx = _nEta * iPhi + iEta;
            if (deltaR<0.4)
              {
               if (trk_pt.at(j1) > ptBkgrThreshold ) _bkgrCones[idx] = 0;   // disable this cone since there is a jet
               if (trk_pt.at(j1) > _bkgrCones_hpT[idx] ) _bkgrCones_hpT[idx]= trk_pt.at(j1);
               }
           }
        }
     }

   _w_ncones = _maxNCones;
   for (int i=0; i<_maxNCones; i++)
     {if (_bkgrCones[i]==0) _w_ncones -= 1;
     }
   //cout << "_w_ncones (track rejection) : " << _w_ncones << endl;
}

void UEEstimator::ExcludeConesByJet(vector<float> &jet_pt,vector<float> &jet_eta,vector<float> &jet_phi)
//@brief: this is to indentify which cones have to be excluded due to possibly containing a jet 
{
   for (int i=0; i<_maxNCones; i++) _bkgrCones[i] = 1;	// any cone can be used unless we prove otherwise

   for (unsigned int j1=0; j1<jet_pt.size(); j1++)
     {
      if (jet_pt.at(j1) < jetptBkgrThreshold ) continue;

      for (int iEta=0; iEta<_nEta; iEta++)
        {for (int iPhi=0; iPhi<_nPhi; iPhi++)
           {Float_t thePhi = (iPhi*2./_nPhi -1 + 1./_nPhi)*TMath::Pi();
            Float_t theEta = (iEta*2./_nEta -1 + 1./_nEta)*2.0;
            Float_t deltaR = DeltaR( jet_phi.at(j1), jet_eta.at(j1), thePhi, theEta );
            Int_t idx = _nEta * iPhi + iEta;
            if (deltaR<m_maxjetdeltaR)
              {
               _bkgrCones[idx] = 0;                           // disable this cone since there is a jet
              }
           }
        }
     }

   _w_ncones = _maxNCones;
   for (int i=0; i<_maxNCones; i++)
     {if (_bkgrCones[i]==0) _w_ncones -= 1;
     }
   //cout << "_w_ncones (jet rejection) : " << _w_ncones << endl;
}

void UEEstimator::ExcludeConesByJetandTrack(vector<float> &trk_pt,vector<float> &trk_eta,vector<float> &trk_phi, vector<float> &jet_pt,vector<float> &jet_eta,vector<float> &jet_phi)
//@brief: this is to indentify which cones have to be excluded due to possibly containing a jet 
{
   for (int i=0; i<_maxNCones; i++) _bkgrCones[i] = 1;	// any cone can be used unless we prove otherwise

   for (unsigned int j1=0; j1<jet_pt.size(); j1++)
     {
      if (jet_pt.at(j1) < jetptBkgrThreshold ) continue;

      for (int iEta=0; iEta<_nEta; iEta++)
        {for (int iPhi=0; iPhi<_nPhi; iPhi++)
           {Float_t thePhi = (iPhi*2./_nPhi -1 + 1./_nPhi)*TMath::Pi();
            Float_t theEta = (iEta*2./_nEta -1 + 1./_nEta)*2.0;
            Float_t deltaR = DeltaR( jet_phi.at(j1), jet_eta.at(j1), thePhi, theEta );
            Int_t idx = _nEta * iPhi + iEta;
            if (deltaR<m_maxjetdeltaR)
              {
               _bkgrCones[idx] = 0;                           // disable this cone since there is a jet
              }
           }
        }
     }
   for (unsigned int j1=0; j1<trk_pt.size(); j1++)
     {
      if (trk_pt.at(j1) < ptBkgrThreshold ) continue;

      for (int iEta=0; iEta<_nEta; iEta++)
        {for (int iPhi=0; iPhi<_nPhi; iPhi++)
           {Float_t thePhi = (iPhi*2./_nPhi -1 + 1./_nPhi)*TMath::Pi();
            Float_t theEta = (iEta*2./_nEta -1 + 1./_nEta)*2.0;
            Float_t deltaR = DeltaR( trk_phi.at(j1), trk_eta.at(j1), thePhi, theEta );
            Int_t idx = _nEta * iPhi + iEta;
            if (deltaR<0.4) 
              {
               _bkgrCones[idx] = 0;                           // disable this cone since there is a jet
               if (trk_pt.at(j1) > _bkgrCones_hpT[idx] ) _bkgrCones_hpT[idx]= trk_pt.at(j1);
              }
           }
        }
     } 
   _w_ncones = _maxNCones;
   for (int i=0; i<_maxNCones; i++)
     {if (_bkgrCones[i]==0) _w_ncones -= 1;
     }
   //cout << "_w_ncones (jet rejection and track) : " << _w_ncones << endl;
}

void UEEstimator::FindCone(float trk_pt,float trk_eta,float trk_phi)
//@brief: find a minimum deltaR between a track and some non-disabled cone, that is find the cone
//@       that will be associated with this background particle, also store this deltaR in _deltaRToConeAxis 
{
   Float_t deltaROrthMin=999.;
   Float_t deltaEtaOrthMin=999.;
   Float_t deltaPhiOrthMin=999.;

   for (int iEta=0; iEta<_nEta; iEta++)
     {for (int iPhi=0; iPhi<_nPhi; iPhi++)
        {if (! _bkgrCones[iPhi*_nEta+iEta]) continue;
         Float_t thePhi = (iPhi*2./_nPhi -1 + 1./_nPhi)*TMath::Pi();
         Float_t theEta = (iEta*2./_nEta -1 + 1./_nEta)*2.0;
         Float_t deltaR = DeltaR( trk_phi, trk_eta, thePhi, theEta );
         if (deltaR < deltaROrthMin)
           {
	     deltaROrthMin = deltaR;
	     deltaEtaOrthMin = trk_eta - theEta; 
	     deltaPhiOrthMin = DeltaPhi(trk_phi, thePhi);
	     
	     _etaOfCone = theEta;
            _phiOfCone = thePhi;
            _maxConePt = _bkgrCones_hpT[iPhi*_nEta+iEta];
            _maxConeIndex = iPhi*_nEta+iEta;
           }
        }
     }

   _deltaRToConeAxis = deltaROrthMin;
   _deltaEtaToConeAxis = deltaEtaOrthMin;
   _deltaPhiToConeAxis = deltaPhiOrthMin;

}


Float_t UEEstimator::CalculateEtaWeight(float trk_pT, float trk_eta, float jet_eta, Int_t icent)
//@brief: calculate a weight that is due to a difference in the yield(eta) between a jet 
//@       position and a position of a given particle that is used to estimate the UE
//@note:  naming of some variable is not optimal as we don't want to deviate from previous codes
{
   
   float pT_temp = trk_pT;
   if (pT_temp < 1.) pT_temp = 1.; // weight in bin bellow 1 GeV is taken as weight at 1 GeV 
   int pT_bin=ptaxis->FindBin(pT_temp);
   //cout << "pt: " << trk_pT <<  " Pt bin" << pT_bin << " eta trk " <<  trk_eta << " jet eta" << jet_eta <<  " cent " << icent <<  " Period " << period << endl;
   Float_t trkEta = trk_eta;
   Float_t nearJetEta = jet_eta;
   Float_t deltaEtaOrthMin = trkEta - _etaOfCone;	// (old version: trkEta - theEta)
   deltaEtaOrthMin *= ((_etaOfCone!=0)&&(nearJetEta!=0))? ( ((_etaOfCone/fabs(_etaOfCone))==(nearJetEta/fabs(nearJetEta)))? 1:-1 ):1;
   //cout << "trk_pT " << trk_pT << " pT_bin " << pT_bin << " cent " << icent << endl; 
    Float_t w_eta = _h_eta_w[pT_bin][icent]->Interpolate( nearJetEta + deltaEtaOrthMin ) / _h_eta_w[pT_bin][icent]->Interpolate( trkEta );
   //Float_t w_eta = _f1_trkEta[icent]->Eval( nearJetEta + deltaEtaOrthMin ) / _f1_trkEta[icent]->Eval( trkEta );

   return w_eta;
}

Float_t UEEstimator::CalculateFlowWeight(float trk_pt,float trk_eta,float trk_phi, float nearJetPhi, float FCalEt)
//@brief: calculate a weight that is due to a difference in the flow between a jet 
//@       position and a position of a given particle that is used to estimate the UE
//@note:  naming of some variable is not optimal as we don't want to deviate from previous codes
{


   Float_t v2 = _h_v2_EP->GetBinContent(_h_v2_EP->GetXaxis()->FindBin(FCalEt),_h_v2_EP->GetYaxis()->FindBin(trk_pt),_h_v2_EP->GetZaxis()->FindBin(trk_eta));

			// calculate the event plane  
   Float_t w_flow  = (v2*cos(2*GetDeltaPsi( trk_phi, Psi))!=-0.5)?
                 ( (1 + 2*v2*cos(2*GetDeltaPsi( nearJetPhi, Psi)) )
                  /(1 + 2*v2*cos(2*GetDeltaPsi( trk_phi, Psi)) ) ):0; 	// modulate/demodulate

   if (w_flow > 1.4) w_flow = 1.4;	// if the correction is too large ...
   if (w_flow < 0.6) w_flow = 0.6;      // ... or too small, then take some maximal/minimal value

   return w_flow;
}

Float_t UEEstimator::GetDeltaPsi(float phi, float psi)
//@brief: distance from the event plane in convention of flow calculations
{
    Float_t diff;
    diff = fabs(phi - psi);
    while (diff > TMath::Pi()/2. ) diff = TMath::Pi() - diff;
    return fabs(diff);
}
