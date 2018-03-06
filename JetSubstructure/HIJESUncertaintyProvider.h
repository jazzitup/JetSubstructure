#ifndef __HIJESUNCERTAINTYPROVIDER_H__
#define __HIJESUNCERTAINTYPROVIDER_H__

#include <string>
#include <vector>
#include <map>
#include <cmath>

#include <TGraph.h>
#include <TAxis.h>
#include <TLorentzVector.h>
#include <TSystem.h>

class HIJESUncertaintyProvider
{
public:
  typedef std::map<std::string,std::vector<TGraph*> > map_t;
  HIJESUncertaintyProvider(){}
  HIJESUncertaintyProvider(std::string fname);
  ~HIJESUncertaintyProvider();

  float GetUncertaintyComponent(std::string comp, float pt, float eta) const;
  void ListUncertaintyComponentKeys() const;
  void GetUncertaintyComponentKeys(std::vector<std::string>& vec) const;
  float GetTotalUncertainty(float pt, float eta) const;


  static std::string s_uncert_graph_prefix;

  inline void UseGeV(bool useGeV=true)
  {
    m_GeV=( useGeV ? 1 : 1e-3);
  }

  inline void UseJESTool(bool useJESTool=true)
  {
    m_use_JES_tool=useJESTool;
  }

  inline void UseAbsEta(bool use_abs_eta=true)
  {
    m_use_abs_eta=use_abs_eta;
  }
  
  ClassDef(HIJESUncertaintyProvider,1);

private:

  float m_GeV;
  TAxis* m_eta_axis;
  bool m_use_abs_eta;
  bool m_use_JES_tool;
  map_t m_uncertainty_graphs;

  unsigned int LookupEtaBin(float eta) const;
  
  inline float GetSign(float eta) const
  {
    return (m_use_abs_eta ?  std::abs(eta) :  eta);
  }





};


#endif
