#include "JetSubstructure/HIJESUncertaintyProvider.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <TFile.h>

#include <TList.h>
#include <set>
#include <cmath>

std::string HIJESUncertaintyProvider::s_uncert_graph_prefix="g_uncert_";

HIJESUncertaintyProvider::HIJESUncertaintyProvider(std::string fname) : m_GeV(1),
									m_eta_axis(0),
									m_use_abs_eta(true),
									m_use_JES_tool(false)
{
  
  TString xfn = gSystem->GetFromPipe("echo $ROOTCOREBIN");
  TFile* fin=TFile::Open(xfn + "/../JetSubstructure/data/" + fname.c_str());
  if(!fin)
  {
    std::cerr << "Instantiating HIJESUncertaintyProvider. Input file does not exist: " << fin << std::endl;
    throw;
  }
  else {
  	std::cout << std::setw(40) << "HIJESUncertaintyProvider : Initialization of HI JES uncertainty provider" << std::endl;
  	std::cout << std::setw(40) << "HIJESUncertaintyProvider : Using file " << xfn << "/../JetSubstructure/data/" << fname.c_str() << std::endl;
  }	
  	
  m_eta_axis=(TAxis*)fin->Get("eta_axis");

  TList* f_keys=fin->GetListOfKeys();
  std::set<std::string> component_keys;
  for(unsigned int ikey=0; ikey<(unsigned int) f_keys->GetSize(); ikey++)
  {
    std::string key_name(f_keys->At(ikey)->GetName());
    size_t pos=key_name.find(s_uncert_graph_prefix);
    if(pos==std::string::npos) continue;
    size_t chop=key_name.find_last_of("_e");
    std::string current_key=key_name.substr(pos+s_uncert_graph_prefix.size(),chop-s_uncert_graph_prefix.size()-1);
    if(current_key.find("total")!=std::string::npos) continue;
    component_keys.insert(current_key);
  }

  std::stringstream ss;	  
  for(std::set<std::string>::const_iterator sItr=component_keys.begin(); sItr!=component_keys.end(); sItr++)
  {
    std::cout << std::setw(40) << "HIJESUncertaintyProvider : Initializing"
	      << std::setw(30) << "adding component " << *sItr
	      << std::endl;
    

    m_uncertainty_graphs.insert(std::pair<std::string,std::vector<TGraph*> >(*sItr,std::vector<TGraph*>(m_eta_axis->GetNbins(),NULL)));
    for(int i=1; i<=m_eta_axis->GetNbins(); i++)
    {
      ss.str("");
      ss << s_uncert_graph_prefix << *sItr << "_e" << i;
      m_uncertainty_graphs[*sItr][i-1]=(TGraph*)fin->Get(ss.str().c_str())->Clone(ss.str().c_str());
    }
  }
  fin->Close();
}

HIJESUncertaintyProvider::~HIJESUncertaintyProvider()
{
  //for(std::vector<TH1D*>::iterator itr=m_uncertainty_graphs.begin(); itr!=m_uncertainty_graphs.end(); itr++) delete (*itr);
  //delete m_eta_axis;
}

float HIJESUncertaintyProvider::GetUncertaintyComponent(std::string comp, float pt, float eta) const
{

  map_t::const_iterator mItr=m_uncertainty_graphs.find(comp);
  if(mItr==m_uncertainty_graphs.end())
  {
    std::cerr << "No component with name " << comp << std::endl;
    throw;
  }
  
  int eta_bin=LookupEtaBin(eta);
  return mItr->second.at(eta_bin)->Eval(pt*m_GeV);
}

void HIJESUncertaintyProvider::ListUncertaintyComponentKeys() const
{
  std::cout << "Listing uncertainty components:" << std::endl;
  for(map_t::const_iterator mItr=m_uncertainty_graphs.begin(); mItr!=m_uncertainty_graphs.end(); mItr++)
  {
    std::cout << std::setw(50) << mItr->first << std::endl;
  }
}

void HIJESUncertaintyProvider::GetUncertaintyComponentKeys(std::vector<std::string>& vec) const
{
  for(map_t::const_iterator mItr=m_uncertainty_graphs.begin(); mItr!=m_uncertainty_graphs.end(); mItr++)
  {
    vec.push_back(mItr->first);
  }
}

float HIJESUncertaintyProvider::GetTotalUncertainty(float pt, float eta) const
{
  int eta_bin=LookupEtaBin(eta);
  float total=0;
  float kine_lim=2760./std::cosh(eta);
  float pt_1=pt;
  if(pt > kine_lim) pt_1=kine_lim;
  for(map_t::const_iterator mItr=m_uncertainty_graphs.begin(); mItr!=m_uncertainty_graphs.end(); mItr++)
  {
    if(mItr->first.compare("baseline")==0 && m_use_JES_tool) continue;
    float uncert=mItr->second.at(eta_bin)->Eval(pt_1*m_GeV);
    total+=uncert*uncert;
  }
  return std::sqrt(total);
}



unsigned int HIJESUncertaintyProvider::LookupEtaBin(float eta) const
{
  int eta_bin=m_eta_axis->FindBin(GetSign(eta));
  if(eta_bin > m_eta_axis->GetNbins()) return m_eta_axis->GetNbins()-1;
  if(eta_bin==0) return 0;
  else return eta_bin-1;
}
