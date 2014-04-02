/////////////////////////////////////////////////////
// FITTING WITH A ROOFIT-USER DEFINED CRYSTAL BALL //
/////////////////////////////////////////////////////
#ifndef DEF_DRAWTURNON
#define DEF_DRAWTURNON

#include "fitEfficiency.h"


void fitEfficiency(float m_norm, float m_alpha, float m_n, float m_mean, float m_sigma,
		   float m_norm, float m_alpha, float m_n, float m_mean, float m_sigma,
		   TString dirResults="/home//llr/cms/ndaci/SKWork/macro/HTauTau/results/TriggerStudies/OtherSresults/",
		   TString lumi="xxx pb", int nCPU=4, 
		   int color1=kBlack, int style1=kFullCircle, int color2=kRed, int style2=kOpenSquare,
		   TString fileIn="*.root", TString image="eff_HLT_others_MuTau")
{

  // STYLE //
  gROOT->Reset();
  loadPresentationStyle();  
  gROOT->ForceStyle();

  // OUTPUT //
  TString name_image = dirResults+"/"+image;

  RooRealVar xaxis("pt","P_{T} [GeV]",0,150) ;

  RooRealVar norm("norm","N",m_norm);
  RooRealVar alpha("alpha","#alpha",m_alpha);
  RooRealVar n("n","n",m_n);
  RooRealVar mean("mean","mean",m_mean);
  RooRealVar sigma("sigma","#sigma",m_sigma);

  FuncCB cb_EB("cb_EB","Fit function EB (cb)",xaxis,mean,sigma,alpha,n,norm) ;

  

}
#endif
