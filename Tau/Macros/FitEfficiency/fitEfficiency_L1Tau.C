/////////////////////////////////////////////////////
// FITTING WITH A ROOFIT-USER DEFINED CRYSTAL BALL //
/////////////////////////////////////////////////////

#ifndef DEF_FASTEFF
#define DEF_FASTEFF

// General C++
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

// RooFit headers
#include "RooGlobalFunc.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooEfficiency.h"
#include "RooDataSet.h"
#include "RooBinning.h"
#include "RooHist.h"
#include "RooWorkspace.h"

// Root headers
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TFrame.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

// Personal headers
#include "FuncCB.h"
#include "../Common/tdrstyle.h"

using namespace RooFit ;

void loadPresentationStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(0.01);
  gStyle->SetLabelOffset(0.005, "XYZ");
  gStyle->SetTitleSize(0.07, "XYZ");
  gStyle->SetTitleFont(22,"X");
  gStyle->SetTitleFont(22,"Y");
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetHistLineWidth(2);

  //gROOT->ProcessLine(".L ../Common/tdrstyle.C");
  //gROOT->ProcessLine("setTDRStyle()");
  setTDRStyle();
}

void fastEfficiencyNadir(unsigned int iEG, int iECAL1, int iColl1, int iECAL2, int iColl2,
			 TString dirIn="/home/llr/cms/ndaci/SKWork/macro/skEfficiency/tagAndProbe/EfficiencyStudy/SingEle-May10ReReco/UPDATE2/Tag80Probe95/", 
			 TString lumi="200 pb", int nCPU=4, 
			 int color1=kBlack, int style1=kFullCircle, int color2=kRed, int style2=kOpenSquare,
			 TString probe="WP80", TString tag="WP80", TString fileIn="tree_effi_TagProbe.root")
{
  // STYLE //
  gROOT->Reset();
  loadPresentationStyle();  
  gROOT->ForceStyle();

  // EG THRESHOLDS //
  const int nEG = 71;
  double thres[nEG];
  for(int i=0 ; i<nEG ; i++) thres[i]=i;

  TString names[nEG];
  ostringstream ossi;
  for(int i=0;i<(int)nEG;i++) {
    ossi.str("");
    ossi << thres[i] ;
    names[i] = ossi.str();
  }

  // NAMES //
  const int nECAL=2;
  const int nColl=2;

  TString name_leg_ecal[nECAL] = {"Barrel","Endcaps"};
  TString name_leg_coll[nColl] = {"Online","Emulation"};  

  TString name_ecal[nECAL] = {"_EB","_EE"};
  TString name_coll[nColl] = {"_N","_M"};

  TString dirResults = dirIn + "/turnons/EG"+names[iEG]+"/" ;
  TString name_image = 
    dirResults + "eff_EG"+names[iEG]+"_tag"+tag+"_probe"+probe+name_ecal[iECAL1]+name_coll[iColl1]+"_vs"+name_ecal[iECAL2]+name_coll[iColl2] ;

  // Output log //
  ofstream fichier(name_image+".txt", ios::out);


  // BINNING //
  const int nbins[nEG] = {29,29,29,29,21,21,21,22,22,21,22,21,22,18,19,18,18,18,18,20,20,20,20,19,20,20,20,20,21,21,
			  21,21,21,21,21,21,21,21,21,21, //EG30
			  22,22,22,22,22,22,22,22,22,22, //EG40
			  29,29,29,29,29,29,29,29,29,29, //EG50
			  29,29,29,29,29,29,29,29,29,29};//EG60

  Double_t bins_0[29] = {1,1.5,1.8,2,2.2,2.4,2.6,2.8, 3, 3.5, 4,4.2,4.5,4.7,5,5.5,6,6.5,7,7.5,8,8.5,9,10,12,15,20,50,150};// EG0
  Double_t bins_1[29] = {1,1.5,1.8,2,2.2,2.4,2.6,2.8, 3, 3.5, 4,4.2,4.5,4.7,5,5.5,6,6.5,7,7.5,8,8.5,9,10,12,15,20,50,150};// EG1
  Double_t bins_2[29] = {1,1.5,1.8,2,2.2,2.4,2.6,2.8, 3, 3.5, 4,4.2,4.5,4.7,5,5.5,6,6.5,7,7.5,8,8.5,9,10,12,15,20,50,150};// EG2 
  Double_t bins_3[29] = {1,1.5,1.8,2,2.2,2.4,2.6,2.8, 3, 3.5, 4,4.2,4.5,4.7,5,5.5,6,6.5,7,7.5,8,8.5,9,10,12,15,20,50,150};// EG3

  Double_t bins_4[21] = {1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 15, 17, 19, 21, 27, 32, 41, 50, 60, 70, 150}; // EG4
  Double_t bins_5[21] = {2, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 24, 31, 40, 50, 60, 70, 150}; // EG5
  Double_t bins_6[21] = {3, 4, 5, 6, 7, 8, 9, 11, 13, 15, 17, 19, 21, 23, 27, 32, 41, 50, 60, 70, 150}; // EG6

  Double_t bins_7[22] = {2, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 24, 26, 31, 40, 50, 60, 70, 150}; // EG7
  Double_t bins_8[22] = {3, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 21, 23, 25, 27, 32, 41, 50, 60, 70, 150}; // EG8
  Double_t bins_9[21] = {4, 6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 20, 22, 24, 26, 31, 40, 50, 60, 70, 150}; // EG9

  Double_t bins_10[22] = {5, 7, 8, 9, 10, 11, 12, 13, 15, 17, 19, 21, 23, 25, 27, 29, 32, 41, 50, 60, 70, 150}; // EG10
  Double_t bins_11[21] = {6, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 24, 26, 28, 31, 40, 50, 60, 70, 150}; // EG11
  Double_t bins_12[22] = {5, 7, 9, 10, 11, 12, 13, 14, 15, 17, 19, 21, 23, 25, 27, 29, 32, 41, 50, 60, 70, 150}; // EG12

  Double_t bins_13[18] = {5, 7, 9, 11, 12, 13, 14, 15, 17, 19, 22, 25, 29, 37, 50, 60, 70, 150}; // EG13
  Double_t bins_14[19] = {6, 8, 10, 12, 13, 14, 15, 16, 18, 20, 22, 25, 30, 35, 40, 50, 60, 70, 150}; // EG14
  Double_t bins_15[18] = {5, 7, 9, 11, 13, 14, 15, 16, 17, 19, 22, 25, 29, 37, 50, 60, 70, 150}; // EG15

  Double_t bins_16[18] = {8, 10, 12, 14, 16, 17, 18, 19, 20, 22, 25, 30, 35, 40, 50, 60, 70, 150}; // EG16
  Double_t bins_17[18] = {9, 11, 13, 15, 16, 17, 18, 19, 21, 23, 25, 30, 35, 40, 50, 60, 70, 150}; // EG17
  Double_t bins_18[18] = {8, 10, 12, 14, 16, 17, 18, 19, 20, 22, 25, 30, 35, 40, 50, 60, 70, 150}; // EG18

  Double_t bins_19[20] = {9, 11, 13, 15, 17, 18, 19, 20, 21, 23, 25, 27, 30, 35, 40, 45, 50, 60, 70, 150}; // EG19
  Double_t bins_20[20] = {8, 10, 12, 14, 16, 18, 19, 20, 21, 22, 24, 26, 30, 35, 40, 45, 50, 60, 70, 100}; // EG20
  Double_t bins_21[20] = {9, 11, 13, 15, 17, 19, 20, 21, 22, 23, 25, 27, 30, 35, 40, 45, 50, 60, 70, 150}; // EG21

  Double_t bins_22[20] = {10, 12, 14, 16, 18, 20, 21, 22, 23, 24, 26, 28, 30, 35, 40, 45, 50, 60, 70, 150}; // EG22
  Double_t bins_23[19] = {11, 13, 15, 17, 19, 21, 22, 23, 24, 25, 27, 30, 35, 40, 45, 50, 60, 70, 150}; // EG23
  Double_t bins_24[20] = {10, 12, 14, 16, 18, 20, 22, 23, 24, 25, 26, 28, 30, 35, 40, 45, 50, 60, 70, 150}; // EG24

  Double_t bins_25[20] = {11, 13, 15, 17, 19, 21, 23, 24, 25, 26, 27, 29, 30, 35, 40, 45, 50, 60, 70, 150}; // EG25
  Double_t bins_26[20] = {10, 12, 14, 16, 18, 20, 22, 24, 25, 26, 27, 28, 30, 35, 40, 45, 50, 60, 70, 150}; // EG26
  Double_t bins_27[20] = {11, 13, 15, 17, 19, 21, 23, 25, 26, 27, 28, 29, 33, 35, 40, 45, 50, 60, 70, 150}; // EG27

  Double_t bins_28[21] = {10, 12, 14, 16, 18, 20, 22, 24, 26, 27, 28, 29, 30, 32, 35, 40, 45, 50, 60, 70, 150}; // EG28
  Double_t bins_29[21] = {11, 13, 15, 17, 19, 21, 23, 25, 27, 28, 29, 30, 31, 33, 35, 40, 45, 50, 60, 70, 150}; // EG29
  Double_t bins_30[21] = {10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 29, 30, 31, 32, 35, 40, 45, 50, 60, 70, 150}; // EG30

  Double_t bins_40[22] = {10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 38, 39, 40, 42, 45, 50, 60, 70, 150}; // EG40
  Double_t bins_50[29] = {10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 48, 50, 55, 60, 70, 90, 110, 130, 150, 170, 190}; // EG50
  Double_t bins_60[29] = {10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 48, 50, 55, 60, 70, 90, 110, 130, 150, 170, 190}; // EG60

  vector< Double_t* > bins;
  bins.push_back( bins_0 ); bins.push_back( bins_1 ); bins.push_back( bins_2 ); bins.push_back( bins_3 ); bins.push_back( bins_4 ); 
  bins.push_back( bins_5 ); bins.push_back( bins_6 ); bins.push_back( bins_7 ); bins.push_back( bins_8 ); bins.push_back( bins_9 ); 
  bins.push_back( bins_10 ); bins.push_back( bins_11 ); bins.push_back( bins_12 ); bins.push_back( bins_13 ); bins.push_back( bins_14 ); 
  bins.push_back( bins_15 ); bins.push_back( bins_16 ); bins.push_back( bins_17 ); bins.push_back( bins_18 ); bins.push_back( bins_19 ); 
  bins.push_back( bins_20 ); bins.push_back( bins_21 ); bins.push_back( bins_22 ); bins.push_back( bins_23 ); bins.push_back( bins_24 ); 
  bins.push_back( bins_25 ); bins.push_back( bins_26 ); bins.push_back( bins_27 ); bins.push_back( bins_28 ); bins.push_back( bins_29 ); 

  for(int iV=0 ; iV<10 ; iV++) bins.push_back( bins_30 );
  for(int iV=0 ; iV<10 ; iV++) bins.push_back( bins_40 );
  for(int iV=0 ; iV<10 ; iV++) bins.push_back( bins_50 );
  for(int iV=0 ; iV<10 ; iV++) bins.push_back( bins_60 );

  RooBinning binning = RooBinning(nbins[iEG]-1, bins[iEG], "binning");


  // INPUT DATA //
  TFile* f1 = TFile::Open(dirIn+"/"+fileIn);

  TTree* treenew;
  TTree* treenew_2;

  treenew = (TTree*) gDirectory->Get( "treenew"+name_ecal[iECAL1]+name_coll[iColl1] ) ;
  treenew_2 = (TTree*) gDirectory->Get( "treenew"+name_ecal[iECAL2]+name_coll[iColl2] ) ;

  TString name_scet[2], name_scdr[2], name_l1bin[2];
  name_scet[0] = "sc_et"+name_ecal[iECAL1]+name_coll[iColl1];
  name_scet[1] = "sc_et"+name_ecal[iECAL2]+name_coll[iColl2];
  name_scdr[0] = "sc_dr"+name_ecal[iECAL1]+name_coll[iColl1];
  name_scdr[1] = "sc_dr"+name_ecal[iECAL2]+name_coll[iColl2];
  
  name_l1bin[0] = "l1_"+names[iEG]+name_ecal[iECAL1]+name_coll[iColl1];
  name_l1bin[1] = "l1_"+names[iEG]+name_ecal[iECAL2]+name_coll[iColl2];

  RooRealVar et_plot(name_scet[0],name_scet[0],0,150) ;
  RooRealVar dr(name_scdr[0],name_scdr[0],0.5,1.5) ; 
  RooRealVar et_plot2(name_scet[1],name_scet[1],0,150) ;
  RooRealVar dr2(name_scdr[1],name_scdr[1],0.5,1.5) ;

  // Acceptance state cut (1 or 0)
  RooCategory cut(name_l1bin[0],name_l1bin[0]) ;
  cut.defineType("accept",1) ;
  cut.defineType("reject",0) ;
  RooCategory cut2(name_l1bin[1],name_l1bin[1]) ;
  cut2.defineType("accept",1) ;
  cut2.defineType("reject",0) ;
  
  // PARAMETRES ROOFIT CRYSTAL BALL
  RooRealVar norm("norm","N",1,0.6,1);
  RooRealVar alpha("alpha","#alpha",0.671034,0.01,8);
  RooRealVar n("n","n",4.07846,1.1,35);
  RooRealVar mean("mean","mean",20.8,0,100);
  //mean.setVal(thres[iEG]);
  RooRealVar sigma("sigma","#sigma",0.972825,0.01,5);
  //RooRealVar pedestal("pedestal","pedestal",0.01,0,0.4);

  RooRealVar norm2("norm2","N",0.999069,0.6,1);
  RooRealVar alpha2("alpha2","#alpha",0.492303,0.01,8);
  RooRealVar n2("n2","n",11.6694,1.1,35);
  RooRealVar mean2("mean2","mean",21.4582,0,100);
  //mean2.setVal(thres[iEG]);
  RooRealVar sigma2("sigma2","#sigma",1.19,0.01,5);
  //RooRealVar pedestal2("pedestal2","pedestal",0.01,0,0.4);

  FuncCB cb("cb","Crystal Ball Integree",et_plot,mean,sigma,alpha,n,norm) ;
  FuncCB cb2("cb2","Crystal Ball Integree",et_plot2,mean2,sigma2,alpha2,n2,norm2) ;
  
  // EFFICIENCY //
  RooEfficiency eff("eff","efficiency",cb,cut,"accept");
  RooEfficiency eff2("eff2","efficiency",cb2,cut2,"accept");

  // DATASETS //
  RooDataSet dataSet("data","data",RooArgSet(et_plot, cut,dr),Import(*treenew)); 
  RooDataSet dataSet2("data2","data2",RooArgSet(et_plot2, cut2,dr2),Import(*treenew_2));

  dataSet.Print();
  dataSet2.Print();
  
  // PLOT //
  RooPlot* frame = et_plot.frame(Bins(18000),Title("Fitted efficiency")) ;
  RooPlot* frame2 = et_plot2.frame(Bins(18000),Title("Fitted efficiency")) ;

  dataSet.plotOn(frame, Binning(binning), Efficiency(cut), MarkerColor(color1), LineColor(color1), MarkerStyle(style1) );
  dataSet2.plotOn(frame2, Binning(binning), Efficiency(cut2), MarkerColor(color2), LineColor(color2), MarkerStyle(style2) );


  /////////////////////// FITTING /////////////////////////////

  double fit_cuts_min = thres[iEG]-1.5 ;
  double fit_cuts_max = 150;

  et_plot.setRange("interesting",fit_cuts_min,fit_cuts_max);
  et_plot2.setRange("interesting",fit_cuts_min,fit_cuts_max);

  RooFitResult* roofitres1 = new RooFitResult("roofitres1","roofitres1");
  RooFitResult* roofitres2 = new RooFitResult("roofitres2","roofitres2");

  fichier << "Fit characteristics :"   << endl ;
  fichier << "EG "     << names[iEG] << endl ;
  fichier << "Fit Range , EB Coll : [" << fit_cuts_min << "," << fit_cuts_max << "]" << endl ;
  fichier << "Fit Range , EE Coll : [" << fit_cuts_min << "," << fit_cuts_max << "]" << endl ;
  fichier << "----------------------"  << endl ;

  // Fit #1 //
  roofitres1 = eff.fitTo(dataSet,ConditionalObservables(et_plot),Range("interesting"),Minos(kTRUE),Warnings(kFALSE),NumCPU(nCPU),Save(kTRUE));
  
  cb.plotOn(frame,LineColor(color1),LineWidth(2));

  double res_norm1  = norm.getVal();
  double err_norm1  = norm.getErrorLo();
  double res_mean1  = mean.getVal();
  double err_mean1  = mean.getError();
  double res_sigma1 = sigma.getVal();
  double err_sigma1 = sigma.getError();
  double res_n1     = n.getVal();
  double err_n1     = n.getError();
  double res_alpha1 = alpha.getVal();
  double err_alpha1 = alpha.getError();

  fichier << "<----------------- EB ----------------->" << endl
	  << "double res_mean="  << res_mean1  << "; "
	  << "double res_sigma=" << res_sigma1 << "; "
          << "double res_alpha=" << res_alpha1 << "; "
          << "double res_n="     << res_n1     << "; "
          << "double res_norm="  << res_norm1  << "; "
	  << endl
	  << "double err_mean="  << err_mean1  << "; "
	  << "double err_sigma=" << err_sigma1 << "; "
          << "double err_alpha=" << err_alpha1 << "; "
          << "double err_n="     << err_n1     << "; "
          << "double err_norm="  << err_norm1  << "; "
	  << endl;

  // Fit #2 //
  roofitres2 = eff2.fitTo(dataSet2,ConditionalObservables(et_plot2),Range("interesting"),Minos(kTRUE),Warnings(kFALSE),NumCPU(nCPU),Save(kTRUE));
 
  cb2.plotOn(frame2,LineColor(color2),LineWidth(2));

  double res_norm2  = norm2.getVal();
  double err_norm2  = norm2.getErrorLo();
  double res_mean2  = mean2.getVal();
  double err_mean2  = mean2.getError();
  double res_sigma2 = sigma2.getVal();
  double err_sigma2 = sigma2.getError();
  double res_n2     = n2.getVal();
  double err_n2     = n2.getError();
  double res_alpha2 = alpha2.getVal();
  double err_alpha2 = alpha2.getError();

  fichier << "<----------------- EE ----------------->" << endl
	  << "double res_mean="  << res_mean2  << "; "
	  << "double res_sigma=" << res_sigma2 << "; "
	  << "double res_alpha=" << res_alpha2 << "; "
	  << "double res_n="     << res_n2     << "; "
	  << "double res_norm="  << res_norm2  << "; "
	  << endl
	  << "double err_mean="  << err_mean2  << "; "
	  << "double err_sigma=" << err_sigma2 << "; "
	  << "double err_alpha=" << err_alpha2 << "; "
	  << "double err_n="     << err_n2     << "; "
	  << "double err_norm="  << err_norm2  << "; "
	  << endl;
    

  ////////////////////////////  DRAWING PLOTS AND LEGENDS /////////////////////////////////
  TCanvas* ca = new TCanvas("ca","Trigger Efficiency") ;

  ca->SetGridx();
  ca->SetGridy();
  ca->cd();
  
  gPad->SetLogx();
  gPad->SetObjectStat(1);

  frame->GetYaxis()->SetRangeUser(0,1.05);
  frame->GetXaxis()->SetRangeUser(1,100.);
  frame->GetYaxis()->SetTitle("Efficiency");
  frame->GetXaxis()->SetTitle("E_{T} [GeV]");
  frame->Draw() ;

  frame2->GetYaxis()->SetRangeUser(0,1.05);
  frame2->GetXaxis()->SetRangeUser(1,100.);
  frame2->GetYaxis()->SetTitle("Efficiency");
  frame2->GetXaxis()->SetTitle("E_{T} [GeV]");
  frame2->Draw("same") ;

  TH1F *SCeta1 = new TH1F("SCeta1","SCeta1",50,-2.5,2.5);
  TH1F *SCeta2 = new TH1F("SCeta2","SCeta2",50,-2.5,2.5);

  SCeta1->SetLineColor(color1) ;
  SCeta1->SetMarkerColor(color1);
  SCeta1->SetMarkerStyle(style1);

  SCeta2->SetLineColor(color2) ;
  SCeta2->SetMarkerColor(color2);
  SCeta2->SetMarkerStyle(style2);

  TLegend *leg = new TLegend(0.246,0.435,0.461,0.560,NULL,"brNDC"); // mid : x=353.5
  leg->SetLineColor(1);
  leg->SetTextColor(1);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetShadowColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->SetMargin(0.25);
  TLegendEntry *entry=leg->AddEntry("NULL","L1_SingleEG"+names[iEG],"h");
//   leg->AddEntry(SCeta1,name_leg_ecal[iECAL1]+" "+name_leg_coll[iColl1],"p");
//   leg->AddEntry(SCeta2,name_leg_ecal[iECAL2]+" "+name_leg_coll[iColl2],"p");
  leg->AddEntry(SCeta1,name_leg_ecal[iECAL1],"p");
  leg->AddEntry(SCeta2,name_leg_ecal[iECAL2],"p");
  leg->Draw();

  leg = new TLegend(0.16,0.725,0.58,0.905,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03);
  leg->SetLineColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry("NULL","CMS Preliminary 2012 pp  #sqrt{s}=8 TeV","h");
  leg->AddEntry("NULL","#int L dt = "+lumi+"^{-1}","h");
  leg->AddEntry("NULL","Threshold : "+names[iEG]+" GeV","h");
  leg->Draw();

  TPaveText *pt2 = new TPaveText(0.220,0.605,0.487,0.685,"brNDC"); // mid : x=353.5                                          
  pt2->SetLineColor(1);
  pt2->SetTextColor(1);
  pt2->SetTextFont(42);
  pt2->SetTextSize(0.03);
  pt2->SetFillColor(kWhite);
  pt2->SetShadowColor(kWhite);
  pt2->AddText("L1 E/Gamma Trigger");
  pt2->AddText("Electrons from Z");
  pt2->Draw();
  
  //TString name_image="eff_EG20_2012_12fb";

  ca->Print(name_image+".cxx","cxx");
  ca->Print(name_image+".png","png");
  ca->Print(name_image+".gif","gif");
  ca->Print(name_image+".pdf","pdf");
  ca->Print(name_image+".ps","ps");
  ca->Print(name_image+".eps","eps");

  /////////////////////////////
  // SAVE THE ROO FIT RESULT //
  /////////////////////////////

  RooWorkspace *w = new RooWorkspace("workspace","workspace") ;

  w->import(dataSet);
  w->import(dataSet2);
  
  w->import(*roofitres1,"roofitres1");
  w->import(*roofitres2,"roofitres2");

  cout << "CREATES WORKSPACE : " << endl;
  w->Print();
  
  w->writeToFile(name_image+"_fitres.root") ;
  //gDirectory->Add(w) ;

  //f1->Close();
}

#endif
