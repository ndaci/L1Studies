#ifndef DEF_FASTEFF_H
#define DEF_FASTEFF_H

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
#include "TChain.h"

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

  //gROOT->ProcessLine(".L /home/llr/cms/ndaci/SKWork/FIT/turnonfit/tdrstyle.C");
  //gROOT->ProcessLine("setTDRStyle()");
  setTDRStyle();
}

#endif
