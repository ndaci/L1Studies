/////////////////////////////////////////////////////
// FITTING WITH A ROOFIT-USER DEFINED CRYSTAL BALL //
/////////////////////////////////////////////////////
#ifndef DEF_FASTEFF
#define DEF_FASTEFF

#include "fitEfficiency.h"

void fitEfficiency(float thres,
		   float mLow=45, float mHigh=70,
		   TString cutdef="pt_HLT>=20", TString cut_choice="HLT",
		   TString dirIn="/data_CMS/cms/ndaci/ndaci_2012/HTauTau/TriggerStudy/SingleMu/MuMu/Run2012A_PRV1/Pairs/", 
		   TString dirResults="/home//llr/cms/ndaci/SKWork/macro/HTauTau/results/TriggerStudies/MuMu/Run2012A_PRV1/",
		   TString lumi="200 pb", int nCPU=4, 
		   int color1=kBlack, int style1=kFullCircle, int color2=kRed, int style2=kOpenSquare,
		   TString fileIn="*.root", TString image="eff_HLT_MuTau")
{
  // STYLE //
  gROOT->Reset();
  loadPresentationStyle();  
  gROOT->ForceStyle();

  // OUTPUT //
  TString name_image = dirResults+"/"+image;
  ofstream fichier(name_image+".txt", ios::out);

  // BINNING //
  const int nbinsEB = 14;
  Double_t binsEB[nbinsEB] = {10., 14., 18., 20., 22., 26., 30., 40., 50., 60., 70., 80, 90, 100};
  RooBinning binningEB = RooBinning(nbinsEB-1, binsEB, "binningEB");

  const int nbinsEE = 14;
  Double_t binsEE[nbinsEE] = {10., 14., 18., 20., 22., 26., 30., 40., 50., 60., 70., 80, 90, 100};
  RooBinning binningEE = RooBinning(nbinsEE-1, binsEE, "binningEE");

  // EB/EE eta cuts //
  TString etacutEB="eta>-1.479 && eta<1.479", etacutEE="eta<=-1.479 || eta>=1.479";
  TString cutdefEB, cutdefEE;
  if(cutdef=="") { 
    cutdefEB=etacutEB;
    cutdefEE=etacutEE;
  }
  else {
    cutdefEB=cutdef+" && "+etacutEB ;
    cutdefEB=cutdef+" && "+etacutEE ;
  }

  // INPUT DATA //
  TChain* treeTnP = new TChain("treeTnP");
  treeTnP->Add(dirIn+"/"+fileIn);

  RooRealVar xaxis("pt","P_{T} [GeV]",0,150) ;
  RooRealVar mass("mass","mass",mLow,mHigh) ;// consider only this mass range when importing data                                      
  RooRealVar eta("eta","eta",-3., 3) ;
  RooRealVar weight("weight","weight",-1,1000) ;
  RooRealVar pt_HLT_tag("pt_HLT_tag","P_{T} [GeV]",0,150) ;
  RooRealVar pt_HLT_tag_sanity("pt_HLT_tag_sanity","P_{T} [GeV]",0,150) ;
  RooRealVar pt_L3("pt_L3","P_{T} [GeV]",0,150) ;
  RooRealVar pt_L25("pt_L25","P_{T} [GeV]",0,150) ;
  RooRealVar pt_L2("pt_L2","P_{T} [GeV]",0,150) ;
  RooRealVar et_L1jet("et_L1_jet","P_{T} [GeV]",0,150) ;
  RooRealVar et_L1tau("et_L1_tau","P_{T} [GeV]",0,150) ;

  // DEFINE EFFICIENCY CUT //
  RooCategory cutHLT("match","") ;
  cutHLT.defineType("accept",1) ;
  cutHLT.defineType("reject",0) ;

  RooCategory cutL1("L1match","") ;
  cutL1.defineType("accept",1) ;
  cutL1.defineType("reject",0) ;

  RooCategory cutL1L2("L1L2match","") ;
  cutL1L2.defineType("accept",1) ;
  cutL1L2.defineType("reject",0) ;

  RooCategory cutL1L2L25("L1L2L25match","") ;
  cutL1L2L25.defineType("accept",1) ;
  cutL1L2L25.defineType("reject",0) ;

  RooCategory cutL1L2L25L3("L1L2L25L3match","") ;
  cutL1L2L25L3.defineType("accept",1) ;
  cutL1L2L25L3.defineType("reject",0) ;

  RooCategory *cut;
  if(cut_choice=="HLT") cut = &cutHLT;
  if(cut_choice=="L1") cut = &cutL1;
  if(cut_choice=="L1L2") cut = &cutL1L2;
  if(cut_choice=="L1L2L25") cut = &cutL1L2L25;
  if(cut_choice=="L1L2L25L3") cut = &cutL1L2L25L3;

  RooDataSet dataSetEB("data","data from tree",
		       RooArgSet(xaxis, *cut,
				 pt_HLT_tag, pt_L3, pt_L25, pt_L2,
				 mass, weight, eta),
		       WeightVar(weight), Import(*treeTnP), Cut(cutdefEB) );
  
  RooDataSet dataSetEE("data","data from tree",
		       RooArgSet(xaxis, *cut,
				 pt_HLT_tag, pt_L3, pt_L25, pt_L2,
				 mass, weight, eta),
		       WeightVar(weight), Import(*treeTnP), Cut(cutdefEE) );
  
  dataSetEB.Print();
  dataSetEE.Print();

  // PLOT //
  RooPlot* frame = xaxis.frame(Bins(18000),Title("Fitted efficiency")) ;
  dataSetEB.plotOn(frame, Binning(binningEB), Efficiency(*cut), MarkerColor(color1), LineColor(color1), MarkerStyle(style1) );
  dataSetEE.plotOn(frame, Binning(binningEE), Efficiency(*cut), MarkerColor(color2), LineColor(color2), MarkerStyle(style2) );

  // PARAMETRES ROOFIT CRYSTAL BALL //
  RooRealVar norm("norm","N",0.95,0.6,1);
  RooRealVar alpha("alpha","#alpha",0.2,0.01,8);
  RooRealVar n("n","n",2,1.1,35);
  RooRealVar mean("mean","mean",10,5,30);
  mean.setVal(thres);
  RooRealVar sigma("sigma","#sigma",0.23,0.01,5);
  RooRealVar pedestal("pedestal","pedestal",0.01,0,0.4);

  FuncCB cb_EB("cb_EB","Fit function EB (cb)",xaxis,mean,sigma,alpha,n,norm) ;
  FuncCB cb_EE("cb_EE","Fit function EE (cb)",xaxis,mean,sigma,alpha,n,norm) ;
  
  RooEfficiency eff_EB("eff_EB","efficiency EB", cb_EB, *cut, "accept");
  RooEfficiency eff_EE("eff_EE","efficiency EE", cb_EE, *cut, "accept");  

  // FIT //
  double fit_cuts_min = thres-1.5 ;
  double fit_cuts_max = 150;
  xaxis.setRange("interesting",fit_cuts_min,fit_cuts_max);

  fichier << "Fit characteristics :" << endl
	  << "Threshold : "          << thres << endl 
	  << "Fit Range : ["         << fit_cuts_min << "," << fit_cuts_max << "]" << endl 
	  << endl << endl;

  RooFitResult* roofitres_EB = new RooFitResult("roofitres_EB","roofitres_EB");
  RooFitResult* roofitres_EE = new RooFitResult("roofitres_EE","roofitres_EE");

  // Fit #1 //
  roofitres_EB = eff_EB.fitTo(dataSetEB,ConditionalObservables(xaxis),Range("interesting"),Minos(kTRUE),Warnings(kFALSE),NumCPU(nCPU),Save(kTRUE),SumW2Error(kTRUE));
  cb_EB.plotOn(frame,LineColor(color1),LineWidth(2));

  fichier << "<----------------- EB ----------------->" << endl
	  << "double res_mean="  << mean.getVal()   << "; "
	  << "double res_sigma=" << sigma.getVal()  << "; "
          << "double res_alpha=" << alpha.getVal()  << "; "
          << "double res_n="     << n.getVal()      << "; "
          << "double res_norm="  << norm.getVal()   << "; "
	  << endl
	  << "double err_mean="  << mean.getError()  << "; "
	  << "double err_sigma=" << sigma.getError() << "; "
          << "double err_alpha=" << alpha.getError() << "; "
          << "double err_n="     << n.getError()     << "; "
          << "double err_norm="  << norm.getErrorLo()<< "; "
	  << endl;

  // Fit #2 //
  roofitres_EE = eff_EE.fitTo(dataSetEE,ConditionalObservables(xaxis),Range("interesting"),Minos(kTRUE),Warnings(kFALSE),NumCPU(nCPU),Save(kTRUE),SumW2Error(kTRUE));
  cb_EE.plotOn(frame,LineColor(color2),LineWidth(2));

  fichier << "<----------------- EE ----------------->" << endl
	  << "double res_mean="  << mean.getVal()   << "; "
	  << "double res_sigma=" << sigma.getVal()  << "; "
          << "double res_alpha=" << alpha.getVal()  << "; "
          << "double res_n="     << n.getVal()      << "; "
          << "double res_norm="  << norm.getVal()   << "; "
	  << endl
	  << "double err_mean="  << mean.getError()  << "; "
	  << "double err_sigma=" << sigma.getError() << "; "
          << "double err_alpha=" << alpha.getError() << "; "
          << "double err_n="     << n.getError()     << "; "
          << "double err_norm="  << norm.getErrorLo()<< "; "
	  << endl;


  ////////////////////////////  DRAWING PLOTS AND LEGENDS /////////////////////////////////
  TCanvas* ca = new TCanvas("ca","Trigger Efficiency") ;

  ca->SetGridx();
  ca->SetGridy();
  ca->cd();
  
  //gPad->SetLogx();
  gPad->SetObjectStat(1);

  frame->GetYaxis()->SetRangeUser(0,1.05);
  frame->GetXaxis()->SetRangeUser(1,150.);
  frame->GetYaxis()->SetTitle("Efficiency");
  frame->GetXaxis()->SetTitle("E_{T} [GeV]");
  frame->Draw() ;

  TH1F *SCeta_EB = new TH1F("SCeta_EB","SCeta_EB",50,-2.5,2.5);
  TH1F *SCeta_EE = new TH1F("SCeta_EE","SCeta_EE",50,-2.5,2.5);

  SCeta_EB->SetLineColor(color1) ;
  SCeta_EB->SetMarkerColor(color1);
  SCeta_EB->SetMarkerStyle(style1);

  SCeta_EE->SetLineColor(color2) ;
  SCeta_EE->SetMarkerColor(color2);
  SCeta_EE->SetMarkerStyle(style2);

  TLegend *leg = new TLegend(0.40, 0.14, 0.63, 0.34, NULL, "brNDC");
  leg->SetLineColor(1);
  leg->SetTextColor(1);
  leg->SetTextFont(42);
  leg->SetTextSize(0.0244755);
  leg->SetShadowColor(kWhite);
  leg->SetFillColor(kWhite);  
  leg->AddEntry("NULL","e #tau electrons efficiency","h");
  //    entry->SetLineColor(1);
  //    entry->SetLineStyle(1);
  //    entry->SetLineWidth(1);
  //    entry->SetMarkerColor(1);
  //    entry->SetMarkerStyle(21);
  //    entry->SetMarkerSize(1);
  //    entry->SetTextFont(62);
  leg->AddEntry(SCeta_EB,"Barrel","p");
  leg->AddEntry(SCeta_EE,"Endcaps","p");
  leg->Draw(); 
  
  ostringstream ossi("");
  ossi << thres;
  TString tossi = ossi.str();

  leg = new TLegend(0.40, 0.30, 0.50, 0.50, NULL, "brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.0297203);
  leg->SetLineColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry("NULL","CMS Preliminary 2012 pp  #sqrt{s}=8 TeV","h");
  leg->AddEntry("NULL","#int L dt = "+lumi+"^{-1}","h");
  leg->AddEntry("NULL","Threshold : "+tossi+" GeV","h");
  leg->Draw();
  
  ca->Print(name_image+".C","C");
  ca->Print(name_image+".cxx","cxx");
  ca->Print(name_image+".png","png");
  ca->Print(name_image+".gif","gif");
  ca->Print(name_image+".pdf","pdf");
  ca->Print(name_image+".ps","ps");

  /////////////////////////////
  // SAVE THE ROO FIT RESULT //
  /////////////////////////////

  RooWorkspace *w = new RooWorkspace("workspace","workspace") ;

  w->import(dataSetEB);
  w->import(dataSetEE);
  
  w->import(*roofitres_EB,"roofitres_EB");
  w->import(*roofitres_EE,"roofitres_EE");

  cout << "CREATES WORKSPACE : " << endl;
  w->Print();
  
  w->writeToFile(name_image+"_fitres.root") ;
  //gDirectory->Add(w) ;

  // DELETE POINTERS
  // int a=0;
//   cin >> a;
//   delete treeTnP; delete cut; delete frame; delete roofitres_EB; delete roofitres_EE; delete ca; delete SCeta_EB; delete SCeta_EE; delete leg; delete w;

}

#endif
