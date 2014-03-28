// FIND eT VALUE FOR WHICH EFFICIENCY HAS A GIVEN VALUE

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

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

#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"

#include "FuncCB_428.h"

using namespace RooFit ;
using namespace std;


double dichotomie(double eff, double a0, double b0, double relErr,
		  FuncCB cb, bool verbose) {
  
  double dicho, effApprox, a, b;
  
  if(a0<b0) {
    a = a0;
    b = b0;
  } else if(a0>b0) {
    a = b0;
    b = a0;
  }
  else {
    cout << "PLEASE CHOOSE DIFFERENT VALUES FOR a AND b" << endl;
    return -999;
  }

  // Test bounds
  if( (cb.valeur(a) > eff) || (cb.valeur(b) < eff) ) {
    cout << "Bounds not large enough : eff(a)=" << cb.valeur(a) 
	 << " ; eff(b)=" << cb.valeur(b) << " ; tested eff=" << eff
	 << endl;
    return -999;
  }

  do {
    dicho = (a+b)/2 ;
    effApprox = cb.valeur(dicho);

    if( effApprox < eff ) {
      a = dicho;
    } else {
      b = dicho;
    }
  }
  while( (fabs(effApprox-eff) / eff) > relErr );

  if(verbose) {
    cout << "relative precision asked (" << relErr*100 << " %) reached !"
	 << endl
	 << "found value of eT : " << dicho << " GeV" 
	 << endl
      //<< "efficiency value : " << 100*efficiency(dicho,mean,sigma,alpha,n,norm) << " %"
	 << endl;
  }

  return dicho;

}

void turnon_sigma(TH1D * h, int nSig, double ref, double& errmin, double& errmax, bool verbose=false)
{
  errmin = 0. ;
  errmax = 0. ;
  double frac = 0. ;

  double nEntries = h->Integral(1,h->GetNbinsX());
  double fracLR = (1-TMath::Erf(nSig/sqrt(2)))/2. ;
  double nEntriesLR = fracLR*nEntries ;
  int binRef = h->FindBin(ref) ;
  if (binRef<=1) return ;

  // left interval
  int binL=-1, binR=-1 ;
  double nEntriesL, nEntriesR ;
  for (int bin=binRef-1 ; bin>=1 ; bin--) {
    nEntriesL = h->Integral(bin,binRef);
    if (nEntriesL >= nEntriesLR) {
      binL = bin ;
      break ;
    }
  }
  if (binL<0) return ;

  // right interval
  int lastFilledBin = 0 ;
  for (int bin=binRef ; bin<= h->GetNbinsX() ; bin++) {
    if (h->GetBinContent(bin)>0) lastFilledBin = bin ;
    nEntriesR = h->Integral(binRef,bin);
    if (nEntriesR >= nEntriesLR) {
      binR = bin ;
      break ;
    }
  }
  if (binR>h->GetNbinsX()) return ;
  if (binR<0) binR = lastFilledBin ;

  //check fraction
  frac = h->Integral(binL, binR)/nEntries ;

  errmin = h->GetBinLowEdge(binL) ;
  errmax = h->GetBinLowEdge(binR) ;

  if (verbose && errmin>ref && errmax<ref) 
    cout << binL << " " << binR << " " << ref << " " 
	 << errmin << " " << errmax << " " << frac << endl;

}

int genHisto(int nbdata, RooDataSet* data, const int nEff, double* eff, double* et, double* et_errmax, double* et_errmin,        
	     double& eff_plateau, double& eff_plateau_errmax, double& eff_plateau_errmin,                                        
	     RooRealVar xaxis, RooRealVar mean, RooRealVar sigma, RooRealVar alpha, RooRealVar n, RooRealVar norm,               
	     bool draw, bool verbose)
{

  FuncCB cb("cb","Crystal Ball Integree",xaxis,mean,sigma,alpha,n,norm) ;

  TString name_h[nEff];
  ostringstream ossi;
  for(int iEff=0;iEff<nEff;iEff++) {
    ossi.str("");
    ossi << iEff;
    name_h[iEff] = ossi.str();
  }
  
  TH1D * hEff = new TH1D("hEff","hEff",11000, 0, 1.1);
  TH1D * hEt[nEff];

  for(int iEff=0 ; iEff<nEff ; iEff++)
    hEt[iEff] = new TH1D("hEt_"+name_h[iEff], "hEt_"+name_h[iEff], 15000, 0, 150);

  double ref=0 ;
  double ref_et[nEff]; for(int iEff=0;iEff<nEff;iEff++) ref_et[iEff]=0;

  double eff_c=0;
  double et_c[nEff], et_p[nEff]; // current and previous values
  for(int iEff=0;iEff<nEff;iEff++) { et_c[iEff]=0; et_p[iEff]=0; }
 
  // Reference values
  cout << "Compute reference values : " << endl
       << "- Plateau efficiency : "; 
  ref = cb.valeur(xaxis.getVal()); 
  cout << ref << " %" << endl;

  for(int iEff=0 ; iEff<nEff ; iEff++) {
    cout << "- eT(" << eff[iEff]*100 << " %) = ";

    ref_et[iEff] = dichotomie(eff[iEff],0,10000,0.0000001,cb,verbose);
    //if(ref_et[iEff]==-999) ref_et[iEff] = ?

    cout << ref_et[iEff] << " GeV" << endl;
  // double dichotomie(double eff, double a0, double b0, double relErr, FuncCB cb, bool verbose)
  }

  // Loop over parameters dataset
  cout << "Loop over parameters dataset" << endl;
  for (int i=0 ; i<nbdata ; i++) {
    if (i%1000==0 && verbose) cout << "=== " << i << "th iteration" << endl;

    const RooArgSet* set = data->get(i) ;
    alpha = set->getRealValue("alpha") ;
    mean = set->getRealValue("mean") ;
    n = set->getRealValue("n") ;
    norm = set->getRealValue("norm") ;
    sigma = set->getRealValue("sigma") ;
    
    if (n.getVal()<=0) { cout << "meaningless n=" << n.getVal() << endl; continue; } 
    // avoid meaningless parameters

    eff_c = cb.valeur(xaxis.getVal());
    hEff->Fill(eff_c);

    cout << "-- iteration : " << i << " : dichotomie" << endl;
    for(int iEff=0 ; iEff<nEff ; iEff++) {
      cout << "--- " << eff[iEff]*100 << " % : ";

      et_c[iEff] = dichotomie(eff[iEff],0,10000,0.0000001,cb,verbose);
      if(et_c[iEff]==-999) et_c[iEff]=et_p[iEff]; // use previous value

      cout << et_c[iEff] << " GeV" << endl;
      hEt[iEff]->Fill(et_c[iEff]);
      
      et_p[iEff]=et_c[iEff]; // save current values for next iteration
    }
  }


  // Extract error bars from histograms //
  cout << "Extract error bars from histograms" << endl;
  double errmin, errmax;
  errmin = errmax = 0;

  // Plateau
  turnon_sigma(hEff, 1, ref, errmin, errmax, verbose);
  eff_plateau = ref;                // eff                                                                         
  eff_plateau_errmax = errmax-ref ; // + err0                      
  eff_plateau_errmin = ref-errmin ; // - err1  

  // et(efficiency)
  for(int iEff=0 ; iEff<nEff ; iEff++) {
    turnon_sigma(hEt[iEff], 1, ref_et[iEff], errmin, errmax, verbose);

    et[iEff] = ref_et[iEff];
    et_errmax[iEff] = errmax - ref_et[iEff] ;
    et_errmin[iEff] = ref_et[iEff] - errmin ;
  }
  /*
  if(verbose)
    cout << "At x=" << xaxis.getVal() << ", eff=" << ref 
	 <<" +" << errmax-ref << " -" << ref-errmin 
	 << " and fraction=" << frac << endl;
  */
  // if(draw) hEff->Draw();

  return 1;

}
int compute(RooFitResult* fit, int nbdata, const int nEff, double* eff, double* et, double* et_errmax, double* et_errmin,
	    double et_plateau, double& eff_plateau, double& eff_plateau_errmax, double& eff_plateau_errmin,                       
	    bool draw, bool verbose)                                                                                              
{
  
  // Extract fit parameters //

  RooArgList param = fit->floatParsFinal() ;
  double err[5] ;
  double mu[5] ;

  for(Int_t i = 0; i < param.getSize(); i++) {

    RooRealVar* var = ( dynamic_cast<RooRealVar*>( param.at(i) ) );
    //var->Print() ;
    mu[i] = var->getVal() ;
    err[i] = var->getError() ;
  }
  
  double min, max;
  
  min = mu[0]-5*err[0] ;
  if (mu[0]-5*err[0]<0) min = 0. ;
  RooRealVar alpha("alpha","#alpha",mu[0],min,mu[0]+5*err[0]);
  
  min = mu[1]-5*err[1] ;
  if (mu[1]-5*err[1]<5) min = 5. ;
  RooRealVar mean("mean","mean",mu[1],min,mu[1]+5*err[1]);
  
  min = mu[2]-5*err[2] ;
  if (mu[2]-5*err[2]<1) min = 1. ;
  RooRealVar n("n","n",mu[2],min,mu[2]+5*err[2]);
  
  min = mu[3]-5*err[3] ;
  if (mu[3]-5*err[3]<0.6) min = 0.6 ; 

  max = mu[3]+5*err[3] ;
  if (mu[3]+5*err[3]>1.) max = 1. ; 
  RooRealVar norm("norm","N",mu[3],min,max);
  
  min = mu[4]-5*err[4] ;
  if (mu[4]-5*err[4]<0.) min = 0. ; 
  RooRealVar sigma("sigma","#sigma",mu[4],min,mu[4]+5*err[4]);
  
  RooRealVar xaxis("x","x",0,150) ;


  // Create PDF and generate nbdata sets of CB parameters
  RooAbsPdf* parabPdf = fit->createHessePdf(RooArgSet(norm,alpha,n,mean,sigma)) ;
  RooDataSet* data = parabPdf->generate(RooArgSet(norm,alpha,n,mean,sigma),nbdata) ;

  // Generate histo to extract error bar on efficiency(xaxis)
  xaxis = et_plateau ;
  cout << "Generate histo to extract error bars" << endl;
  genHisto(nbdata, data, nEff, eff, et, et_errmax, et_errmin, eff_plateau, eff_plateau_errmax, eff_plateau_errmin,
	   xaxis, mean, sigma, alpha, n, norm, draw, verbose);

  /* int genHisto(int nbdata, RooDataSet* data, const int nEff, double* eff, double* et, double* et_errmax, double* et_errmin,
                  double& eff_plateau, double& eff_plateau_errmax, double& eff_plateau_errmin, 
		  RooRealVar xaxis, RooRealVar mean, RooRealVar sigma, RooRealVar alpha, RooRealVar n, RooRealVar norm,
		  bool draw, bool verbose)

     int compute(RooFitResult* fit, int nbdata, const int nEff, double* eff, double* et, double* et_errmax, double* et_errmin,
                 double et_plateau, double& eff_plateau, double& eff_plateau_errmax, double& eff_plateau_errmin,
                 bool draw, bool verbose)
  */

  return 1;
}

int prodResults(TString dir, TString file, TString namews="workspace", int nbdata=10000, double et_plateau=100, 
		bool draw=false, bool verbose=false) 
{

  // Define efficiency/et values and errors
  const int nEff=9; // # valeurs d'efficacite
  const int nE=2;   // EB / EE

  double eff[nEff] = {0.25, 0.33, 0.5, 0.66, 0.75, 0.9, 0.95, 0.98, 0.99}; // efficacites a dichotomiser
  double et[nE][nEff];        // valeurs eT(eff)
  double et_errmax[nE][nEff]; // erreur sup sur eT(eff)
  double et_errmin[nE][nEff]; // erreur inf sur eT(eff)

  double eff_plateau[nE]={0,0};        // eff(et_plateau)
  double eff_plateau_errmax[nE]={0,0}; // erreur sup
  double eff_plateau_errmin[nE]={0,0}; // erreur sup

  for(int iE=0 ; iE<nE ; iE++) {
    for(int iEff=0 ; iEff<nEff ; iEff++) {
      et[iE][iEff]=0;
      et_errmax[iE][iEff]=0;
      et_errmin[iE][iEff]=0;
    }
  }

  // Open file and retrieve workspace
  cout << "Open file and retrieve workspace" << endl
       << "dir : "  << dir << endl
       << "file : " << file << endl;

  TFile *f = new TFile(dir+"/"file+".root") ;
  RooWorkspace* w = (RooWorkspace*) f->Get(namews) ;

  // Retrieve fit results
  TString namefit[nE]= {"roofitres1", "roofitres2"};
  RooFitResult* fitres[nE];

  for(int iE=0 ; iE<nE ; iE++)
    fitres[iE]= (RooFitResult*) w->obj(namefit[iE]) ;
  

  // Compute efficiency @ 100 GeV and et values (with error bars)
  for(int iE=0 ; iE<nE ; iE++) {
    cout << "Compute efficiency @ 100 GeV and et values (with error bars) : #" << iE << endl;
    compute(fitres[iE], nbdata, nEff, eff, et[iE], et_errmax[iE], et_errmin[iE],
	    et_plateau, eff_plateau[iE], eff_plateau_errmax[iE], eff_plateau_errmin[iE], draw, verbose);
  }
  /*
    int compute(RooFitResult* fit, int nbdata, const int nEff, double* eff, 
                double* et, double* et_errmax, double* et_errmin, 
		double et_plateau, double& eff_plateau, double& eff_plateau_errmax, double& eff_plateau_errmin, 
		bool draw, bool verbose)
  */

  // Output results
  TString nameE[nE] = {"BARREL","ENDCAPS"};
  for(int iE=0 ; iE<nE ; iE++) {
    cout << "<------ "+nameE[iE]+" RESULTS ------>" << endl << endl;

    cout << "Fit Results : " << endl;
    fitres[iE]->Print();
    cout << endl << endl;

    cout << "Transverse energy :" << endl;
    for(int iEff=0 ; iEff<nEff ; iEff++)
      cout << 100*eff[iEff] << " % : " << et[iE][iEff] 
	   << " + " << et_errmax[iE][iEff] 
	   << " - " << et_errmin[iE][iEff] 
	   << " GeV" << endl;
    cout << endl;

    cout << "Plateau @ " << et_plateau << " GeV : " << eff_plateau[iE] 
	 << " + " << eff_plateau_errmax[iE]
	 << " _ " << eff_plateau_errmin[iE]
	 << " %" << endl;

    cout << "__________________________________________________________________________" << endl;
  }

  return 1;
}
