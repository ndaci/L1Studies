#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveStats.h"
//#include "treesList.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TFrame.h"

using namespace std;

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

  gROOT->ProcessLine(".L /home/llr/cms/ndaci/SKWork/FIT/turnonfit/tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
}

float errPoissCount(float n, float ntot) {
  
  if(ntot==0) {
    cout << "errPoissCount ERROR : ntot=0" << endl;
    return 0;
  }

  float x = n*(1-n/ntot);
  
  if( x<0 ) {
    cout << "errPoissCount ERROR : sqrt(<0)" << endl;
    return 0;
  }

  return (1/ntot)*sqrt(n*(1-n/ntot)) ;
}

float errPoissCount_Alex(float n, float ntot) {

  if(ntot==0) {
    cout << "errPoissCount ERROR : ntot=0" << endl;
    return 0;
  }

  float x = n*(1-n)/ntot ;
  
  if( x<0 ) {
    cout << "errPoissCount ERROR : sqrt(<0)" <<endl;
    return 0;
  }

  return sqrt(x);
}

void spike_summarize(TString dirOut="./HPU2012/2specRuns/", int factor=3, TString graph_style="A*") {
  
  loadPresentationStyle();
  /*
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetCanvasDefH(600); //Height of canvas                                                     
    gStyle->SetCanvasDefW(600); //Width of canvas                                                                   
    gStyle->SetCanvasDefX(0);   //POsition on screen                                                                 
    gStyle->SetCanvasDefY(0);

    gStyle->SetTitleFont(42);
    gStyle->SetTitleColor(1);
    gStyle->SetTitleTextColor(1);
  */
  //gStyle->SetTitleFillColor(10);
  /*
    gStyle->SetTitleFontSize(0.05);
    gStyle->SetTitleColor(1, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
  */
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetTitleYOffset(1.2);
  /*
    gStyle->SetLabelColor(1, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetLabelOffset(0.007, "XYZ");
    gStyle->SetLabelSize(0.05, "XYZ");

    gStyle->SetAxisColor(1, "XYZ");
    gStyle->SetStripDecimals(kTRUE);
    gStyle->SetTickLength(0.03, "XYZ");
    gStyle->SetNdivisions(510, "XYZ");
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame                                                  
    gStyle->SetPadTickY(1);
  */
  //const int nSrc = 3; // 2011A, 2011B, highPU
  const int nSrc=6; // data 2011B, newkill 2011B, data highPU, emul highPU, 2012 online HPU run1, run2
  const int nStrict = 2;
  const int nEG = 8;
  const int nVtx = 40;
  const int nVTX[nSrc]={40,40,40,40,60};
  const int nEG_g=4; // nEG for graph : EG12,15,20,30

  //TString src_name[nSrc]={"2011A","2011B","highPU"};
  TString src_name[nSrc]={"data_2011B","newsetting_2011B", "data_highPU", "d2012_HPU_1", "d2012_HPU_2"};
  TString strict_name[nStrict] = {"large","strict"};

  TString trigname[nEG]={"2","5","8","10","12","15","20","30"};
  int trigthresh[nEG]={2,5,8,10,12,15,20,30};

  TH1F * h_evts_trigBy_any[nEG][nSrc];
  TH1F * h_evts_trigBy_spikes[nEG][nStrict][nSrc];

  cout << "objects declared and counters initialized" << endl;

  TString name_histo;

  for( int iSrc=0 ; iSrc<nSrc ; iSrc++ ) {
    for( int iEG=0 ; iEG<nEG ; iEG++ ) {
      name_histo = "h_evts_trigBy_any_EG"+trigname[iEG]+"_"+src_name[iSrc];
      h_evts_trigBy_any[iEG][iSrc] = new TH1F(name_histo,name_histo,nVTX[iSrc],0,nVTX[iSrc]);

      for( int iStrict=0 ; iStrict<nStrict ; iStrict++ ) {
	name_histo = "h_evts_trigBy_spikes_EG"+trigname[iEG]+"_"+strict_name[iStrict]+"_"+src_name[iSrc];
	h_evts_trigBy_spikes[iEG][iStrict][iSrc] = new TH1F(name_histo,name_histo,nVTX[iSrc],0,nVTX[iSrc]);
      }
    }
  }

  TString file[nSrc];

  // FILES : HLTEG12
  file[0]="/home/llr/cms/ndaci/SKWork/macro/skEfficiency/tagAndProbe/Spike2011A/Commi_2011B_HLTEG12_fromfilter/spike_plots_half_0_1.root";

  file[1]="/home/llr/cms/ndaci/SKWork/macro/skEfficiency/tagAndProbe/Spike2011A/CommiEmulReco_Run2011B/GetContamNewkill/spike_plots_half_0_1.root";

  file[2]="/home/llr/cms/ndaci/SKWork/macro/skEfficiency/tagAndProbe/Spike2011A/Commi_2011B_highPU_HLTEG_repro/spike_plots_half_0_1.root" ;

  file[3]="/home/llr/cms/ndaci/SKWork/macro/skEfficiency/Spikes/HighPuContamination/try2/spike_plots.root" ;

  //file[4]="/home/llr/cms/ndaci/SKWork/macro/skEfficiency/Spikes2012/HPU/spike_plots_half_0_1.root" ;

  file[4]="/home/llr/cms/ndaci/SKWork/macro/skEfficiency/Spikes2012/HPU/Run198603/spike_plots_half_0_1.root" ;
  file[5]="/home/llr/cms/ndaci/SKWork/macro/skEfficiency/Spikes2012/HPU/Run198609/spike_plots_half_0_1.root" ;

  TFile * f;

  cout << "going to open files" << endl;

  for(int iSrc=0 ; iSrc<nSrc ; iSrc++) {

    if(iSrc==3) continue;

    f = TFile::Open(file[iSrc]);
    cout << "file #" << iSrc << " opened" << endl;
    for( int iEG=0 ; iEG<nEG ; iEG++ ) {
      name_histo = "h_evts_trigBy_any_EG"+trigname[iEG];
      h_evts_trigBy_any[iEG][iSrc] -> Add( (TH1F*)gDirectory->Get(name_histo) );
  
      if(iSrc==4) h_evts_trigBy_any[iEG][iSrc] -> Rebin(2);

      cout << " got histogram " + name_histo << endl;

      for( int iStrict=0 ; iStrict<nStrict ; iStrict++ ) {
        name_histo = "h_evts_trigBy_spikes_EG"+trigname[iEG]+"_"+strict_name[iStrict];
        h_evts_trigBy_spikes[iEG][iStrict][iSrc] -> Add( (TH1F*)gDirectory->Get(name_histo) );

	if(iSrc==4) h_evts_trigBy_spikes[iEG][iStrict][iSrc] -> Rebin(2);

	cout << "got histogram "+ name_histo <<endl;
      }
    }
    f->Close();
  } 
  
  // highPU emul is special case
  f = TFile::Open(file[3]);
  cout << "file #" << 3 << " opened" << endl;
  for( int iEG=0 ; iEG<nEG ; iEG++ ) {
    name_histo = "h_evts_trigBy_any_EG"+trigname[iEG]+"_M";
    h_evts_trigBy_any[iEG][3] -> Add( (TH1F*)gDirectory->Get(name_histo) );
    //
    name_histo = "h_evts_trigBy_spikes_EG"+trigname[iEG]+"_M";
    h_evts_trigBy_spikes[iEG][0][3] -> Add( (TH1F*)gDirectory->Get(name_histo) );
    h_evts_trigBy_spikes[iEG][1][3] -> Add( (TH1F*)gDirectory->Get(name_histo) );
  }

  const int nFr=2; // low/up frontier of iBin
  const int nBins_src[nSrc][nEG_g]={ {21,15,19,16} , {16,15,13,10} , {28,25,24,21} , {28,25,24,21} , {1,1,1,1} , {1,1,1,1}};

  int idxBinsFr[nSrc][nEG_g][nFr]={ { {2,22},{1,15},{2,20},{3,18} },
				    { {3,18},{1,15},{4,16},{5,14} }, 
				    { {12,39},{14,38},{14,37},{14,34} }, 
				    { {12,39},{14,38},{14,37},{14,34} },
				    { {44,44},{44,44},{44,44},{44,44} },
				    { {66,66},{66,66},{66,66},{66,66} } };


  int fact[nSrc]={1,1,1,1,factor,factor};

  vector< float > n_evts_trigBy_any[nEG][nSrc]; // double tableau de vecteurs de floats
  vector< float > n_evts_trigBy_spikes[nEG][nStrict][nSrc]; // triple 
  vector< float > n_contam[nEG][nStrict][nSrc]; 
  vector< float > err_contam[nEG][nStrict][nSrc];

  // --> chaque element du tableau est un vecteur contenant les valeurs des bins des histos

  for(int iEG=0 ; iEG<nEG_g ; iEG++) {
    for(int iSrc=0 ; iSrc<nSrc ; iSrc++) {

      const int nBins = nBins_src[iSrc][iEG] ;
	       
      for(int iBin=0 ; iBin<nBins ; iBin++) {
	n_evts_trigBy_any[iEG][iSrc].push_back( 
					       (float)( h_evts_trigBy_any[iEG+4][iSrc] -> GetBinContent(idxBinsFr[iSrc][iEG][0]/fact[iSrc]+iBin+1) ) );
	for(int iStrict=0 ; iStrict<nStrict ; iStrict++ ) {
	  n_evts_trigBy_spikes[iEG][iStrict][iSrc].push_back( 
							     (float)( h_evts_trigBy_spikes[iEG+4][iStrict][iSrc] 
								      -> GetBinContent(idxBinsFr[iSrc][iEG][0]/fact[iSrc]+iBin+1) ) );

	  if( n_evts_trigBy_any[iEG][iSrc][iBin] !=0 ) {
	    n_contam[iEG][iStrict][iSrc].push_back(
						   n_evts_trigBy_spikes[iEG][iStrict][iSrc][iBin]/n_evts_trigBy_any[iEG][iSrc][iBin]);
	    err_contam[iEG][iStrict][iSrc].push_back( 
						     errPoissCount( n_evts_trigBy_spikes[iEG][iStrict][iSrc][iBin], n_evts_trigBy_any[iEG][iSrc][iBin]) );
	  }

	  else {
	    n_contam[iEG][iStrict][iSrc].push_back(-0.01);
	    err_contam[iEG][iStrict][iSrc].push_back(0);
	  }
	} // loop over iStrict
      } // loop over iBin
    } // loop over iSrc
  } // loop over iEG


  /////////////////////////////////////////////////////////
  // GET MEAN VALUE FOR 2 SPECIAL RUNS 2012              //
  const int nSpecRuns=2;
  const int numDen=3; // 0:numerator ; 1:denominator ; 2:ratio
  float meanValSpecRuns[nSpecRuns][nEG_g][nStrict][numDen];
  float errValSpecRuns[nSpecRuns][nEG_g][nStrict];

  for(int iS=0 ; iS<nStrict ; iS++) {
    for(int iR=0 ; iR<nSpecRuns ; iR++) { 
      for(int iEG=0 ; iEG<nEG_g ; iEG++) {
	meanValSpecRuns[iR][iEG][iS][0] = h_evts_trigBy_spikes[iEG+4][iR+4][iS]->GetMean();
	meanValSpecRuns[iR][iEG][iS][1] = h_evts_trigBy_any[iEG+4][iR+4]->GetMean();

	meanValSpecRuns[iR][iEG][iS][1] != 0  ? 
	  meanValSpecRuns[iR][iEG][iS][2] = meanValSpecRuns[iR][iEG][iS][0] / meanValSpecRuns[iR][iEG][iS][1] : -1.0 ;

	errValSpecRuns[iR][iEG][iS] = errPoissCount( meanValSpecRuns[iR][iEG][iS][0] , meanValSpecRuns[iR][iEG][iS][1] );

      }
    }
  }
  /////////////////////////////////////////////////////////


  for(int iEG=0 ; iEG<nEG_g ; iEG++) {
    for(int iSrc=0 ; iSrc<nSrc ; iSrc++) {
      int nBins = n_evts_trigBy_any[iEG][iSrc].size();
      cout << "EG" << trigname[iEG+4] << " src=" << src_name[iSrc] << endl;
      for(int iBin=0 ; iBin<nBins ; iBin++)
	cout << n_evts_trigBy_any[iEG][iSrc][iBin] << "   " ;
      cout << endl;
      
      for(int iStrict=0 ; iStrict<nStrict ; iStrict++ ) {
	cout << "strict="<< strict_name[iStrict]<< endl;

	int nBins_s = n_evts_trigBy_spikes[iEG][iStrict][iSrc].size();
	for(int iBin=0 ; iBin<nBins_s ; iBin++)
	  cout << n_evts_trigBy_spikes[iEG][iStrict][iSrc][iBin] << "   ";
	cout << endl;

	nBins_s = n_contam[iEG][iStrict][iSrc].size();
	for(int iBin=0 ; iBin<nBins_s ; iBin++)
	  cout << n_contam[iEG][iStrict][iSrc][iBin] << "   ";
	cout << endl;
      }
      cout << endl;
    }
    cout << endl;
  }

  //int mark_style[nSrc] = {22,20,21}; // kFullTriangleUp, kFullCircle, kFullSquare
  //int mark_col[nSrc] = {4,2,1}; // bleu, rouge, noir

  int mark_style[nSrc] = {kFullTriangleUp,kFullTriangleUp,kFullSquare,kFullSquare,kFullCircle,kFullCircle};
  int mark_col[nSrc] = {kRed,kGreen+2,kRed,kGreen+2,kBlue,kViolet};

  TGraph * g_evt_contam[nEG][nStrict][nSrc];
  TGraphErrors * g_evt_contam_err[nEG][nStrict][nSrc];
				    
  for(int iSrc=0 ; iSrc<nSrc ; iSrc++) {
    for(int iStrict=0 ; iStrict<nStrict ; iStrict++) {
      for(int iEG=0 ; iEG<nEG_g ; iEG++) {

	const int nPts = nBins_src[iSrc][iEG];
	float x_vtx[nPts];
	float err_vtx[nPts];
	for(int iVtx=0 ; iVtx<nPts ; iVtx++) {
	  x_vtx[iVtx] = (float)( idxBinsFr[iSrc][iEG][0] + iVtx*fact[iSrc] ) ;
	  cout << x_vtx[iVtx] << " | " ;
	  err_vtx[iVtx] = 0.49*fact[iSrc];
	}
	cout << endl;

	const int nDots = n_contam[iEG][iStrict][iSrc].size();
	if( nDots != nPts ) {
	  cout << "ERROR : nDots != nPts ! exit..." << endl;
	  return;
	}

	float contam_vtx[nPts];
	float err_contam_vtx[nPts];

	for(int iVtx=0 ; iVtx<nPts ; iVtx++) {
	  contam_vtx[iVtx] = n_contam[iEG][iStrict][iSrc][iVtx];
	  err_contam_vtx[iVtx] = err_contam[iEG][iStrict][iSrc][iVtx];

	  if(iSrc>3) {
	    contam_vtx[iVtx]     = meanValSpecRuns[iSrc-4][iEG][iStrict][2];
	    err_contam_vtx[iVtx] = errValSpecRuns[iSrc-4][iEG][iStrict];
	  }
	  cout << n_contam[iEG][iStrict][iSrc][iVtx] << " |||" ;
	  cout << contam_vtx[iVtx] << " || ";

	}
	cout << endl;

        g_evt_contam[iEG][iStrict][iSrc] = new TGraph( nPts, x_vtx, contam_vtx );
	g_evt_contam_err[iEG][iStrict][iSrc] = new TGraphErrors( nPts, x_vtx, contam_vtx, err_vtx, err_contam_vtx );

        g_evt_contam[iEG][iStrict][iSrc]->SetName(
						  "Event_"+strict_name[iStrict]+"_contamination_EG"+trigname[iEG+4]+"_"+src_name[iSrc] );
        g_evt_contam[iEG][iStrict][iSrc]->SetTitle(
						   "Event "+strict_name[iStrict]+" contamination (EG"+trigname[iEG+4]+") "+src_name[iSrc] );
	g_evt_contam[iEG][iStrict][iSrc]->SetMarkerStyle(mark_style[iSrc]);
	g_evt_contam[iEG][iStrict][iSrc]->SetMarkerColor(mark_col[iSrc]);
	g_evt_contam[iEG][iStrict][iSrc]->SetFillColor(kWhite);
	cout << "N POINTS = " << g_evt_contam[iEG][iStrict][iSrc]->GetN() << endl;

        g_evt_contam_err[iEG][iStrict][iSrc]->SetName(
						      "Event_"+strict_name[iStrict]+"_contamination_error_EG"+trigname[iEG+4]+"_"+src_name[iSrc] );
        g_evt_contam_err[iEG][iStrict][iSrc]->SetTitle(
						       "Event "+strict_name[iStrict]+" contamination (EG"+trigname[iEG+4]+") "+src_name[iSrc] );
        g_evt_contam_err[iEG][iStrict][iSrc]->SetMarkerStyle(mark_style[iSrc]);
        g_evt_contam_err[iEG][iStrict][iSrc]->SetMarkerColor(mark_col[iSrc]);
        g_evt_contam_err[iEG][iStrict][iSrc]->SetFillColor(kWhite);



      }
    }
  }
  cout << endl << endl;

  for(int iSrc=0 ; iSrc<nSrc ; iSrc++)
    for(int iStrict=0 ; iStrict<nStrict ; iStrict++)
      for(int iEG=0 ; iEG<nEG_g ; iEG++)
	cout << "N POINTS = " << g_evt_contam[iEG][iStrict][iSrc]->GetN() << endl;
  //for(int iVtx=0 ; iVtx<nPts ; iVtx++)
  //cout << n_contam[iEG][iStrict][iSrc][iVtx] << " || " <<
	    
  //////////////////////////////////////////////////////////
  // PLOTTING //////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  TCanvas * c_evt_contam[nEG_g][nStrict];
  TCanvas * c_evt_contam_err[nEG_g][nStrict];

  float g_max[nEG_g] = {0.3,0.4,0.6,0.35};

  //for(int iEG=0 ; iEG<nEG_g ; iEG++) {
  /*
    for(int iEG=0 ; iEG<2 ; iEG++) {
    for(int iStrict=0 ; iStrict<nStrict ; iStrict++) {
      
    c_evt_contam[iEG][iStrict] = new TCanvas("c_evt_contam_EG"+trigname[iEG+4]+"_"+strict_name[iStrict],
    "Spike contamination "+strict_name[iStrict]+" (EG"+trigname[iEG+4]+")",
    0,0,800,600);
      
    c_evt_contam[iEG][iStrict]->SetFillColor(kWhite);
    c_evt_contam[iEG][iStrict]->GetFrame()->SetFillColor(kWhite);

    g_evt_contam[iEG][iStrict][0]->SetMinimum(0);
    g_evt_contam[iEG][iStrict][0]->SetMaximum(g_max[iEG]);      
    g_evt_contam[iEG][iStrict][0]->GetXaxis()->Set(40,0,40);
    g_evt_contam[iEG][iStrict][0]->GetXaxis()->SetTitle("Number of vertices");
    g_evt_contam[iEG][iStrict][0]->GetYaxis()->SetTitle("Fraction of EG events trig. by spikes");

    g_evt_contam[iEG][iStrict][0]->SetTitle( "Spike contamination (EG"+trigname[iEG+4]+")" );

    g_evt_contam[iEG][iStrict][0]->Draw("AP");
    g_evt_contam[iEG][iStrict][1]->Draw("P");
    g_evt_contam[iEG][iStrict][2]->Draw("P");
    g_evt_contam[iEG][iStrict][3]->Draw("P");

    TLegend * leg = new TLegend(0.5,0.70,0.885,0.90,"CMS Preliminary 2011 pp #sqrt{s}=7 TeV","brNDC"); 
    //(x1,y1),(x2,y2)->bottom left , top right
    leg->SetLineColor(1);
    leg->SetTextColor(1);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetMargin(0.2);
    leg->SetShadowColor(kWhite);
    leg->SetFillColor(kWhite);
    //
    leg->AddEntry(g_evt_contam[iEG][iStrict][0],"Run 2011B Data","P");
    leg->AddEntry(g_evt_contam[iEG][iStrict][1],"Run 2011B Emul new setting (data)","P");
    leg->AddEntry(g_evt_contam[iEG][iStrict][2],"High PU runs Data","P");
    leg->AddEntry(g_evt_contam[iEG][iStrict][3],"High PU runs Emul new setting (data)","P");
    leg->Draw();

    c_evt_contam[iEG][iStrict]->Print(dirOut+"evt_contam_EG"+trigname[iEG+4]+"_"+strict_name[iStrict]+".C");
    c_evt_contam[iEG][iStrict]->Print(dirOut+"evt_contam_EG"+trigname[iEG+4]+"_"+strict_name[iStrict]+".gif");
    c_evt_contam[iEG][iStrict]->Print(dirOut+"evt_contam_EG"+trigname[iEG+4]+"_"+strict_name[iStrict]+".png");
    c_evt_contam[iEG][iStrict]->Print(dirOut+"evt_contam_EG"+trigname[iEG+4]+"_"+strict_name[iStrict]+".eps");
    c_evt_contam[iEG][iStrict]->Print(dirOut+"evt_contam_EG"+trigname[iEG+4]+"_"+strict_name[iStrict]+".ps");
    c_evt_contam[iEG][iStrict]->Print(dirOut+"evt_contam_EG"+trigname[iEG+4]+"_"+strict_name[iStrict]+".pdf");
    }
    }
  */
  cout << "<--- FINE UNTIL HERE --->" << endl;

  //for(int iEG=0 ; iEG<nEG_g ; iEG++) {
  for(int iEG=1 ; iEG<3 ; iEG++) {
    for(int iStrict=0 ; iStrict<nStrict ; iStrict++) {

      c_evt_contam_err[iEG][iStrict] = new TCanvas("c_evt_contam_err_EG"+trigname[iEG+4]+"_"+strict_name[iStrict],
						   "Spike contamination "+strict_name[iStrict]+" (EG"+trigname[iEG+4]+")",
						   0,0,800,600);

      c_evt_contam_err[iEG][iStrict]->SetFillColor(kWhite);
      c_evt_contam_err[iEG][iStrict]->GetFrame()->SetFillColor(kWhite);

      g_evt_contam_err[iEG][iStrict][0]->SetMinimum(0);
      g_evt_contam_err[iEG][iStrict][0]->SetMaximum(g_max[iEG]);
      g_evt_contam_err[iEG][iStrict][0]->GetXaxis()->Set(60,0,60);

      g_evt_contam_err[iEG][iStrict][0]->SetTitle( "Spike contamination (EG"+trigname[iEG+4]+")" );
      g_evt_contam_err[iEG][iStrict][0]->GetXaxis()->SetTitle("Number of vertices");
      g_evt_contam_err[iEG][iStrict][0]->GetYaxis()->SetTitle("Fraction of EG events trig. by spikes");

      g_evt_contam_err[iEG][iStrict][0]->Draw("AP");
      g_evt_contam_err[iEG][iStrict][1]->Draw("P");
      g_evt_contam_err[iEG][iStrict][2]->Draw("P");
      g_evt_contam_err[iEG][iStrict][3]->Draw("P");
      g_evt_contam_err[iEG][iStrict][4]->Draw("P");
      g_evt_contam_err[iEG][iStrict][5]->Draw("P");

      //TLegend * leg = new TLegend(0.17,0.695,0.58,0.90,"CMS Preliminary 2012 pp  #sqrt{s}=7+8 TeV","brNDC");
      TLegend * leg = new TLegend(0.17,0.675,0.58,0.88,"CMS Preliminary 2012 pp  #sqrt{s}=7+8 TeV","brNDC");
      //(x1,y1),(x2,y2)->bottom left , top right
      leg->SetMargin(0.15);
      leg->SetLineColor(1);
      leg->SetTextColor(1);
      leg->SetTextFont(42);
      leg->SetTextSize(0.03);
      leg->SetShadowColor(kWhite);
      leg->SetFillColor(kWhite);
      //
      leg->AddEntry(g_evt_contam[iEG][iStrict][0],"Run 2011B Data","P");
      leg->AddEntry(g_evt_contam[iEG][iStrict][1],"Run 2011B Emul new setting (data)","P");
      leg->AddEntry(g_evt_contam[iEG][iStrict][2],"High PU runs Data","P");
      leg->AddEntry(g_evt_contam[iEG][iStrict][3],"High PU runs Emul new setting (data)","P");
      leg->AddEntry(g_evt_contam[iEG][iStrict][4],"High PU 2012 run 198603 (data)","P");
      leg->AddEntry(g_evt_contam[iEG][iStrict][5],"High PU 2012 run 198609 (data)","P");
      leg->Draw();

      /*
	TLegend * leg2 = new TLegend(0.4,0.78,0.885,0.88,"","brNDC"); //(x1,y1),(x2,y2):bottom left,top right
	leg2->SetLineColor(1);
	leg2->SetTextColor(1);
	leg2->SetTextFont(42);
	leg2->SetTextSize(0.03);
	leg2->SetShadowColor(kWhite);
	leg2->SetFillColor(kWhite);
      
	leg2->AddEntry("NULL","Contamination = #frac{N(events triggered by a spike)}{N(events triggered by any object)}","h");
      
	leg2->Draw();
      */

      c_evt_contam_err[iEG][iStrict]->Print(dirOut+"evt_contam_err_EG"+trigname[iEG+4]+"_"+strict_name[iStrict]+".C");
      c_evt_contam_err[iEG][iStrict]->Print(dirOut+"evt_contam_err_EG"+trigname[iEG+4]+"_"+strict_name[iStrict]+".gif");
      c_evt_contam_err[iEG][iStrict]->Print(dirOut+"evt_contam_err_EG"+trigname[iEG+4]+"_"+strict_name[iStrict]+".png");
      c_evt_contam_err[iEG][iStrict]->Print(dirOut+"evt_contam_err_EG"+trigname[iEG+4]+"_"+strict_name[iStrict]+".eps");
      c_evt_contam_err[iEG][iStrict]->Print(dirOut+"evt_contam_err_EG"+trigname[iEG+4]+"_"+strict_name[iStrict]+".ps");
      c_evt_contam_err[iEG][iStrict]->Print(dirOut+"evt_contam_err_EG"+trigname[iEG+4]+"_"+strict_name[iStrict]+".pdf");

    }
  }

}
