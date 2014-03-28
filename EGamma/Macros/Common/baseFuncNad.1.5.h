/////////////////////////////////////////////////////////////////////////////
/////////////////////////////    baseFuncNad.h    ///////////////////////////
//////////////    Contains basic functions for my macros    /////////////////
/////////////////////////////////////////////////////////////////////////////
#ifndef BASEFUNC
#define BASEFUNC

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <utility>

#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "THStack.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
//#include "treesList.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "../Common/readJSONFile.cc"

using namespace std;

// To map iEntry to (nRun,nEvent) --> help reading trigOnly data
typedef map< pair<int,int> , int > MAPTrigOnly;
typedef map<int, vector<pair<int, int> > > MAPJson;
typedef map< pair<int,int> , pair<int,int> > MAPTT;
typedef map< pair<int,int> , pair< int , pair< vector<int> , vector<int> > > > MAPPY ;
// { <RCT ieta,RCT iphi> , < < iEle, <[eg1,eg2,eg5,eg8,eg10,eg12]N ,[...]M>  > }

typedef map< pair<int,int> , pair< pair<int,int>,pair<int,int> > > MAPTTS;
//           TT(ieta,iphi)           (adc,adcM)  (sFGVB,sFGVB_M)
// adc = mapTT[coords].first.first
// adcM = mapTT[coords].first.second
// sFGVB = maTT[coords].second.first
// sFGVB_M = maTT[coords].second.second

int detJson(int nRun) {

  // 0:2010B ; 1:May10 ; 2:Aug05 ; 3:Prompt

  if(nRun<160404) // 2010B
    return 0;
  if(nRun>=160404 && nRun<=163869) // May10
    return 1;
  if(nRun>=165071 && nRun<=168437) // PRV4
    return 3;
  if(nRun>=170053 && nRun<=172619) // Aug05 (PRV5)
    return 2;
  if(nRun>=172620 && nRun<=175770) // Oct03 (PRV6)
    return 3;
  if(nRun>=175832 && nRun<=180296) // 2011B (PRV1)
     return 3;

  return -1;
}
// check EB/EE localization
bool isInEB( int iRCTeta, int ietaTT ) {
  
  if( iRCTeta>=7 && iRCTeta<=14 ) return true;
  else if( iRCTeta==6 && ietaTT==-17 ) return true;
  else if( iRCTeta==15 && ietaTT==17 ) return true;
  else return false;

  return false;

}

//GETTING RCT/TT coordinates
// ====================================================================================
int getGCTRegionPhi(int ttphi)
// ====================================================================================
{
  int gctphi=0;
  gctphi = (ttphi+1)/4;
  if(ttphi<=2) gctphi=0;
  if(ttphi>=71) gctphi=0;
	
  return gctphi;
}

// ====================================================================================
int getGCTRegionEta(int tteta)
// ====================================================================================
{
  int gcteta = 0;
	
  if(tteta>0) gcteta = (tteta-1)/4 + 11;
  else if(tteta<0) gcteta = (tteta+1)/4 + 10;
	
  return gcteta;
}

vector<int> getPhiTT(int gctphi){
  vector<int> ttphi;

  if(gctphi==0){
    ttphi.push_back(71);
    ttphi.push_back(72);
    ttphi.push_back(1);
    ttphi.push_back(2);
  }
  else
    for(int i=0;i!=4;++i) ttphi.push_back(gctphi*4-1+i);

  return ttphi;
}

vector<int> getEtaTT(int gcteta){
  vector<int> tteta; 

  if(gcteta>=11)
    for(int i=1;i<=4;++i)
      tteta.push_back( (gcteta-11)*4+i );

  else
    for(int i=0;i!=4;++i)
      tteta.push_back( (gcteta-11)*4+i);

  return tteta;
}

void fireL1(int L1noniso, int L1iso, vector<int> & firedEG, vector<int> EG) {
  
  const int n=EG.size();
  
  int fired = -1;

  for(int i=0 ; i<n ; i++) {
    if( L1iso >= EG[n-1-i] || L1noniso >= EG[n-1-i] ) {
    //if(i>3) {
      fired = n-i-1;
      //cout << "L1iso=" << L1iso << " | L1noniso=" << L1noniso << " | fired EG : " << EG[n-1-i] << endl;
      break;
    }
  }
  
  if(fired>=0)
    for(int i=0 ; i<fired+1 ; i++)
      firedEG[i]=1;

}

void globalFireL1(double eleRCT_eT, int noniso, int iso, int nonisoM, int isoM, 
		  vector<int> & firedEG_N, vector<int> & firedEG_M, vector<int> menu) {

  if( eleRCT_eT > 2 ) {// interesting regions <=> et > 2GeV
    fireL1( noniso , iso , firedEG_N, menu );
    fireL1( nonisoM , isoM , firedEG_M, menu );    
  }
}

void globalFireL1_Normal(double eleRCT_eT, int noniso, int iso,
			 vector<int> & firedEG_N, vector<int> menu) {

  if( eleRCT_eT > 2 ) {// interesting regions <=> et > 2GeV
    fireL1( noniso , iso , firedEG_N, menu );
  }
}


int threshCateg(int nRun) {

  if(nRun>=146240 && nRun<=147042) return 7;
  if(nRun==147043) return 9;
  if(nRun==147048) return 0;
  if(nRun>=147114 && nRun<=147116) return 4;
  if(nRun>=147195 && nRun<=147222) return 6;
  if(nRun==147284) return 3;
  if(nRun==147390) return 5;
  if(nRun>=147450 && nRun<=147454 && nRun!=147451) return 8;
  if(nRun>=147485 && nRun<148524) return 2;
  if(nRun>=148524 && nRun<149378) return 1;
  if(nRun>=149378 &&  nRun<=149510) return 9;

  return -1 ;
}

void formatHisto(TH1F* histo, string title, string xtitle, string ytitle, int color) {
  histo->SetLineColor(color) ;
  histo->GetXaxis()->SetTitle( xtitle.c_str() ) ;
  histo->GetYaxis()->SetTitle( ytitle.c_str() ) ;
  histo->SetMarkerColor(color);
  histo->SetFillColor(color);
  //histo->Rebin(4);
  histo->SetMarkerStyle(21);
  histo->SetMarkerSize(1.);
  histo->SetLineWidth(2);
  histo->SetFillStyle (3018) ;
  histo->SetTitle( title.c_str() );
}

vector<int> det4MaxTab(vector<double> tab, int n) {

  //cout << "n=" << n << endl;
  if(n<=4) {
    vector<int> iMax;
    for(int i=0 ; i<n ; i++)
      iMax.push_back(i);
    //cout << "designed 1 size idxEle = " << iMax.size() << endl;
    return iMax;
  }
  
  else { 
    vector<int> iMax(4,0);
    // 0/1 <-> dont-look/look tab[i] 
    vector<int> idxToLook(tab.size(),1);
     
    // cherche l'indice des 4 maxima
    for(int j=0 ; j<4 ; j++) {
      for(int i=1 ; i<(int)tab.size() ; i++) {
	if(tab[i]>=tab[iMax[j]] && idxToLook[i]==1) iMax[j] = i;
      }
      idxToLook[iMax[j]] = 0 ; // ne plus regarder tab[iMax[j]]
    }    
    //cout << "designed 2 size idxEle = " << iMax.size() << endl;
    return iMax;
  }
}

int max(int val1, int val2) {
  if(val1>val2) return val1;
  else return val2;
}

vector<int> detPT(double pt) {

  vector<int> indices;
  indices.push_back(0);

  if(pt>=0 && pt<10) indices.push_back(1);
  else if(pt>=10 && pt<15) indices.push_back(2);
  else {
    if(pt>=15) indices.push_back(3);
    if(pt>=20) indices.push_back(4);
    if(pt>=50) indices.push_back(5);
    if(pt>=100) indices.push_back(6);
  }
  return indices ;
}

vector<int> detISP(int severity, int outoftime) {
  vector<int> indices;
  indices.clear();
  indices.push_back(0);
  if(severity==3) {
    if(outoftime==0) {
      indices.push_back(1);
      indices.push_back(2);
      indices.push_back(7);
      indices.push_back(8);
    }
    else if(outoftime==1) {
      indices.push_back(1);
      indices.push_back(2);
      indices.push_back(4);
      indices.push_back(7);
    }
  }
  else if(severity==4) {
      indices.push_back(1);
      indices.push_back(2);
      indices.push_back(3);
      indices.push_back(4);
  }  
  else {
    if(outoftime==0) {
      indices.push_back(5);
      indices.push_back(6);
      indices.push_back(7);
      indices.push_back(8);
    }
    else if(outoftime==1) {
      indices.push_back(5);
      indices.push_back(7);
    }
  }
  return indices;
}


bool VBTFcuts(TString asked, TString RunPhase,
	      double pt, double eT, double eta, double etaEle, double ele_tkSumPt_dr03, double ele_ecalRecHitSumEt_dr03, 
	      double ele_hcalDepth1TowerSumEt_dr03, double ele_hcalDepth2TowerSumEt_dr03, double ele_expected_inner_hits,
	      double ele_deltaphiin, double ele_deltaetain, double ele_he, double ele_sigmaietaieta,
	      double ele_conv_dist, double ele_conv_dcot, double ele_fbrem, int ele_isConversion)
{
  //isolation variables
  double trackIsoRel03 = ele_tkSumPt_dr03 / pt;
  double ecalIsoRel03 = ele_ecalRecHitSumEt_dr03 / pt;
  double hcalIsoRel03 = (ele_hcalDepth1TowerSumEt_dr03+ele_hcalDepth2TowerSumEt_dr03) / pt;

  // define cuts
  // expectedInnerHits ; EB:|DPhiIn|,|DEtaIn|,H/E,|tkIso|,|hIso|,|eIso| ; EE:idem
  double cuts[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};

  double cuts95[13] = {1 , 0.8 , 0.007 , 0.15 , 0.15 , 0.12 , 2.00 , 0.7 , 0.01 , 0.07 , 0.08 , 0.05 , 0.06} ;
  double cuts80[13] = {0 , 0.06 , 0.004 , 0.04 , 0.09 , 0.10 , 0.07 , 0.03 , 0.007 , 0.025 , 0.04 , 0.025 , 0.05} ;
  double cuts60[13] = {0 , 0.025 , 0.004 , 0.025 , 0.04 , 0.03 , 0.04 , 0.02 , 0.005 , 0.025 , 0.025 , 0.02 , 0.02} ;

  if(RunPhase=="2011A") {
    cuts95[9] = cuts80[9] = cuts60[9] = 0.15 ;
  }

  // define HLT_Ele27 cuts
  //double trackIsoRel03_H = ele_tkSumPt_dr03 / eT;
  //double ecalIsoRel03_H = ele_ecalRecHitSumEt_dr03 / eT;
  //double hcalIsoRel03_H = (ele_hcalDepth1TowerSumEt_dr03+ele_hcalDepth2TowerSumEt_dr03) / eT;

  // HoE , sigmaietaieta , ecaliso/et , hcaliso/et , deta, dphi, tkiso/et
  double cutsHLT_Ele27_IdIso_EB[7] = {0.05 , 0.011 , 0.15 , 0.15 , 0.008 , 0.07 , 0.15};
  double cutsHLT_Ele27_IdIso_EE[7] = {0.05 , 0.031 , 0.075 , 0.075 , 0.008 , 0.05 , 0.075};

  double values[7] = { ele_he , ele_sigmaietaieta , ecalIsoRel03 , hcalIsoRel03 , fabs(ele_deltaetain) , fabs(ele_deltaphiin) , trackIsoRel03 } ;

  // "conversion" selection
  if(asked=="conv") {
    // ID
    if ( fabs(ele_deltaphiin) >= 0.1 ) return false;
    if ( ele_he >= 0.1 ) return false;
    if ( ele_fbrem <= -0.1 ) return false;
    if ( fabs(etaEle)<1.479 ) {
      if ( ele_sigmaietaieta <= 0.008 ) return false;
    }
    else {
      return false;
      if ( ele_sigmaietaieta <= 0.02 ) return false;
    }
    if( fabs(ele_deltaetain) >= 0.01 ) return false;

    // ISO
    if ( fabs(ele_tkSumPt_dr03) >= 0.5 ) return false;
    if ( fabs(ele_hcalDepth1TowerSumEt_dr03) >= 0.5 ) return false;

    // "conversion ID"
    if( ele_isConversion != 1 ) return false;

    // if passed all cuts...
    return true;
  }

  // HLT-like selection
  else if(asked=="HLT_Ele27") {
    for(int i=0 ; i<7 ; i++) {
      if( etaEle<1.479 ) {
	if( values[i] >= cutsHLT_Ele27_IdIso_EB[i] ) return false;
      }
      else {
	if( values[i] >= cutsHLT_Ele27_IdIso_EE[i] ) return false;
      }
    }
    return true;
  }
  
  // VBTF selection
  else if(asked=="WP95" || asked=="WP80" || asked=="WP60") {

    if(asked=="WP95") {
      for(int i=0 ; i<13 ; i++)
	cuts[i] = cuts95[i] ;
    }
    else if(asked=="WP80") {
      for(int i=0 ; i<13 ; i++)
	cuts[i] = cuts80[i] ;
    }
    else if(asked=="WP60") {
      for(int i=0 ; i<13 ; i++)
	cuts[i] = cuts60[i] ;
    }

    // missing hits
    if ( ele_expected_inner_hits > cuts[0] ) return false ;
	  
    if( fabs(etaEle)<1.479 ){ // barrel
      // idiEletification
      if ( fabs(ele_deltaphiin) >= cuts[1] )   return false ;
      if ( fabs(ele_deltaetain) >= cuts[2] ) return false ;
      if ( ele_he >= cuts[3] )                 return false ;
      if ( ele_sigmaietaieta >= 0.01 ) return false ;
	    
      // Relative isolation
      if ( fabs(trackIsoRel03)>= cuts[4] ) return false ;
      if ( fabs(hcalIsoRel03) >= cuts[5] ) return false ;
      if ( fabs(ecalIsoRel03) >= cuts[6] ) return false;

      // Combined Isolation
      
    } else { // endcap
      //return false; // remove endcap electrons
      // identification
      if ( fabs(ele_deltaphiin) >= cuts[7] )   return false ;
      if ( fabs(ele_deltaetain) >= cuts[8] ) return false ;
      if ( ele_he >= cuts[9] )                 return false ;
      if ( ele_sigmaietaieta >= 0.03 ) return false ;

      // Relative isolation
      if ( fabs(trackIsoRel03)>= cuts[10] ) return false ;
      if ( fabs(hcalIsoRel03) >= cuts[11] ) return false ;
      if ( fabs(ecalIsoRel03) >= cuts[12] ) return false;
    }

    //fiducial isolation
    //if ( fabs(eta) >= 2.5 || (fabs(eta)< 1.566 && fabs(eta)>1.4442)) return false ;
 
    // conversion rejection
    if(asked=="WP80")
      if ( fabs(ele_conv_dist) <= 0.02 && fabs(ele_conv_dcot) <=0.02 ) return false;
    if(asked=="WP60")
      if ( fabs(ele_conv_dist) <= 0.02 && fabs(ele_conv_dcot) <=0.02 ) return false;

    // if passed all cuts...
    return true;
  }

  else {
    cout << "PLEASE CHOOSE CORRECT ASKED STRING FOR CUTS" << endl;
    return false;
  }

}

vector<int> passedVBTFcuts(TString asked, 
			   double pt, double eta, double etaEle, double ele_tkSumPt_dr03, double ele_ecalRecHitSumEt_dr03, 
			   double ele_hcalDepth1TowerSumEt_dr03, double ele_hcalDepth2TowerSumEt_dr03, double ele_expected_inner_hits,
			   double ele_deltaphiin, double ele_deltaetain, double ele_he, double ele_sigmaietaieta,
			   double ele_conv_dist, double ele_conv_dcot, double ele_fbrem, int ele_isConversion)
{

  // output vector
  vector<int> passed; 
  passed.push_back(0);
  /*
    0 : all
    1 to 15 : common cuts
    16 : additionnal cut for WP80
    17 : EB / 18 : EE
    19 : EB passing all cuts / 20 : EE passing all cuts
    --> 21 indices
  */

  //isolation variables
  double trackIsoRel03 = ele_tkSumPt_dr03 / pt;
  double ecalIsoRel03 = ele_ecalRecHitSumEt_dr03 / pt;
  double hcalIsoRel03 = (ele_hcalDepth1TowerSumEt_dr03+ele_hcalDepth2TowerSumEt_dr03) / pt;

  // define cuts
  const int nCuts = 15 ; // + dernier cut pour wp80
  double cuts[nCuts] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double cuts95[nCuts] = {1 , 0.8 , 0.007 , 0.15 , 0.15 , 0.12 , 2.00 , 0.7 , 0.01 , 0.07 , 0.08 , 0.05 , 0.06 , 0.01 , 0.03} ;
  double cuts80[nCuts] = {0 , 0.06 , 0.004 , 0.04 , 0.09 , 0.10 , 0.07 , 0.03 , 0.007 , 0.025 , 0.04 , 0.025 , 0.05 , 0.01 , 0.03} ;

  // for selection "conversion"

  // determine which cuts to use
  if(asked=="WP95")
    for(int i=0 ; i<nCuts ; i++)
      cuts[i] = cuts95[i] ;
  else if(asked=="WP80")
    for(int i=0 ; i<nCuts ; i++)
      cuts[i] = cuts80[i] ;
  else {
    return passed ;
  }

  bool failed;
  failed = false;

  // missing hits
  if ( ele_expected_inner_hits <= cuts[0] ) passed.push_back(16)  ;
  else failed = true;

  if( fabs(etaEle)<1.479 ){ // barrel : dernier bit
    passed.push_back(17); // 17<->EB
    // idiEletification
    if ( fabs(ele_deltaphiin) < cuts[1] ) passed.push_back(1)  ;
    else failed = true;
    if ( fabs(ele_deltaetain) < cuts[2] )  passed.push_back(2) ;
    else failed = true;
    if ( ele_he < cuts[3] )                  passed.push_back(3) ;
    else failed = true;
    if ( ele_sigmaietaieta < cuts[13] )  passed.push_back(13) ;
    else failed = true;
	    
    // Relative isolation
    if ( fabs(trackIsoRel03)< cuts[4] )  passed.push_back(4) ;
    else failed = true;
    if ( fabs(hcalIsoRel03) < cuts[5] )  passed.push_back(5) ;
    else failed = true;
    if ( fabs(ecalIsoRel03) < cuts[6] )  passed.push_back(6) ;
    else failed = true;
      
  } else { // endcap
    passed.push_back(18); // 18<->EE
    // identification
    if ( fabs(ele_deltaphiin) < cuts[7] )    passed.push_back(7) ;
    else failed = true;
    if ( fabs(ele_deltaetain) < cuts[8] )  passed.push_back(8) ;
    else failed = true;
    if ( ele_he < cuts[9] )                  passed.push_back(9) ;
    else failed = true;
    if ( ele_sigmaietaieta < cuts[14] )  passed.push_back(14) ;
    else failed = true;

    // Relative isolation
    if ( fabs(trackIsoRel03)< cuts[10] )  passed.push_back(10) ;
    else failed = true;
    if ( fabs(hcalIsoRel03) < cuts[11] )  passed.push_back(11) ;
    else failed = true;
    if ( fabs(ecalIsoRel03) < cuts[12] )  passed.push_back(12) ;
    else failed = true;
  }

  //fiducial isolation
  //if ( fabs(eta) >= 2.5 || (fabs(eta)< 1.566 && fabs(eta)>1.4442)) return false ;
 
  // conversion rejection : cut #15
  if(asked=="WP80") {
    if ( fabs(ele_conv_dist) > 0.02 || fabs(ele_conv_dcot) > 0.02 ) passed.push_back(15);
    else failed = true;
  }

  if( !failed ) {
    if( fabs(etaEle)<1.479 )
      passed.push_back(19);
    else passed.push_back(20);
  }

  return passed ;

}

void writeTab(ofstream& outtab, int sfgvb, int eleMatched, int eleKilled, 
	      int eleMatchedSevNo34, int eleMatchedSev4, int eleMatchedSev34,
	      int eleKilledSevNo34, int eleKilledSev4, int eleKilledSev34)
{

  outtab << sfgvb << " | "
	 << eleMatched << " | "
	 << eleKilled << " | "
	 << (double)eleKilled/((double)eleMatched) << " | "
	 << eleMatchedSevNo34 << " | "
	 << eleKilledSevNo34 << " | "
	 << (double)eleKilledSevNo34/((double)eleMatchedSevNo34) << " | "
	 << eleMatchedSev4 << " | "
	 << eleKilledSev4 << " | "
	 << (double)eleKilledSev4/((double)eleMatchedSev4) << " | "
	 << eleMatchedSev34 << " | "
	 << eleKilledSev34 << " | "
	 << (double)eleKilledSev34/((double)eleMatchedSev34) << " | "
	 << endl;
}

int whichCuts(TString RunPhase, double pt, double et, double eta, double etaEle, double ele_tkSumPt_dr03, double ele_ecalRecHitSumEt_dr03, 
	      double ele_hcalDepth1TowerSumEt_dr03, double ele_hcalDepth2TowerSumEt_dr03, double ele_expected_inner_hits,
	      double ele_deltaphiin, double ele_deltaetain, double ele_he, double ele_sigmaietaieta,
	      double ele_conv_dist, double ele_conv_dcot, double ele_fbrem, int ele_isConversion)
{
  int res = 0;

  bool wp95, wp80, wp60, conv;
  wp95 = wp80 = wp60 = conv = false ;

  wp95 = VBTFcuts( "WP95", RunPhase,
		   pt,  et, eta,  etaEle,  ele_tkSumPt_dr03,  ele_ecalRecHitSumEt_dr03, 
		   ele_hcalDepth1TowerSumEt_dr03,  ele_hcalDepth2TowerSumEt_dr03,  ele_expected_inner_hits,
		   ele_deltaphiin,  ele_deltaetain,  ele_he,  ele_sigmaietaieta,
		   ele_conv_dist,  ele_conv_dcot,  ele_fbrem, ele_isConversion ) ;

  wp80 = VBTFcuts( "WP80", RunPhase,
		   pt,  et, eta,  etaEle,  ele_tkSumPt_dr03,  ele_ecalRecHitSumEt_dr03, 
		   ele_hcalDepth1TowerSumEt_dr03,  ele_hcalDepth2TowerSumEt_dr03,  ele_expected_inner_hits,
		   ele_deltaphiin,  ele_deltaetain,  ele_he,  ele_sigmaietaieta,
		   ele_conv_dist,  ele_conv_dcot,  ele_fbrem, ele_isConversion ) ;

  wp60 = VBTFcuts( "WP60", RunPhase,
		   pt,  et, eta,  etaEle,  ele_tkSumPt_dr03,  ele_ecalRecHitSumEt_dr03, 
		   ele_hcalDepth1TowerSumEt_dr03,  ele_hcalDepth2TowerSumEt_dr03,  ele_expected_inner_hits,
		   ele_deltaphiin,  ele_deltaetain,  ele_he,  ele_sigmaietaieta,
		   ele_conv_dist,  ele_conv_dcot,  ele_fbrem, ele_isConversion ) ;

  if(wp95) res = 1;
  if(wp80) res = 2;
  if(wp60) res = 3;

  return res ;

}

#endif
