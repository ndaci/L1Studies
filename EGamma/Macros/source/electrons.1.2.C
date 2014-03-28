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
//#include "treesList.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "/home/llr/cms/broutin/cmssw/CMSSW_3_8_4/src/EGamma/ECGelec/test/macro/Charge/new38x/readJSONFile.cc"

using namespace std;
typedef map< pair<int,int> , pair<int,int> > MAPTT;
typedef map< pair<int,int> , pair< int , pair< vector<int> , vector<int> > > > MAPPY ;
// { <RCT ieta,RCT iphi> , < < iEle, <[eg1,eg2,eg5,eg8,eg10,eg12]N ,[...]M>  > }

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

void fireL1(int L1noniso, int L1iso, std::vector<int> & firedEG) {
  
  const int n=6;
  int EG[n] = {2,4,10,16,20,24};  
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
  if(val1>=val2) return val1;
  if(val1<=val2) return val2;
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

bool VBTFcuts(TString asked, 
	      double pt, double eta, double etaEle, double ele_tkSumPt_dr03, double ele_ecalRecHitSumEt_dr03, 
	      double ele_hcalDepth1TowerSumEt_dr03, double ele_hcalDepth2TowerSumEt_dr03, double ele_expected_inner_hits,
	      double ele_deltaphiin, double ele_deltaetain, double ele_he, double ele_sigmaietaieta,
	      double ele_conv_dist, double ele_conv_dcot, double ele_fbrem, int ele_isConversion)
{
  //isolation variables
  double trackIsoRel03 = ele_tkSumPt_dr03 / pt;
  double ecalIsoRel03 = ele_ecalRecHitSumEt_dr03 / pt;
  double hcalIsoRel03 = (ele_hcalDepth1TowerSumEt_dr03+ele_hcalDepth2TowerSumEt_dr03) / pt;

  // define cuts
  double cuts[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
  double cuts95[13] = {1 , 0.8 , 0.007 , 0.15 , 0.15 , 0.12 , 2.00 , 0.7 , 0.01 , 0.07 , 0.08 , 0.05 , 0.06} ;
  double cuts80[13] = {0 , 0.06 , 0.004 , 0.04 , 0.09 , 0.10 , 0.07 , 0.03 , 0.007 , 0.025 , 0.04 , 0.025 , 0.05} ;

  // for selection "conversion"
  bool wait = false;

  // determine which cuts to use
  if(asked=="WP95")
    for(int i=0 ; i<13 ; i++)
      cuts[i] = cuts95[i] ;
  else if(asked=="WP80")
    for(int i=0 ; i<13 ; i++)
      cuts[i] = cuts80[i] ;
  else if(asked=="conv")
    wait = true;
  else
    return true;

  // "conversion" selection
  if(wait) {
    // ID
    if ( fabs(ele_deltaphiin) >= 0.1 ) return false;
    if ( ele_he >= 0.1 ) return false;
    if ( ele_fbrem <= -0.1 ) return false;
    if ( etaEle<1.479 ) {
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
  
  // VBTF selection
  else {
    // missing hits
    if ( ele_expected_inner_hits > cuts[0] ) return false ;
	  
    if( etaEle<1.479 ){ // barrel
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
      return false; // remove endcap electrons
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

    // if passed all cuts...
    return true;
  }

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

void TreeToHistos(TChain* myChain, TString dirOut, TString asked, int nEntries, ofstream& outlog, bool debug=false)  {
  
  // Output trees
  if(debug) cout << "declare output trees" << endl;
  Int_t l1_1;
  Int_t l1_2;
  Int_t l1_5;
  Int_t l1_8;
  Int_t l1_10;
  Int_t l1_12;
  float sc_et;
  float sc_dr;

  Int_t l1_1_2;
  Int_t l1_2_2;
  Int_t l1_5_2;
  Int_t l1_8_2;
  Int_t l1_10_2;
  Int_t l1_12_2;
  float sc_et_2;
  float sc_dr_2;

  TTree* treenew[11][9][3];
  TTree* treenew_2[11][9][3];

  int thresh[11] = {10,15,17,19,21,23,30,35,40,50,-1} ;
  TString baseCategory[11]={"t10","t15","t17","t19","t21","t23","t30","t35","t40","t50","tAll"} ;
  TString spName[9] = {"_allEle","_sev34","_sev34_or_oot1","_sev4","_sev4_or_oot1",
		       "_sevNo34","_sevNo34_and_oot0","_sevNo4","_sevNo4_and_oot0"};
  TString ptBin[3] = {"_AllPt","_Under20GeV","_Above20GeV"};

  int EG[7] = {0,2,4,10,16,20,24} ;
  TString EGval[7] = {"0","1","2","5","8","10","12"} ;

  int nCat = 11;
  int nSp = 9;
  int nPt = 3;

  TString name1,name2;

  if(debug) cout << "output trees declared" << endl << "define these trees" << endl;
  
  for(int iCat=0 ; iCat<nCat ; iCat++) {
    for(int iSp=0 ; iSp<nSp ; iSp++) {
      for(int iPt=0 ; iPt<nPt ; iPt++ ) {

	name1 = "treenew_"+baseCategory[iCat]+spName[iSp]+ptBin[iPt];
	name2 = "treenew_2_"+baseCategory[iCat]+spName[iSp]+ptBin[iPt];

	treenew[iCat][iSp][iPt] = new TTree(name1,name1);
	treenew[iCat][iSp][iPt]->Branch("l1_1" ,&l1_1,  "l1_1/I");
	treenew[iCat][iSp][iPt]->Branch("l1_2" ,&l1_2,  "l1_2/I");
	treenew[iCat][iSp][iPt]->Branch("l1_5" ,&l1_5,  "l1_5/I");
	treenew[iCat][iSp][iPt]->Branch("l1_8" ,&l1_8,  "l1_8/I");
	treenew[iCat][iSp][iPt]->Branch("l1_10" ,&l1_10,  "l1_10/I");
	treenew[iCat][iSp][iPt]->Branch("l1_12" ,&l1_12,  "l1_12/I");
	treenew[iCat][iSp][iPt]->Branch("sc_et",&sc_et, "sc_et/F");
	treenew[iCat][iSp][iPt]->Branch("sc_dr",&sc_dr, "sc_dr/F");
	
	treenew_2[iCat][iSp][iPt]=new TTree(name2,name2);
	treenew_2[iCat][iSp][iPt]->Branch("l1_1_2" ,&l1_1_2,  "l1_1_2/I");
	treenew_2[iCat][iSp][iPt]->Branch("l1_2_2" ,&l1_2_2,  "l1_2_2/I");
	treenew_2[iCat][iSp][iPt]->Branch("l1_5_2" ,&l1_5_2,  "l1_5_2/I");
	treenew_2[iCat][iSp][iPt]->Branch("l1_8_2" ,&l1_8_2,  "l1_8_2/I");
	treenew_2[iCat][iSp][iPt]->Branch("l1_10_2" ,&l1_10_2,  "l1_10_2/I");
	treenew_2[iCat][iSp][iPt]->Branch("l1_12_2" ,&l1_12_2,  "l1_12_2/I");
	treenew_2[iCat][iSp][iPt]->Branch("sc_et_2",&sc_et_2, "sc_et_2/F"); 
	treenew_2[iCat][iSp][iPt]->Branch("sc_dr_2",&sc_dr_2, "sc_dr_2/F");  
      }
    }
  }

  if(debug) cout << "get the input tree" << endl;

  // get the tree
  int nEvent, nRun, nLumi, trig_isUnbiased;   
  // Electrons
  TClonesArray * electrons = new TClonesArray ("TLorentzVector");
  int ele_N, sc_hybrid_N; 
  int ele_outOfTimeSeed[25],ele_severityLevelSeed[25];
  double ele_e1[25], ele_e33[25];
  double ele_he[25], ele_sigmaietaieta[25];
  double ele_hcalDepth1TowerSumEt_dr03[25], ele_hcalDepth2TowerSumEt_dr03[25];
  double ele_ecalRecHitSumEt_dr03[25], ele_tkSumPt_dr03[25];
  double ele_sclEta[25];
  double ecalIsoRel03,hcalIsoRel03,trackIsoRel03;
  double ele_deltaphiin[25], ele_deltaetain[25];
  double ele_conv_dist[25], ele_conv_dcot[25];
  double ele_fbrem[25];
  int ele_expected_inner_hits[25];
  int ele_ambiguousGsfTracks[25];
  int ele_isConversion[25];
  //
  int ele_RCTeta[25], ele_RCTphi[25], ele_RCTL1iso[25], ele_RCTL1noniso[25], ele_RCTL1iso_modif[25], 
    ele_RCTL1noniso_modif[25];
//   int ele_TTetaVect[25][50], ele_TTphiVect[25][50];
//   double ele_TTetVect[25][50];
  int ele_RCTetaVect[25][10], ele_RCTphiVect[25][10], ele_RCTL1isoVect[25][10], 
    ele_RCTL1nonisoVect[25][10],ele_RCTL1isoVect_modif[25][10], ele_RCTL1nonisoVect_modif[25][10];
  double ele_RCTetVect[25][10];
  // TP info
  int trig_tower_N,trig_tower_ieta[4032],trig_tower_iphi[4032],trig_tower_adc[4032]; 
  int trig_tower_N_modif,trig_tower_ieta_modif[4032],trig_tower_iphi_modif[4032],trig_tower_adc_modif[4032]; 
  int trig_L1emIso_N, trig_L1emNonIso_N, trig_L1emIso_N_modif, trig_L1emNonIso_N_modif;

  // L1 candidates info
  int trig_L1emIso_ieta[4], trig_L1emIso_iphi[4], trig_L1emIso_rank[4];
  int trig_L1emNonIso_ieta[4], trig_L1emNonIso_iphi[4], trig_L1emNonIso_rank[4];
  int trig_L1emIso_ieta_modif[4], trig_L1emIso_iphi_modif[4], trig_L1emIso_rank_modif[4];
  int trig_L1emNonIso_ieta_modif[4], trig_L1emNonIso_iphi_modif[4], trig_L1emNonIso_rank_modif[4];

  // Global
  myChain->SetBranchAddress("nEvent",&nEvent);
  myChain->SetBranchAddress("nRun",&nRun);
  myChain->SetBranchAddress("nLumi",&nLumi);
  myChain->SetBranchAddress("trig_isUnbiased",&trig_isUnbiased);

  // SC
  myChain->SetBranchAddress("sc_hybrid_N",   &sc_hybrid_N);
  
  // Electrons
  myChain->SetBranchAddress("ele_N",   &ele_N);
  myChain->SetBranchAddress("electrons",&electrons);
  myChain->SetBranchAddress("ele_outOfTimeSeed",&ele_outOfTimeSeed);
  myChain->SetBranchAddress("ele_severityLevelSeed", &ele_severityLevelSeed);
  myChain->SetBranchAddress("ele_e1", &ele_e1);
  myChain->SetBranchAddress("ele_e33", &ele_e33);
  myChain->SetBranchAddress("ele_he",&ele_he);
  myChain->SetBranchAddress("ele_sigmaietaieta",&ele_sigmaietaieta);
  myChain->SetBranchAddress("ele_hcalDepth1TowerSumEt_dr03", &ele_hcalDepth1TowerSumEt_dr03);
  myChain->SetBranchAddress("ele_hcalDepth2TowerSumEt_dr03", &ele_hcalDepth2TowerSumEt_dr03);
  myChain->SetBranchAddress("ele_ecalRecHitSumEt_dr03", &ele_ecalRecHitSumEt_dr03);
  myChain->SetBranchAddress("ele_tkSumPt_dr03",&ele_tkSumPt_dr03);
  myChain->SetBranchAddress("ele_sclEta",&ele_sclEta);
  myChain->SetBranchAddress("ele_expected_inner_hits",&ele_expected_inner_hits);
  myChain->SetBranchAddress("ele_deltaphiin",&ele_deltaphiin);
  myChain->SetBranchAddress("ele_deltaetain",&ele_deltaetain);
  myChain->SetBranchAddress("ele_conv_dist",&ele_conv_dist);
  myChain->SetBranchAddress("ele_conv_dcot",&ele_conv_dcot);
  myChain->SetBranchAddress("ele_fbrem",&ele_fbrem);
  myChain->SetBranchAddress("ele_ambiguousGsfTracks", &ele_ambiguousGsfTracks);
  myChain->SetBranchAddress("ele_isConversion",&ele_isConversion);
 
  // L1 stuff
//   myChain->SetBranchAddress("ele_TTetaVect", &ele_TTetaVect);
//   myChain->SetBranchAddress("ele_TTphiVect", &ele_TTphiVect);
//   myChain->SetBranchAddress("ele_TTetVect", &ele_TTetVect);
//   myChain->SetBranchAddress("ele_RCTetaVect", &ele_RCTetaVect);
//   myChain->SetBranchAddress("ele_RCTphiVect", &ele_RCTphiVect);
//   myChain->SetBranchAddress("ele_RCTetVect", &ele_RCTetVect);
//   myChain->SetBranchAddress("ele_RCTL1isoVect", &ele_RCTL1isoVect);
//   myChain->SetBranchAddress("ele_RCTL1nonisoVect", &ele_RCTL1nonisoVect);
  //
  myChain->SetBranchAddress("ele_RCTeta", &ele_RCTeta);
  myChain->SetBranchAddress("ele_RCTphi", &ele_RCTphi);
  myChain->SetBranchAddress("ele_RCTL1iso", &ele_RCTL1iso);
  myChain->SetBranchAddress("ele_RCTL1noniso", &ele_RCTL1noniso);
  myChain->SetBranchAddress("ele_RCTL1iso_modif", &ele_RCTL1iso_modif);
  myChain->SetBranchAddress("ele_RCTL1noniso_modif", &ele_RCTL1noniso_modif);

  myChain->SetBranchAddress("ele_RCTetaVect", &ele_RCTetaVect);
  myChain->SetBranchAddress("ele_RCTphiVect", &ele_RCTphiVect);
  myChain->SetBranchAddress("ele_RCTetVect", &ele_RCTetVect);
  myChain->SetBranchAddress("ele_RCTL1isoVect", &ele_RCTL1isoVect);
  myChain->SetBranchAddress("ele_RCTL1nonisoVect", &ele_RCTL1nonisoVect);
  myChain->SetBranchAddress("ele_RCTL1isoVect_modif", &ele_RCTL1isoVect_modif);
  myChain->SetBranchAddress("ele_RCTL1nonisoVect_modif", &ele_RCTL1nonisoVect_modif);


  // L1 candidates
  myChain->SetBranchAddress("trig_L1emIso_N", &trig_L1emIso_N);
  myChain->SetBranchAddress("trig_L1emIso_ieta", &trig_L1emIso_ieta);
  myChain->SetBranchAddress("trig_L1emIso_iphi", &trig_L1emIso_iphi);
  myChain->SetBranchAddress("trig_L1emIso_rank", &trig_L1emIso_rank);

  myChain->SetBranchAddress("trig_L1emNonIso_N", &trig_L1emNonIso_N);
  myChain->SetBranchAddress("trig_L1emNonIso_ieta", &trig_L1emNonIso_ieta);
  myChain->SetBranchAddress("trig_L1emNonIso_iphi", &trig_L1emNonIso_iphi); 
  myChain->SetBranchAddress("trig_L1emNonIso_rank", &trig_L1emNonIso_rank);

  myChain->SetBranchAddress("trig_L1emIso_N_modif", &trig_L1emIso_N_modif);
  myChain->SetBranchAddress("trig_L1emIso_ieta_modif", &trig_L1emIso_ieta_modif);
  myChain->SetBranchAddress("trig_L1emIso_iphi_modif", &trig_L1emIso_iphi_modif);
  myChain->SetBranchAddress("trig_L1emIso_rank_modif", &trig_L1emIso_rank_modif);

  myChain->SetBranchAddress("trig_L1emNonIso_N_modif", &trig_L1emNonIso_N_modif);
  myChain->SetBranchAddress("trig_L1emNonIso_ieta_modif", &trig_L1emNonIso_ieta_modif);
  myChain->SetBranchAddress("trig_L1emNonIso_iphi_modif", &trig_L1emNonIso_iphi_modif);
  myChain->SetBranchAddress("trig_L1emNonIso_rank_modif", &trig_L1emNonIso_rank_modif);

  // Trigger Towers
  // normal collection
  myChain->SetBranchAddress("trig_tower_N", &trig_tower_N);
  myChain->SetBranchAddress("trig_tower_ieta",  &trig_tower_ieta);
  myChain->SetBranchAddress("trig_tower_iphi",  &trig_tower_iphi);
  myChain->SetBranchAddress("trig_tower_adc",  &trig_tower_adc);
 
  // modified collection
  myChain->SetBranchAddress("trig_tower_N_modif", &trig_tower_N_modif);
  myChain->SetBranchAddress("trig_tower_ieta_modif",  &trig_tower_ieta_modif);
  myChain->SetBranchAddress("trig_tower_iphi_modif",  &trig_tower_iphi_modif);
  myChain->SetBranchAddress("trig_tower_adc_modif",  &trig_tower_adc_modif);

  if(debug) cout << "got the input tree" << endl;

  TLorentzVector* theEle;
  bool isGoodRun,check_barrel,selected,output,outputRCT,outputRegion,weird,outOfTime;
  vector<int> idxEle;
  string flag;
  int iEle, iEleInMap;
  int iso, noniso, isoM, nonisoM, isop, nonisop, isopM, nonisopM, adc, adc_modif;
  double eta, phi, etaEle;
  float dr,drM;
  vector<double> ele_Et;
  vector<int> ietaTTele, iphiTTele, ietaTT, iphiTT;
  vector<int> firedEG_N; // EG 1/2/5/8/10/12 
  vector<int> firedEG_M; 

  // Counters
  int nEle[11][9][3];    // [sFGVB][spikyness][ptbin]
  int nMatched[11][9][3][7]; // ... + [EGX]
  int nMatchedM[11][9][3][7];
  int nLost[11][9][3][7];
  int nAppeared[11][9][3][7];
  int nSwitch[11][9][3]; 
  int nLower[11][9][3];  
  int nHigher[11][9][3];
  int nSeverEleSameReg[11][9][3];

  int nLostL1NoSpiky = 0;

  for(int iCat=0 ; iCat<nCat ; iCat++) {
    for(int iSp=0 ; iSp<nSp ; iSp++) {
      for(int iPt=0 ; iPt<nPt ; iPt++) {
	nEle[iCat][iSp][iPt] = 0 ;
	for(int iEG=0 ; iEG<7 ; iEG++) {
	  nMatched[iCat][iSp][iPt][iEG] = 0 ;
	  nMatchedM[iCat][iSp][iPt][iEG] = 0 ;
	  nLost[iCat][iSp][iPt][iEG] = 0 ;
	  nAppeared[iCat][iSp][iPt][iEG] = 0 ;
	}
	nSwitch[iCat][iSp][iPt] = 0 ; 
	nLower[iCat][iSp][iPt] = 0 ;  
	nHigher[iCat][iSp][iPt] = 0 ;
	nSeverEleSameReg[iCat][iSp][iPt] = 0 ;
      }
    }
  }
  
  MAPTT adcTT;
  MAPTT::iterator iterTT;
  pair<int,int> coords;

  MAPPY mapEt;
  MAPPY::iterator iter;

  TString filename;

  int nLostCheck[10], nMatchedCheck[10], nAppearedCheck[10], nMatchedCheckM[10], nEvents[11][6], nEventsM[11][6];
  for(int i=0 ; i<10 ; i++) {
    nLostCheck[i]=nMatchedCheck[i]=nAppearedCheck[i]=nMatchedCheckM[i]= 0;
  }
  for(int i=0 ; i<11 ; i++)
    for(int j=0 ; j<6 ; j++)
      nEventsM[i][j] = nEvents[i][j] = 0;

  bool eventN[6], eventM[6];
  bool matcheck, matcheckM;  

  vector<int> IdxCat;
  vector<int> IdxSp;
  vector<int> IdxPt;
  
  // -------------------------------------------------------------------------------
  // JSON FILE READER
  // -------------------------------------------------------------------------------
  // define map of run/LS
  string jsonFile = "/home/llr/cms/ochando/Analysis/Ana_Multi/json/goodrunlist_json.txt";
  map<int, vector<pair<int, int> > > jsonMap = readJSONFile(jsonFile);   
  
  if(debug) cout << "gonna loop over events" << endl;

  // loop over events
  int numEntries = myChain->GetEntries () ;
  int nProcess = numEntries;
  if(nEntries>=0 && nEntries<numEntries)
    nProcess = nEntries;

  int nCurrentRun;
  outlog << "will process " << nProcess << "/" << numEntries << "entries" << endl;
  for (int iEvent = 0 ; iEvent < nProcess ; iEvent++ )
    { 

      myChain->GetEntry (iEvent) ;
     
      // run selection
      isGoodRun = AcceptEventByRunAndLumiSection(nRun, nLumi, jsonMap);
      if(!isGoodRun) continue;

      // HLT selection
      //if(trig_isUnbiased!=1) continue;

      // show processed file
      if(iEvent==0) {
        filename = myChain->GetFile()->GetName() ;
        outlog << "File : " << filename << endl << endl;
      }
      else if( filename != myChain->GetFile()->GetName() ) {
        filename = myChain->GetFile()->GetName() ;
        outlog << "File : " << myChain->GetFile()->GetName() << endl << endl;
      }

      // det sFGVB threshold category
      int iCat = threshCateg(nRun);
      if(iCat<0 || iCat>9) continue;
      //if(iCat != 9) continue;
      
      // show current run/iCat processed
      if(iEvent==0) {
	nCurrentRun = nRun ;
	outlog << "nRun=" << nRun << " | iCat=" << iCat << endl;
      }
      else if(nRun!=nCurrentRun) {
	nCurrentRun=nRun ;
	outlog << "nRun=" << nRun << " | iCat=" << iCat << endl;
      }

      // check position
      check_barrel = false;
      ele_Et.clear();
      for(int i=0 ; i<ele_N ; i++) {
	theEle = (TLorentzVector*)(electrons->At(i));
	ele_Et.push_back( theEle->Et() );
	if( fabs(theEle->Eta()) < 1.479 ) check_barrel = true;
      }
      if(!check_barrel) continue;

      //cout << "iCat=" << iCat << endl;      
      //cout << "nEle=" << nEle << endl;

      // looking every TP
      adcTT.clear();
      for(int t=0 ; t<trig_tower_N ; t++) {
	coords = make_pair( trig_tower_ieta[t] , trig_tower_iphi[t] );
	adcTT[coords].first = trig_tower_adc[t];
	//h_adc->Fill(trig_tower_adc[t]);
	//cout << "RCTeta=" << getGCTRegionEta(trig_tower_ieta[t]) 
	//   << " | RCTphi=" << getGCTRegionPhi(trig_tower_iphi[t])
	//   << endl;
      }
      for(int t=0 ; t<trig_tower_N_modif ; t++) {
	coords = make_pair( trig_tower_ieta_modif[t] , trig_tower_iphi_modif[t] );
	iterTT = adcTT.find( coords );
	if( iterTT != adcTT.end() )
	  adcTT[coords].second = trig_tower_adc_modif[t];
	else {
	  outlog << "mapping problem" << endl;
	  adcTT[coords].first = -777;
	  adcTT[coords].second = trig_tower_adc_modif[t];
	}
	//h_adc_modif->Fill(trig_tower_adc_modif[t]);
	//cout << "RCTetaM=" << getGCTRegionEta(trig_tower_ieta_modif[t]) 
	//   << " | RCTphiM=" << getGCTRegionPhi(trig_tower_iphi_modif[t])
	//   << endl;
      }
      // looking every TP

      // det indices of electrons to process (the 4 electrons of highest Et)
      //idxEle.clear();
      //idxEle = det4MaxTab(ele_Et,ele_N);

      // map of electrons firing L1 EGs
      mapEt.clear();

      // loop over barrel ele

      outputRCT = false;
      for(int iEG=0 ; iEG<6 ; iEG++)
	eventN[iEG] = eventM[iEG] = false; // to calculate the rate
      
      for (iEle=0 ; iEle<ele_N ; iEle++) 
	{
	  theEle = (TLorentzVector*) (electrons->At (iEle)) ;

	  // SELECTION
	  selected = true;
	  if (ele_severityLevelSeed[iEle]==5 || ele_severityLevelSeed[iEle]==5) selected = false;
	  bool cutselect = VBTFcuts(asked, 
				    theEle->Pt(), ele_sclEta[iEle], theEle->Eta(), ele_tkSumPt_dr03[iEle], ele_ecalRecHitSumEt_dr03[iEle], 
				    ele_hcalDepth1TowerSumEt_dr03[iEle], ele_hcalDepth2TowerSumEt_dr03[iEle], ele_expected_inner_hits[iEle],
				    ele_deltaphiin[iEle], ele_deltaetain[iEle], ele_he[iEle], ele_sigmaietaieta[iEle],
				    ele_conv_dist[iEle], ele_conv_dcot[iEle], ele_fbrem[iEle], ele_isConversion[iEle]) ;
	  if( (!selected) || (!cutselect) ) continue;
	  // SELECTION

	  IdxCat.clear();
	  IdxSp.clear();
	  IdxPt.clear();

	  // sFGVB threshold category
	  IdxCat.push_back(iCat);
	  IdxCat.push_back(10);

	  //spike-ness of the electron
	  IdxSp = detISP(ele_severityLevelSeed[iEle] , ele_outOfTimeSeed[iEle]);
	  
	  // pt bin 
	  IdxPt.push_back(0);
	  if(theEle->Pt()<20.) IdxPt.push_back(1);
	  else IdxPt.push_back(2);
	  
	  if(debug) cout << "retrieve information" << endl;

	  eta = theEle->Eta();
	  phi = theEle->Phi();
	  iso     = ele_RCTL1iso[iEle];
	  noniso  = ele_RCTL1noniso[iEle];
	  isoM    = ele_RCTL1iso_modif[iEle];
	  nonisoM = ele_RCTL1noniso_modif[iEle];

	  if(debug) cout << "mapping" << endl;
	  
	  // MAPPING

	  // Det L1 EG fired by the ele
	  firedEG_N.clear();
	  firedEG_N.resize(6,0);	  
	  firedEG_M.clear();
	  firedEG_M.resize(6,0);	  

	  matcheck = matcheckM = false;

	  for(int iR=0 ; iR < 10 ; iR++) {
	    if(ele_RCTetVect[iEle][iR]>2) {// interesting regions <=> et > 2GeV
	      
	      fireL1(ele_RCTL1nonisoVect[iEle][iR] , ele_RCTL1isoVect[iEle][iR] ,
		     firedEG_N);
	      
	      fireL1(ele_RCTL1nonisoVect_modif[iEle][iR] , ele_RCTL1isoVect_modif[iEle][iR] ,
		     firedEG_M);
	      
	      // to check
	      if(ele_RCTL1nonisoVect[iEle][iR]>=16 || ele_RCTL1isoVect[iEle][iR]>=16) {
		matcheck=true;
	      }
	      if(ele_RCTL1nonisoVect_modif[iEle][iR]>=16 || ele_RCTL1isoVect_modif[iEle][iR]>=16) {
		matcheckM=true;
	      }
	     
	    }
	  }

	  if(matcheck) {
	    nMatchedCheck[iCat]++ ;
	    if(!matcheckM) nLostCheck[iCat]++ ;
	  }
	  if(matcheckM) {
	    nMatchedCheckM[iCat]++ ;
	    if(!matcheck) nAppearedCheck[iCat]++ ;
	  }

	  // MAP : < <etaR,phiR> ; < iEle , <firedEG_N,firedEG_M> > >
	  // for >1 ele whose main region is (etaR,phiR) : keep only highest Et ele
	  coords = make_pair( ele_RCTeta[iEle] , ele_RCTphi[iEle] );
	  bool severalEleSameReg = false ; 	  
	  iter = mapEt.find( coords );
	  if( iter != mapEt.end() ) {
	    severalEleSameReg = true ;      
	    iEleInMap = (iter->second).first;
	    if(ele_Et[iEle] > ele_Et[iEleInMap])
	      mapEt[coords] = make_pair( iEle , make_pair(firedEG_N,firedEG_M) ) ;
	  }
	  else
	    mapEt[coords] = make_pair( iEle , make_pair(firedEG_N,firedEG_M) ) ;
	  

	  // COUNTING
	  if(debug) cout << "counting" << endl;

	  flag = "";
	  output = false;
	  bool lost[6], appeared[6];
	  bool several = false;
	  for(int iEG=0 ; iEG<6 ; iEG++)
	    lost[iEG]=appeared[iEG]=false;
	  
	  if(debug) cout << "gonna loop over counters" << endl;
	  /*
	    int nEle[11][9][3];    // [sFGVB][spikyness][ptbin]
	    int nMatched[11][9][3];
	    int nMatchedM[11][9][3];
	    int nLost[11][9][3];
	    int nAppeared[11][9][3];
	    int nSwitch[11][9][3]; 
	    int nLower[11][9][3];  
	    int nHigher[11][9][3];
	    int nSeverEleSameReg[11][9][3];
	  */
	  for(int i=0 ; i<IdxCat.size() ; i++) {
	    for(int j=0 ; j<IdxSp.size() ; j++) {
	      for(int k=0 ; k<IdxPt.size() ; k++) {
		if(debug) cout << i << " " << j << " " << k << endl;
		nEle[IdxCat[i]][IdxSp[j]][IdxPt[k]]++ ;
		//cout << nEle[IdxCat[i]][IdxSp[j]][IdxPt[k]] << endl;

		for(int iEG=0 ; iEG<6 ; iEG++) {
		  //cout << iEG << endl;
		  if(firedEG_N[iEG]==1) {
		    nMatched[IdxCat[i]][IdxSp[j]][IdxPt[k]][iEG+1]++;
		    eventN[iEG] = true;
		  }

		  if(firedEG_M[iEG]==1) {
		    nMatchedM[IdxCat[i]][IdxSp[j]][IdxPt[k]][iEG+1]++;
		    eventM[iEG] = true;
		  }
		  
		  if(firedEG_M[iEG]==0 && firedEG_N[iEG]==1) {
		    nLost[IdxCat[i]][IdxSp[j]][IdxPt[k]][iEG+1]++ ;
		    lost[iEG]=true;
		    output=true;
		    outputRCT=true;
		  }
		  if(firedEG_M[i]==1 && firedEG_N[i]==0) {
		    nAppeared[IdxCat[i]][IdxSp[j]][IdxPt[k]][iEG+1]++ ;
		    appeared[iEG]=true;
		    output=true;
		    outputRCT=true;
		  }
		}
		if(iso!=-999 || noniso!=-999) {
		  nMatched[IdxCat[i]][IdxSp[j]][IdxPt[k]][0]++ ;
		  if(isoM==-999 && nonisoM==-999) {
		    nLost[IdxCat[i]][IdxSp[j]][IdxPt[k]][0]++ ;
		    if(ele_severityLevelSeed[iEle]<3) nLostL1NoSpiky++ ;
		  }
		}
		if(isoM!=-999 || nonisoM!=-999) {
		  nMatchedM[IdxCat[i]][IdxSp[j]][IdxPt[k]][0]++ ;
		  if(iso==-999 && noniso==-999)
		    nAppeared[IdxCat[i]][IdxSp[j]][IdxPt[k]][0]++ ;
		}
		if(iso!=-999 && isoM==-999 && noniso==-999 && nonisoM!=-999)
		  nSwitch[IdxCat[i]][IdxSp[j]][IdxPt[k]]++ ;
		else if(iso==-999 && isoM!=-999 && noniso!=-999 && nonisoM==-999)
		  nSwitch[IdxCat[i]][IdxSp[j]][IdxPt[k]]++ ;
		
		if(debug) cout << "ok gonna max" << endl;

		int L1max = max(iso,noniso);
		int L1maxM = max(isoM,nonisoM);
		if(L1max>L1maxM) nLower[IdxCat[i]][IdxSp[j]][IdxPt[k]]++ ;
		else if(L1max<L1maxM) nHigher[IdxCat[i]][IdxSp[j]][IdxPt[k]]++ ;

		if(severalEleSameReg) {
		  nSeverEleSameReg[IdxCat[i]][IdxSp[j]][IdxPt[k]]++ ;
		  flag += " several";
		}
	      }
	    }
	  }
	  if(debug) cout << "has looped" << endl;

	  for(int iEG=0 ; iEG<6 ; iEG++) {
	    if(lost[iEG]) flag += " lostEG"+EGval[iEG+1];
	    if(appeared[iEG]) flag += " appearedEG"+EGval[iEG+1];
	  }

	  // decide output
	  if(debug) cout << "decide output" << endl;
	  if( (iso!=isoM) || (noniso!=nonisoM) ) {
	    output = true;
	    outputRCT = true;
	    //flag += " different";
	  }

	  if( ele_severityLevelSeed[iEle]==4 && (isoM!=-999 || nonisoM!=-999) ) {
	    output = true;
	    outputRCT = true;
	    flag += " spiky-non-killed!";
	  }

	  if(debug) cout << "gonna output if true" << endl;

	  if(output) {
	    
	    // OUTPUT : ELECTRON + ASSOCIATED L1
	    outlog<< "Run #"    << nRun
		  << " Evt #"   << nEvent
		  << "("        << iEvent
		  << ") Ele #"  << iEle
		  << " eta="    << eta
		  << " phi="    << phi
		  << " Et="     << theEle->Et()
		  << " RCTeta=" << ele_RCTeta[iEle]
		  << " RCTphi=" << ele_RCTphi[iEle]
		  << endl
		  << "outOfTime="   << ele_outOfTimeSeed[iEle]
		  << "   severity=" << ele_severityLevelSeed[iEle]
		  << endl
		  << "normal L1 : "
		  << " RCTL1iso="          << iso
		  << " RCTL1noniso="       << noniso
		  << endl
		  << "modif L1 : "
		  << " RCTL1iso="          << isoM
		  << " RCTL1noniso="       << nonisoM
		  << endl 
		  << "firedEG : ";

	    for(int iEG=0 ; iEG<6 ; iEG++)
	      outlog << firedEG_N[iEG] << "   " ; 
		
	    outlog << "firedEGM : ";
	    for(int iEG=0 ; iEG<6 ; iEG++)
	      outlog << firedEG_M[iEG] << "   " ; 
	    
	    outlog << flag 
		 << endl;

	    // OUTPUT : L1 CANDIDATES OF THE EVENT
	    outlog << "L1 Candidates of the event" << endl
		   << "Normal collection" << endl;

 	    outlog << "trig_L1emIso_N="             << trig_L1emIso_N
		   << " | trig_L1emNonIso_N="       << trig_L1emNonIso_N
		   << " | trig_L1emIso_N_modif="    << trig_L1emIso_N_modif
		   << " | trig_L1emNonIso_N_modif=" << trig_L1emNonIso_N_modif
		   << endl;
 	    
	    //for(int i=0 ; i<trig_L1emIso_N ; i++)
	    for(int i=0 ; i<4 ; i++)
	      //if(trig_L1emIso_rank[i] != -999)
		outlog << "iso #"   << i << " : "
		       << "iEta="   << trig_L1emIso_ieta[i]
		       << "   iPhi="<< trig_L1emIso_iphi[i]
		       << "   rank="<< trig_L1emIso_rank[i]
		       << endl;

	    //for(int i=0 ; i<trig_L1emNonIso_N ; i++)
	    for(int i=0 ; i<4 ; i++)
	      //if(trig_L1emNonIso_rank[i] != -999)
		outlog << "noniso #"   << i << " : "
		       << "iEta="   << trig_L1emNonIso_ieta[i]
		       << "   iPhi="<< trig_L1emNonIso_iphi[i]
		       << "   rank="<< trig_L1emNonIso_rank[i]
		       << endl;
	    
	    outlog << "Modif collection" << endl;

	    //for(int i=0 ; i<trig_L1emIso_N_modif ; i++)
	    for(int i=0 ; i<4 ; i++)
	      //if(trig_L1emIso_rank_modif[i] != -999)
		outlog << "isoM #"   << i << " : "
		       << "iEta="   << trig_L1emIso_ieta_modif[i]
		       << "   iPhi="<< trig_L1emIso_iphi_modif[i]
		       << "   rank="<< trig_L1emIso_rank_modif[i]
		       << endl;

	    //for(int i=0 ; i<trig_L1emNonIso_N_modif ; i++)
	    for(int i=0 ; i<4 ; i++)
	      //if(trig_L1emNonIso_rank_modif[i] != -999)
		outlog << "nonisoM #"   << i << " : "
		       << "iEta="   << trig_L1emNonIso_ieta_modif[i]
		       << "   iPhi="<< trig_L1emNonIso_iphi_modif[i]
		       << "   rank="<< trig_L1emNonIso_rank_modif[i]
		       << endl;
	  }
	  //ietaTTele = getEtaTT( ele_RCTeta[iEle] );
	  //iphiTTele = getPhiTT( ele_RCTphi[iEle] );
	    
	}//loop over barrel Eles
	
      // To calc rate
      for(int iEG=0 ; iEG<6 ; iEG++) {
	if(eventN[iEG]) nEvents[iCat][iEG]++ ;
	if(eventM[iEG]) nEventsM[iCat][iEG]++ ;
      }

      if(debug) cout << "has ended the loop over barrel ele" << endl;
      
      // loop over every RCT regions
      if(outputRCT) {
	if(debug) cout << "looping over RCT" << endl;
	for(int iRCTeta=4 ; iRCTeta<=17 ; iRCTeta++) {
	  for(int iRCTphi=0 ; iRCTphi<=17 ; iRCTphi++) {
	    ietaTT = getEtaTT( iRCTeta );
	    iphiTT = getPhiTT( iRCTphi );
	    // loop over TTs of the region
	    outputRegion = false;
	    for (int i=0 ; i<(int)ietaTT.size() ; i++) {
	      for (int j=0 ; j<(int)iphiTT.size() ; j++) {
		coords = make_pair(ietaTT[i],iphiTT[j]);
		adc = adcTT[coords].first;
		adc_modif = adcTT[coords].second;
		
		// OUTPUT : TT
		if(adc!=0 || adc_modif!=0) {
		  outlog << "Reta="    << iRCTeta 
		       << " Rphi="   << iRCTphi   
		       << " Teta=" << ietaTT[i] 
		       << " Tphi=" << iphiTT[j] 
		       << " adc="       << adc       
		       << " adc_modif=" << adc_modif ;
		  if(iRCTeta==ele_RCTeta[iEle] && iRCTphi==ele_RCTphi[iEle])
		    outlog << "   ele-region";
		  if(adc!=adc_modif) outlog << "   !! DIFF !!" << endl;
		  else outlog << endl;
		  outputRegion = true;
		}//endif		
	      }
	    } // end loop over TTs of the region
	    if(outputRegion) outlog << endl;
	  }
	}// end loop over RCT regions
	outlog << "----------------------------------------------------------------------------" << endl;
      } // endif
      if(debug) cout << "ended eventual loop over RCT" << endl
		     << "reading the map" << endl;
      
      // READING THE MAP
      for(iter=mapEt.begin(); iter!=mapEt.end(); iter++) {
	iEle = (iter->second).first ;
	// cout << "iEle=" << iEle << endl;
	firedEG_N = ((iter->second).second).first ;
	firedEG_M = ((iter->second).second).second ;
	
	theEle = (TLorentzVector*) (electrons->At (iEle)) ;
	etaEle = theEle->Eta();
	dr = drM = 0;
	if(ele_RCTL1iso[iEle]>0 || ele_RCTL1noniso[iEle]>0) dr=1;
	if(ele_RCTL1iso_modif[iEle]>0 || ele_RCTL1noniso_modif[iEle]>0) drM=1;

	IdxCat.clear();
	IdxSp.clear();
	IdxPt.clear();
	
	// sFGVB threshold category
	IdxCat.push_back(iCat);
	IdxCat.push_back(10);
	
	//spike-ness of the electron
	IdxSp = detISP(ele_severityLevelSeed[iEle] , ele_outOfTimeSeed[iEle]);
	
	// pt bin 
	IdxPt.push_back(0);
	if(theEle->Pt()<20.) IdxPt.push_back(1);
	else IdxPt.push_back(2);
		
	if(etaEle<1.479) {

	  l1_1 = firedEG_N[0]; 
	  l1_2 = firedEG_N[1];
	  l1_5 = firedEG_N[2];
	  l1_8 = firedEG_N[3];
	  l1_10 = firedEG_N[4];
	  l1_12 = firedEG_N[5];
	  //sc_et = ele_Et[iEle];
	  sc_et = theEle->Et();
	  sc_dr = dr;
	  
	  l1_1_2 = firedEG_M[0]; 
	  l1_2_2 = firedEG_M[1];
	  l1_5_2 = firedEG_M[2];
	  l1_8_2 = firedEG_M[3];
	  l1_10_2 = firedEG_M[4];
	  l1_12_2 = firedEG_M[5];
	  sc_et_2 = theEle->Et();
	  sc_dr_2 = drM;
	  
	  // loop over indices of trees to fill
	  for(int i=0 ; i<IdxCat.size() ; i++) {
	    for(int j=0 ; j<IdxSp.size() ; j++) {
	      for(int k=0 ; k<IdxPt.size() ; k++) {		  
		treenew[IdxCat[i]][IdxSp[j]][IdxPt[k]]->Fill();
		treenew_2[IdxCat[i]][IdxSp[j]][IdxPt[k]]->Fill();
	      }
	    }
	  }// loop indices
	    
	  //nEleFired++ ;
	} 
      }// loop over the map
      if(debug) cout << "looped map" << endl;
      
    }//loop over events

  outlog << "nLostL1NoSpiky=" << nLostL1NoSpiky << endl; 
	    
  /*
    int nEle[11][9][3];    // [sFGVB][spikyness][ptbin]
    int nMatched[11][9][3][7]; // ... + [EGX]
    int nMatchedM[11][9][3][7];
    int nLost[11][9][3][7];
    int nAppeared[11][9][3][7];
    int nSwitch[11][9][3]; 
    int nLower[11][9][3];  results/goodReprocess3/
    int nHigher[11][9][3];
    int nSeverEleSameReg[11][9][3];
  */

  for(int iEG=0 ; iEG<6 ; iEG++) {
    for(int iCat=0 ; iCat<10 ; iCat++) {
      nEvents[10][iEG] += nEvents[iCat][iEG];
      nEventsM[10][iEG] += nEventsM[iCat][iEG];
    }
  }
  
  for(int iPt=0 ; iPt<nPt ; iPt++) {
    ofstream outtab(dirOut+"counters"+ptBin[iPt]+".csv", ios::out);
    outtab << "sFGVB thresh | # ele matched | # Killed | Fraction | "
	   << "# ele matched sev no 3,4 | # Killed | Fraction | "
	   << "# ele matched sev=4 | # Killed | Fraction | "
	   << "# ele matched sev=3,4 | # Killed | Fraction | "
	   << endl;
    for(int iCat=0 ; iCat<nCat ; iCat++) {
	writeTab(outtab , thresh[iCat] , nMatched[iCat][0][iPt][4] , nLost[iCat][0][iPt][4], 
	       nMatched[iCat][5][iPt][4] , nMatched[iCat][3][iPt][4] , nMatched[iCat][1][iPt][4] ,
	       nLost[iCat][5][iPt][4] ,nLost[iCat][3][iPt][4] , nLost[iCat][1][iPt][4]);
    }
  }
  // iSp : All=0  /  No34=5  /  4=3  /  34=1
  /*
    writeTab(ofstream& outtab, int sfgvb, int eleMatched, int eleKilled, 
	     int eleMatchedSevNo34, int eleMatchedSev4, int eleMatchedSev34,
	     int eleKilledSevNo34, int eleKilledSev4, int eleKilledSev34)
  */

  ofstream output1(dirOut+"counters.txt", ios::out);

  output1 << "RATE STUFF" << endl;
  for(int iCat=0 ; iCat<11 ; iCat++) {
    output1 << "Category " << iCat << "(" << thresh[iCat] << ")"
	    << endl;
    output1 << "Normal Collection : ";
    for(int iEG=0 ; iEG<6 ; iEG++) {
      output1 << nEvents[iCat][iEG] << "   " ;
    }
    output1 << endl << "Modified Collection : ";
    for(int iEG=0 ; iEG<6 ; iEG++) {
      output1 << nEventsM[iCat][iEG] << "   " ;
    }
    output1 << endl;
  }
  output1 << endl;
  
  for(int iPt=0 ; iPt<nPt ; iPt++) {
    for(int iSp=0 ; iSp<nSp ; iSp++) {
 
      output1 << endl
	      << "Pt Bin : " << ptBin[iPt] << "   "
	      << "spikyness : " << spName[iSp] << endl;
      
      for(int iCat=0 ; iCat<nCat ; iCat++) {
    
	output1 << "------------------------- CATEGORY " << iCat 
		<< " : sFGVB threshold = "               << thresh[iCat]
		<<" adc -------------------------"       << endl
		<< "nEle="            << nEle[iCat][iSp][iPt]
		<< " | nMatchedEG8="  << nMatched[iCat][iSp][iPt][4]
		<< " | nMatchedMEG8=" << nMatchedM[iCat][iSp][iPt][4]
		<< " | nLostEG8="     << nLost[iCat][iSp][iPt][4]
		<< " | nAppearedEG8=" << nAppeared[iCat][iSp][iPt][4]
		<< endl
		<< "nMatchedCheck="   << nMatchedCheck[iCat]
		<< " | nMatchedCheckM="  << nMatchedCheckM[iCat]
		<< " | nLostCheck="      << nLostCheck[iCat]
		<< " | nAppearedCheck="  << nAppearedCheck[iCat]
		<< endl
		<< "nMatchedEG :";

	for(int i=0 ; i<6 ; i++)
	  output1 << nMatched[iCat][iSp][iPt][i] << " ";
	
	output1 << endl
		<< "nMatchedEGM : " ;
	for(int i=0 ; i<6 ; i++)
	  output1 << nMatchedM[iCat][iSp][iPt][i] << " ";
	
	output1 << endl
		<< "nLostEG : ";
	for(int i=0 ; i<6 ; i++)
	  output1 << nLost[iCat][iSp][iPt][i] << " ";

	output1 << endl
		<< "nAppearedEG : ";
	for(int i=0 ; i<6 ; i++)
	  output1 << nAppeared[iCat][iSp][iPt][i] << " ";

	output1 << endl
		<< "nHigher="         << nHigher[iCat][iSp][iPt]
		<< " | nLower="       << nLower[iCat][iSp][iPt]
		<< " | nSwitch="      << nSwitch[iCat][iSp][iPt]
		<< "  nSeverEleSameReg=" << nSeverEleSameReg[iCat][iSp][iPt]
		<< endl;
      }
      output1 << endl;

    }
  }
  
  // Output file
  TFile *outfile;  
  for(int iCat=0 ; iCat<nCat ; iCat++) {
    for(int iPt=0 ; iPt<nPt ; iPt++) {
      outfile = 
	new TFile( dirOut + baseCategory[iCat]+"/effiReducTree_"+baseCategory[iCat]+ptBin[iPt]+".root" 
		   , "RECREATE");
      
      for(int iSp=0 ; iSp<nSp ; iSp++) {
	treenew[iCat][iSp][iPt]->Write();
	treenew_2[iCat][iSp][iPt]->Write();
      }
      outfile->Close();
    }
  }
}

  

void electrons(int nEntries=-1,
	       TString dirOut="/home/llr/cms/ndaci/SKWork/macro/skEfficiency/18GeV/electrons/results/forgot/",
	       TString dirIn="/data_CMS/cms/ndaci/ndaci_2010B/EGMonitorModifComplete/Skimmed/",
	       TString cuts="WP95",bool canvas=false,bool debug=false)
{
  // Output Log
  ofstream outlog(dirOut+"log.txt",ios::out);

  // Input trees
  TChain * myChain = new TChain ("eIDSimpleTree");
  TString file = "*.root";
  myChain->Add(dirIn+file);

  // Process the tree
  if(debug) cout << "process the tree" << endl;

  TreeToHistos(myChain,dirOut,cuts,nEntries,outlog,debug);
  
}
