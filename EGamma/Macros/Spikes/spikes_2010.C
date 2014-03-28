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
#include "TPaveStats.h"
#include "TStyle.h"
#include "TSystem.h"
//#include "treesList.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "../Common/readJSONFile.cc"

// DATA TREE LISTS
// include here your files treelist.h
// #include "[path]/treelist.h"                                                                                                          

using namespace std;
typedef map< pair<int,int> , pair<int,int> > MAPTT;

typedef map< pair<int,int> , pair< pair<int,int>,pair<int,int> > > MAPTTS;
//           TT(ieta,iphi)           (adc,adcM)  (sFGVB,sFGVB_M)
// adc = mapTT[coords].first.first
// adcM = mapTT[coords].first.second
// sFGVB = maTT[coords].second.first
// sFGVB_M = maTT[coords].second.second

typedef map< pair<int,int> , pair< int , pair< vector<int> , vector<int> > > > MAPPY ;
// { <RCT ieta,RCT iphi> , < iEle , <[eg1,eg2,eg5,eg8,eg10,eg12]N ,[...]M>  > }

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

void TreeToHistos(TChain* myChain, TString zeroThresh, TString dirOut, ofstream& outlog, int nEntries, TString halftag, bool spikeCut, bool debug=false)  {
  
  ostringstream ossi;

  int nCat = 13;
  int iCat;
  TString category[13]={"_t10","_t15","_t17","_t19","_t21","_t23","_t30","_t35","_t40","_t50","_tAll","_t148524","_tRest"} ;
  //int thresh[10] = {10,15,17,19,21,23,30,35,40,50} ;
  //string EGval[6] = {"1","2","5","8","10","12"} ;

  // get the tree
  int nEvent, nRun, nLumi, trig_isUnbiased, trig_isL1SingleEG2, trig_isL1SingleEG5, trig_isL1SingleEG8;   

  // Masked TT
  int trig_nMaskedCh;
  int trig_iMaskedTTeta[4032], trig_iMaskedTTphi[4032];
 
  // TP info
  int trig_tower_N,trig_tower_ieta[4032],trig_tower_iphi[4032],trig_tower_adc[4032],trig_tower_sFGVB[4032]; 
  int trig_tower_N_modif,trig_tower_ieta_modif[4032],trig_tower_iphi_modif[4032],trig_tower_adc_modif[4032],trig_tower_sFGVB_modif[4032]; 
  int trig_L1emIso_N, trig_L1emNonIso_N, trig_L1emIso_N_modif, trig_L1emNonIso_N_modif;
  int trig_tower_N_emul,trig_tower_ieta_emul[4032],trig_tower_iphi_emul[4032],trig_tower_adc_emul[4032][5],trig_tower_sFGVB_emul[4032][5];

  // L1 candidates info
  int trig_L1emIso_ieta[4], trig_L1emIso_iphi[4], trig_L1emIso_rank[4];
  int trig_L1emNonIso_ieta[4], trig_L1emNonIso_iphi[4], trig_L1emNonIso_rank[4];
  int trig_L1emIso_ieta_modif[4], trig_L1emIso_iphi_modif[4], trig_L1emIso_rank_modif[4];
  int trig_L1emNonIso_ieta_modif[4], trig_L1emNonIso_iphi_modif[4], trig_L1emNonIso_rank_modif[4];

  // Spikes
  int spike_N,spike_TTieta[5000], spike_TTiphi[5000], spike_Rieta[5000], spike_Riphi[5000], spike_severityLevel[5000], spike_outOfTime[5000];
  double  spike_Et[5000], spike_eta[5000], spike_phi[5000], spike_theta[5000];

  // Global
  myChain->SetBranchAddress("nEvent",&nEvent);
  myChain->SetBranchAddress("nRun",&nRun);
  myChain->SetBranchAddress("nLumi",&nLumi);
  myChain->SetBranchAddress("trig_isUnbiased",&trig_isUnbiased);
  myChain->SetBranchAddress("trig_isL1SingleEG2",&trig_isL1SingleEG2);
  myChain->SetBranchAddress("trig_isL1SingleEG5",&trig_isL1SingleEG5);
  myChain->SetBranchAddress("trig_isL1SingleEG8",&trig_isL1SingleEG8);

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

  // masked
//   myChain->SetBranchAddress("trig_nMaskedRCT", &trig_nMaskedRCT);
//   myChain->SetBranchAddress("trig_iMaskedRCTeta", &trig_iMaskedRCTeta);
//   myChain->SetBranchAddress("trig_iMaskedRCTcrate", &trig_iMaskedRCTcrate);
//   myChain->SetBranchAddress("trig_iMaskedRCTphi", &trig_iMaskedRCTphi);
  myChain->SetBranchAddress("trig_nMaskedCh", &trig_nMaskedCh);
  myChain->SetBranchAddress("trig_iMaskedTTeta", &trig_iMaskedTTeta);
  myChain->SetBranchAddress("trig_iMaskedTTphi", &trig_iMaskedTTphi);

  // normal collection
  myChain->SetBranchAddress("trig_tower_N", &trig_tower_N);
  myChain->SetBranchAddress("trig_tower_ieta",  &trig_tower_ieta);
  myChain->SetBranchAddress("trig_tower_iphi",  &trig_tower_iphi);
  myChain->SetBranchAddress("trig_tower_adc",  &trig_tower_adc);
  myChain->SetBranchAddress("trig_tower_sFGVB",  &trig_tower_sFGVB);  
  
  // modified collection
  myChain->SetBranchAddress("trig_tower_N_modif", &trig_tower_N_modif);
  myChain->SetBranchAddress("trig_tower_ieta_modif",  &trig_tower_ieta_modif);
  myChain->SetBranchAddress("trig_tower_iphi_modif",  &trig_tower_iphi_modif);
  myChain->SetBranchAddress("trig_tower_adc_modif",  &trig_tower_adc_modif);
  myChain->SetBranchAddress("trig_tower_sFGVB_modif",  &trig_tower_sFGVB_modif);  

  myChain->SetBranchAddress("trig_tower_N_emul", &trig_tower_N_emul);
  myChain->SetBranchAddress("trig_tower_ieta_emul",  &trig_tower_ieta_emul);
  myChain->SetBranchAddress("trig_tower_iphi_emul",  &trig_tower_iphi_emul);
  myChain->SetBranchAddress("trig_tower_adc_emul",  &trig_tower_adc_emul);
  myChain->SetBranchAddress("trig_tower_sFGVB_emul",  &trig_tower_sFGVB_emul);

  // Spikes
  myChain->SetBranchAddress("spike_N",&spike_N);
  myChain->SetBranchAddress("spike_TTieta",&spike_TTieta);
  myChain->SetBranchAddress("spike_TTiphi",&spike_TTiphi);
  myChain->SetBranchAddress("spike_Rieta",&spike_Rieta);
  myChain->SetBranchAddress("spike_Riphi",&spike_Riphi);
  myChain->SetBranchAddress("spike_severityLevel",&spike_severityLevel);
  myChain->SetBranchAddress("spike_outOfTime",&spike_outOfTime);
  myChain->SetBranchAddress("spike_Et",&spike_Et);
  myChain->SetBranchAddress("spike_eta",&spike_eta);
  myChain->SetBranchAddress("spike_phi",&spike_phi);
  myChain->SetBranchAddress("spike_theta",&spike_theta);

  TString n_coll[2] = {"Normal","Modif."};
  TString n_coll_1[2] = {"_N","_M"};
  TString n_match[3] = {"_all","_L1","_L1_EG8"};

  TH1F* h_adc[6][13];
  TH1F* h_sp[6][13];
  TH2F* h_sptp[6][13];
  TString n_histo[6][13];

  int nHistos = 6;

  for(int i=0 ; i<nHistos ; i++) {
    for(iCat=0 ; iCat<nCat ; iCat++) {
      n_histo[i][iCat] = zeroThresh + n_coll_1[i/3] + n_match[i%3] + category[iCat];

      h_adc[i][iCat] = new TH1F("h_adc_"+n_histo[i][iCat] , "h_adc_"+n_histo[i][iCat] , 128 , 0 , 64);
      h_adc[i][iCat]->SetXTitle("TP value (GeV)");
      h_adc[i][iCat]->SetYTitle("N_{TP}");
      h_adc[i][iCat]->SetTitle(n_coll[i/3]+" collection TP value");
     
      h_sp[i][iCat] = new TH1F("h_sp_"+n_histo[i][iCat] , "h_sp_"+n_histo[i][iCat] , 800 , 0 , 1600);
      h_sp[i][iCat]->SetXTitle("Spike E_{T} (GeV)");
      h_sp[i][iCat]->SetYTitle("N_{TP}");
      h_sp[i][iCat]->SetTitle(n_coll[i/3]+" collection TP value");
     
      h_sptp[i][iCat] = new TH2F("h_sptp_"+n_histo[i][iCat] , "h_sptp_"+n_histo[i][iCat] , 128 , 0 , 64 , 800 , 0 , 1600);
      h_sptp[i][iCat]->SetXTitle("TP value (GeV)");
      h_sptp[i][iCat]->SetYTitle("Spike E_{T} (GeV)");
      h_sptp[i][iCat]->SetTitle("Spike E_{T} vs TP value (" + n_coll[i/3] + " coll.)");
    }
  }

  bool isGoodRun,output,saturates,wellId,missed,eliminated,zeroedTT,trigEG8,trigEG8M;
  string flag;
  vector<int> firedEG_N; // EG 1/2/5/8/10/12 
  vector<int> firedEG_M; 
  // loop over events

  vector<int> matchIso, matchNonIso, matchIsoM, matchNonIsoM;
  int nSpikyL1[13][7], nLostSpikyL1[13][7], nSpikes[13], nSpikyL1sat[13], nSpikesWellId[13], nDiffSFGVB[13], nMatchedTT[13], nLostTT[13];

  for(int i=0 ; i<nCat ; i++) {
    nSpikes[i] = nSpikyL1sat[i] = nSpikesWellId[i] = nDiffSFGVB[i] 
      = nMatchedTT[i] = nLostTT[i] = 0;
    for(int itg=0 ; itg<7 ; itg++) {
      nSpikyL1[i][itg] = nLostSpikyL1[i][itg] = 0 ;
    }
  }
  
  vector<int> IdxCat;
  vector< pair<int,int> > maskedTT;

  MAPTTS adcTT;
  MAPTTS::iterator iterTT;
  pair<int,int> coords;

  TString trigname[7] = {"2","5","8","10","12","15","20"} ;
  int trigthresh[7] = {4,10,16,20,24,30,40} ;
  bool triggEG[7], triggEGM[7];
  for(int i=0 ; i<7 ; i++) {
    triggEG[i] = false ;
    triggEGM[i] = false ;
  }

  //ofstream ttfile("ttfile.txt",ios::out);

  // -------------------------------------------------------------------------------
  // JSON FILE READER
  // -------------------------------------------------------------------------------
  // define map of run/LS
  string jsonFile = "/home/llr/cms/ochando/Analysis/Ana_Multi/json/goodrunlist_json.txt";
  map<int, vector<pair<int, int> > > jsonMap = readJSONFile(jsonFile);   
  
  if(debug) cout << "gonna loop over events" << endl;
  TString filename;

  int nCurrentRun;

  int numEntries = myChain->GetEntries () ;
  outlog << "numEntries=" << numEntries << endl;

  int nProcess = numEntries;
  if(nEntries>=0 && nEntries<numEntries)
    nProcess = nEntries;

  // loop over events
  for (int iEvent = 0 ; iEvent < nProcess ; iEvent++ )
    { 
      //if(iEvent%500 != 0) continue;
      myChain->GetEntry (iEvent) ;

      // run selection
      isGoodRun = AcceptEventByRunAndLumiSection(nRun, nLumi, jsonMap);
      if(!isGoodRun) continue;
      
      // HLT selection
      //if(trig_isL1SingleEG8 != 1)

      if(debug) cout << "passed json" << endl;

      // show which file is being processed
      if(iEvent==0) {
        filename = myChain->GetFile()->GetName() ;
        outlog << "File : " << filename << endl;
      }
      else if( filename != myChain->GetFile()->GetName() ) {
        filename = myChain->GetFile()->GetName() ;
        outlog << "File : " << myChain->GetFile()->GetName() << endl;
      }

      // sFGVB  threshold category
      iCat = threshCateg(nRun);
      if(iCat<0 || iCat>9) continue;

      // indices
      IdxCat.clear();
      IdxCat.push_back(iCat);
      IdxCat.push_back(10);
      if(iCat==1) {
	if(nRun==148524) IdxCat.push_back(11);
	else IdxCat.push_back(12);
      }

      // show which run/category is being processed
      if(iEvent==0) {
	nCurrentRun = nRun ;
	outlog << "nRun=" << nRun << " | iCat=" << iCat << endl;
      }
      else if(nRun!=nCurrentRun) {
	nCurrentRun=nRun ;
	outlog << "nRun=" << nRun << " | iCat=" << iCat << endl;
      }
      //cout << nRun << " " << nEvent << endl;

      // map the towers
      adcTT.clear();
      for(int t=0 ; t<trig_tower_N ; t++) {
	coords = make_pair( trig_tower_ieta[t] , trig_tower_iphi[t] );
	adcTT[coords].first.first = trig_tower_adc[t];
	adcTT[coords].second.first = trig_tower_sFGVB[t];
	if(trig_tower_ieta[t]>-17 && trig_tower_ieta[t]<18)
	  for(int i=0 ; i<IdxCat.size() ; i++)
	    h_adc[0][IdxCat[i]]->Fill(0.25*trig_tower_adc[t]);
	/*ttfile << "eta="    << trig_tower_ieta[t]
	       << "   phi=" << trig_tower_iphi[t]
	  //<< "   adc=" <<
	       << "   sFGVB=" << trig_tower_sFGVB[t]
	       << endl;*/
      }
      for(int t=0 ; t<trig_tower_N_emul ; t++) {
	coords = make_pair( trig_tower_ieta_emul[t] , trig_tower_iphi_emul[t] );
	iterTT = adcTT.find( coords );
	if( iterTT != adcTT.end() ) {
	  adcTT[coords].first.second = trig_tower_adc_emul[t][2];
	  adcTT[coords].second.second = trig_tower_sFGVB_emul[t][2];
	  if(trig_tower_ieta_emul[t]>-17 && trig_tower_ieta_emul[t]<18) 
	    for(int i=0 ; i<IdxCat.size() ; i++)
	      h_adc[3][iCat]->Fill(0.25*trig_tower_adc_emul[t][2]);
	  /*ttfile << "etaM="    << trig_tower_ieta_modif[t]
		 << "   phiM=" << trig_tower_iphi_modif[t]
		 << "   sFGVBM=" << trig_tower_sFGVB_modif[t]
		 << endl;*/
	}
	else {
	  outlog << "mapping problem" << endl;
	  adcTT[coords].first.first = -777;
	  adcTT[coords].first.second = trig_tower_adc_emul[t][2];
	  adcTT[coords].second.first = -777;
	  adcTT[coords].second.second = trig_tower_sFGVB_emul[t][2];
	  
	}
      }
//       for(int t=0 ; t<trig_tower_N_modif ; t++) {
// 	coords = make_pair( trig_tower_ieta_modif[t] , trig_tower_iphi_modif[t] );
// 	iterTT = adcTT.find( coords );
// 	if( iterTT != adcTT.end() ) {
// 	  adcTT[coords].first.second = trig_tower_adc_modif[t];
// 	  adcTT[coords].second.second = trig_tower_sFGVB_modif[t];
// 	  if(trig_tower_ieta_modif[t]>-17 && trig_tower_ieta_modif[t]<18) 
// 	    for(int i=0 ; i<IdxCat.size() ; i++)
// 	      h_adc[3][iCat]->Fill(0.25*trig_tower_adc_modif[t]);
// 	  /*ttfile << "etaM="    << trig_tower_ieta_modif[t]
// 		 << "   phiM=" << trig_tower_iphi_modif[t]
// 		 << "   sFGVBM=" << trig_tower_sFGVB_modif[t]
// 		 << endl;*/
// 	}
// 	else {
// 	  outlog << "mapping problem" << endl;
// 	  adcTT[coords].first.first = -777;
// 	  adcTT[coords].first.second = trig_tower_adc_modif[t];
// 	  adcTT[coords].second.first = -777;
// 	  adcTT[coords].second.second = trig_tower_sFGVB_modif[t];
	  
// 	}
//       }

      // get masked towers
      maskedTT.clear();
      for(int iTT=0 ; iTT<trig_nMaskedCh ; iTT++) {
	maskedTT.push_back( make_pair( trig_iMaskedTTeta[iTT] , trig_iMaskedTTphi[iTT] ) );
      }

      // loop over evil rechits     
      if(debug) cout << "loop over evil rechits" << endl;
      saturates = false;
      for(int irh=0 ; irh<spike_N ; irh++) {
	if(spike_severityLevel[irh]==4 || spike_severityLevel[irh]==4) { // (kTime) + kWeird

	  // selection on the energy
	  if(spikeCut && spike_Et[irh]<8) continue;
	  
// 	  if(iCat==1 || iCat==4 || iCat==7 || iCat==9) 
// 	    if(spike_Et[irh]<11) continue;
	  
	  for(int i=0 ; i<IdxCat.size() ; i++) {
	    nSpikes[IdxCat[i]] ++ ;
	  }
	  
	  matchIso.clear();
	  matchNonIso.clear();
	  matchIsoM.clear();
	  matchNonIsoM.clear();
	  
	  if(spike_Et[irh]>0) {

	    output=false;
	    saturates = false;
	    for(int itg=0 ; itg<7 ; itg++) triggEG[itg] = triggEGM[itg] = false;
	    
	    for(int ic=0 ; ic<4 ; ic++) {
	      if(trig_L1emIso_ieta[ic]==spike_Rieta[irh] && trig_L1emIso_iphi[ic]==spike_Riphi[irh]) {
		output = true;
		matchIso.push_back(ic);
		if(trig_L1emIso_rank[ic]>=63) saturates = true;
		for(int itg=0 ; itg<7 ; itg++)
		  if(trig_L1emIso_rank[ic]>=trigthresh[itg]) triggEG[itg] = true;		
	      }
	      if(trig_L1emNonIso_ieta[ic]==spike_Rieta[irh] && trig_L1emNonIso_iphi[ic]==spike_Riphi[irh]) {
		output = true;
		matchNonIso.push_back(ic);
		if(trig_L1emNonIso_rank[ic]>=63) saturates = true;
		for(int itg=0 ; itg<7 ; itg++)
		  if(trig_L1emNonIso_rank[ic]>=trigthresh[itg]) triggEG[itg] = true;		
	      }
	      if(trig_L1emIso_ieta_modif[ic]==spike_Rieta[irh] && trig_L1emIso_iphi_modif[ic]==spike_Riphi[irh]) {
		output = true;
		matchIsoM.push_back(ic);
		if(trig_L1emIso_rank_modif[ic]>=63) saturates = true;
		for(int itg=0 ; itg<7 ; itg++)
		  if(trig_L1emIso_rank_modif[ic]>=trigthresh[itg]) triggEGM[itg] = true;		
	      }
	      if(trig_L1emNonIso_ieta_modif[ic]==spike_Rieta[irh] && trig_L1emNonIso_iphi_modif[ic]==spike_Riphi[irh]) {
		output = true;
		matchNonIsoM.push_back(ic);
		if(trig_L1emNonIso_rank_modif[ic]>=63) saturates = true;
		for(int itg=0 ; itg<7 ; itg++)
		  if(trig_L1emNonIso_rank_modif[ic]>=trigthresh[itg]) triggEGM[itg] = true;		
	      }
	    }

	    coords = make_pair(spike_TTieta[irh] , spike_TTiphi[irh]);

	    bool b_nuage = false;
	    b_nuage = (spike_Et[irh]>600)&&(adcTT[coords].first.first>50)&&(adcTT[coords].first.first<150) ;
	    bool b_zero = false;
	    b_zero =  (adcTT[coords].first.first < 10) && (spike_Et[irh] > 5.) ;
	    
	    bool b_mask = false;
	    for(int iTT=0 ; iTT<maskedTT.size() ; iTT++) {
	      if( maskedTT[iTT] == coords )
		b_mask = true;
	    }


	    if(output) {
	      outlog << "Spike   nRun=" << nRun << "   iCat=" << iCat << "   nEvent=" << nEvent << "("<<iEvent<<")" << endl
		     << "Rieta="         << spike_Rieta[irh]
		     << "   Riphi="      << spike_Riphi[irh]
		     << "   TTieta="     << spike_TTieta[irh]
		     << "   TTiphi="     << spike_TTiphi[irh]
		//<< "   outOfTime="  << spike_outOfTime[irh]
		     << "   severity="   << spike_severityLevel[irh]
		     << "   Et="         << spike_Et[irh];
	      
	      if(b_nuage) outlog << " | nuage ";
	      if(b_zero) outlog << " | zeroblem ";
	      if(b_mask) outlog << " | masked" ;
	      outlog << endl;
	      
	      if(matchIso.size()>0)
		for(int ic=0 ; ic<matchIso.size() ; ic++)
		  outlog << "Iso rank=" << trig_L1emIso_rank[matchIso[ic]] << "      ";
	      if(matchNonIso.size()>0)
		for(int ic=0 ; ic<matchNonIso.size() ; ic++)
		  outlog << "NonIso rank=" << trig_L1emNonIso_rank[matchNonIso[ic]] << "      ";
	      if(matchIsoM.size()>0)
		for(int ic=0 ; ic<matchIsoM.size() ; ic++)
		  outlog << "IsoM rank=" << trig_L1emIso_rank_modif[matchIsoM[ic]] << "      ";
	      if(matchNonIsoM.size()>0)
		for(int ic=0 ; ic<matchNonIsoM.size() ; ic++)
		  outlog << "NonIsoM rank=" << trig_L1emNonIso_rank_modif[matchNonIsoM[ic]] ;
	      
	      outlog<<endl;
	      
	      outlog << "Trigger Tower : adc=" << adcTT[coords].first.first
		     << "   adcM="             << adcTT[coords].first.second
		     << "   sFGVB="            << adcTT[coords].second.first
		     << "   sFGVB_M="          << adcTT[coords].second.second ;
	      if(adcTT[coords].second.first != adcTT[coords].second.second) {
		for(int i=0 ; i<IdxCat.size() ; i++)
		  nDiffSFGVB[IdxCat[i]]++ ;
		outlog << "   !! sFGVB diff !!" ;
	      }
	      outlog << endl;
	      
	      wellId = zeroedTT = eliminated = missed = false;

	      // loop over categories
	      for(int i=0 ; i<IdxCat.size() ; i++) {
		
		// non L1-matching-dependant properties
		if(adcTT[coords].second.second==0) {
		  wellId = true;
		  nSpikesWellId[IdxCat[i]]++ ;
		}
		if(adcTT[coords].first.first>0) {
		  nMatchedTT[IdxCat[i]]++;
		  if(adcTT[coords].first.second==0) {
		    nLostTT[IdxCat[i]]++;
		    zeroedTT = true;
		  }
		}		  
		h_sp[0][IdxCat[i]]->Fill(spike_Et[irh]);
		h_sp[3][IdxCat[i]]->Fill(spike_Et[irh]);
		h_sptp[0][IdxCat[i]]->Fill(0.25*adcTT[coords].first.first , spike_Et[irh]); // spEt vs tp (N_all)
		h_sptp[3][IdxCat[i]]->Fill(0.25*adcTT[coords].first.second , spike_Et[irh]); // spEt vs tp (M_all)
		// non L1-matching
		
		// L1-matching-dependant properties
		if( matchIso.size()>0 || matchNonIso.size()>0 ) {
		  
		  h_sp[1][IdxCat[i]]->Fill(spike_Et[irh]);
		  h_adc[1][IdxCat[i]]->Fill(0.25*adcTT[coords].first.first); // tp(spike) N_L1
		  h_sptp[1][IdxCat[i]]->Fill(0.25*adcTT[coords].first.first , spike_Et[irh]); // spEt vs tp (N_L1)
		  
		  if(saturates) {
		    //saturates = true;
		    nSpikyL1sat[IdxCat[i]]++ ;
		  }
		  
		} // L1-matching
		
		// EGX
		for(int itg=0 ; itg<7 ; itg++) {
		  if( triggEG[itg] ) {
		    nSpikyL1[IdxCat[i]][itg]++ ;
		    
		    if(itg==2) {
		      h_sp[2][IdxCat[i]]->Fill(spike_Et[irh]);
		      h_adc[2][IdxCat[i]]->Fill(0.25*adcTT[coords].first.first); // tp(spike) N_L1_EG8
		      h_sptp[2][IdxCat[i]]->Fill(0.25*adcTT[coords].first.first , spike_Et[irh]); // spEt vs tp (N_L1_EG8)
		    }
		    if( !triggEGM[itg] ) {
		      nLostSpikyL1[IdxCat[i]][itg]++ ;
		      eliminated = true;
		    }
		  } // EGX-L1-matching
		}
		// Modified L1-matching
		if( matchIsoM.size()>0 || matchNonIsoM.size()>0 ) {
		  h_sp[4][IdxCat[i]]->Fill(spike_Et[irh]);
		  h_adc[4][IdxCat[i]]->Fill(0.25*adcTT[coords].first.second); // tp(spike) M_L1
		  h_sptp[4][IdxCat[i]]->Fill(0.25*adcTT[coords].first.second , spike_Et[irh]); // spEt vs tp (M_L1)
		}
		
		// EG8
		if(triggEGM[2]) {
		  h_sp[5][IdxCat[i]]->Fill(spike_Et[irh]);
		  h_adc[5][IdxCat[i]]->Fill(0.25*adcTT[coords].first.second); // tp(spike) M_L1_EG8
		  h_sptp[5][IdxCat[i]]->Fill(0.25*adcTT[coords].first.second , spike_Et[irh]); // spEt vs tp (M_L1_EG8)
		  missed = true;
		}

	      } // end loop over categories

	      if(wellId) outlog << " wellId";
	      if(zeroedTT) outlog << " zeroedTT" ;
	      if(eliminated) outlog << " eliminatedEG8";
	      if(missed) outlog << " missedEG8";
	      if(saturates) outlog << " saturates";
	      outlog << endl ;
	      
	      // show all the trigger towers of the region
	      if(b_zero || b_nuage || b_mask) {
		outlog << "Trigger towers of the region" << endl;
		vector<int> ietaTT = getEtaTT( spike_Rieta[irh] );
		vector<int> iphiTT = getPhiTT( spike_Riphi[irh] );
		pair<int,int> etaphi;
		for(int i=0 ; i<ietaTT.size() ; i++) {
		  for(int j=0 ; j<iphiTT.size() ; j++) {
		    etaphi = make_pair( ietaTT[i] , iphiTT[j] );
		    if( adcTT[etaphi].first.first != 0 ) {
		      outlog << "ietaTT="    << ietaTT[i]
			     << "   iphiTT=" << iphiTT[j]
			     << "   adc="    << adcTT[etaphi].first.first
			     << "   adcM="   << adcTT[etaphi].first.second
			     << endl;
		    }
		  }
		}
	      }
	      outlog << endl;

	    } // endif output
	  } //endif
	}//endif
      } // loop over rechits 
    
    }//loop over events  


  // show counters
  outlog << "COUNTING RESULTS" << endl;
  for(int i=0 ; i<nCat ; i++) {
    outlog << category[i] 
	   << " : nSpikes="      << nSpikes[i] 
	   << " | nSpikesWellId="<< nSpikesWellId[i]
	   << " | nSpikyL1EG8="  << nSpikyL1[i][2]
	   << " | nSpikyL1sat="  << nSpikyL1sat[i]
	   << endl
	   << " | nLostSpikyL1EG8=" << nLostSpikyL1[i][2] 
	   << " | nMatchedTT="      << nMatchedTT[i]
	   << " | nLostTT="         << nLostTT[i]
	   << " | nDiffSFGVB="      << nDiffSFGVB[i]
	   << endl;
  }

  ofstream outcount(dirOut+"counters"+halftag+".csv",ios::out);
  for(int itg=0 ; itg<7 ; itg++) {
    for(int i=0 ; i<nCat ; i++) {
      outcount << category[i] 
	       << " | " << nSpikes[i] 
	       << " | " << nSpikesWellId[i]
	       << " | " << nSpikyL1[i][itg]
	       << " | " << nLostSpikyL1[i][itg]
	       << " | " << nMatchedTT[i]
	       << " | " << nLostTT[i]
	       << " | " << nDiffSFGVB[i]
	       << " | " << nSpikyL1sat[i] 
	       << endl;
    }
    outcount << endl << endl;
  }

  // plots
  if(debug) cout << "is gonna plot in file" << endl;

  TFile *outplot = new TFile(dirOut+"spike_plots"+halftag+".root","RECREATE");
  for(iCat=0 ; iCat<nCat ; iCat++) {
    for(int i=0 ; i<nHistos ; i++) {
      h_adc[i][iCat]->Write();
      h_sp[i][iCat]->Write();     
      h_sptp[i][iCat]->Write();
    }
  }
  outplot->Close();

  if(debug) cout << "plotted in file" << endl;

  gStyle->SetPalette(1);

  if(debug) cout << "is gonna plot" << endl;

  TCanvas *c_sptp[6][13];
  TCanvas *c_spikeplots[6][13];
  TCanvas *c_adcplots[4][13];
  int cIdx[8] = { 0,3,0,3,1,4,2,5 };  
  int i1,i2;

//   for(iCat=0 ; iCat<nCat ; iCat++) {

//     // SPTP PLOTS
//     for(int i=0 ; i<6 ; i++) {
//       c_sptp[i][iCat]= new TCanvas("c_sptp_"+n_histo[i][iCat] , "c_sptp_"+n_histo[i][iCat] , 0,0,1200,800);
//       c_sptp[i][iCat]->SetLogy();
//       h_sptp[i][iCat]->Draw("colz");
//       c_sptp[i][iCat]->Update();
//       TPaveStats *st = (TPaveStats*)h_sptp[i][iCat]->FindObject("stats");
//       st->SetX1NDC(0.75); //new x start position.
//       st->SetX2NDC(0.9); //new x end position
//       st->SetY1NDC(0.85); //new y start position
//       st->SetY2NDC(1.); //new y end position
//       st->SetOptStat(111);
//       c_sptp[i][iCat]->Modified();
//       c_sptp[i][iCat]->Update();
//       c_sptp[i][iCat]->Print(dirOut+"sptp"+n_histo[i][iCat]+halftag+".gif","gif");
//     }

//     // ADC PLOTS    
//     for(int i=0 ; i<4 ; i++) {
//       i1 = cIdx[i];
//       i2 = cIdx[i+4];

//       c_adcplots[i][iCat] = new TCanvas( "c_adc_"+n_histo[i2][iCat] , "c_adc_"+n_histo[i2][iCat] ,0,0,1200,800);
//       c_adcplots[i][iCat]->SetLogy();
      
//       h_adc[i1][iCat]->SetLineColor(kBlack);
//       h_adc[i1][iCat]->SetMarkerColor(kBlack);
//       h_adc[i1][iCat]->SetMarkerStyle(21);
//       h_adc[i1][iCat]->SetMinimum(1);
//       h_adc[i1][iCat]->Draw();
      
//       h_adc[i2][iCat]->SetLineColor(kRed);
//       h_adc[i2][iCat]->SetMarkerColor(kRed);
//       h_adc[i2][iCat]->SetMarkerStyle(22);
//       h_adc[i2][iCat]->SetMinimum(1);
//       h_adc[i2][iCat]->Draw("same");
      
//       c_adcplots[i][iCat]->Print(dirOut+"adcplots_"+n_histo[i2][iCat]+halftag+".gif","gif");  
//     }
    
//     // SPIKE PLOTS    
//     for(int i=0 ; i<4 ; i++) {
//       i1 = cIdx[i];
//       i2 = cIdx[i+4];

//       c_spikeplots[i][iCat] = new TCanvas( "c_spike_"+n_histo[i2][iCat] , "c_spike_"+n_histo[i2][iCat] ,0,0,1200,800);
//       c_spikeplots[i][iCat]->SetLogy();
      
//       h_sp[i1][iCat]->SetLineColor(kBlack);
//       h_sp[i1][iCat]->SetMarkerColor(kBlack);
//       h_sp[i1][iCat]->SetMarkerStyle(21);
//       h_sp[i1][iCat]->SetMinimum(1);
//       h_sp[i1][iCat]->Draw();
      
//       h_sp[i2][iCat]->SetLineColor(kRed);
//       h_sp[i2][iCat]->SetMarkerColor(kRed);
//       h_sp[i2][iCat]->SetMarkerStyle(22);
//       h_sp[i2][iCat]->SetMinimum(1);
//       h_sp[i2][iCat]->Draw("same");
      
//       c_spikeplots[i][iCat]->Print(dirOut+"spikesET_"+n_histo[i2][iCat]+halftag+".gif","gif");  
//     }
//   }
  
    }
  

int spikes(int iHalf=0, int nHalf=1, TString zeroThresh="18GeV", int nEntries=-1, 
	   TString dirOut="/home/llr/cms/ndaci/SKWork/macro/skEfficiency/18GeV/spikes/results/testeur/",
	   bool spikeCut=true, bool debug=false)
{

  // def output tag according to half
  ostringstream ossi;
  TString halftag = "_half_";
  ossi.str("");
  ossi << iHalf;
  halftag += ossi.str();
  ossi.str("");
  ossi << nHalf;
  halftag += "_" + ossi.str();
 
  // Output Log
  ofstream outlog(dirOut+"log"+halftag+".txt",ios::out);

  // Input trees
  TString dirIn;
  vector<TString> listTrees;

  if(zeroThresh=="18GeV") {
    dirIn = "/data_CMS/cms/ndaci/ndaci_2010B/EGMonitorModifSpikes2/";
    listTrees = listFilesSpikes2(dirIn);
  }
  else if(zeroThresh=="8GeV") {    
    dirIn = "/data_CMS/cms/ndaci/ndaci_2010B/EGMonitorModif8GeV/";
    listTrees = listFiles8GeV(dirIn);
  }
  else if(zeroThresh=="12GeV") {
    dirIn = "/data_CMS/cms/ndaci/ndaci_2010B/EGMonitorModif12GeV3/";
    listTrees = listFiles12GeV(dirIn);
  }
  else if(zeroThresh=="8G17a") {
    listTrees = list_8G17a();
  }
  else if(zeroThresh=="8G35a") {
    listTrees = list_8G35a();
  }
  else if(zeroThresh=="3517") {
    listTrees = list_8G3517();
  }
  else if(zeroThresh=="1517") {
    listTrees = list_1517();
  }
  else if(zeroThresh=="2117") {
    listTrees = list_2117();
  }
  else if(zeroThresh=="5017") {
    listTrees = list_5017();
  }
  else {
    cout << "problem of dataset tag" << endl;
    return 0;
  }
  
  if(listTrees.size()==0) {
    cout << "tree list empty" << endl;
    return 0;
  }
  
  if(iHalf<0 || nHalf<1 || iHalf>=nHalf) {
    cout << "problem with iHalf and nHalf" << endl;
    return 0;
  }

  if(nHalf>listTrees.size()) {
    cout << "requested more halfes than available trees ; set nHalf to listTrees.size()=" << listTrees.size() << endl;
    nHalf = listTrees.size();
  }

  int nTrees, nStart, nEnd;
  nTrees = listTrees.size() / nHalf;
  if(iHalf < nHalf-1) {
    nStart = iHalf*nTrees;
    nEnd = (iHalf + 1)*nTrees;
  }
  else if(iHalf==nHalf-1) {
    nStart = iHalf*nTrees;
    nEnd = listTrees.size();
  }

  if(debug) cout << "gonna add the trees to the chain" << endl;

  TChain * myChain = new TChain ("produceNtuple/eIDSimpleTree");
  for(int i=nStart ; i<nEnd ; i++)
    myChain->Add(listTrees[i]);

  if(debug) cout << "added" << endl;

  cout << "tree list size : " << listTrees.size() << endl
	 << "nStart=" << nStart << " nEnd=" << nEnd << " nTrees=" << nTrees << endl; 
  //for(int i=nStart ; i<nEnd ; i++)
  //outlog << listTrees[i] << endl;

  if(nEntries==-999) {
    cout << "no tree processing requested ; exiting" << endl;
    return 0;
  }

  // Process the tree
  if(debug) cout << "process the tree" << endl;
  TreeToHistos(myChain, zeroThresh, dirOut, outlog, nEntries, halftag, spikeCut, debug);

  return 1;
}
