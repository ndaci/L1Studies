// ancestor : spikes_2011A_EmulReco_GetContam.C
// adding N(trig by spikes)_M

#include "../Common/baseFuncNad.h"

// tree lists
#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/Spikes-May10/treelist.h"
#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/PromptRecoV4/treelist.h"
#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/treelist_Spikes_2011A_1fb.h"
#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/treelist_Spikes_2011A_1fb_try2.h"
//#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/HLTEG5/treelist.h"

// Run 2011A
#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/treelist_Run2011A.h"
#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/PromptRecoV6/treelist.h"

// Run 2011AB HLTEG5
//#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/HLTEG5/treelist_Run2011A.h" // 2011A
//#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/HLTEG5/Run2011B_PRV1/treelist.h" // 2011B
//#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/HLTEG5/treelist_Run2011AB.h" // 2011A + 2011B

#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/HLTEG/treelist_2011A_HLTEG.h"
#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/HLTEG/Run2011B_PRV1/treelist.h"
#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/HLTEG/HighPU/treelist.h"

#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/CommissioningEmulReco/REPRO/treelist_EmulRecoNoKill_Run2011A.h"

#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/CommissioningEmulReco/NewKill/Run2011A/treelist.h"
#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/CommissioningEmulReco/NewKill/Run2011B/treelist.h"

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

int detIdx_summary(TString data_sum, int iVtx) {

  int tab_2011A[3] = {4,7,14};
  int tab_2011B[3] = {6,14,19};
  int tab_high_PU[3] = {6,14,19};

  for(int i=0 ; i<3 ; i++) {
    if(data_sum=="2011A") {
      if(iVtx<=tab_2011A[i]) return i;
    }
    else if(data_sum=="2011B") {
      if(iVtx<=tab_2011B[i]) return i;
    }
    else if(data_sum=="high_PU") {
      if(iVtx<=tab_high_PU[i]) return i;
    }
    else break;
  }

  return -1;

}

void TreeToHistos(TChain* myChain, TString data_sum, TString dirOut, ofstream& outlog, int nEntries, 
		  TString halftag, float spikeCut, bool debug, TString graph_style, int sev_user, TString hlt_sel)  {
  
  ostringstream ossi;

  //int nCat = 13;
  //int iCat;
  //TString category[13]={"_t10","_t15","_t17","_t19","_t21","_t23","_t30","_t35","_t40","_t50","_tAll","_t148524","_tRest"} ;
  //int thresh[10] = {10,15,17,19,21,23,30,35,40,50} ;
  //string EGval[6] = {"1","2","5","8","10","12"} ;

  // get the tree
  int nEvent, nRun, nLumi, trig_isUnbiased, trig_isL1SingleEG2, trig_isL1SingleEG5, trig_isL1SingleEG8, vtx_N;
  int trig_HLT_path[4];

  // Masked TT
  int trig_nMaskedCh;
  int trig_iMaskedTTeta[4032], trig_iMaskedTTphi[4032];
 
  // TP info
  int trig_tower_N,trig_tower_ieta[4032],trig_tower_iphi[4032],trig_tower_adc[4032],trig_tower_sFGVB[4032]; 
  int trig_tower_N_modif,trig_tower_ieta_modif[4032],trig_tower_iphi_modif[4032],trig_tower_adc_modif[4032],trig_tower_sFGVB_modif[4032]; 
  int trig_tower_N_emul,trig_tower_ieta_emul[4032],trig_tower_iphi_emul[4032],trig_tower_adc_emul[4032][5],trig_tower_sFGVB_emul[4032][5];

  // L1 candidates info
  int trig_L1emIso_N, trig_L1emNonIso_N, trig_L1emIso_N_M, trig_L1emNonIso_N_M;
  int trig_L1emIso_ieta[4], trig_L1emIso_iphi[4], trig_L1emIso_rank[4];
  int trig_L1emNonIso_ieta[4], trig_L1emNonIso_iphi[4], trig_L1emNonIso_rank[4];
  int trig_L1emIso_ieta_M[4], trig_L1emIso_iphi_M[4], trig_L1emIso_rank_M[4];
  int trig_L1emNonIso_ieta_M[4], trig_L1emNonIso_iphi_M[4], trig_L1emNonIso_rank_M[4];

  // Spikes
  int spike_N,spike_TTieta[5000], spike_TTiphi[5000], spike_Rieta[5000], spike_Riphi[5000], spike_severityLevel[5000], spike_outOfTime[5000];
  double  spike_Et[5000], spike_eta[5000], spike_phi[5000], spike_theta[5000];

  // Global
  myChain->SetBranchAddress("nEvent",&nEvent);
  myChain->SetBranchAddress("nRun",&nRun);
  myChain->SetBranchAddress("nLumi",&nLumi);

  // Vertices
  myChain->SetBranchAddress("vtx_N",&vtx_N);

  // Trigger
  myChain->SetBranchAddress("trig_HLT_path",&trig_HLT_path);
  //  myChain->SetBranchAddress("trig_isUnbiased",&trig_isUnbiased);
  //  myChain->SetBranchAddress("trig_isL1SingleEG2",&trig_isL1SingleEG2);
  //  myChain->SetBranchAddress("trig_isL1SingleEG5",&trig_isL1SingleEG5);
  //  myChain->SetBranchAddress("trig_isL1SingleEG8",&trig_isL1SingleEG8);

  // L1 candidates
  myChain->SetBranchAddress("trig_L1emIso_N", &trig_L1emIso_N);
  myChain->SetBranchAddress("trig_L1emIso_ieta", &trig_L1emIso_ieta);
  myChain->SetBranchAddress("trig_L1emIso_iphi", &trig_L1emIso_iphi);
  myChain->SetBranchAddress("trig_L1emIso_rank", &trig_L1emIso_rank);

  myChain->SetBranchAddress("trig_L1emNonIso_N", &trig_L1emNonIso_N);
  myChain->SetBranchAddress("trig_L1emNonIso_ieta", &trig_L1emNonIso_ieta);
  myChain->SetBranchAddress("trig_L1emNonIso_iphi", &trig_L1emNonIso_iphi); 
  myChain->SetBranchAddress("trig_L1emNonIso_rank", &trig_L1emNonIso_rank);

  myChain->SetBranchAddress("trig_L1emIso_N_M", &trig_L1emIso_N_M);
  myChain->SetBranchAddress("trig_L1emIso_ieta_M", &trig_L1emIso_ieta_M);
  myChain->SetBranchAddress("trig_L1emIso_iphi_M", &trig_L1emIso_iphi_M);
  myChain->SetBranchAddress("trig_L1emIso_rank_M", &trig_L1emIso_rank_M);

  myChain->SetBranchAddress("trig_L1emNonIso_N_M", &trig_L1emNonIso_N_M);
  myChain->SetBranchAddress("trig_L1emNonIso_ieta_M", &trig_L1emNonIso_ieta_M);
  myChain->SetBranchAddress("trig_L1emNonIso_iphi_M", &trig_L1emNonIso_iphi_M);
  myChain->SetBranchAddress("trig_L1emNonIso_rank_M", &trig_L1emNonIso_rank_M);

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


  // VARIOUS //
  bool isGoodRun,output,saturates,wellId,missed,eliminated,zeroedTT,trigEG8,trigEG8M;
  string flag;
  vector<int> firedEG_N; // EG 1/2/5/8/10/12 
  vector<int> firedEG_M; 

  // TRIGGERING //
  const int nEG = 8;
  const int nCand = 8;
  TString trigname[nEG] = {"2","5","8","10","12","15","20","30"} ;
  //int trigthresh[7] = {4,10,16,20,24,30,40} ;
  int trigthresh[nEG] = {2,5,8,10,12,15,20,30} ;

  // BOOLE //
  bool triggEG[nEG], triggEGM[nEG], anyTrigEG[nEG], anyTrigEGM[nEG];
  bool isL1Spiky[nCand][nEG], isL1SpikyM[nCand][nEG], isnotL1Spiky[nCand][nEG], isnotL1SpikyM[nCand][nEG];
  for(int iEG=0 ; iEG<nEG ; iEG++) {
    triggEG[iEG] = false ;
    triggEGM[iEG] = false ;
    anyTrigEG[iEG] = false ;
    anyTrigEGM[iEG] = false ;
    
    for(int iCand=0; iCand<nCand ; iCand++) {
      isL1Spiky[iCand][iEG] = false;
      isL1SpikyM[iCand][iEG] = false;
      isnotL1Spiky[iCand][iEG] = false;
      isnotL1SpikyM[iCand][iEG] = false;
    }
  }

  int nBiased=0;

  // HISTOGRAMS //
  TH1F * h_spikes;
  TH1F * h_rejected_spikes[ nEG ];
  TH1F * h_triggering_spikes[ nEG ];

  const int nStrict = 2; // 0: #evts(spike trig EGn) ; 1: #evts(!spike trig EGn)
  TString strict_name[nStrict] = {"large","strict"};

  const int nM = 2;
  TString name_M[nM] = {"_N","_M"};

  TH1F * h_evts_trigBy_any[nEG][nM];
  TH1F * h_evts_trigBy_spikes[nEG][nStrict][nM];
  //TH1F * h_contamination[nEG][nStrict]; // #evts[iStrict] / #evts(EGn triggered)

  const int nVtx = 40;
  TString name_histo;

  h_spikes = new TH1F("h_spikes","h_spikes",nVtx,0,nVtx); // abscisse = vtx_N

  for( int i=0 ; i<nEG ; i++ ) {
    h_rejected_spikes[i] = new TH1F("h_rejected_spikes_EG"+trigname[i],"h_rejected_spikes_EG"+trigname[i],nVtx,0,nVtx);
    h_triggering_spikes[i] = new TH1F("h_triggering_spikes_EG"+trigname[i],"h_triggering_spikes_EG"+trigname[i],nVtx,0,nVtx);
    
    for( int iM=0 ; iM<nM ; iM++) {
      name_histo = "h_evts_trigBy_any_EG"+trigname[i]+name_M[iM];
      h_evts_trigBy_any[i][iM] = new TH1F(name_histo,name_histo,nVtx,0,nVtx);
    
      for( int iStrict=0 ; iStrict<nStrict ; iStrict++ ) {
	name_histo = "h_evts_trigBy_spikes_EG"+trigname[i]+"_"+strict_name[iStrict]+name_M[iM];
	h_evts_trigBy_spikes[i][iStrict][iM] = new TH1F(name_histo,name_histo,nVtx,0,nVtx);
      }
    }
  }

  // SUMMARY
  int iVtx_sum = 3;
  TH1F * h_evts_trigBy_any_sum[nEG];
  TH1F * h_evts_trigBy_spikes_sum[nEG][nStrict];

  for( int i=0 ; i<nEG ; i++ ) {
    name_histo = "h_evts_trigBy_any_sum_EG"+trigname[i];
    h_evts_trigBy_any_sum[i] = new TH1F(name_histo,name_histo,3,0,3);
    
    for( int iStrict=0 ; iStrict<nStrict ; iStrict++ ) {
      name_histo = "h_evts_trigBy_spikes_sum_EG"+trigname[i]+"_"+strict_name[iStrict];
      h_evts_trigBy_spikes_sum[i][iStrict] = new TH1F(name_histo,name_histo,3,0,3);
    }
  }
  
  // detIdx_summary(TString data_sum, int iVtx)

  // COUNTERS //
  vector<int> matchIso, matchNonIso, matchIsoM, matchNonIsoM;
  int nSpikyL1[nEG], nLostSpikyL1[nEG], nSpikes, nSpikyL1sat, nSpikesWellId, nDiffSFGVB, nMatchedTT, nLostTT;

  nSpikes = nSpikyL1sat = nSpikesWellId = nDiffSFGVB = nMatchedTT = nLostTT = 0;
  for(int itg=0 ; itg<nEG ; itg++) {
    nSpikyL1[itg] = nLostSpikyL1[itg] = 0 ;
  }
  
  vector<int> IdxCat;
  vector< pair<int,int> > maskedTT;

  MAPTTS adcTT;
  MAPTTS::iterator iterTT;
  pair<int,int> coords;


  //ofstream ttfile("ttfile.txt",ios::out);

  // -------------------------------------------------------------------------------
  // JSON FILE READER
  // -------------------------------------------------------------------------------
  // define map of run/LS
  //string jsonFile = "/home/llr/cms/ochando/Analysis/Ana_Multi/json/goodrunlist_json.txt";
  //string jsonFile = "/data_CMS/cms/ndaci/ndaci_2011A/JSON/Cert_160404-163869_May10ReReco_163870-167913_7TeV_PromptReco_Collisions11_JSON.txt";
  //string jsonFile1 = "/data_CMS/cms/ndaci/ndaci_2011A/JSON/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v2.txt";
  //string jsonFile2 = "/data_CMS/cms/ndaci/ndaci_2011A/JSON/Cert_160404-167913_7TeV_PromptReco_Collisions11_JSON.txt"; 

  string jsonFile1,jsonFile2;
  jsonFile1 = jsonFile2 = "/data_CMS/cms/ndaci/ndaci_2011A/JSON/Cert_160404-178078_7TeV_PromptReco_Collisions11_JSON.txt" ;

  map<int, vector<pair<int, int> > > jsonMap1 = readJSONFile(jsonFile1);   
  //map<int, vector<pair<int, int> > > jsonMap2 = readJSONFile(jsonFile2);   
  
  if(debug) cout << "gonna loop over events" << endl;
  TString filename;

  int nCurrentRun;

  int numEntries = myChain->GetEntries () ;
  outlog << "numEntries=" << numEntries << endl;

  int nProcess = numEntries;
  if(nEntries>=0 && nEntries<numEntries)
    nProcess = nEntries;

  //int nVtx_max = 1;

  // loop over events
  for (int iEvent = 0 ; iEvent < nProcess ; iEvent++ )
    { 
      //if(iEvent%500 != 0) continue;
      myChain->GetEntry (iEvent) ;

      // run selection
      /*
      if( nRun>=160404 && nRun <= 163869 ) 
 	isGoodRun = AcceptEventByRunAndLumiSection(nRun, nLumi, jsonMap1);
      else if( nRun>=163870 && nRun<=167913 )
 	isGoodRun = AcceptEventByRunAndLumiSection(nRun, nLumi, jsonMap2) ;
      else isGoodRun = false;
      */
      //isGoodRun = AcceptEventByRunAndLumiSection(nRun, nLumi, jsonMap1);
      //if(!isGoodRun) continue;
      
      // HLT selection
      //if(trig_isL1SingleEG8 != 1)
      if(!trig_HLT_path[0]) {
	nBiased++ ;
	if(hlt_sel=="unbias")
	  continue;
      }
    

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


      // show which run is being processed
      if(iEvent==0) {
	nCurrentRun = nRun ;
	outlog << "nRun=" << nRun  << endl;
      }
      else if(nRun!=nCurrentRun) {
	nCurrentRun=nRun ;
	outlog << "nRun=" << nRun  << endl;
      }
      //cout << nRun << " " << nEvent << endl;

      // map the towers
      /*
      adcTT.clear();
      for(int t=0 ; t<trig_tower_N ; t++) {
	coords = make_pair( trig_tower_ieta[t] , trig_tower_iphi[t] );
	adcTT[coords].first.first = trig_tower_adc[t];
	adcTT[coords].second.first = trig_tower_sFGVB[t];
	if(trig_tower_ieta[t]>-17 && trig_tower_ieta[t]<18)
	  for(int i=0 ; i<IdxCat.size() ; i++)
	    h_adc[0][IdxCat[i]]->Fill(0.25*trig_tower_adc[t]);
	    //ttfile << "eta="    << trig_tower_ieta[t]
	    //<< "   phi=" << trig_tower_iphi[t]
	  //<< "   adc=" <<
	  //<< "   sFGVB=" << trig_tower_sFGVB[t]
	  //<< endl;
      }
      for(int t=0 ; t<trig_tower_N_emul ; t++) {
	coords = make_pair( trig_tower_ieta_emul[t] , trig_tower_iphi_emul[t] );
	iterTT = adcTT.find( coords );
	if( iterTT != adcTT.end() ) {
	  adcTT[coords].first.second = trig_tower_adc_emul[t][2];
	  adcTT[coords].second.second = trig_tower_sFGVB_emul[t][2];
	  if(trig_tower_ieta_emul[t]>-17 && trig_tower_ieta_emul[t]<18) 
	    for(int i=0 ; i<IdxCat.size() ; i++)
	      h_adc[3]->Fill(0.25*trig_tower_adc_emul[t][2]);
	      //ttfile << "etaM="    << trig_tower_ieta_modif[t]
	      // << "   phiM=" << trig_tower_iphi_modif[t]
	      //<< "   sFGVBM=" << trig_tower_sFGVB_modif[t]
	      //<< endl;
	}
	else {
	  outlog << "mapping problem" << endl;
	  adcTT[coords].first.first = -777;
	  adcTT[coords].first.second = trig_tower_adc_emul[t][2];
	  adcTT[coords].second.first = -777;
	  adcTT[coords].second.second = trig_tower_sFGVB_emul[t][2];
	  
	}
      }
      */

//       for(int t=0 ; t<trig_tower_N_modif ; t++) {
// 	coords = make_pair( trig_tower_ieta_modif[t] , trig_tower_iphi_modif[t] );
// 	iterTT = adcTT.find( coords );
// 	if( iterTT != adcTT.end() ) {
// 	  adcTT[coords].first.second = trig_tower_adc_modif[t];
// 	  adcTT[coords].second.second = trig_tower_sFGVB_modif[t];
// 	  if(trig_tower_ieta_modif[t]>-17 && trig_tower_ieta_modif[t]<18) 
// 	    for(int i=0 ; i<IdxCat.size() ; i++)
// 	      h_adc[3]->Fill(0.25*trig_tower_adc_modif[t]);
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


      // initialize event triggering booleans
      for(int iEG=0 ; iEG<nEG ; iEG++) {
	anyTrigEG[iEG] = anyTrigEGM[iEG] = false; // any object in the event triggers iEG
	for(int iCand=0 ; iCand<nCand ; iCand++)
	  isL1Spiky[iCand][iEG] = isL1SpikyM[iCand][iEG] = // candidate iCand is spiky and triggers iEG
	    isnotL1Spiky[iCand][iEG] = isnotL1SpikyM[iCand][iEG] = false; // candidate iCand is not spiky, and triggers iEG
      }

      // loop over evil rechits     
      if(debug) cout << "loop over evil rechits" << endl;
      saturates = false;

      for(int irh=0 ; irh<spike_N ; irh++) {
	if(spike_severityLevel[irh]>= sev_user) { // (kTime) + kWeird

	  // selection on the energy
	  if( spike_Et[irh] < spikeCut ) continue;
	  	  
	  nSpikes ++ ;
	  h_spikes->Fill(vtx_N);
	  
	  //if( vtx_N > nVtx_max ) nVtx_max = vtx_N ;
	  
	  matchIso.clear();
	  matchNonIso.clear();
	  matchIsoM.clear();
	  matchNonIsoM.clear();
	  
	  if(spike_Et[irh]>0) {

	    output=false;
	    saturates = false;
	    for(int itg=0 ; itg<nEG ; itg++)
	      triggEG[itg] = triggEGM[itg] = false;
	    
	    for(int ic=0 ; ic<4 ; ic++) {
	      if(trig_L1emIso_ieta_M[ic]==spike_Rieta[irh] && trig_L1emIso_iphi_M[ic]==spike_Riphi[irh]) {
		output = true;
		matchIso.push_back(ic);
		if(trig_L1emIso_rank_M[ic]>=63) saturates = true;
		for(int itg=0 ; itg<nEG ; itg++) {
		  if(trig_L1emIso_rank_M[ic]>=trigthresh[itg]) {
		    triggEG[itg] = true; // spike irh triggers EG[itg]
		    isL1Spiky[ic][itg] = true; // candidate ic is spiky and triggers EG[itg]
		  }
		}
	      }
	      if(trig_L1emNonIso_ieta_M[ic]==spike_Rieta[irh] && trig_L1emNonIso_iphi_M[ic]==spike_Riphi[irh]) {
		output = true;
		matchNonIso.push_back(ic);
		if(trig_L1emNonIso_rank_M[ic]>=63) saturates = true;
		for(int itg=0 ; itg<nEG ; itg++) {
		  if(trig_L1emNonIso_rank_M[ic]>=trigthresh[itg]) {
		    triggEG[itg] = true; // spike irh triggers EG[itg]
		    isL1Spiky[ic+4][itg] = true; // candidate ic+4 is spiky and triggers EG[itg]
		  }
		}
	      }
	      if(trig_L1emIso_ieta[ic]==spike_Rieta[irh] && trig_L1emIso_iphi[ic]==spike_Riphi[irh]) {
		output = true;
		matchIsoM.push_back(ic);
		if(trig_L1emIso_rank[ic]>=63) saturates = true;
		for(int itg=0 ; itg<nEG ; itg++) {
		  if(trig_L1emIso_rank[ic]>=trigthresh[itg]) {
		    triggEGM[itg] = true; // spike irh triggers EG[itg] in modified collection		
		    isL1SpikyM[ic][itg] = true; // candidate ic is spiky and triggers EG[itg] in modif coll
		  }
		}
	      }
	      if(trig_L1emNonIso_ieta[ic]==spike_Rieta[irh] && trig_L1emNonIso_iphi[ic]==spike_Riphi[irh]) {
		output = true;
		matchNonIsoM.push_back(ic);
		if(trig_L1emNonIso_rank[ic]>=63) saturates = true;
		for(int itg=0 ; itg<nEG ; itg++) {
		  if(trig_L1emNonIso_rank[ic]>=trigthresh[itg]) {
		    triggEGM[itg] = true; // spike irh triggers EG[itg] in modified collection
		    isL1SpikyM[ic+4][itg] = true; // candidate ic+4 is spiky and triggers EG[itg] in modif coll
		  }
		}
	      }
	    }

	    coords = make_pair(spike_TTieta[irh] , spike_TTiphi[irh]);
	    
	    bool b_mask = false;
	    for(int iTT=0 ; iTT<maskedTT.size() ; iTT++) {
	      if( maskedTT[iTT] == coords )
		b_mask = true;
	    }

	    if(output || true) {
	      outlog << "Spike   nRun=" << nRun << "   nEvent=" << nEvent << "("<<iEvent<<")" << endl
		     << "Rieta="         << spike_Rieta[irh]
		     << "   Riphi="      << spike_Riphi[irh]
		     << "   TTieta="     << spike_TTieta[irh]
		     << "   TTiphi="     << spike_TTiphi[irh]
		//<< "   outOfTime="  << spike_outOfTime[irh]
		     << "   severity="   << spike_severityLevel[irh]
		     << "   Et="         << spike_Et[irh];
	      
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
		  outlog << "IsoM rank=" << trig_L1emIso_rank_M[matchIsoM[ic]] << "      ";
	      if(matchNonIsoM.size()>0)
		for(int ic=0 ; ic<matchNonIsoM.size() ; ic++)
		  outlog << "NonIsoM rank=" << trig_L1emNonIso_rank_M[matchNonIsoM[ic]] ;
	      
	      outlog<<endl;
	      
// 	      outlog << "Trigger Tower : adc=" << adcTT[coords].first.first
// 		     << "   adcM="             << adcTT[coords].first.second
// 		     << "   sFGVB="            << adcTT[coords].second.first
// 		     << "   sFGVB_M="          << adcTT[coords].second.second ;
// 	      if(adcTT[coords].second.first != adcTT[coords].second.second) {
// 		for(int i=0 ; i<IdxCat.size() ; i++)
// 		  nDiffSFGVB++ ;
// 		outlog << "   !! sFGVB diff !!" ;
// 	      }
// 	      outlog << endl;
	      
	      wellId = zeroedTT = eliminated = missed = false;

		
	      // non L1-matching-dependant properties
	      if(adcTT[coords].second.second==0) {
		wellId = true;
		nSpikesWellId++ ;
	      }
	      if(adcTT[coords].first.first>0) {
		nMatchedTT++;
		if(adcTT[coords].first.second==0) {
		  nLostTT++;
		  zeroedTT = true;
		}
	      }		  
	      // non L1-matching
		
	      // L1-matching-dependant properties
	      if( matchIso.size()>0 || matchNonIso.size()>0 ) {
		  		  
		if(saturates) {
		  //saturates = true;
		  nSpikyL1sat++ ;
		}
		  
	      } // L1-matching
		
	      // EGX
	      for(int itg=0 ; itg<nEG ; itg++) {
		if( triggEG[itg] ) {
		  nSpikyL1[itg]++ ;
		  h_triggering_spikes[itg]->Fill(vtx_N);

		  //if(itg==2) {
		  //}
		  if( !triggEGM[itg] ) {
		    nLostSpikyL1[itg]++ ;
		    eliminated = true;
		  }
		} // EGX-L1-matching
		else {
		  h_rejected_spikes[itg]->Fill(vtx_N);
		}
	      }
	      // Modified L1-matching
// 	      if( matchIsoM.size()>0 || matchNonIsoM.size()>0 ) {
// 	      }
		
	      // EG8
	      if(triggEGM[2]) {
		missed = true;
	      }

	      if(wellId) outlog << " wellId";
	      if(zeroedTT) outlog << " zeroedTT" ;
	      if(eliminated) outlog << " eliminatedEG8";
	      if(missed) outlog << " missedEG8";
	      if(saturates) outlog << " saturates";
	      outlog << endl ;
	      
	      // show all the trigger towers of the region
	      /*
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
	      */
	      //outlog << endl;

	    } // endif output
	  } //endif
	}//endif
      } // loop over rechits 
    
      // loop over candidates to determinate the contamination
      for(int ic=0 ; ic<4 ; ic++) {
	for(int itg=0 ; itg<nEG ; itg++) {
	  if(trig_L1emIso_rank[ic]>=trigthresh[itg]) {
	    anyTrigEG[itg] = true;
	    if(!isL1Spiky[ic][itg]) isnotL1Spiky[ic][itg] = true;
	  }
	  if(trig_L1emNonIso_rank[ic]>=trigthresh[itg]) {
	    anyTrigEG[itg] = true;
	    if(!isL1Spiky[ic+4][itg]) isnotL1Spiky[ic+4][itg] = true;
	  }
	  if(trig_L1emIso_rank_M[ic]>=trigthresh[itg]) {
	    anyTrigEGM[itg] = true;
	    if(!isL1SpikyM[ic][itg]) isnotL1SpikyM[ic][itg] = true;
	  }
	  if(trig_L1emNonIso_rank_M[ic]>=trigthresh[itg]) {
	    anyTrigEGM[itg] = true;
	    if(!isL1SpikyM[ic+4][itg]) isnotL1SpikyM[ic+4][itg] = true;
	  }
	}
      }
      
      bool spikyL1 = false;
      bool notOnlySpikyL1 = false;

      bool spikyL1M = false;
      bool notOnlySpikyL1M = false;

      for(int iEG=0 ; iEG<nEG ; iEG++) {
	
	if( anyTrigEG[iEG] ) {
	  h_evts_trigBy_any[iEG][0]->Fill(vtx_N);

	  iVtx_sum = detIdx_summary(data_sum,vtx_N);
	  if(iVtx_sum > -1) h_evts_trigBy_any_sum[iEG]->Fill(vtx_N);
	}
	if( anyTrigEGM[iEG] ) {
	  h_evts_trigBy_any[iEG][1]->Fill(vtx_N);
	}
	// determine nature of triggering objects
	spikyL1 = notOnlySpikyL1 = spikyL1M = notOnlySpikyL1M = false;

	for(int ic=0 ; ic<8 ; ic++) {
	  if( isL1Spiky[ic][iEG] ) spikyL1 = true; // there is a spiky L1 cand triggering iEG
	  if( isnotL1Spiky[ic][iEG] ) notOnlySpikyL1 = true; // there is a non-spiky L1 cand triggering iEG 
	  if( isL1SpikyM[ic][iEG] ) spikyL1M = true; // there is a spiky L1 cand triggering iEG
	  if( isnotL1SpikyM[ic][iEG] ) notOnlySpikyL1M = true; // there is a non-spiky L1 cand triggering iEG 
	}

	if( spikyL1 ) {
	  h_evts_trigBy_spikes[iEG][0][0]->Fill(vtx_N);

	  iVtx_sum = detIdx_summary(data_sum,vtx_N);
	  if(iVtx_sum > -1) h_evts_trigBy_spikes_sum[iEG][0]->Fill(vtx_N);

	  // there was no non-spiky L1 cand triggering iEG : 
	  // the event was "strictly" triggered by a spike
	  if( !notOnlySpikyL1 ) {
	    h_evts_trigBy_spikes[iEG][1][0]->Fill(vtx_N);
	    
	    iVtx_sum = detIdx_summary(data_sum,vtx_N);
	    if(iVtx_sum > -1) h_evts_trigBy_spikes_sum[iEG][1]->Fill(vtx_N);

	  }
	}
	if( spikyL1M ) {
	  h_evts_trigBy_spikes[iEG][0][1]->Fill(vtx_N);

	  iVtx_sum = detIdx_summary(data_sum,vtx_N);
	  if(iVtx_sum > -1) h_evts_trigBy_spikes_sum[iEG][0]->Fill(vtx_N);

	  // there was no non-spiky L1 cand triggering iEG : 
	  // the event was "strictly" triggered by a spike
	  if( !notOnlySpikyL1M ) {
	    h_evts_trigBy_spikes[iEG][1][1]->Fill(vtx_N);
	    
	    iVtx_sum = detIdx_summary(data_sum,vtx_N);
	    if(iVtx_sum > -1) h_evts_trigBy_spikes_sum[iEG][1]->Fill(vtx_N);

	  }
	}

      }
      


    }//loop over events  


  // show counters
  outlog << "COUNTING RESULTS" << endl;
  outlog << " : nSpikes="      << nSpikes 
	 << " | nSpikesWellId="<< nSpikesWellId
	 << " | nSpikyL1EG8="  << nSpikyL1[2]
	 << " | nSpikyL1sat="  << nSpikyL1sat
	 << endl
	 << " | nLostSpikyL1EG8=" << nLostSpikyL1[2] 
	 << " | nMatchedTT="      << nMatchedTT
	 << " | nLostTT="         << nLostTT
	 << " | nDiffSFGVB="      << nDiffSFGVB
	 << endl;
  

  ofstream outcount(dirOut+"counters"+halftag+".csv",ios::out);
  for(int itg=0 ; itg<nEG ; itg++) {
    outcount << " | " << nSpikes 
	     << " | " << nSpikesWellId
	     << " | " << nSpikyL1[itg]
	     << " | " << nLostSpikyL1[itg]
	     << " | " << nMatchedTT
	     << " | " << nLostTT
	     << " | " << nDiffSFGVB
	     << " | " << nSpikyL1sat 
	     << endl;
  
    outcount << endl << endl;
  }

  // PLOTS //

  // Graph //
  
  float n_spikes[nVtx+1]; // last bin = overall N_vtx

  float n_rejected_spikes[nEG][nVtx+1];
  float n_triggering_spikes[nEG][nVtx+1];

  float n_rejection[nEG][nVtx+1];
  float n_unreject[nEG][nVtx+1];

  float n_evts_trigBy_any[nEG][nVtx+1][nM];
  float n_evts_trigBy_spikes[nEG][nStrict][nVtx+1][nM];
  float n_event_contamination[nEG][nStrict][nVtx+1][nM];

  for(int iBin=0 ; iBin<nVtx+1 ; iBin++) {
    n_spikes[iBin] = 0;
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      n_rejected_spikes[iEG][iBin] = 0;
      n_triggering_spikes[iEG][iBin] = 0;

      for(int iM=0 ; iM<nM ; iM++) {
	n_evts_trigBy_any[iEG][iBin][iM] = 0;
	for(int iStrict=0 ; iStrict<nStrict ; iStrict++) {
	  n_evts_trigBy_spikes[iEG][iStrict][iBin][iM] = 0;
	  n_event_contamination[iEG][iStrict][iBin][iM] = 0;
	}
      }
    }
  }
  
  for(int iBin=0 ; iBin<nVtx ; iBin++) {
    n_spikes[iBin] = (float)(h_spikes->GetBinContent(iBin+1));
    n_spikes[nVtx] += n_spikes[iBin];

    for(int iEG=0 ; iEG<nEG ; iEG++) {

      // Event contamination
      for(int iM=0; iM<nM ; iM++) {
	n_evts_trigBy_any[iEG][iBin][iM] = (float)(h_evts_trigBy_any[iEG][iM]->GetBinContent(iBin+1));
        n_evts_trigBy_any[iEG][nVtx][iM] += n_evts_trigBy_any[iEG][iBin][iM];

	for(int iStrict=0 ; iStrict<nStrict ; iStrict++) {

	  n_evts_trigBy_spikes[iEG][iStrict][iBin][iM] = (float)(h_evts_trigBy_spikes[iEG][iStrict][iM]->GetBinContent(iBin+1));
	  n_evts_trigBy_spikes[iEG][iStrict][nVtx][iM] += n_evts_trigBy_spikes[iEG][iStrict][iBin][iM] ;

	  if( n_evts_trigBy_any[iEG][iBin][iM] != 0 )
	    n_event_contamination[iEG][iStrict][iBin][iM] = 
	      n_evts_trigBy_spikes[iEG][iStrict][iBin][iM] / n_evts_trigBy_any[iEG][iBin][iM] ;
	  else
	    n_event_contamination[iEG][iStrict][iBin][iM] = -0.01;
	}
      }
      
      // Spike rejection/contamination
      n_rejected_spikes[iEG][iBin] = (float)(h_rejected_spikes[iEG]->GetBinContent(iBin+1));
      n_triggering_spikes[iEG][iBin] = (float)(h_triggering_spikes[iEG]->GetBinContent(iBin+1));

      n_rejected_spikes[iEG][nVtx] += n_rejected_spikes[iEG][iBin];
      n_triggering_spikes[iEG][nVtx] += n_triggering_spikes[iEG][iBin];      

      if(n_spikes[iBin]!=0) {
	n_rejection[iEG][iBin] = n_rejected_spikes[iEG][iBin] / n_spikes[iBin] ;
	n_unreject[iEG][iBin] = n_triggering_spikes[iEG][iBin] / n_spikes[iBin] ;
      } else {
	n_rejection[iEG][iBin] = -0.01 ;
	n_unreject[iEG][iBin] = -0.01 ;
      }
    }
  }
  // last bin for TGraph
  for(int iEG=0 ; iEG<nEG ; iEG++) {

    // Event contamination
    for(int iStrict=0 ; iStrict<nStrict ; iStrict++) {
	for(int iM=0; iM<nM ; iM++) {
	  if(n_evts_trigBy_any[iEG][nVtx][iM] != 0 )
	    n_event_contamination[iEG][iStrict][nVtx][iM] = n_evts_trigBy_spikes[iEG][iStrict][nVtx][iM] / n_evts_trigBy_any[iEG][nVtx][iM] ;
	  else
	    n_event_contamination[iEG][iStrict][nVtx][iM] = -0.01 ;
	}
    }

    // Spike rejection/contamination
    if(n_spikes[nVtx]!=0) {
      n_rejection[iEG][nVtx] = n_rejected_spikes[iEG][nVtx] / n_spikes[nVtx] ;
      n_unreject[iEG][nVtx] = n_triggering_spikes[iEG][nVtx] / n_spikes[nVtx] ;
    } else {
      n_rejection[iEG][nVtx] = -0.01 ;
      n_unreject[iEG][nVtx] = -0.01 ;
    }
  }

  // put it in the log

  // Spike rejection/contamination
  outlog << "New results" << endl
	 << "N_Spikes[ iVtx ] = " ;
  for(int iBin=0 ; iBin<nVtx+1 ; iBin++)
    outlog << n_spikes[iBin] << " | " ;

  outlog << endl << endl << "N_Rejected_Spikes : " << endl;
  for(int iBin=0 ; iBin<nVtx+1 ; iBin++) {
    outlog << " iVtx=" << iBin << " : ";
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outlog << n_rejected_spikes[iEG][iBin] << " | " ;
    }
    outlog << endl;
  }
  
  outlog << endl << endl << "N_Unreject_Spikes : " << endl;
  for(int iBin=0 ; iBin<nVtx+1 ; iBin++) {
    outlog << " iVtx=" << iBin << " : ";
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outlog << n_triggering_spikes[iEG][iBin] << " | " ;
    }
    outlog << endl;
  }

  outlog << endl << endl << "rejected fraction : " << endl;
  for(int iBin=0 ; iBin<nVtx+1 ; iBin++) {
    outlog << " iVtx=" << iBin << " : ";
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outlog << n_rejection[iEG][iBin] << " | " ;
    }
    outlog << endl;
  }

  outlog << endl << endl << "unreject fraction : " << endl;
  for(int iBin=0 ; iBin<nVtx+1 ; iBin++) {
    outlog << " iVtx=" << iBin << " : ";
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outlog << n_unreject[iEG][iBin] << " | " ;
    }
    outlog << endl;
  }

  // Event contamination
  outlog << endl << endl << "n_evts_trigBy_any N : " << endl;
  for(int iBin=0 ; iBin<nVtx+1 ; iBin++) {
    outlog << " iVtx=" << iBin << " : ";
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outlog << n_evts_trigBy_any[iEG][iBin][0] << " | " ;
    }
    outlog << endl;
  }

  outlog << endl << endl << "n_evts_trigBy_any M : " << endl;
  for(int iBin=0 ; iBin<nVtx+1 ; iBin++) {
    outlog << " iVtx=" << iBin << " : ";
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outlog << n_evts_trigBy_any[iEG][iBin][1] << " | " ;
    }
    outlog << endl;
  }


  outlog << endl << endl << "n_evts_trigBy_spikes LARGE : " << endl;
  for(int iBin=0 ; iBin<nVtx+1 ; iBin++) {
    outlog << " iVtx=" << iBin << " : ";
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outlog << n_evts_trigBy_spikes[iEG][0][iBin][0] << " | " ;
    }
    outlog << endl;
  }

  outlog << endl << endl << "n_evts_trigBy_spikes LARGE M : " << endl;
  for(int iBin=0 ; iBin<nVtx+1 ; iBin++) {
    outlog << " iVtx=" << iBin << " : ";
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outlog << n_evts_trigBy_spikes[iEG][0][iBin][1] << " | " ;
    }
    outlog << endl;
  }

  outlog << endl << endl << "fraction evts_trigBy_spikes LARGE : " << endl;
  for(int iBin=0 ; iBin<nVtx+1 ; iBin++) {
    outlog << " iVtx=" << iBin << " : ";
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outlog << n_event_contamination[iEG][0][iBin][0] << " | " ;
    }
    outlog << endl;
  }  

  outlog << endl << endl << "fraction evts_trigBy_spikes LARGE M: " << endl;
  for(int iBin=0 ; iBin<nVtx+1 ; iBin++) {
    outlog << " iVtx=" << iBin << " : ";
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outlog << n_event_contamination[iEG][0][iBin][1] << " | " ;
    }
    outlog << endl;
  }  

  outlog << endl << endl << "n_evts_trigBy_spikes STRICT : " << endl;
  for(int iBin=0 ; iBin<nVtx+1 ; iBin++) {
    outlog << " iVtx=" << iBin << " : ";
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outlog << n_evts_trigBy_spikes[iEG][1][iBin][0] << " | " ;
    }
    outlog << endl;
  }

  outlog << endl << endl << "fraction evts_trigBy_spikes STRICT : " << endl;
  for(int iBin=0 ; iBin<nVtx+1 ; iBin++) {
    outlog << " iVtx=" << iBin << " : ";
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outlog << n_event_contamination[iEG][1][iBin][0] << " | " ;
    }
    outlog << endl;
  }  

  outlog << endl << endl << "n_evts_trigBy_spikes STRICT M : " << endl;
  for(int iBin=0 ; iBin<nVtx+1 ; iBin++) {
    outlog << " iVtx=" << iBin << " : ";
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outlog << n_evts_trigBy_spikes[iEG][1][iBin][1] << " | " ;
    }
    outlog << endl;
  }

  outlog << endl << endl << "fraction evts_trigBy_spikes STRICT M : " << endl;
  for(int iBin=0 ; iBin<nVtx+1 ; iBin++) {
    outlog << " iVtx=" << iBin << " : ";
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outlog << n_event_contamination[iEG][1][iBin][1] << " | " ;
    }
    outlog << endl;
  }  

  // abscisse
  float x_vtx[nVtx+1];
  for(int iVtx=0 ; iVtx<nVtx+1 ; iVtx++)
    x_vtx[iVtx] = iVtx;
  
  // graphs
  TGraph * g_rejection[nEG];
  for(int iEG=0 ; iEG<nEG ; iEG++) {
    g_rejection[iEG] = new TGraph( nVtx+1, x_vtx, n_rejection[iEG] );
    //g_rejection[iEG] = new TGraph( h_rejected_spikes[iEG] );
    g_rejection[iEG]->SetName("Spike L1 rejection (EG"+trigname[iEG]+")");
    g_rejection[iEG]->SetTitle("Spike L1 rejection (EG"+trigname[iEG]+")");
  }

  TGraph * g_unreject[nEG];
  for(int iEG=0 ; iEG<nEG ; iEG++) {
    g_unreject[iEG] = new TGraph( nVtx+1, x_vtx, n_unreject[iEG] );
    g_unreject[iEG]->SetName("Spike L1 unrejection (EG"+trigname[iEG]+")");
    g_unreject[iEG]->SetTitle("Spike L1 unrejection (EG"+trigname[iEG]+")");
  }

  TGraph * g_evt_contam[nEG][nStrict][nM];
  for(int iStrict=0 ; iStrict<nStrict ; iStrict++) {
    for(int iM=0; iM<nM ; iM++) {
      for(int iEG=0 ; iEG<nEG ; iEG++) {
	g_evt_contam[iEG][iStrict][iM] = new TGraph( nVtx+1, x_vtx, n_event_contamination[iEG][iStrict][iM] );
	g_evt_contam[iEG][iStrict][iM]->SetName("Event "+strict_name[iStrict]+" contamination (EG"+trigname[iEG]+")");
	g_evt_contam[iEG][iStrict][iM]->SetTitle("Event "+strict_name[iStrict]+" contamination (EG"+trigname[iEG]+")");
      }
    }
  }

  /*
  // canvas
  TCanvas * c_spikes;
  TCanvas * c_rejected_spikes[ nEG ];
  TCanvas * c_triggering_spikes[ nEG ];
  TCanvas * c_rejection[ nEG ];
  TCanvas * c_unreject[ nEG ];  

  TCanvas * c_evts_trigBy_any[nEG];
  TCanvas * c_evts_trigBy_spikes[nEG][nStrict];
  TCanvas * c_evt_contam[nEG][nStrict];

  c_spikes = new TCanvas("h_N_Spikes","N_Spikes",0,0,800,600);
  h_spikes->Draw();

  for(int iEG=0 ; iEG<nEG ; iEG++) {
    c_rejected_spikes[iEG] = new TCanvas("h_N_Rejected_Spikes_EG"+trigname[iEG],"N_Rejected_Spikes_EG"+trigname[iEG],0,0,800,600);
    h_rejected_spikes[iEG]->Draw();
    c_rejected_spikes[iEG]->Print(dirOut+"rejected_spikes_EG"+trigname[iEG]+".gif");

    c_triggering_spikes[iEG] = new TCanvas("h_N_Unreject_Spikes_EG"+trigname[iEG],"N_Unreject_Spikes_EG"+trigname[iEG],0,0,800,600);
    h_triggering_spikes[iEG]->Draw();
    c_triggering_spikes[iEG]->Print(dirOut+"triggering_spikes_EG"+trigname[iEG]+".gif");

    c_rejection[iEG] = new TCanvas("c_spike_L1_rejection_EG"+trigname[iEG],"Spike L1 rejection (EG"+trigname[iEG]+")",0,0,800,600);
    g_rejection[iEG]->Draw(graph_style);
    c_rejection[iEG]->Print(dirOut+"rejection_EG"+trigname[iEG]+".gif");
 
    c_unreject[iEG] = new TCanvas("c_spike_L1_unreject_EG"+trigname[iEG],"Spike L1 unrejection (EG"+trigname[iEG]+")",0,0,800,600);
    g_unreject[iEG]->Draw(graph_style);
    c_unreject[iEG]->Print(dirOut+"unreject_EG"+trigname[iEG]+".gif");

    c_evts_trigBy_any[iEG] = new TCanvas("c_trigBy_any_EG"+trigname[iEG],"Events triggered by any object (EG"+trigname[iEG]+")",0,0,800,600);
    h_evts_trigBy_any[iEG]->Draw();
    c_evts_trigBy_any[iEG]->Print(dirOut+"events_trigByAny_EG"+trigname[iEG]+".gif");

    for(int iStrict=0 ; iStrict<nStrict ; iStrict++) {

      c_evts_trigBy_spikes[iEG][iStrict] = new TCanvas("c_evt_trigBySpikes_EG"+trigname[iEG]+"_"+strict_name[iStrict],
					       "Events triggered by spikes ("+strict_name[iStrict]+") : EG"+trigname[iEG],
					       0,0,800,600);
      h_evts_trigBy_spikes[iEG][iStrict]->Draw();
      c_evts_trigBy_spikes[iEG][iStrict]->Print(dirOut+"events_trigBySpikes_"+strict_name[iStrict]+"_EG"+trigname[iEG]+".gif");


      c_evt_contam[iEG][iStrict] = new TCanvas("c_evt_contam_EG"+trigname[iEG]+"_"+strict_name[iStrict],
					       "Event contamination "+strict_name[iStrict]+" (EG"+trigname[iEG]+")",
					       0,0,800,600);
      g_evt_contam[iEG][iStrict]->Draw(graph_style);
      c_evt_contam[iEG][iStrict]->Print(dirOut+"evt_contam_EG"+trigname[iEG]+"_"+strict_name[iStrict]+".gif");    
    }
  }
  */

  if(debug) cout << "is gonna plot in file" << endl;
  TFile *outplot = new TFile(dirOut+"spike_plots"+halftag+".root","RECREATE");
  h_spikes->Write();

  for(int iEG=0 ; iEG<nEG ; iEG++) {
    h_rejected_spikes[iEG]->Write();
    h_triggering_spikes[iEG]->Write();
    g_rejection[iEG]->Write();
    g_unreject[iEG]->Write();

    for(int iM=0; iM<nM ; iM++) {
      h_evts_trigBy_any[iEG][iM]->Write();
    }
    h_evts_trigBy_any_sum[iEG]->Write();

    for(int iStrict=0 ; iStrict<nStrict ; iStrict++) {
      h_evts_trigBy_spikes_sum[iEG][iStrict]->Write();
      for(int iM=0; iM<nM ; iM++) {
	h_evts_trigBy_spikes[iEG][iStrict][iM]->Write();
	g_evt_contam[iEG][iStrict][iM]->Write();
      }
    }
  }
  outplot->Close();
  if(debug) cout << "plotted in file" << endl;
  
}
  

int spikes(int nEntries=-1, int iHalf=0, int nHalf=1, int sev_user=3, TString hlt_sel="",
	   TString dirOut="/home/llr/cms/ndaci/SKWork/macro/skEfficiency/tagAndProbe/Spike2011A/CommissioningRun2011A/",
	   TString graph_style="A*", int marker_style=1, TString data="Spikes_Run2011A", TString data_sum="2011A", float spikeCut=8., bool debug=false)
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
  
  if( data=="Run2011A_EmulReco" ) {
    listTrees = list_Spikes_EmulRecoNoKill_Run2011A();
  }
  else if( data=="Run2011A_EmulReco_Newkill") {
    listTrees = list_Spikes_Run2011A_Newkill();
  }
  else if( data=="Run2011B_EmulReco_Newkill") {
    listTrees = list_Spikes_Run2011B_Newkill();
  }
  /*
  if(data=="Spikes_Run2011A_HLTEG") {
    listTrees = list_Spikes_Run2011A_HLTEG();
  }
  else if(data=="Spikes_Run2011B_HLTEG") {
    listTrees = list_Spikes_Run2011B_HLTEG();
  }
  else if(data=="Spikes_Run2011B_highPU_HLTEG") {
    listTrees = list_Spikes_highPU_HLTEG();
  }
  else if(data=="Spikes_2011A_HLTEG5") {
    listTrees = list_Spikes_2011A_HLTEG5();
  }
  else if(data=="Spikes_2011B_HLTEG5") {
    listTrees = list_Spikes_2011B_PRV1_HLTEG5();
  }
  else if(data=="Spikes_2011AB_HLTEG5") {
    listTrees = list_Spikes_2011AB();
  }
  else if(data=="Spikes_Run2011A") {
    listTrees = list_Spikes_Run2011A();
  }
  else if(data=="Spikes2011A_PRV6") {
    listTrees = list_Spikes_2011A_PromptRecoV6();
  }
  else if(data=="Spikes2011A") {
    listTrees = list_Spikes_2011A();
  }
  else if(data=="Spikes2011A_PromptRecoV4") {
    listTrees = list_Spikes_2011A_PromptRecoV4();
  }
  else if(data=="Spikes2011A_1fb") {
    listTrees = list_Spikes_2011A_1fb();
  }
  else if(data=="Spikes2011A_1fb_try2") {
    listTrees = list_Spikes_2011A_1fb_try2();
  }
  */
  /*
  else if(data=="Spikes2011A_HLTEG5") {
    listTrees = list_Spikes_2011A_HLTEG5();
  }
  */
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
  TreeToHistos(myChain, data_sum, dirOut, outlog, nEntries, halftag, spikeCut, debug, graph_style,sev_user,hlt_sel);

  return 1;
}
