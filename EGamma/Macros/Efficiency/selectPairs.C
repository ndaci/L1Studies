//////////////////////////////////////////////////////////////////////////////
///////////////////////////    selectPairs.C    //////////////////////////////
////                                                                      ////
////    Input : makePairs.C output trees containing pre-selected pairs    ////
////                                                                      ////
////    - further selection of the pairs                                  ////
////      --> tag : VBTF WP80 + eT > 20 GeV                               ////
////      + trig EG12/15 (HLT selection's L1 seed early 2011 SingleEle)   ////
////      --> probe : VBTF WP95 + eT > 5 GeV                              ////
////      --> mass : between 60 GeV and 120 GeV                           ////
////    - if pair selected : check if probe trig EGn                      ////
////                                                                      ////
////    Output : - trees readable by fastEfficiencyNadir.C to             ////
////               plot the efficiency turn-on curve                      ////
////             - counters of electrons / matched / matchedM etc         ////
////                                                                      ////
////   Designed to printout all the information for inefficient electrons ////
////                                                                      ////
////   Adding EG18 and EG22 (HLT menu 3e33)                               ////
////   + map efficiency eta/phi                                           ////
///////////////////////////    !!! ENJOY IT ;-) !!!    ///////////////////////
//////////////////////////////////////////////////////////////////////////////

// v2.6 : correcting Tag&Probe method --> if both ele of the pair are able to be tag => consider both possibilities
// v2.7 : better treenew and counters definition (removing pT + arrays of "l1_" branches etc)
// v2.7.1 : menu = every integer in [0;30]

#include "../Common/baseFuncNad.1.5.h"

void process(TChain* myChain, TString dirOut, int iEff, int iEff_M,
	     TString tagCuts, double tagEt, TString probeCuts, double probeEt, double massCut1, double massCut2,
	     int nEntries, ofstream& outlog, TString tag_hlt , bool checkCharge, bool data2011, 
	     TString eg_menu, bool debug)  {

  // OUTPUT FILE //
  //TFile *outfile = new TFile(dirOut+"tree_effi_TagProbe.root","RECREATE");
  
  // LOG FOR PROBES RUN:EVENT //
  ofstream out_elepairs_json(dirOut+"elepairs.json",ios::out);
  ofstream outNotToZero(dirOut+"notToZero.txt", ios::out);
  ofstream out_probes_trigEG15(dirOut+"probes_trigEG15.txt", ios::out);
  ofstream out_probes_notrigEG15(dirOut+"probes_notrigEG15.txt", ios::out);
  ofstream out_notrigEG15_RLE(dirOut+"probes_notrigEG15_RLE.txt",ios::out);
  ofstream out_notrigEG15_R(dirOut+"probes_notrigEG15_R.txt",ios::out);
  ofstream out_custom_json(dirOut+"inefficiency.json",ios::out);

  out_custom_json << "{" << endl;
  out_elepairs_json << "{" << endl;

  // OUTPUT TREES //
  if(debug) cout << "declare output trees" << endl;

  // L1 EG menu
  vector<int> menu;
  if(eg_menu=="L1_EG_MENU_2010") {
    int EG_2010[10] = {0,2,4,10,16,20,24,30,40,60};
    for(int i=0 ; i<10 ; i++)
      menu.push_back(EG_2010[i]);
  }
  else if(eg_menu=="L1_EG_MENU_2011") {
    for(int i=0; i<31 ; i++)
      menu.push_back(i);
  }
  const int nEG=menu.size();

  TString EGval[nEG];
  ostringstream ossi;
  for(int iEG=0;iEG<nEG;iEG++) {
    ossi.str("");
    ossi << menu[iEG] ;
    EGval[iEG] = ossi.str();
  }
  ////

  const int nECAL = 2; // EB / EE
  const int nColl = 2; // _N / _M

  TString name_ecal[nECAL] = {"_EB","_EE"};
  TString name_coll[nColl] = {"_N","_M"};  

  TTree* treenew[nECAL][nColl];
  float sc_et[nECAL][nColl];
  float sc_dr[nECAL][nColl];
  Int_t l1_bin[nECAL][nColl][nEG];

  int nHLT;
  vector<TString> HLT_names;
  TString HLT_names_2011[10] = {"unbias","EG2","EG5","EG8","EG12","EG12old","Ele15","Ele17","Ele8calIso","Ele8tkIso"} ;
  TString HLT_names_2010[4] = {"unbias","EG2","EG5","EG8"} ;

  HLT_names.clear();
  if(data2011) {
    nHLT = 10;
    for(int i=0 ; i<nHLT ; i++) 
      HLT_names.push_back(HLT_names_2011[i]);
  }
  else {
    nHLT = 4;
    for(int i=0 ; i<nHLT ; i++) 
      HLT_names.push_back(HLT_names_2010[i]);
  }

  TString name, name_scet, name_scdr, name_l1bin;
  for(int iECAL=0 ; iECAL<nECAL ; iECAL++ ) {
    for(int iColl=0 ; iColl<nColl ; iColl++) {

      name = "treenew" + name_ecal[iECAL] + name_coll[iColl];
      treenew[iECAL][iColl] = new TTree(name,name);

      name_scet = "sc_et"+name_ecal[iECAL]+name_coll[iColl];
      name_scdr = "sc_dr"+name_ecal[iECAL]+name_coll[iColl]; 

      treenew[iECAL][iColl]->Branch(name_scet, &sc_et[iECAL][iColl], name_scet+"/F");
      treenew[iECAL][iColl]->Branch(name_scdr, &sc_dr[iECAL][iColl], name_scdr+"/F");	

      for(int iEG=0 ; iEG<nEG ; iEG++) {
	name_l1bin = "l1_"+EGval[iEG]+name_ecal[iECAL]+name_coll[iColl];
	treenew[iECAL][iColl]->Branch(name_l1bin , &l1_bin[iECAL][iColl][iEG] , name_l1bin+"/I");	
      }
    }
  }

  if(debug) cout << "get the input tree" << endl;

  // INPUT TREE //
  int nEvent, nRun, nLumi ;
  int trig_isUnbiased, trig_isL1SingleEG2, trig_isL1SingleEG5, trig_isL1SingleEG8, 
    trig_isL1SingleEG12, trig_isL1SingleEG12_old, trig_is_HLT_Ele15, trig_is_HLT_Ele17, trig_is_HLT_Ele8_calIso, trig_is_HLT_Ele8_tkIso ;   
  
  trig_isUnbiased =  trig_isL1SingleEG2 =  trig_isL1SingleEG5 =  trig_isL1SingleEG8 =  
    trig_isL1SingleEG12 =  trig_isL1SingleEG12_old =  trig_is_HLT_Ele15 =  trig_is_HLT_Ele17 =  trig_is_HLT_Ele8_calIso =  trig_is_HLT_Ele8_tkIso = 0;
  
  // TP info
  int trig_tower_N,trig_tower_ieta[4032],trig_tower_iphi[4032],trig_tower_adc[4032],trig_tower_sFGVB[4032];
  int trig_tower_N_modif,trig_tower_ieta_modif[4032],trig_tower_iphi_modif[4032],trig_tower_adc_modif[4032],trig_tower_sFGVB_modif[4032]; 
  int trig_tower_N_emul,trig_tower_ieta_emul[4032],trig_tower_iphi_emul[4032],trig_tower_adc_emul[4032][5],trig_tower_sFGVB_emul[4032][5]; 

  trig_tower_N = trig_tower_N_modif = trig_tower_N_emul = 0;
  for(int iT=0;iT<4032;iT++) {
    trig_tower_ieta[iT]=trig_tower_iphi[iT]=trig_tower_adc[iT]=trig_tower_sFGVB[iT]=-999;
    trig_tower_ieta_modif[iT]=trig_tower_iphi_modif[iT]=trig_tower_adc_modif[iT]=trig_tower_sFGVB_modif[iT]=-999;
    trig_tower_ieta_emul[iT]=trig_tower_iphi_emul[iT]=-999;
    for(int iEm=0;iEm<5;iEm++)
      trig_tower_adc_emul[iT][iEm]=trig_tower_sFGVB_emul[iT][iEm]=-999;
  }

  // L1 candidates info
  int trig_L1emIso_N, trig_L1emNonIso_N, trig_L1emIso_N_M, trig_L1emNonIso_N_M;
  int trig_L1emIso_ieta[4], trig_L1emIso_iphi[4], trig_L1emIso_rank[4];
  int trig_L1emNonIso_ieta[4], trig_L1emNonIso_iphi[4], trig_L1emNonIso_rank[4];
  int trig_L1emIso_ieta_M[4], trig_L1emIso_iphi_M[4], trig_L1emIso_rank_M[4];
  int trig_L1emNonIso_ieta_M[4], trig_L1emNonIso_iphi_M[4], trig_L1emNonIso_rank_M[4];

  // L1 prefiring
  int trig_preL1emIso_N; 
  int trig_preL1emNonIso_N;
  int trig_preL1emIso_ieta[4], trig_preL1emIso_iphi[4], trig_preL1emIso_rank[4]; 
  int trig_preL1emNonIso_ieta[4], trig_preL1emNonIso_iphi[4],trig_preL1emNonIso_rank[4];
  // L1 postfiring
  int trig_postL1emIso_N; 
  int trig_postL1emNonIso_N;
  int trig_postL1emIso_ieta[4], trig_postL1emIso_iphi[4], trig_postL1emIso_rank[4]; 
  int trig_postL1emNonIso_ieta[4], trig_postL1emNonIso_iphi[4],trig_postL1emNonIso_rank[4];
  
  // Masking
  int trig_nMaskedRCT, trig_nMaskedCh;
  int trig_iMaskedRCTeta[100], trig_iMaskedRCTphi[100], trig_iMaskedRCTcrate[100], trig_iMaskedTTeta[100], trig_iMaskedTTphi[100];


  // Pairs info
  double pair_M;
  double pair_eta[2], pair_sclEta[2], pair_sclEt[2], pair_phi[2], pair_pT[2], pair_eT[2], pair_E[2];
  int pair_cuts[2], pair_HLT_Ele27_cut[2], pair_fidu[2], pair_charge[2], pair_RCTeta[2], pair_RCTphi[2], 
    pair_L1iso[2], pair_L1noniso[2], pair_L1iso_M[2], pair_L1noniso_M[2];
  int pair_RCTetaVect[2][10], pair_RCTphiVect[2][10], 
    pair_L1isoVect[2][10], pair_L1nonisoVect[2][10], 
    pair_L1isoVect_M[2][10], pair_L1nonisoVect_M[2][10];
  int pair_TTetaVect[2][50], pair_TTphiVect[2][50];
  double pair_TTetVect[2][50];
  double pair_RCTetVect[2][10];

  // Global
  myChain->SetBranchAddress("nEvent",&nEvent);
  myChain->SetBranchAddress("nRun",&nRun);
  myChain->SetBranchAddress("nLumi",&nLumi);

  // Trigger
  /*
  myChain->SetBranchAddress("trig_isUnbiased",&trig_isUnbiased);
  myChain->SetBranchAddress("trig_isL1SingleEG2",&trig_isL1SingleEG2);
  myChain->SetBranchAddress("trig_isL1SingleEG5",&trig_isL1SingleEG5);
  myChain->SetBranchAddress("trig_isL1SingleEG8",&trig_isL1SingleEG8);
  */
//   if( data2011 ) {
//     myChain->SetBranchAddress("trig_isL1SingleEG12",&trig_isL1SingleEG12);
//     myChain->SetBranchAddress("trig_isL1SingleEG12_old",&trig_isL1SingleEG12_old);
//     myChain->SetBranchAddress("trig_is_HLT_Ele15",&trig_is_HLT_Ele15);
//     myChain->SetBranchAddress("trig_is_HLT_Ele17",&trig_is_HLT_Ele17);
//     myChain->SetBranchAddress("trig_is_HLT_Ele8_calIso",&trig_is_HLT_Ele8_calIso);
//     myChain->SetBranchAddress("trig_is_HLT_Ele8_tkIso",&trig_is_HLT_Ele8_tkIso);
//   }
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

  // pre/post-firing L1 candidates
  myChain->SetBranchAddress("trig_preL1emIso_N",     &trig_preL1emIso_N);
  myChain->SetBranchAddress("trig_preL1emIso_ieta",  &trig_preL1emIso_ieta);
  myChain->SetBranchAddress("trig_preL1emIso_iphi",  &trig_preL1emIso_iphi);
  myChain->SetBranchAddress("trig_preL1emIso_rank",  &trig_preL1emIso_rank);
  //
  myChain->SetBranchAddress("trig_preL1emNonIso_N",     &trig_preL1emNonIso_N);
  myChain->SetBranchAddress("trig_preL1emNonIso_ieta",  &trig_preL1emNonIso_ieta);
  myChain->SetBranchAddress("trig_preL1emNonIso_iphi",  &trig_preL1emNonIso_iphi);
  myChain->SetBranchAddress("trig_preL1emNonIso_rank",  &trig_preL1emNonIso_rank);
  //
  myChain->SetBranchAddress("trig_postL1emIso_N",     &trig_postL1emIso_N);
  myChain->SetBranchAddress("trig_postL1emIso_ieta",  &trig_postL1emIso_ieta);
  myChain->SetBranchAddress("trig_postL1emIso_iphi",  &trig_postL1emIso_iphi);
  myChain->SetBranchAddress("trig_postL1emIso_rank",  &trig_postL1emIso_rank);
  //
  myChain->SetBranchAddress("trig_postL1emNonIso_N",     &trig_postL1emNonIso_N);
  myChain->SetBranchAddress("trig_postL1emNonIso_ieta",  &trig_postL1emNonIso_ieta);
  myChain->SetBranchAddress("trig_postL1emNonIso_iphi",  &trig_postL1emNonIso_iphi);
  myChain->SetBranchAddress("trig_postL1emNonIso_rank",  &trig_postL1emNonIso_rank);


  // Trigger Towers
  // normal collection
  /*
  myChain->SetBranchAddress("trig_tower_N", &trig_tower_N);
  myChain->SetBranchAddress("trig_tower_ieta",  &trig_tower_ieta);
  myChain->SetBranchAddress("trig_tower_iphi",  &trig_tower_iphi);
  myChain->SetBranchAddress("trig_tower_adc",  &trig_tower_adc);
  myChain->SetBranchAddress("trig_tower_sFGVB",  &trig_tower_sFGVB);
 
  // modified collection
  myChain->SetBranchAddress("trig_tower_N_emul", &trig_tower_N_emul);
  myChain->SetBranchAddress("trig_tower_ieta_emul",  &trig_tower_ieta_emul);
  myChain->SetBranchAddress("trig_tower_iphi_emul",  &trig_tower_iphi_emul);
  myChain->SetBranchAddress("trig_tower_adc_emul",  &trig_tower_adc_emul);
  myChain->SetBranchAddress("trig_tower_sFGVB_emul",  &trig_tower_sFGVB_emul);

  // emulated collection
  myChain->SetBranchAddress("trig_tower_N_emul", &trig_tower_N_emul);
  myChain->SetBranchAddress("trig_tower_ieta_emul",  &trig_tower_ieta_emul);
  myChain->SetBranchAddress("trig_tower_iphi_emul",  &trig_tower_iphi_emul);
  myChain->SetBranchAddress("trig_tower_adc_emul",  &trig_tower_adc_emul);
  myChain->SetBranchAddress("trig_tower_sFGVB_emul",  &trig_tower_sFGVB_emul);
  */

  // Masking
  myChain->SetBranchAddress("trig_nMaskedRCT",      &trig_nMaskedRCT);      
  myChain->SetBranchAddress("trig_iMaskedRCTeta",   &trig_iMaskedRCTeta);                                          
  myChain->SetBranchAddress("trig_iMaskedRCTcrate", &trig_iMaskedRCTcrate);
  myChain->SetBranchAddress("trig_iMaskedRCTphi",   &trig_iMaskedRCTphi);
  myChain->SetBranchAddress("trig_nMaskedCh",       &trig_nMaskedCh);    
  myChain->SetBranchAddress("trig_iMaskedTTeta",    &trig_iMaskedTTeta);   
  myChain->SetBranchAddress("trig_iMaskedTTphi",    &trig_iMaskedTTphi);      	


  // Pair Info
  myChain->SetBranchAddress("pair_M",&pair_M);
  //
  myChain->SetBranchAddress("pair_cuts",&pair_cuts); 
  // 0 : noCut | 1 : VBTF 95 | 2 : VBTF 80 | 3 : VBTF 60
  myChain->SetBranchAddress("pair_HLT_Ele27_cut",&pair_HLT_Ele27_cut); 
  myChain->SetBranchAddress("pair_fidu",&pair_fidu);
  //
  myChain->SetBranchAddress("pair_eta",&pair_eta);
  myChain->SetBranchAddress("pair_sclEta",&pair_sclEta);
  myChain->SetBranchAddress("pair_phi",&pair_phi);
  myChain->SetBranchAddress("pair_RCTeta",&pair_RCTeta);
  myChain->SetBranchAddress("pair_RCTphi",&pair_RCTphi);
  //
  myChain->SetBranchAddress("pair_charge",&pair_charge);
  myChain->SetBranchAddress("pair_pT",&pair_pT);
  myChain->SetBranchAddress("pair_eT",&pair_eT);
  myChain->SetBranchAddress("pair_sclEt",&pair_sclEt);
  myChain->SetBranchAddress("pair_E",&pair_E);

  myChain->SetBranchAddress("pair_TTetaVect",&pair_TTetaVect);
  myChain->SetBranchAddress("pair_TTphiVect",&pair_TTphiVect);
  myChain->SetBranchAddress("pair_TTetVect",&pair_TTetVect);  

  myChain->SetBranchAddress("pair_L1iso",&pair_L1iso);
  myChain->SetBranchAddress("pair_L1noniso",&pair_L1noniso);
  myChain->SetBranchAddress("pair_L1iso_M",&pair_L1iso_M);
  myChain->SetBranchAddress("pair_L1noniso_M",&pair_L1noniso_M);
  //
  myChain->SetBranchAddress("pair_RCTetVect",&pair_RCTetVect);
  myChain->SetBranchAddress("pair_RCTetaVect",&pair_RCTetaVect);
  myChain->SetBranchAddress("pair_RCTphiVect",&pair_RCTphiVect);
  myChain->SetBranchAddress("pair_L1isoVect",&pair_L1isoVect);
  myChain->SetBranchAddress("pair_L1nonisoVect",&pair_L1nonisoVect);
  myChain->SetBranchAddress("pair_L1isoVect_M",&pair_L1isoVect_M);
  myChain->SetBranchAddress("pair_L1nonisoVect_M",&pair_L1nonisoVect_M);

  // DEFINITION OF THE TAG AND THE PROBE //
  int iTag, iProbe;
  vector<int> goodTags, goodProbes;
  int cutsDef[2]; // 0:tag ; 1:probe
  bool passingCuts[2][2]; // [iEle][iTPDefCuts] iTPDefCuts=0(1) <-> Tag(Probe)
  TString cutsName[5] = {"noCuts","WP95","WP80","WP60","HLT_Ele27"};
  double etDef[2] = {tagEt,probeEt};
  TString asked[2] = {tagCuts,probeCuts};
  for(int iEle=0 ; iEle<2 ; iEle++) {
    for(int iCut=0 ; iCut<5 ; iCut++) {
      if(asked[iEle]==cutsName[iCut]) cutsDef[iEle] = iCut ;
    }
  }


  // USEFUL VARIABLES //
  int b_hlt, nCurrentRun,nCurrentEvt;
  nCurrentRun = nCurrentEvt = 0;
  //bool isGoodRun,isGoodEvent,b_notToZero;
  TString filename = "";
  TString flag = "";
  vector<int> IdxPt;
  TString nameTagProbe[2] = {"TAG","PROBE"} ;

  MAPTT adcTT;
  MAPTT::iterator iterTT;
  pair<int,int> coords;
  vector<int> ietaTT, iphiTT;
  //int tOri, tM, adc, sFGVB;
  //int sFGVB_M[5], adc_emul[5];

  // TRIGGERING PROPERTIES //
  vector<int> trigHLT;
  vector<int> firedEG[2][2]; // [N/M][tag/probe]
  int iTrig[2][2]; // [N/M][tag/probe]
  int tagEG;

  // COUNTERS //
  int nCases[3] = {0,0,0}; // all ; same ; different
  int nZeroStatus[3] = {0,0,0}; // ToBeZeroed ; Zeroed ; NotZeroed
  int nEB_NotToZero_Zeroed = 0;
  int nEE_NotToZero_Zeroed = 0;
  int n_fails_fidu = 0;

  int nEle[nECAL][nColl], nMatched[nECAL][nColl][nEG], nLost[nECAL][nColl][nEG], nAppeared[nECAL][nColl][nEG];
  for(int iECAL=0 ; iECAL<nECAL ; iECAL++ ) {
    for(int iColl=0 ; iColl<nColl ; iColl++) {
      nEle[iECAL][iColl] = 0;
      for(int iEG=0 ; iEG<nEG ; iEG++) {
	nMatched[iECAL][iColl][iEG] = nLost[iECAL][iColl][iEG] = nAppeared[iECAL][iColl][iEG] = 0;
      }
    }
  }

  int nProbesTrig15, nProbesNoTrig15, nCurrentRun_ineff, nCurrentLumi_ineff,nCurrentRun_elepairs, nCurrentLumi_elepairs;
  nProbesTrig15 = nProbesNoTrig15 = 0;
  nCurrentRun_ineff = nCurrentLumi_ineff = -1;
  nCurrentRun_elepairs = nCurrentLumi_elepairs = -1;

  // Z MASS PEAK PLOT //
  TH1F* h_Mee[3];
  TString nameZ[3] = {"h_EBEB","h_EEEE","h_EBEE"};
  for(int i=0;i<3;i++) {
    h_Mee[i] = new TH1F(nameZ[i],nameZ[i],80,50,130);
  }

  // ETA/PHI EFFICIENCY MAP
  TH2F* h_map_trig[nEG];
  TH2F* h_map_all;

  for(int iEG=0 ; iEG<nEG ; iEG++) {
    h_map_trig[iEG] = new TH2F("h_map_trig_EG"+EGval[iEG],"eta/phi map ele trig EG"+EGval[iEG],100,-2.5,2.5,100,-3.14,3.14);
  }
  h_map_all = new TH2F("h_map_all","eta/phi map all ele",100,-2.5,2.5,100,-3.14,3.14);

  // ETA/ET PHI/ET ETA/PHI plots
  TH2F* h_Map[3];
  TString name_Map[3] = {"h_EtaEt","h_PhiEt","h_EtaPhi"};
  h_Map[0] = new TH2F(name_Map[0],name_Map[0],50,-2.5,2.5,50,0,100);
  h_Map[1] = new TH2F(name_Map[1],name_Map[1],63,-3.14,3.14,50,0,100);
  h_Map[2] = new TH2F(name_Map[2],name_Map[2],50,-2.5,2.5,63,-3.14,3.14);

  TH1F* h_Map1D[3];
  TString name_Map1D[3] = {"h_Eta","h_Phi","h_Et"};
  h_Map1D[0] = new TH1F(name_Map1D[0],name_Map1D[0],50,-2.5,2.5); // eta
  h_Map1D[1] = new TH1F(name_Map1D[1],name_Map1D[1],63,-3.14,3.14); // phi
  h_Map1D[2] = new TH1F(name_Map1D[2],name_Map1D[2],50,0,100); // eT
 

  // JSON FILE READER //
  
  string jsonDir = "/data_CMS/cms/ndaci/" ;
  const int nJson = 9;
  string jsonFile[nJson]; // 2010B, May10, Aug05, Prompt
  map<int, vector<pair<int, int> > > jsonMap[nJson] ;

  jsonFile[0] = jsonDir + "ndaci_2011A/JSON/goodrunlist_json.txt" ;
  jsonFile[1] = jsonDir + "ndaci_2011A/JSON/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.txt" ;
  jsonFile[2] = jsonDir + "ndaci_2011A/JSON/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v3.txt" ;
  jsonFile[3] = jsonDir + "ndaci_2011A/JSON/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt" ;
  jsonFile[4] = jsonDir + "ndaci_2011A/JSON/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt" ;

  // 2012
  jsonFile[5] = jsonDir + "ndaci_2012/JSON/Cert_190456-203002_8TeV_PromptReco_Collisions12_JSON.txt";
  jsonFile[6] = jsonDir + "ndaci_2012/JSON/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt";
  jsonFile[7] = jsonDir + "ndaci_2012/JSON/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt";
  jsonFile[8] = jsonDir + "ndaci_2012/JSON/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt";

  for(int i=0 ; i<nJson ; i++)
    jsonMap[i] = readJSONFile(jsonFile[i]);


  // LOOP OVER PAIRS //
  
  int numEntries = myChain->GetEntries () ;
  int nProcess = numEntries;
  if(nEntries>=0 && nEntries<numEntries)
    nProcess = nEntries;
  outlog << "will process " << nProcess << "/" << numEntries << " pairs" << endl;

  for (int iPair = 0 ; iPair < nProcess ; iPair++ ) { 
    
    if(debug) cout << "is getting the entry #" << iPair << endl;
    myChain->GetEntry (iPair) ;
    if(debug) cout << "got it" << endl;

    // Get HLT information //
    int trigHLT_2011[10] = {trig_isUnbiased, trig_isL1SingleEG2, trig_isL1SingleEG5, trig_isL1SingleEG8,
			    trig_isL1SingleEG12, trig_isL1SingleEG12_old, trig_is_HLT_Ele15, trig_is_HLT_Ele17,
			    trig_is_HLT_Ele8_calIso, trig_is_HLT_Ele8_tkIso} ;
    int trigHLT_2010[4] = {trig_isUnbiased, trig_isL1SingleEG2, trig_isL1SingleEG5, trig_isL1SingleEG8} ;

    trigHLT.clear();
    for(int iHLT=0 ; iHLT<nHLT ; iHLT++) {
      if(data2011) 
	trigHLT.push_back( trigHLT_2011[iHLT] ) ;
      else 
	trigHLT.push_back( trigHLT_2010[iHLT] ) ;
    }

    // LOG : FILE / RUN / EVENT / HLT / L1 //
    if( iPair==0 || filename != myChain->GetFile()->GetName() ) {
      filename = myChain->GetFile()->GetName() ;
      outlog << "<------ File : " << filename << " ------>" << endl;
    }

    if( iPair==0 || nRun!=nCurrentRun ) {
      nCurrentRun = nRun ;
      outlog << endl << "---- nRun=" << nRun << " ----" << endl ;
    }
    

    // HLT SELECTION //
    b_hlt=1;
    if(tag_hlt=="HLT_EG2")
      b_hlt = trig_isL1SingleEG2 ;
    if(tag_hlt=="HLT_EG5")
      b_hlt = trig_isL1SingleEG5 ;
    if(tag_hlt=="HLT_EG8")
      b_hlt = trig_isL1SingleEG8 ;
    if(tag_hlt=="HLT_unbias")
      b_hlt = trig_isUnbiased ;
    if(tag_hlt=="noHLT")
      b_hlt = 1;
    if(b_hlt != 1) {
      //if( iPair==0 || nEvent!=nCurrentEvt ) { 
      //outlog << "-- nEvent=" << nEvent << " failed HLT selection" << endl;
      //}
      continue;
    }

    if(debug) cout << "passed HLT selection" << endl;

    // MASS CUT //
    if( massCut1 <= massCut2 ) {
      if( pair_M < massCut1 || pair_M > massCut2 ) continue;
    }
    else {
      if( pair_M < massCut2 || pair_M > massCut1 ) continue;
    }

    if(debug) cout << "passed mass cut" << endl;

    // CHARGE OPPOSITION //
    if( checkCharge )
      if( pair_charge[0]*pair_charge[1] != -1 ) continue;

    if(debug) cout << "passed opposite-charge selection" << endl;
    
    // CHECK TRIGGERING PROPERTIES //
    for(int iEle=0 ; iEle<2 ; iEle++) {
      for(int iColl=0 ; iColl<nColl ; iColl++) {
	firedEG[iColl][iEle].clear();
	firedEG[iColl][iEle].resize(nEG,0);
      }
    }

    //pair<int,int> pair_trig_coords_N[2], pair_trig_coords_M[2];
    //int iL1cand_trig_N, iL1cand_trig_M;

    for(int iEle=0 ; iEle<2 ; iEle++) {
      for(int iR=0 ; iR < 10 ; iR++) {

	globalFireL1_Normal( pair_RCTetVect[iEle][iR] , pair_L1nonisoVect[iEle][iR],
			     pair_L1isoVect[iEle][iR], firedEG[0][iEle], menu); 
	
	globalFireL1_Normal( pair_RCTetVect[iEle][iR] , pair_L1nonisoVect_M[iEle][iR],
			     pair_L1isoVect_M[iEle][iR], firedEG[1][iEle], menu); 
	
	if(debug) cout << pair_RCTetVect[iEle][iR] << " " << pair_L1nonisoVect[iEle][iR] << " " << pair_L1isoVect[iEle][iR] << endl;

      }
    }

    if(debug) {
      for(int iColl=0 ; iColl<nColl ; iColl++) {
	for(int iEle=0 ; iEle<2 ; iEle++) {
	  cout << "Ele" << iEle << " :: " ;
	  for(int i=0 ; i<(int)firedEG[iColl][iEle].size() ; i++) {
	    cout << "EG" << i << " : " << firedEG[iColl][iEle][i] << " | ";
	  }
	  cout << endl;
	}      
      }
    }

    for(int iColl=0 ; iColl<nColl ; iColl++)
      for(int iEle=0 ; iEle<2 ; iEle++)
	iTrig[iColl][iEle] = 0;
    
    for(int iColl=0 ; iColl<nColl ; iColl++)    
      for(int iEle=0 ; iEle<2 ; iEle++)
	for(int iEG=0 ; iEG<nEG ; iEG++)
	  if(firedEG[iColl][iEle][iEG]==1) iTrig[iColl][iEle]=iEG;

    if(debug) cout << "TRIG : ele0=" << iTrig[0][0] << " | ele1=" << iTrig[0][1] << endl;

    // EG that the tag have to trigger
    /*
    if(nRun==161217) tagEG = 6; // min L1 seed in HLT menu : EG12
    else tagEG = 7;             // ... : EG15
    */
    tagEG=15; // EG15

    // DECIDE WHO IS THE TAG/PROBE //
    iTag = iProbe = -1 ;
    goodTags.clear();
    goodProbes.clear();
    
    //bool passingCuts[2][2]; // [iEle][iTPDefCuts] iTPDefCuts=0(1) <-> Tag(Probe)
    for(int i=0 ; i<2 ; i++)
      for(int j=0 ; j<2 ; j++)
	passingCuts[i][j] = false;

    for(int iEle=0 ; iEle<2 ; iEle++) {
      for(int iTP=0 ; iTP<2 ; iTP++) {
	if( cutsDef[iTP]==4 ) {
	  if( pair_HLT_Ele27_cut[iEle]==1 ) { 
	    if( pair_fidu[iEle]==1 ) {
	      passingCuts[iEle][iTP]=true; 
	      if(debug) cout<<"pass cuts" << endl; 
	    }
	  }
	}
	else {
	  if( pair_cuts[iEle] >= cutsDef[iTP] ) { 
	    if( pair_fidu[iEle]==1 ) {
	      passingCuts[iEle][iTP]=true; 
	      if(debug) cout<<"pass cuts"<<endl; 
	    }
	  }
	}
      }
    }

//     for( int iEle=0 ; iEle<2 ; iEle++ ) {
//       if ( passingCuts[iEle][0] && pair_sclEt[iEle] >= etDef[0] && iTrig[iEle] >= tagEG ) goodTags.push_back(iEle);
//       if ( passingCuts[iEle][1] && pair_sclEt[iEle] >= etDef[1] ) goodProbes.push_back(iEle);
//     }

    for( int iEle=0 ; iEle<2 ; iEle++ ) {
      if ( passingCuts[iEle][0] ) {
	if(debug) cout << "ele pass cuts tag" << endl;
	if(pair_sclEt[iEle] >= etDef[0]) {
	  if(debug) cout << "et is ok : " << pair_sclEt[iEle] << endl; 
	  if(iTrig[0][iEle] >= tagEG ) {
	    if(debug) cout << "triggering " << iTrig[0][iEle] << endl;
	    goodTags.push_back(iEle);
	  }
	}
      }
      
      if ( passingCuts[iEle][1] ) {
	if(debug) cout << "ele pass cuts probe" << endl;
	if(pair_sclEt[iEle] >= etDef[1] ) {
	  if(debug) cout << "et is ok : " << pair_sclEt[iEle] << endl;
	  goodProbes.push_back(iEle);
	}
      }
    }
   
    if(debug) cout << "goodTags.size()=" << goodTags.size() << " | goodProbes.size()=" << goodProbes.size() << endl;

    if( goodTags.size()==0 || goodProbes.size()==0 ) continue; // no tag or no probe

    // Look every cases
    for(int iCaseT=0 ; iCaseT<(int)goodTags.size() ; iCaseT++) {

      iTag = goodTags[iCaseT];

      for(int iCaseP=0 ; iCaseP<(int)goodProbes.size() ; iCaseP++) {
	iProbe = goodProbes[iCaseP];
	if(iProbe==iTag) continue;
	
	// CUSTOM ELEPAIRS JSON
	if( nRun!=nCurrentRun_elepairs ) {
	  if(nCurrentRun_elepairs!=-1) {
	    out_elepairs_json << "]," << endl;
	  }
	  nCurrentLumi_elepairs=-1;
	  out_elepairs_json << '"' << nRun << '"' << ": [" ;
	  nCurrentRun_elepairs = nRun;
	}
	if( nLumi!=nCurrentLumi_elepairs ) {
	  if( nCurrentLumi_elepairs!=-1 ) {
	    out_elepairs_json << "," ;
	  }
	  out_elepairs_json << "[" << nLumi << "," << nLumi << "]" ;
	  nCurrentLumi_elepairs=nLumi;
	}

    
	// FILL TRIGGER TREE FOR EFFICIENCY ///////////////////////
	int idxECAL=-1;
	if( pair_eta[iProbe]<1.479 ) idxECAL=0;
	else idxECAL = 1;
	
	for(int iColl=0 ; iColl<nColl ; iColl++) {
	  sc_et[idxECAL][iColl] = pair_sclEt[iProbe];
	  sc_dr[idxECAL][iColl] = 0;
	  
	  if(pair_L1iso[iProbe]>0 || pair_L1noniso[iProbe]>0)
	    sc_dr[idxECAL][iColl] = 1;
	  
	  for(int iEG=0 ; iEG<nEG ; iEG++)
	    l1_bin[idxECAL][iColl][iEG]  = firedEG[iColl][iProbe][iEG];
	}
	
	for(int iColl=0 ; iColl<nColl ; iColl++)
	  treenew[idxECAL][iColl]->Fill();
	////////////////////////////////////////////////////////////

	// FILL EFFICIENCY MAP TOOLS
	h_map_all->Fill( pair_eta[iProbe] , pair_phi[iProbe] );
	for(int iEG=0 ; iEG<nEG ; iEG++) {
	  if( firedEG[0][iProbe][iEG]==1 ) {
	    h_map_trig[iEG]->Fill( pair_eta[iProbe] , pair_phi[iProbe] );
	  }
	}

	// Z MASS PEAK PLOT // EBEB EEEE EBEE
	if( abs(pair_eta[iProbe])<1.479 && abs(pair_eta[iTag])<1.479 ) h_Mee[0]->Fill(pair_M);
	if( abs(pair_eta[iProbe])>1.479 && abs(pair_eta[iTag])>1.479 ) h_Mee[1]->Fill(pair_M);
	else h_Mee[2]->Fill(pair_M);

	// OUTPUT JSON //
	if( nRun!=nCurrentRun_ineff ) {
	  if(nCurrentRun_ineff!=-1) {
	    out_custom_json << "]," << endl;
	  }
	  nCurrentLumi_ineff=-1;
	  out_custom_json << '"' << nRun << '"' << ": [" ;
	  nCurrentRun_ineff = nRun;
	}
	if( nLumi!=nCurrentLumi_ineff ) {
	  if( nCurrentLumi_ineff!=-1 ) {
	    out_custom_json << "," ;
	  }
	  out_custom_json << "[" << nLumi << "," << nLumi << "]" ;
	  nCurrentLumi_ineff=nLumi;
	}


	// LOG PROBES RUN:EVENT //
	if( pair_sclEt[iProbe]>15. && firedEG[0][iProbe][7]==iEff ) {
	  out_probes_trigEG15 << "'" << nRun << ":" << nEvent << "'" ;
	  //out_probes_trigEG15 << "'" << nRun ;
	  if(iPair < nProcess-1 ) out_probes_trigEG15 << "," ;
	  nProbesTrig15++ ;
    
	  //else if( pair_sclEt[iProbe]>15. && firedEG_N[iProbe][7]==0) {
      
	  //       if( nRun!=nCurrentRun_ineff ) {
	  // 	if(nCurrentRun_ineff!=-1) {
	  // 	  out_custom_json << "]," << endl;
	  // 	}
	  // 	nCurrentLumi_ineff=-1;
	  // 	out_custom_json << '"' << nRun << '"' << ": [" ;
	  // 	nCurrentRun_ineff = nRun;
	  //       }
	  //       if( nLumi!=nCurrentLumi_ineff ) {
	  // 	if( nCurrentLumi_ineff!=-1 ) {
	  // 	  out_custom_json << "," ;
	  // 	}
	  // 	out_custom_json << "[" << nLumi << "," << nLumi << "]" ;
	  // 	nCurrentLumi_ineff=nLumi;
	  //       }

	  out_probes_notrigEG15 << "'" << nRun << ":" << nEvent << "'" ;
	  out_notrigEG15_RLE << "'" << nRun << ":" << nLumi << ":" << nEvent << "'" ;
	  out_notrigEG15_R  << nRun ;
	  //out_probes_notrigEG15 << "'" << nRun ;
	  if(iPair < nProcess-1 ) {
	    out_probes_notrigEG15 << "," ;
	    out_notrigEG15_RLE << "," ;
	    out_notrigEG15_R << "," ;
	  }
	  nProbesNoTrig15++ ;
	}

	// COUNTING //
	flag = "" ;
	for(int iECAL=0 ; iECAL<nECAL ; iECAL++ ) {
	  for(int iColl=0 ; iColl<nColl ; iColl++) {
	    nEle[iECAL][iColl]++ ;
	    for(int iEG=0 ; iEG<nEG ; iEG++) {
	      if(firedEG[0][iProbe][iEG]==1)  nMatched[iECAL][0][iEG]++ ;
	      /*
	      if(firedEG_M[iProbe][iEG]==1) nMatched_M[IdxPt[k]][iEG]++ ;
	      if(firedEG_N[iProbe][iEG]==1 && firedEG_M[iProbe][iEG]==0) {
	      nLost[IdxPt[k]][iEG]++ ;
	      if(k==0)
	      flag += " lostEG"+EGval[iEG];
	      }
	      if(firedEG_N[iProbe][iEG]==0 && firedEG_M[iProbe][iEG]==1) {
	      nAppeared[IdxPt[k]][iEG]++ ;
	      if(k==0)
	      flag += " appearedEG"+EGval[iEG];
	      }
	    */
	    }
	  }
	}
  
	// MAP THE TP //
	adcTT.clear();
	for(int t=0 ; t<trig_tower_N ; t++) {
	  coords = make_pair( trig_tower_ieta[t] , trig_tower_iphi[t] );
	  adcTT[coords].first = t ;
	}
	for(int t=0 ; t<trig_tower_N_emul ; t++) {
	  coords = make_pair( trig_tower_ieta_emul[t] , trig_tower_iphi_emul[t] );
	  iterTT = adcTT.find( coords );
	  if( iterTT != adcTT.end() )
	    adcTT[coords].second = t;
	  else {
	    outlog << "mapping problem" << endl;
	    adcTT[coords].first = -777;
	    adcTT[coords].second = t;
	  }
	}


	// PUT DETAILS AND COUNTERS IN LOGS //    
	if( pair_sclEt[iProbe]>15. && firedEG[0][iProbe][7]==iEff && firedEG[1][iProbe][7]==iEff_M) {

	  if( iPair==0 || nEvent!=nCurrentEvt ) {

	    outlog << "-----------------------------------------------------------------------------------------"
		   << endl;

	    nCurrentEvt = nEvent ;
	    outlog << endl
		   << "-- nEvent=" << nEvent << " -- nLumi=" << nLumi << " --" << endl
		   << "HLT : " ;
	    for(int iHLT=0 ; iHLT<nHLT ; iHLT++) {
	      outlog << HLT_names[iHLT] << ":" << trigHLT[iHLT] << " | " ;
	    }


	    outlog << endl
		   << "ISO            NONISO        ISO_M         NONISO_M         preISO       preNonIso       postIso       postNonIso" 
		   << endl;
	    for(int iL1=0 ; iL1<4 ; iL1++) {
	      outlog << "(" << trig_L1emIso_ieta[iL1]      << ";" << trig_L1emIso_iphi[iL1]      << ")=" << trig_L1emIso_rank[iL1]      << "      "
		     << "(" << trig_L1emNonIso_ieta[iL1]   << ";" << trig_L1emNonIso_iphi[iL1]   << ")=" << trig_L1emNonIso_rank[iL1]   << "      "
		     << "(" << trig_L1emIso_ieta_M[iL1]    << ";" << trig_L1emIso_iphi_M[iL1]    << ")=" << trig_L1emIso_rank_M[iL1]    << "      "
		     << "(" << trig_L1emNonIso_ieta_M[iL1] << ";" << trig_L1emNonIso_iphi_M[iL1] << ")=" << trig_L1emNonIso_rank_M[iL1] << "      ";

	      if(trig_preL1emIso_N>0)
		outlog << "(" << trig_preL1emIso_ieta[iL1]   << ";" << trig_preL1emIso_iphi[iL1]   << ")=" << trig_preL1emIso_rank[iL1]   << "      ";

	      if(trig_preL1emNonIso_N>0)
		outlog << "(" << trig_preL1emNonIso_ieta[iL1]<< ";" << trig_preL1emNonIso_iphi[iL1]<< ")=" << trig_preL1emNonIso_rank[iL1]<< "      ";

	      if(trig_postL1emIso_N>0)
		outlog << "(" << trig_postL1emIso_ieta[iL1]   << ";" << trig_postL1emIso_iphi[iL1]   << ")=" << trig_postL1emIso_rank[iL1]<< "      ";

	      if(trig_postL1emNonIso_N>0)
		outlog << "(" << trig_postL1emNonIso_ieta[iL1]<< ";" << trig_postL1emNonIso_iphi[iL1]<< ")=" << trig_postL1emNonIso_rank[iL1] << "      ";
	    
	      outlog  << endl;
	    }
	
	    if( trig_preL1emIso_N>0 || trig_preL1emNonIso_N>0 || trig_postL1emIso_N>0 || trig_postL1emNonIso_N>0 )
	      outlog << "PRE-POST-FIRING" << endl;
	    outlog << endl;
	  }
           
	  // PAIR INFORMATION
	  //outlog << "PPPPPPPPPPPPPPPPP" << endl;
	  outlog << "Pair M=" << pair_M << endl;
      
	  // COMPARISON ONLINE/EMULATOR
	  nCases[0]++ ;
	  if(iTrig[0][iProbe]==iTrig[1][iProbe]) {
	    nCases[1]++ ;
	    outlog << "ONLINE/EMULATOR = SAME" << endl;
	  }
	  else {
	    nCases[2]++ ;
	    outlog << "ONLINE/EMULATOR = DIFFERENT" << endl;
	  }
      
	  int IdxTagProbe[2] = {iTag , iProbe} ;
	  //TString nameTagProbe[2] = {"TAG","PROBE"} ;

	  // FILLING THE MAPS
	  h_Map[0]->Fill( pair_eta[IdxTagProbe[1]] , pair_eT[IdxTagProbe[1]] ); // eta/et
	  h_Map[1]->Fill( pair_phi[IdxTagProbe[1]] , pair_eT[IdxTagProbe[1]] ); // phi/et
	  h_Map[2]->Fill( pair_eta[IdxTagProbe[1]] , pair_phi[IdxTagProbe[1]] ); // eta/phi

	  h_Map1D[0]->Fill( pair_eta[IdxTagProbe[1]] ); // eta
	  h_Map1D[1]->Fill( pair_phi[IdxTagProbe[1]] ); // phi
	  h_Map1D[2]->Fill( pair_eT[IdxTagProbe[1]] ); // et

	  for(int i=0 ; i<2 ; i++) {
	    outlog << nameTagProbe[i] << " : eta=" << pair_eta[IdxTagProbe[i]] << " phi=" << pair_phi[IdxTagProbe[i]] 
		   << " pT=" << pair_pT[IdxTagProbe[i]] << " eT=" << pair_eT[IdxTagProbe[i]]
		   << " RCTeta=" << pair_RCTeta[IdxTagProbe[i]] << " RCTphi=" << pair_RCTphi[IdxTagProbe[i]]
		   << " cuts=" << cutsName[ pair_cuts[IdxTagProbe[i]] ]
		   << endl
		   << "Fired EG : " ;
	    for(int iEG=0 ; iEG<(int)firedEG[0][IdxTagProbe[i]].size() ; iEG++) {
	      if(firedEG[0][IdxTagProbe[i]][iEG]==1) outlog << "EG" << EGval[iEG] << " " ;
	    }
	    outlog << endl;
	
	    outlog << "Fired EG M : " ;
	    for(int iEG=0 ; iEG<(int)firedEG[1][IdxTagProbe[i]].size() ; iEG++) {
	      if(firedEG[1][IdxTagProbe[i]][iEG]==1) outlog << "EG" << EGval[iEG] << " " ;
	    }
	    outlog << endl;
	
	  }

	  // LOCALIZATION
	  if( pair_eta[IdxTagProbe[1]] < 1.479 ) outlog << "PROBE IN EB" << endl;
	  else outlog << "PROBE IN EE" << endl;
	  outlog << endl;

	  // FIDUCIAL CUT
	  if ( fabs(pair_eta[IdxTagProbe[1]]) >= 2.5 || (fabs(pair_eta[IdxTagProbe[1]])< 1.566 && fabs(pair_eta[IdxTagProbe[1]])>1.4442)) {
	    outlog << "FAILS THE FIDUCIAL CUT" << endl;
	    n_fails_fidu++ ;
	  }

	  // Masking Information
	  outlog << "MASKING" << endl;
	  for(int iM=0 ; iM<trig_nMaskedCh ; iM++) {
	
	    outlog << "Reta="  << getGCTRegionEta(trig_iMaskedTTeta[iM])
		   << " Rphi=" << getGCTRegionPhi(trig_iMaskedTTphi[iM])
		   << " Teta=" << trig_iMaskedTTeta[iM]
		   << " Tphi=" << trig_iMaskedTTphi[iM]
		   << endl;
	  }
	  outlog << endl;

	  outlog << endl 
		 << "------------------------------------------------------------------------------------------------"
		 << endl;
	} //endif probe inefficient

      } // loop over goodProbes
    } // loop over goodTags
    
  } // loop over pairs

  // OUTPUT FILE //
  TFile *outfile = new TFile(dirOut+"tree_effi_TagProbe.root","RECREATE");
  for(int iECAL=0 ; iECAL<nECAL ; iECAL++ ) {
    for(int iColl=0 ; iColl<nColl ; iColl++) {
      treenew[iECAL][iColl]->Write();
    }
  }
  outfile->Close();
  cout << "...filled !" << endl;

  // PLOTTING Z MASS
  cout << "Plotting Z mass peaks..." << endl;
  TCanvas* c_Mee[3];
  TString nameC[3] = {"c_ECAL","c_EB","c_EE"};
  for(int i=0 ; i<3 ; i++) {
    c_Mee[i] = new TCanvas(nameC[i],nameC[i],0,0,800,600);
    if(debug) cout << "TCanvas created" << endl;
    h_Mee[i]->Draw();
    if(debug) cout << "h_Mee drawn" << endl;
    c_Mee[i]->Print(dirOut+"mass_Zee_"+nameC[i]+".gif","gif");
  }
  cout << "...plotted !" << endl;

  /*
  TFile *outplot = new TFile(dirOut+"mass_Z_plots.root","RECREATE");
  for(int i=0 ; i<3 ; i++)
    h_Mee[i]->Write();
  outplot->Close();
  */

  // PLOTTING MAPS

  cout << "Plotting maps" << endl;
  TCanvas* c_Maps[3];
  TString name_c_Maps[3] = {"_EtaEt","_PhiEt","_EtaPhi"};
  gStyle->SetPalette(1);
  for(int iC=0 ; iC<3 ; iC++) {
    c_Maps[iC] = new TCanvas("c"+name_c_Maps[iC],"c"+name_c_Maps[iC],0,0,800,600);
    h_Map[iC]->Draw("colz");
    c_Maps[iC]->Print(dirOut+"maps"+name_c_Maps[iC]+".gif","gif");
  }

  TCanvas* c_Maps1D[3];
  TString name_c_Maps1D[3] = {"_Eta","_Phi","_eT"};
  gStyle->SetPalette(1);
  for(int iC=0 ; iC<3 ; iC++) {
    c_Maps1D[iC] = new TCanvas("c"+name_c_Maps1D[iC],"c"+name_c_Maps1D[iC],0,0,800,600);
    h_Map1D[iC]->Draw("");
    c_Maps1D[iC]->Print(dirOut+"maps"+name_c_Maps1D[iC]+".gif","gif");
  }
  cout << "plotted" << endl;


  // Efficiency eta/phi map
  double n_ele_all[10000];
  double n_ele_trig[10000][nEG];
  double bin_eta[10000], bin_phi[10000], effi_val[10000];

  for(int iBinEta=0 ; iBinEta<100 ; iBinEta++) {
    for(int iBinPhi=0 ; iBinPhi<100 ; iBinPhi++) {

      bin_eta[ iBinEta*100 + iBinPhi ] = -2.5  + iBinEta*(5./100.);
      bin_phi[ iBinEta*100 + iBinPhi ] = -3.14 + iBinPhi*(6.28/100.);

      n_ele_all[ iBinEta*100 + iBinPhi ] = (double)( h_map_all->GetBinContent(iBinEta+1,iBinPhi+1) );

      for(int iEG=0 ; iEG<nEG ; iEG++) {
	n_ele_trig[ iBinEta*100 + iBinPhi ][iEG] = (double)( h_map_trig[iEG]->GetBinContent(iBinEta+1,iBinPhi+1) );

	if( n_ele_all[ iBinEta*100 + iBinPhi ] != 0 ) {
	  effi_val[ iBinEta*100 + iBinPhi ] = n_ele_trig[ iBinEta*100 + iBinPhi ][iEG] / n_ele_all[ iBinEta*100 + iBinPhi ];
	}
	else {
	  effi_val[ iBinEta*100 + iBinPhi ] = -0.1;
	}

      }

    }
  }


  TGraph2D* g_effi_map[nEG];
  TCanvas*  c_effi_map[nEG];

  //for(int iEG=0 ; iEG<nEG ; iEG++) {
  int jEG=15;
  g_effi_map[jEG] = new TGraph2D("g_effi_map_EG"+EGval[jEG],"Efficiency Eta/Phi Map : EG"+EGval[jEG],10000,bin_eta,bin_phi,effi_val);
  c_effi_map[jEG] = new TCanvas("c_effi_map_EG"+EGval[jEG],"Efficiency Eta/Phi Map : EG"+EGval[jEG],0,0,800,600);
  gStyle->SetPalette(1);
  g_effi_map[jEG]->Draw("colz");
  c_effi_map[jEG]->Print(dirOut+"g_effi_map_EG"+EGval[jEG]+".gif");
  c_effi_map[jEG]->Print(dirOut+"g_effi_map_EG"+EGval[jEG]+".C");
    //}

  // OUTPUT LOGS //
  if(debug) cout << "Filling output logs..." << endl;
  ofstream outCounters(dirOut+"counters_selectPairs.csv", ios::out);

  for(int iECAL=0 ; iECAL<nECAL ; iECAL++ ) {
    for(int iColl=0 ; iColl<nColl ; iColl++) {
      outCounters << nEle[iECAL][iColl] << " | " ;
      for(int iEG=0 ; iEG<nEG ; iEG++) {
	outCounters << nMatched[iECAL][iColl][iEG] << "," ;
      }
    }

    outCounters << endl;
  }    

  outCounters << endl << "nProbesTrig15=" << nProbesTrig15 
	      << " | nProbesNoTrig15=" << nProbesNoTrig15 << endl;

  // MASS PLOTS output
  TFile *outplot = new TFile(dirOut+"mass_Z_plots.root","RECREATE");
  for(int i=0 ; i<3 ; i++)
    h_Mee[i]->Write();
  outplot->Close();
  
  //int pause=0;
  //cin >> pause;  

  //outfile->Close();

  out_custom_json << ']' << endl << '}' ;
  out_elepairs_json << ']' << endl << '}' ;

  outlog << endl << endl;
  for(int i=0 ; i<3 ; i++)
    outlog << "-----" ;
  outlog << " RESULTS " ;
  for(int i=0 ; i<4 ; i++)
    outlog << "-----" ;
  outlog << endl
	 << "nCases = " << nCases[0] 
	 << " | nCases SAME = " << nCases[1]
	 << " | nCases DIFF = " << nCases[2]
	 << endl
	 << "nToBeZeroed = " << nZeroStatus[0]
	 << " | nZeroed = " << nZeroStatus[1]
	 << " | nNotZeroed = " << nZeroStatus[2]
	 << endl
	 << "nEB_NotToZero_Zeroed = " << nEB_NotToZero_Zeroed
	 << " | nEE_NotToZero_Zeroed = " << nEE_NotToZero_Zeroed
	 << endl
	 << "n_fails_fidu = " << n_fails_fidu
	 << endl;
  for(int i=0 ; i<9 ; i++)
    outlog << "-----" ;

  
}

void selectPairs(int nEntries=-1, int iEff=0, int iEff_M=0, TString eg_menu="L1_EG_MENU_2011",
		 TString dirOut="",
		 TString dirIn="",
		 TString tagCuts="WP80", double tagEt=20, TString probeCuts="WP95", double probeEt=5, 
		 double massCut1=60, double massCut2=120, TString tag_hlt="noHLT", TString nameChain="ElePairs",
		 bool checkCharge=true, bool data2011=true, bool debug=false)
{
  // Output Log
  ofstream outlog(dirOut+"log_selectPairs.txt",ios::out);

  // Input trees
  TChain * myChain = new TChain(nameChain);
  TString file = "*.root";

  /*
  ostringstream ossi("");
  int nDir=3;
  for(int iDir=0 ; iDir<nDir ; iDir++) {
    ossi.str("");
    ossi << iDir+4 ;
    myChain->Add(dirIn+"PromptRecoV"+ossi.str()+"/Elepairs/"+file);
  }
  */

  myChain->Add(dirIn+file);

  if(debug) cout << myChain->GetEntries() << " entries" << endl;

  // Process the tree
  if(debug) cout << "process the tree" << endl;

  process(myChain, dirOut, iEff, iEff_M, tagCuts, tagEt, probeCuts, probeEt, massCut1, massCut2,
	  nEntries, outlog, tag_hlt, checkCharge, data2011, eg_menu, debug);
  
  cout << "DONE" << endl;

}
