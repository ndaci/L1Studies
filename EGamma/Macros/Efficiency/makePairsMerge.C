///////////////////////////////////////////////////////////////////////////////
///////////////////////////////    makePairs.C    /////////////////////////////
////                                                                       ////
////    Input : Nadir-SimpleNtple output trees                             ////
////            OR (>0 VBTF WP95 ele / evt)-skimmed ones                   ////
////                                                                       ////
////    - select events passing JSON file (2010B / 2011A) with ele_N > 1   ////
////    - make candidate pairs where tag = [VBTF WP95 + eT >= 5 GeV]       ////
////      and probe = any other ele in event                               ////
////      --> pre-selection : M > 30 GeV                                   ////
////      --> no double counting of the pairs !                            ////
////                                                                       ////
////    Output : - tree dirIn/elepairs.root "ElePairs" (1entry=1pair)      ////
////             - dirOut/log.txt                                          ////
////                                                                       ////
////    ==> further selection of the pairs in selectPairs.C                ////
////                                                                       ////
/////////////////////////////    !!! ENJOY IT ;-)    !!!///////////////////////
///////////////////////////////////////////////////////////////////////////////

// makePairs.2.2.C : to read outputs from cleaned MinCode Analyzer (ElectronL1Study)
// makePairsMerge  : will also get the trigger emulation informations (TrigOnly)
// makePairsMerge.1.4.C : choose between both
// makePairsMerge.1.5.C : adding ele_RCTL1iso_To --> 3 sets of ranks : data, data To, M (M To)

#include "../Common/baseFuncNad.1.5.h"

// DATA TREE LISTS
// include here your files treelist.h
// #include "[path]/treelist.h"

//////////////////////////////////////
// MATCH L1 CANDIDATES TO ELECTRONS //
//////////////////////////////////////
int matchL1toEle(int trig_L1emIso_N, int* trig_L1emIso_ieta, int* trig_L1emIso_iphi, int* trig_L1emIso_rank,
		 int ele_RCTeta, int ele_RCTphi, int &ele_RCTL1iso,
		 int* ele_RCTetaVect, int* ele_RCTphiVect, int* ele_RCTL1isoVect, const int nReg, ofstream& outcheckL1) {

  int nC = trig_L1emIso_N;
  if( nC>4 || nC<0 ) nC=4;

  outcheckL1 << "ELE : (" << ele_RCTeta << ";" << ele_RCTphi << ")=" << ele_RCTL1iso << endl;
  
  for(int iC=0 ; iC<nC ; iC++) {
    if( trig_L1emIso_ieta[iC]== ele_RCTeta && trig_L1emIso_iphi[iC]== ele_RCTphi )
      ele_RCTL1iso = trig_L1emIso_rank[iC];

    for(int iR=0 ; iR<nReg ; iR++)
      if( trig_L1emIso_ieta[iC]== ele_RCTetaVect[iR] && trig_L1emIso_iphi[iC]== ele_RCTphiVect[iR] )
	ele_RCTL1isoVect[iR] = trig_L1emIso_rank[iC];
  }

  return 1;
}
////////////////////////////////////////////////////////////////////

//////////////////////////////
// MAP THE TrigOnly ENTRIES //
//////////////////////////////
MAPTrigOnly MapTrigOnlyEntries(TChain *trigOnlyChain, int nEntriesTo, ofstream& outcheck, 
			       MAPJson *jsonMap, bool dojson, bool debug, bool b_mapcheck) {
  
  // Get trigOnlyChain information //
  if(debug) cout << "--- define trigOnlyChain branches" << endl;
  int nEvent_To, nRun_To, nLumi_To;
  //
  trigOnlyChain->SetBranchStatus("trig*",0);
  //
  trigOnlyChain->SetBranchAddress("nRun", &nRun_To);
  trigOnlyChain->SetBranchAddress("nLumi", &nLumi_To);
  trigOnlyChain->SetBranchAddress("nEvent", &nEvent_To);

  // Build the map //
  MAPTrigOnly mapTo; // map< pair<int,int> , int >
  MAPTrigOnly::iterator iterMapTo; //  < <nRun,nEvt> , iEntry >
  pair<int,int> runevtTo;

  // Fill the map //
  if(debug) cout << "--- fill the map" << endl;
  int nEnt = trigOnlyChain->GetEntries();
  if( nEntriesTo>0 && nEntriesTo<nEnt ) nEnt = nEntriesTo;

  for(int iEntry=0 ; iEntry<nEnt ; iEntry++) {
    trigOnlyChain->GetEntry(iEntry);
    if(nEvent_To<0) continue;

    runevtTo = make_pair(nRun_To,nEvent_To);
    iterMapTo = mapTo.find( runevtTo );

    if( iterMapTo == mapTo.end() )
      mapTo[ runevtTo ] = iEntry;

    if(b_mapcheck)
      outcheck << "nRun=" << nRun_To << " | nEvent=" << nEvent_To << " | nEntry=" << iEntry << endl;
  }

  if(debug) cout << "--- return the map" << endl;
  return mapTo;
  
}
//////////////////////////////////////////////////////////////////////////

int process(TChain *trigOnlyChain, MAPJson *jsonMap, MAPTrigOnly mapTo, TString nameChain, TString file, int iFile, int nEntries, 
	    TString dirIn, TString dirOut, bool dojson, TString RunPhase, bool GetTrigOnly, bool debug)  {


  // TrigOnly MAP //
  MAPTrigOnly::iterator iterMapTo;
  pair<int,int> runevtTo;
  int iEntryTo=-1;

  // OUTPUT FILE //
  if(debug) cout << "--- define output file : " ;
  ostringstream ossi("");
  ossi << iFile ;
  //
  TFile *outfile = new TFile(dirIn+"elepairs_"+ossi.str()+".root","RECREATE");
  ofstream outlog(dirOut+"/logs/log_"+ossi.str()+".txt",ios::out);
  ofstream outcheckL1(dirOut+"/checkL1_"+ossi.str()+".txt",ios::out);
  //
  if(debug) cout << "elepairs_"+ossi.str()+".root" << endl;
  ossi.str("");
  
  // INPUT TREE //
  if(debug) cout << "--- add data trees to myChain" << endl;
  TChain * myChain = new TChain(nameChain);
  myChain->Add(file);
  
  // VARIABLES //
  int nEvent, nRun, nLumi ; // data
  int nEvent_To, nRun_To, nLumi_To; // trigOnly

  // Vertices //
  int _vtx_N;
  double _vtx_x[200], _vtx_y[200], _vtx_z[200];
  double _vtx_normalizedChi2[200], _vtx_ndof[200], _vtx_nTracks[200], _vtx_d0[200];

  // Trigger Paths //
  int trig_hltInfo[250];
  //int _trig_isEleHLTpath;
  int trig_HLT_path[4]; // unbias, EG5, EG8, EG12
  char trig_fired_names[5000];
  //
  vector<string> m_HLT_pathsV;
  vector<string> m_HLT_triggered;
  vector<int> m_HLT_pathsV_check;

  // Electrons
  TClonesArray * electrons = new TClonesArray ("TLorentzVector");
  int ele_N; //, sc_hybrid_N; 
  int ele_severityLevelSeed[10];
  double ele_he[10], ele_sigmaietaieta[10];
  double ele_hcalDepth1TowerSumEt_dr03[10], ele_hcalDepth2TowerSumEt_dr03[10];
  double ele_ecalRecHitSumEt_dr03[10], ele_tkSumPt_dr03[10];
  double ele_sclEta[10], ele_sclEt[10];
  //double ecalIsoRel03,hcalIsoRel03,trackIsoRel03;
  double ele_deltaphiin[10], ele_deltaetain[10];
  double ele_conv_dist[10], ele_conv_dcot[10];
  double ele_fbrem[10];
  int ele_expected_inner_hits[10];
  
  //int ele_ambiguousGsfTracks[10];
  int ele_isConversion[10];
  int ele_echarge[10];
  //
  int ele_RCTeta[10], ele_RCTphi[10], ele_RCTL1iso[10], ele_RCTL1noniso[10], ele_RCTL1iso_M[10], ele_RCTL1noniso_M[10], ele_RCTL1iso_To[10], ele_RCTL1noniso_To[10];
  int ele_TTetaVect[10][50], ele_TTphiVect[10][50];
  double ele_TTetVect[10][50];
  int ele_RCTetaVect[10][10], ele_RCTphiVect[10][10], ele_RCTL1isoVect[10][10], 
    ele_RCTL1nonisoVect[10][10],ele_RCTL1isoVect_M[10][10], ele_RCTL1nonisoVect_M[10][10],ele_RCTL1isoVect_To[10][10], ele_RCTL1nonisoVect_To[10][10];
  double ele_RCTetVect[10][10];

  // TP info
  const int nTow = 4032;
  int trig_tower_N,trig_tower_ieta[nTow],trig_tower_iphi[nTow],trig_tower_adc[nTow],trig_tower_sFGVB[nTow],trig_tower_FG[nTow]; 
  int trig_tower_N_M,trig_tower_ieta_M[nTow],trig_tower_iphi_M[nTow],trig_tower_adc_M[nTow],trig_tower_sFGVB_M[nTow],trig_tower_FG_M[nTow]; 
  int trig_tower_N_E,trig_tower_ieta_E[nTow],trig_tower_iphi_E[nTow],trig_tower_adc_E[nTow][5],trig_tower_sFGVB_E[nTow][5],trig_tower_FG_E[nTow][5];

  // HCAL TP
  int trig_tower_hcal_N, trig_tower_hcal_ieta[4032], trig_tower_hcal_iphi[4032], trig_tower_hcal_FG[4032],trig_tower_hcal_et[4032];

  // L1 Candidates //
  int trig_L1emIso_N, trig_L1emNonIso_N; //trig_L1emIso_N_M, trig_L1emNonIso_N_M;
  int trig_L1emIso_N_To, trig_L1emNonIso_N_To, trig_L1emIso_N_M_To, trig_L1emNonIso_N_M_To;

  // L1 candidates info
  int trig_L1emIso_ieta[4], trig_L1emIso_iphi[4], trig_L1emIso_rank[4];
  int trig_L1emNonIso_ieta[4], trig_L1emNonIso_iphi[4], trig_L1emNonIso_rank[4];
  // int trig_L1emIso_ieta_M[4], trig_L1emIso_iphi_M[4], trig_L1emIso_rank_M[4];
  // int trig_L1emNonIso_ieta_M[4], trig_L1emNonIso_iphi_M[4], trig_L1emNonIso_rank_M[4];

  int trig_L1emIso_ieta_To[4], trig_L1emIso_iphi_To[4], trig_L1emIso_rank_To[4];
  int trig_L1emNonIso_ieta_To[4], trig_L1emNonIso_iphi_To[4], trig_L1emNonIso_rank_To[4];
  int trig_L1emIso_ieta_M_To[4], trig_L1emIso_iphi_M_To[4], trig_L1emIso_rank_M_To[4];
  int trig_L1emNonIso_ieta_M_To[4], trig_L1emNonIso_iphi_M_To[4], trig_L1emNonIso_rank_M_To[4];

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

  int trig_strip_mask_N;
  int trig_strip_mask_TTieta[1000], trig_strip_mask_TTiphi[1000], trig_strip_mask_status[1000],
    trig_strip_mask_StripID[1000], trig_strip_mask_PseudoStripID[1000], trig_strip_mask_TccID[1000], trig_strip_mask_CCU[1000],
    trig_strip_mask_xtal_ix[1000][5], trig_strip_mask_xtal_iy[1000][5], trig_strip_mask_xtal_iz[1000][5];

  int trig_xtal_mask_N; // [EB+EE]
  int trig_xtal_mask_ieta[1000],trig_xtal_mask_iphi[1000], // for EE : xtal ieta->ix ; iphi -> iy
    trig_xtal_mask_TTieta[1000],trig_xtal_mask_TTiphi[1000], // but for EE towers, still ieta, iphi...
    trig_xtal_mask_Rieta[1000],trig_xtal_mask_Riphi[1000],
    trig_xtal_mask_status[1000], trig_xtal_mask_EBEE[1000]; // EBEE = {0,1} => 0=EB ; 1=EE
  //double trig_xtal_mask_eT[1000];  


  // INITIALIZATION //
  //
  // Global
  nEvent = 0;
  nRun = 0;
  nLumi = 0;
  //
  // Vertices
  _vtx_N = 0; 
  for(int iv=0;iv<200;iv++) {
    _vtx_normalizedChi2[iv] = 0.;
    _vtx_ndof[iv] = 0.;
    _vtx_nTracks[iv] = 0.;
    _vtx_d0[iv] = 0.;
    _vtx_x[iv] = 0.;
    _vtx_y[iv] = 0.;
    _vtx_z[iv] = 0.;
  }
  //
  // L1 candidates
  trig_L1emIso_N    = 0; 
  trig_L1emNonIso_N = 0;
  trig_L1emIso_N_M_To  = 0; 
  trig_L1emNonIso_N_M_To = 0;
  trig_preL1emIso_N   = 0; 
  trig_preL1emNonIso_N  = 0;
  trig_postL1emIso_N    = 0; 
  trig_postL1emNonIso_N = 0;
  //
  for(int il1=0 ; il1<4 ; il1++) {
    trig_L1emIso_ieta[il1] = 0; 
    trig_L1emIso_iphi[il1] = 0; 
    trig_L1emIso_rank[il1] = 0; 
    trig_L1emNonIso_ieta[il1] = 0; 
    trig_L1emNonIso_iphi[il1] = 0; 
    trig_L1emNonIso_rank[il1] = 0; 

    trig_L1emIso_ieta_To[il1] = 0; 
    trig_L1emIso_iphi_To[il1] = 0; 
    trig_L1emIso_rank_To[il1] = 0; 
    trig_L1emNonIso_ieta_To[il1] = 0; 
    trig_L1emNonIso_iphi_To[il1] = 0; 
    trig_L1emNonIso_rank_To[il1] = 0; 

    trig_L1emIso_ieta_M_To[il1] = 0; 
    trig_L1emIso_iphi_M_To[il1] = 0; 
    trig_L1emIso_rank_M_To[il1] = 0; 
    trig_L1emNonIso_ieta_M_To[il1] = 0; 
    trig_L1emNonIso_iphi_M_To[il1] = 0; 
    trig_L1emNonIso_rank_M_To[il1] = 0; 
		
    trig_preL1emIso_ieta[il1] = 0; 
    trig_preL1emIso_iphi[il1] = 0; 
    trig_preL1emIso_rank[il1] = 0;
    trig_preL1emNonIso_ieta[il1] = 0; 
    trig_preL1emNonIso_iphi[il1] = 0; 
    trig_preL1emNonIso_rank[il1] = 0; 
		
    trig_postL1emIso_ieta[il1] = 0; 
    trig_postL1emIso_iphi[il1] = 0; 
    trig_postL1emIso_rank[il1] = 0;
    trig_postL1emNonIso_ieta[il1] = 0; 
    trig_postL1emNonIso_iphi[il1] = 0; 
    trig_postL1emNonIso_rank[il1] = 0;  
  }
  // 
  // Trigger towers
  trig_tower_N = 0;
  for(int iTow=0 ; iTow<nTow ; iTow++) {
    trig_tower_ieta[iTow] = trig_tower_iphi[iTow]  = -999;
    trig_tower_adc[iTow]  = trig_tower_sFGVB[iTow] = trig_tower_FG[iTow] = -999;
  }
  trig_tower_N_M = 0;
  for(int iTow=0 ; iTow<nTow ; iTow++) {
    trig_tower_ieta_M[iTow] = trig_tower_iphi_M[iTow]  = -999;
    trig_tower_adc_M[iTow]  = trig_tower_sFGVB_M[iTow] = trig_tower_FG_M[iTow] = -999;
  }
  trig_tower_N_E = 0;
  for(int iTow=0 ; iTow<nTow ; iTow++) {
    trig_tower_ieta_E[iTow] = trig_tower_iphi_E[iTow] = -999;
    for(int i=0 ; i<5 ; i++)
      trig_tower_adc_E[iTow][i] = trig_tower_sFGVB_E[iTow][i] = trig_tower_FG_E[iTow][i] = -999;
  }
  trig_tower_hcal_N = 0;
  for(int iTow=0 ; iTow<nTow ; iTow++) {
    trig_tower_hcal_ieta[iTow] = trig_tower_hcal_iphi[iTow]  = -999;
    trig_tower_hcal_FG[iTow]  = trig_tower_hcal_et[iTow] = -999;
  }
  //
  // Masked Towers
  trig_nMaskedRCT=0;
  trig_nMaskedCh=0;
  //
  for (int ii=0;ii<100;ii++) {
    trig_iMaskedRCTeta[ii]   = -999;
    trig_iMaskedRCTphi[ii]   = -999;
    trig_iMaskedRCTcrate[ii] = -999;
    trig_iMaskedTTeta[ii]    = -999;
    trig_iMaskedTTphi[ii]    = -999;
  }
  //
  // Masked strip/xtals
  trig_strip_mask_N = 0;
  trig_xtal_mask_N = 0;
  //
  for(int i=0 ; i<1000 ; i++) {
    trig_strip_mask_TTieta[i] = -999;
    trig_strip_mask_TTiphi[i] = -999;
    trig_strip_mask_status[i] = -999;
    trig_strip_mask_StripID[i] = -999;
    trig_strip_mask_PseudoStripID[i] = -999;
    trig_strip_mask_TccID[i] = -999;
    trig_strip_mask_CCU[i] = -999;
    //
    for(int j=0 ; j<5 ; j++) {
      trig_strip_mask_xtal_ix[i][j] = -999;
      trig_strip_mask_xtal_iy[i][j] = -999;
      trig_strip_mask_xtal_iz[i][j] = -999;
    }
    trig_xtal_mask_ieta[i] = -999;
    trig_xtal_mask_iphi[i] = -999;
    trig_xtal_mask_TTieta[i] = -999;
    trig_xtal_mask_TTiphi[i] = -999;
    trig_xtal_mask_Rieta[i] = -999;
    trig_xtal_mask_Riphi[i] = -999;
    trig_xtal_mask_status[i] = -999;
    trig_xtal_mask_EBEE[i] = -999;
  }
  //

  // Disable useless branches
  if(debug) cout << "--- disable useless branches (sc*)" << endl;
  //myChain->SetBranchStatus("spike_*",0);
  //myChain->SetBranchStatus("vtx_*",0);
  //myChain->SetBranchStatus("skim_*",0);
  //myChain->SetBranchStatus("trig_pre*",0);
  //myChain->SetBranchStatus("trig_post*",0);
  //myChain->SetBranchStatus("trig_HLT*",0);
  //myChain->SetBranchStatus("BS*",0);
  //myChain->SetBranchStatus("MC_*",0);
  //myChain->SetBranchStatus("ele_MC*",0);
  ////myChain->SetBranchStatus("ele_eid*",0);
  ////myChain->SetBranchStatus("ele_Seed*",0);
  //myChain->SetBranchStatus("ele_charge*",0);
  //myChain->SetBranchStatus("met_*",0);
  //myChain->SetBranchStatus("muons*",0);
  //myChain->SetBranchStatus("jets*",0);
  myChain->SetBranchStatus("sc*",0);
  //myChain->SetBranchStatus("sc_hybrid_N",1);
  //myChain->SetBranchStatus("",0);

  if(debug) cout << "--- define myChain branches" << endl;
  // Global
  myChain->SetBranchAddress("nEvent",&nEvent);
  myChain->SetBranchAddress("nRun",&nRun);
  myChain->SetBranchAddress("nLumi",&nLumi);

  // Trigger
//   myChain->SetBranchAddress ("trig_HLT_triggered", &m_HLT_triggered);
//   myChain->SetBranchAddress ("trig_HLT_pathsV", &m_HLT_pathsV);
//   myChain->SetBranchAddress ("trig_HLT_pathsV_check", &m_HLT_pathsV_check);
  //
  myChain->SetBranchAddress("trig_HLT_path",&trig_HLT_path);
  // unbias, EG5, EG8, EG12
  //
  myChain->SetBranchAddress("trig_fired_names",&trig_fired_names);
  myChain->SetBranchAddress("trig_hltInfo",&trig_hltInfo);

  // SC
  //myChain->SetBranchAddress("sc_hybrid_N",   &sc_hybrid_N);
  
  // Electrons
  myChain->SetBranchAddress("ele_N",   &ele_N);
  myChain->SetBranchAddress("electrons",&electrons);
  myChain->SetBranchAddress("ele_severityLevelSeed", &ele_severityLevelSeed);
  myChain->SetBranchAddress("ele_he",&ele_he);
  myChain->SetBranchAddress("ele_sigmaietaieta",&ele_sigmaietaieta);
  myChain->SetBranchAddress("ele_hcalDepth1TowerSumEt_dr03", &ele_hcalDepth1TowerSumEt_dr03);
  myChain->SetBranchAddress("ele_hcalDepth2TowerSumEt_dr03", &ele_hcalDepth2TowerSumEt_dr03);
  myChain->SetBranchAddress("ele_ecalRecHitSumEt_dr03", &ele_ecalRecHitSumEt_dr03);
  myChain->SetBranchAddress("ele_tkSumPt_dr03",&ele_tkSumPt_dr03);
  myChain->SetBranchAddress("ele_sclEta",&ele_sclEta);
  myChain->SetBranchAddress("ele_sclEt",&ele_sclEt);
  myChain->SetBranchAddress("ele_expected_inner_hits_aod",&ele_expected_inner_hits);
  myChain->SetBranchAddress("ele_deltaphiin",&ele_deltaphiin);
  myChain->SetBranchAddress("ele_deltaetain",&ele_deltaetain);
  myChain->SetBranchAddress("ele_conv_dist",&ele_conv_dist);
  myChain->SetBranchAddress("ele_conv_dcot",&ele_conv_dcot);
  myChain->SetBranchAddress("ele_fbrem",&ele_fbrem);
  //myChain->SetBranchAddress("ele_ambiguousGsfTracks", &ele_ambiguousGsfTracks);
  myChain->SetBranchAddress("ele_isConversion",&ele_isConversion);
  myChain->SetBranchAddress("ele_echarge",&ele_echarge);

  // L1 electron informations
  myChain->SetBranchAddress("ele_TTetaVect", &ele_TTetaVect);
  myChain->SetBranchAddress("ele_TTphiVect", &ele_TTphiVect);
  myChain->SetBranchAddress("ele_TTetVect", &ele_TTetVect);
  //
  myChain->SetBranchAddress("ele_RCTeta", &ele_RCTeta);
  myChain->SetBranchAddress("ele_RCTphi", &ele_RCTphi);
  myChain->SetBranchAddress("ele_RCTL1iso", &ele_RCTL1iso);
  myChain->SetBranchAddress("ele_RCTL1noniso", &ele_RCTL1noniso);
  myChain->SetBranchAddress("ele_RCTL1iso_M", &ele_RCTL1iso_M);
  myChain->SetBranchAddress("ele_RCTL1noniso_M", &ele_RCTL1noniso_M);

  myChain->SetBranchAddress("ele_RCTetaVect", &ele_RCTetaVect);
  myChain->SetBranchAddress("ele_RCTphiVect", &ele_RCTphiVect);
  myChain->SetBranchAddress("ele_RCTetVect", &ele_RCTetVect);
  myChain->SetBranchAddress("ele_RCTL1isoVect", &ele_RCTL1isoVect);
  myChain->SetBranchAddress("ele_RCTL1nonisoVect", &ele_RCTL1nonisoVect);
  myChain->SetBranchAddress("ele_RCTL1isoVect_M", &ele_RCTL1isoVect_M);
  myChain->SetBranchAddress("ele_RCTL1nonisoVect_M", &ele_RCTL1nonisoVect_M);

  // L1 candidates
  myChain->SetBranchAddress("trig_L1emIso_N", &trig_L1emIso_N);
  myChain->SetBranchAddress("trig_L1emIso_ieta", &trig_L1emIso_ieta);
  myChain->SetBranchAddress("trig_L1emIso_iphi", &trig_L1emIso_iphi);
  myChain->SetBranchAddress("trig_L1emIso_rank", &trig_L1emIso_rank);

  myChain->SetBranchAddress("trig_L1emNonIso_N", &trig_L1emNonIso_N);
  myChain->SetBranchAddress("trig_L1emNonIso_ieta", &trig_L1emNonIso_ieta);
  myChain->SetBranchAddress("trig_L1emNonIso_iphi", &trig_L1emNonIso_iphi); 
  myChain->SetBranchAddress("trig_L1emNonIso_rank", &trig_L1emNonIso_rank);

  // Pre/post - firing L1 candidates
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

  // Strip masking
  myChain->SetBranchAddress("trig_strip_mask_N", &trig_strip_mask_N);
  //myChain->SetBranchAddress("trig_strip_mask_status", &trig_strip_mask_status);
  myChain->SetBranchAddress("trig_strip_mask_TTieta", &trig_strip_mask_TTieta);
  myChain->SetBranchAddress("trig_strip_mask_TTiphi", &trig_strip_mask_TTiphi);
  myChain->SetBranchAddress("trig_strip_mask_StripID", &trig_strip_mask_StripID);
  myChain->SetBranchAddress("trig_strip_mask_PseudoStripID", &trig_strip_mask_PseudoStripID);
  myChain->SetBranchAddress("trig_strip_mask_TccID", &trig_strip_mask_TccID);
  myChain->SetBranchAddress("trig_strip_mask_CCU", &trig_strip_mask_CCU);
  myChain->SetBranchAddress("trig_strip_mask_xtal_ix", &trig_strip_mask_xtal_ix);
  myChain->SetBranchAddress("trig_strip_mask_xtal_iy", &trig_strip_mask_xtal_iy);
  myChain->SetBranchAddress("trig_strip_mask_xtal_iz", &trig_strip_mask_xtal_iz);
  //
  // Crystal masking
  myChain->SetBranchAddress("trig_xtal_mask_N", &trig_xtal_mask_N);
  myChain->SetBranchAddress("trig_xtal_mask_ieta", &trig_xtal_mask_ieta);
  myChain->SetBranchAddress("trig_xtal_mask_iphi", &trig_xtal_mask_iphi);
  myChain->SetBranchAddress("trig_xtal_mask_TTieta", &trig_xtal_mask_TTieta);
  myChain->SetBranchAddress("trig_xtal_mask_TTiphi", &trig_xtal_mask_TTiphi);
  myChain->SetBranchAddress("trig_xtal_mask_Rieta", &trig_xtal_mask_Rieta);
  myChain->SetBranchAddress("trig_xtal_mask_Riphi", &trig_xtal_mask_Riphi);
  myChain->SetBranchAddress("trig_xtal_mask_status", &trig_xtal_mask_status);
  myChain->SetBranchAddress("trig_xtal_mask_EBEE", &trig_xtal_mask_EBEE);

  // Masking
  myChain->SetBranchAddress("trig_nMaskedRCT",      &trig_nMaskedRCT);      
  myChain->SetBranchAddress("trig_iMaskedRCTeta",   &trig_iMaskedRCTeta);                                          
  myChain->SetBranchAddress("trig_iMaskedRCTcrate", &trig_iMaskedRCTcrate);
  myChain->SetBranchAddress("trig_iMaskedRCTphi",   &trig_iMaskedRCTphi);
  myChain->SetBranchAddress("trig_nMaskedCh",       &trig_nMaskedCh);    
  myChain->SetBranchAddress("trig_iMaskedTTeta",    &trig_iMaskedTTeta);   
  myChain->SetBranchAddress("trig_iMaskedTTphi",    &trig_iMaskedTTphi);      	


  if(debug) cout << "--- define trigOnlyChain branches" << endl;

  // INPUT TrigOnly (emulated trigger information) DATA //
  trigOnlyChain->SetBranchAddress("nRun", &nRun_To);
  trigOnlyChain->SetBranchAddress("nLumi", &nLumi_To);
  trigOnlyChain->SetBranchAddress("nEvent", &nEvent_To);

  // Trigger Towers
  // normal collection
  trigOnlyChain->SetBranchAddress("trig_tower_N", &trig_tower_N);
  trigOnlyChain->SetBranchAddress("trig_tower_ieta",  &trig_tower_ieta);
  trigOnlyChain->SetBranchAddress("trig_tower_iphi",  &trig_tower_iphi);
  trigOnlyChain->SetBranchAddress("trig_tower_adc",  &trig_tower_adc);
  trigOnlyChain->SetBranchAddress("trig_tower_sFGVB",  &trig_tower_sFGVB);
  trigOnlyChain->SetBranchAddress("trig_tower_FG",  &trig_tower_FG);
 
  // modified collections
  trigOnlyChain->SetBranchAddress("trig_tower_N_M", &trig_tower_N_M);
  trigOnlyChain->SetBranchAddress("trig_tower_ieta_M",  &trig_tower_ieta_M);
  trigOnlyChain->SetBranchAddress("trig_tower_iphi_M",  &trig_tower_iphi_M);
  trigOnlyChain->SetBranchAddress("trig_tower_adc_M",  &trig_tower_adc_M);
  trigOnlyChain->SetBranchAddress("trig_tower_sFGVB_M",  &trig_tower_sFGVB_M);
  trigOnlyChain->SetBranchAddress("trig_tower_FG_M",  &trig_tower_FG_M);

  trigOnlyChain->SetBranchAddress("trig_tower_N_E",     &trig_tower_N_E);
  trigOnlyChain->SetBranchAddress("trig_tower_ieta_E",  &trig_tower_ieta_E);
  trigOnlyChain->SetBranchAddress("trig_tower_iphi_E",  &trig_tower_iphi_E);
  trigOnlyChain->SetBranchAddress("trig_tower_adc_E",   &trig_tower_adc_E);
  trigOnlyChain->SetBranchAddress("trig_tower_sFGVB_E",  &trig_tower_sFGVB_E);
  trigOnlyChain->SetBranchAddress("trig_tower_FG_E", &trig_tower_FG_E);
  
  // HCAL TP
  trigOnlyChain->SetBranchAddress("trig_tower_hcal_N", &trig_tower_hcal_N);
  trigOnlyChain->SetBranchAddress("trig_tower_hcal_ieta",  &trig_tower_hcal_ieta);
  trigOnlyChain->SetBranchAddress("trig_tower_hcal_iphi",  &trig_tower_hcal_iphi);
  trigOnlyChain->SetBranchAddress("trig_tower_hcal_et",  &trig_tower_hcal_et);
  trigOnlyChain->SetBranchAddress("trig_tower_hcal_FG",  &trig_tower_hcal_FG);

  // L1 candidates collections
  trigOnlyChain->SetBranchAddress("trig_L1emIso_N", &trig_L1emIso_N_To);
  trigOnlyChain->SetBranchAddress("trig_L1emIso_ieta", &trig_L1emIso_ieta_To);
  trigOnlyChain->SetBranchAddress("trig_L1emIso_iphi", &trig_L1emIso_iphi_To);
  trigOnlyChain->SetBranchAddress("trig_L1emIso_rank", &trig_L1emIso_rank_To);

  trigOnlyChain->SetBranchAddress("trig_L1emNonIso_N", &trig_L1emNonIso_N_To);
  trigOnlyChain->SetBranchAddress("trig_L1emNonIso_ieta", &trig_L1emNonIso_ieta_To);
  trigOnlyChain->SetBranchAddress("trig_L1emNonIso_iphi", &trig_L1emNonIso_iphi_To);
  trigOnlyChain->SetBranchAddress("trig_L1emNonIso_rank", &trig_L1emNonIso_rank_To);

  trigOnlyChain->SetBranchAddress("trig_L1emIso_N_M", &trig_L1emIso_N_M_To);
  trigOnlyChain->SetBranchAddress("trig_L1emIso_ieta_M", &trig_L1emIso_ieta_M_To);
  trigOnlyChain->SetBranchAddress("trig_L1emIso_iphi_M", &trig_L1emIso_iphi_M_To);
  trigOnlyChain->SetBranchAddress("trig_L1emIso_rank_M", &trig_L1emIso_rank_M_To);
  
  trigOnlyChain->SetBranchAddress("trig_L1emNonIso_N_M", &trig_L1emNonIso_N_M_To);
  trigOnlyChain->SetBranchAddress("trig_L1emNonIso_ieta_M", &trig_L1emNonIso_ieta_M_To);
  trigOnlyChain->SetBranchAddress("trig_L1emNonIso_iphi_M", &trig_L1emNonIso_iphi_M_To);
  trigOnlyChain->SetBranchAddress("trig_L1emNonIso_rank_M", &trig_L1emNonIso_rank_M_To);
  
  // OUTPUT TREE //
  if(debug) cout << "--- define output tree ElePairs" << endl;

  TTree * outtree = new TTree("ElePairs","ElePairs");

  // General informations
  outtree->Branch("nRun",&nRun,"nRun/I");
  outtree->Branch("nLumi",&nLumi,"nLumi/I");
  outtree->Branch("nEvent",&nEvent,"nEvent/I");

  // Vertices
  outtree->Branch("vtx_N",&_vtx_N,"vtx_N/I");
  outtree->Branch("vtx_normalizedChi2",&_vtx_normalizedChi2,"vtx_normalizedChi2[200]/D");
  outtree->Branch("vtx_ndof",&_vtx_ndof,"vtx_ndof[200]/D");
  outtree->Branch("vtx_nTracks",&_vtx_nTracks,"vtx_nTracks[200]/D");
  outtree->Branch("vtx_d0",&_vtx_d0,"vtx_d0[200]/D");
  outtree->Branch("vtx_x",&_vtx_x,"vtx_x[200]/D");
  outtree->Branch("vtx_y",&_vtx_y,"vtx_y[200]/D");
  outtree->Branch("vtx_z",&_vtx_z,"vtx_z[200]/D");

  // HLT informations
  outtree->Branch ("trig_HLT_triggered", &m_HLT_triggered, 256000,0);
  outtree->Branch ("trig_HLT_pathsV", &m_HLT_pathsV, 256000,0);
  outtree->Branch ("trig_HLT_pathsV_check", &m_HLT_pathsV_check, 256000,0);
  //
  outtree->Branch("trig_HLT_path",&trig_HLT_path,"trig_HLT_path[4]/I");
  // unbias, EG5, EG8, EG12
  //
  outtree->Branch("trig_fired_names",&trig_fired_names,"trig_fired_names[5000]/C");
  outtree->Branch("trig_hltInfo",&trig_hltInfo,"trig_hltInfo[250]/I");  

  // Trigger towers
  outtree->Branch("trig_tower_N",&trig_tower_N,"trig_tower_N/I");
  outtree->Branch("trig_tower_ieta",&trig_tower_ieta,"trig_tower_ieta[4032]/I");
  outtree->Branch("trig_tower_iphi",&trig_tower_iphi,"trig_tower_iphi[4032]/I");
  outtree->Branch("trig_tower_adc",&trig_tower_adc,"trig_tower_adc[4032]/I");
  outtree->Branch("trig_tower_sFGVB",&trig_tower_sFGVB,"trig_tower_sFGVB[4032]/I");
  outtree->Branch("trig_tower_FG",&trig_tower_FG,"trig_tower_FG[4032]/I");
  //
  outtree->Branch("trig_tower_N_M",&trig_tower_N_M,"trig_tower_N_M/I");
  outtree->Branch("trig_tower_ieta_M",&trig_tower_ieta_M,"trig_tower_ieta_M[4032]/I");
  outtree->Branch("trig_tower_iphi_M",&trig_tower_iphi_M,"trig_tower_iphi_M[4032]/I");
  outtree->Branch("trig_tower_adc_M",&trig_tower_adc_M,"trig_tower_adc_M[4032]/I");
  outtree->Branch("trig_tower_sFGVB_M",&trig_tower_sFGVB_M,"trig_tower_sFGVB_M[4032]/I");
  outtree->Branch("trig_tower_FG_M",&trig_tower_FG_M,"trig_tower_FG_M[4032]/I");
  //
  outtree->Branch("trig_tower_N_E",&trig_tower_N_E,"trig_tower_N_E/I");
  outtree->Branch("trig_tower_ieta_E",&trig_tower_ieta_E,"trig_tower_ieta_E[4032]/I");
  outtree->Branch("trig_tower_iphi_E",&trig_tower_iphi_E,"trig_tower_iphi_E[4032]/I");
  outtree->Branch("trig_tower_adc_E",&trig_tower_adc_E,"trig_tower_adc_E[4032][5]/I");
  outtree->Branch("trig_tower_sFGVB_E",&trig_tower_sFGVB_E,"trig_tower_sFGVB_E[4032][5]/I");
  outtree->Branch("trig_tower_FG_E",&trig_tower_FG_E,"trig_tower_FG_E[4032][5]/I");

  // HCAL TP
  outtree->Branch("trig_tower_hcal_N", &trig_tower_hcal_N, "trig_tower_hcal_N/I");
  outtree->Branch("trig_tower_hcal_ieta",  &trig_tower_hcal_ieta,  "trig_tower_hcal_ieta[trig_tower_N]/I");
  outtree->Branch("trig_tower_hcal_iphi",  &trig_tower_hcal_iphi,  "trig_tower_hcal_iphi[trig_tower_N]/I");
  outtree->Branch("trig_tower_hcal_et",  &trig_tower_hcal_et,  "trig_tower_hcal_et[trig_tower_N]/I");
  outtree->Branch("trig_tower_hcal_FG",  &trig_tower_hcal_FG,  "trig_tower_hcal_FG[trig_tower_N]/I");

  // Strip masking
  outtree->Branch("trig_strip_mask_N", &trig_strip_mask_N, "trig_strip_mask_N/I");
  //outtree->Branch("trig_strip_mask_status", &trig_strip_mask_status, "trig_strip_mask_status[trig_strip_mask_N]/I");
  outtree->Branch("trig_strip_mask_TTieta", &trig_strip_mask_TTieta, "trig_strip_mask_TTieta[trig_strip_mask_N]/I");
  outtree->Branch("trig_strip_mask_TTiphi", &trig_strip_mask_TTiphi, "trig_strip_mask_TTiphi[trig_strip_mask_N]/I");
  outtree->Branch("trig_strip_mask_StripID", &trig_strip_mask_StripID, "trig_strip_mask_StripID[trig_strip_mask_N]/I");
  outtree->Branch("trig_strip_mask_PseudoStripID", &trig_strip_mask_PseudoStripID, "trig_strip_mask_PseudoStripID[trig_strip_mask_N]/I");
  outtree->Branch("trig_strip_mask_TccID", &trig_strip_mask_TccID, "trig_strip_mask_TccID[trig_strip_mask_N]/I");
  outtree->Branch("trig_strip_mask_CCU", &trig_strip_mask_CCU, "trig_strip_mask_CCU[trig_strip_mask_N]/I");
  outtree->Branch("trig_strip_mask_xtal_ix", &trig_strip_mask_xtal_ix, "trig_strip_mask_xtal_ix[trig_strip_mask_N][5]/I");
  outtree->Branch("trig_strip_mask_xtal_iy", &trig_strip_mask_xtal_iy, "trig_strip_mask_xtal_iy[trig_strip_mask_N][5]/I");
  outtree->Branch("trig_strip_mask_xtal_iz", &trig_strip_mask_xtal_iz, "trig_strip_mask_xtal_iz[trig_strip_mask_N][5]/I");
  //
  // Crystal masking
  outtree->Branch("trig_xtal_mask_N", &trig_xtal_mask_N, "trig_xtal_mask_N/I");
  outtree->Branch("trig_xtal_mask_ieta", &trig_xtal_mask_ieta, "trig_xtal_mask_ieta[trig_xtal_mask_N]/I");
  outtree->Branch("trig_xtal_mask_iphi", &trig_xtal_mask_iphi, "trig_xtal_mask_iphi[trig_xtal_mask_N]/I");
  outtree->Branch("trig_xtal_mask_TTieta", &trig_xtal_mask_TTieta, "trig_xtal_mask_TTieta[trig_xtal_mask_N]/I");
  outtree->Branch("trig_xtal_mask_TTiphi", &trig_xtal_mask_TTiphi, "trig_xtal_mask_TTiphi[trig_xtal_mask_N]/I");
  outtree->Branch("trig_xtal_mask_Rieta", &trig_xtal_mask_Rieta, "trig_xtal_mask_Rieta[trig_xtal_mask_N]/I");
  outtree->Branch("trig_xtal_mask_Riphi", &trig_xtal_mask_Riphi, "trig_xtal_mask_Riphi[trig_xtal_mask_N]/I");
  outtree->Branch("trig_xtal_mask_status", &trig_xtal_mask_status, "trig_xtal_mask_status[trig_xtal_mask_N]/I");
  outtree->Branch("trig_xtal_mask_EBEE", &trig_xtal_mask_EBEE, "trig_xtal_mask_EBEE[trig_xtal_mask_N]/I");

  // L1 candidates from Data
  outtree->Branch("trig_L1emIso_N",     &trig_L1emIso_N,     "trig_L1emIso_N/I");
  outtree->Branch("trig_L1emIso_ieta",  &trig_L1emIso_ieta,  "trig_L1emIso_ieta[4]/I");
  outtree->Branch("trig_L1emIso_iphi",  &trig_L1emIso_iphi,  "trig_L1emIso_iphi[4]/I");
  outtree->Branch("trig_L1emIso_rank",  &trig_L1emIso_rank,  "trig_L1emIso_rank[4]/I");
  //   outtree->Branch("trig_L1emIso_et",    &trig_L1emIso_et,    "trig_L1emIso_et[4]/D");
  //
  outtree->Branch("trig_L1emNonIso_N",     &trig_L1emNonIso_N,     "trig_L1emNonIso_N/I");
  outtree->Branch("trig_L1emNonIso_ieta",  &trig_L1emNonIso_ieta,  "trig_L1emNonIso_ieta[4]/I");
  outtree->Branch("trig_L1emNonIso_iphi",  &trig_L1emNonIso_iphi,  "trig_L1emNonIso_iphi[4]/I");
  outtree->Branch("trig_L1emNonIso_rank",  &trig_L1emNonIso_rank,  "trig_L1emNonIso_rank[4]/I");
  //   outtree->Branch("trig_L1emNonIso_et",    &trig_L1emNonIso_et,    "trig_L1emNonIso_et[4]/D");

  // L1 candidates from TrigOnly
  outtree->Branch("trig_L1emIso_N_To",     &trig_L1emIso_N_To,     "trig_L1emIso_N_To/I");
  outtree->Branch("trig_L1emIso_ieta_To",  &trig_L1emIso_ieta_To,  "trig_L1emIso_ieta_To[4]/I");
  outtree->Branch("trig_L1emIso_iphi_To",  &trig_L1emIso_iphi_To,  "trig_L1emIso_iphi_To[4]/I");
  outtree->Branch("trig_L1emIso_rank_To",  &trig_L1emIso_rank_To,  "trig_L1emIso_rank_To[4]/I");
  //   outtree->Branch("trig_L1emIso_et_To",    &trig_L1emIso_et_To,    "trig_L1emIso_et_To[4]/D");
  //
  outtree->Branch("trig_L1emNonIso_N_To", &trig_L1emNonIso_N_To, "trig_L1emNonIso_N_To/I");
  outtree->Branch("trig_L1emNonIso_ieta_To", &trig_L1emNonIso_ieta_To, "trig_L1emNonIso_ieta_To[4]/I");
  outtree->Branch("trig_L1emNonIso_iphi_To", &trig_L1emNonIso_iphi_To, "trig_L1emNonIso_iphi_To[4]/I");
  outtree->Branch("trig_L1emNonIso_rank_To", &trig_L1emNonIso_rank_To, "trig_L1emNonIso_rank_To[4]/I");
  //   outtree->Branch("trig_L1emNonIso_et_To", &trig_L1emNonIso_et_To, "trig_L1emNonIso_et_To[4]/D");
  //
  outtree->Branch("trig_L1emIso_N_M_To",     &trig_L1emIso_N_M_To,     "trig_L1emIso_N_M_To/I");
  outtree->Branch("trig_L1emIso_ieta_M_To",  &trig_L1emIso_ieta_M_To,  "trig_L1emIso_ieta_M_To[4]/I");
  outtree->Branch("trig_L1emIso_iphi_M_To",  &trig_L1emIso_iphi_M_To,  "trig_L1emIso_iphi_M_To[4]/I");
  outtree->Branch("trig_L1emIso_rank_M_To",  &trig_L1emIso_rank_M_To,  "trig_L1emIso_rank_M_To[4]/I");
  //   outtree->Branch("trig_L1emIso_et_M_To",    &trig_L1emIso_et_M_To,    "trig_L1emIso_et_M_To[4]/D");
  //
  outtree->Branch("trig_L1emNonIso_N_M_To", &trig_L1emNonIso_N_M_To, "trig_L1emNonIso_N_M_To/I");
  outtree->Branch("trig_L1emNonIso_ieta_M_To", &trig_L1emNonIso_ieta_M_To, "trig_L1emNonIso_ieta_M_To[4]/I");
  outtree->Branch("trig_L1emNonIso_iphi_M_To", &trig_L1emNonIso_iphi_M_To, "trig_L1emNonIso_iphi_M_To[4]/I");
  outtree->Branch("trig_L1emNonIso_rank_M_To", &trig_L1emNonIso_rank_M_To, "trig_L1emNonIso_rank_M_To[4]/I");
  //   outtree->Branch("trig_L1emNonIso_et_M_To", &trig_L1emNonIso_et_M_To, "trig_L1emNonIso_et_M_To[4]/D");

  // pre-post firing L1 candidates
  outtree->Branch("trig_preL1emIso_N",     &trig_preL1emIso_N,     "trig_preL1emIso_N/I");
  outtree->Branch("trig_preL1emIso_ieta",  &trig_preL1emIso_ieta,  "trig_preL1emIso_ieta[4]/I");
  outtree->Branch("trig_preL1emIso_iphi",  &trig_preL1emIso_iphi,  "trig_preL1emIso_iphi[4]/I");
  outtree->Branch("trig_preL1emIso_rank",  &trig_preL1emIso_rank,  "trig_preL1emIso_rank[4]/I");
  //
  outtree->Branch("trig_preL1emNonIso_N",     &trig_preL1emNonIso_N,     "trig_preL1emNonIso_N/I");
  outtree->Branch("trig_preL1emNonIso_ieta",  &trig_preL1emNonIso_ieta,  "trig_preL1emNonIso_ieta[4]/I");
  outtree->Branch("trig_preL1emNonIso_iphi",  &trig_preL1emNonIso_iphi,  "trig_preL1emNonIso_iphi[4]/I");
  outtree->Branch("trig_preL1emNonIso_rank",  &trig_preL1emNonIso_rank,  "trig_preL1emNonIso_rank[4]/I");
  //
  outtree->Branch("trig_postL1emIso_N",     &trig_postL1emIso_N,     "trig_postL1emIso_N/I");
  outtree->Branch("trig_postL1emIso_ieta",  &trig_postL1emIso_ieta,  "trig_postL1emIso_ieta[4]/I");
  outtree->Branch("trig_postL1emIso_iphi",  &trig_postL1emIso_iphi,  "trig_postL1emIso_iphi[4]/I");
  outtree->Branch("trig_postL1emIso_rank",  &trig_postL1emIso_rank,  "trig_postL1emIso_rank[4]/I");
  //
  outtree->Branch("trig_postL1emNonIso_N",     &trig_postL1emNonIso_N,     "trig_postL1emNonIso_N/I");
  outtree->Branch("trig_postL1emNonIso_ieta",  &trig_postL1emNonIso_ieta,  "trig_postL1emNonIso_ieta[4]/I");
  outtree->Branch("trig_postL1emNonIso_iphi",  &trig_postL1emNonIso_iphi,  "trig_postL1emNonIso_iphi[4]/I");
  outtree->Branch("trig_postL1emNonIso_rank",  &trig_postL1emNonIso_rank,  "trig_postL1emNonIso_rank[4]/I");
  //
  outtree->Branch("trig_nMaskedRCT",      &trig_nMaskedRCT,     "trig_nMaskedRCT/I");      
  outtree->Branch("trig_iMaskedRCTeta",   &trig_iMaskedRCTeta,  "trig_iMaskedRCTeta[trig_nMaskedRCT]/I");                                          
  outtree->Branch("trig_iMaskedRCTcrate", &trig_iMaskedRCTcrate,"trig_iMaskedRCTcrate[trig_nMaskedRCT]/I");                                                    
  outtree->Branch("trig_iMaskedRCTphi",   &trig_iMaskedRCTphi,  "trig_iMaskedRCTphi[trig_nMaskedRCT]/I");
  outtree->Branch("trig_nMaskedCh",       &trig_nMaskedCh,      "trig_nMaskedCh/I");    
  outtree->Branch("trig_iMaskedTTeta",    &trig_iMaskedTTeta,   "trig_iMaskedTTeta[trig_nMaskedCh]/I");   
  outtree->Branch("trig_iMaskedTTphi",    &trig_iMaskedTTphi,   "trig_iMaskedTTphi[trig_nMaskedCh]/I");      	


  // Pairs informations
  if(debug) cout << "----> adding pair informations branches" << endl;
  double pair_M;
  double pair_eta[2], pair_sclEta[2], pair_sclEt[2], pair_phi[2], pair_pT[2], pair_eT[2], pair_E[2];
  int pair_cuts[2], pair_HLT_Ele27_cut[2], pair_fidu[2], pair_charge[2], pair_RCTeta[2], pair_RCTphi[2], 
    pair_L1iso[2], pair_L1noniso[2], pair_L1iso_M[2], pair_L1noniso_M[2], pair_L1iso_To[2], pair_L1noniso_To[2];
  int pair_RCTetaVect[2][10], pair_RCTphiVect[2][10], 
    pair_L1isoVect[2][10], pair_L1nonisoVect[2][10],pair_L1isoVect_M[2][10], pair_L1nonisoVect_M[2][10],
    pair_L1isoVect_To[2][10], pair_L1nonisoVect_To[2][10];
  double pair_RCTetVect[2][10];
  int pair_TTetaVect[2][50], pair_TTphiVect[2][50];
  double pair_TTetVect[2][50];

  /*
  for(int i=0 ; i<2 ; i++) {
    pair_L1iso[i]=0;
    pair_L1noniso[i]=0;
    pair_L1iso_M[i]=0;
    pair_L1noniso_M[i]=0;
  }
  */

  //
  outtree->Branch("pair_M",&pair_M,"pair_M/D");
  //
  outtree->Branch("pair_cuts",&pair_cuts,"pair_cuts[2]/I");
  // 0 : noCut | 1 : VBTF 95 | 2 : VBTF 80 | 3 : VBTF 60
  outtree->Branch("pair_HLT_Ele27_cut",&pair_HLT_Ele27_cut,"pair_HLT_Ele27_cut[2]/I");
  outtree->Branch("pair_fidu",&pair_fidu,"pair_fidu[2]/I");
  //
  outtree->Branch("pair_eta",&pair_eta,"pair_eta[2]/D");
  outtree->Branch("pair_sclEta",&pair_sclEta,"pair_sclEta[2]/D");
  outtree->Branch("pair_phi",&pair_phi,"pair_phi[2]/D");
  outtree->Branch("pair_RCTeta",&pair_RCTeta,"pair_RCTeta[2]/I");
  outtree->Branch("pair_RCTphi",&pair_RCTphi,"pair_RCTphi[2]/I");
  //
  outtree->Branch("pair_charge",&pair_charge,"pair_charge[2]/I");
  outtree->Branch("pair_pT",&pair_pT,"pair_pT[2]/D");
  outtree->Branch("pair_eT",&pair_eT,"pair_eT[2]/D");
  outtree->Branch("pair_sclEt",&pair_sclEt,"pair_sclEt[2]/D");
  outtree->Branch("pair_E",&pair_E,"pair_E[2]/D");

  outtree->Branch("pair_TTetaVect", &pair_TTetaVect,"pair_TTetaVect[2][50]/I");
  outtree->Branch("pair_TTphiVect", &pair_TTphiVect,"pair_TTphiVect[2][50]/I");
  outtree->Branch("pair_TTetVect", &pair_TTetVect,"pair_TTetVect[2][50]/D");

  outtree->Branch("pair_L1iso",&pair_L1iso,"pair_L1iso[2]/I");
  outtree->Branch("pair_L1noniso",&pair_L1noniso,"pair_L1noniso[2]/I");
  outtree->Branch("pair_L1iso_M",&pair_L1iso_M,"pair_L1iso_M[2]/I");
  outtree->Branch("pair_L1noniso_M",&pair_L1noniso_M,"pair_L1noniso_M[2]/I");
  outtree->Branch("pair_L1iso_To",&pair_L1iso_To,"pair_L1iso_To[2]/I");
  outtree->Branch("pair_L1noniso_To",&pair_L1noniso_To,"pair_L1noniso_To[2]/I");
  //
  outtree->Branch("pair_RCTetVect",&pair_RCTetVect,"pair_RCTetVect[2][10]/D");
  outtree->Branch("pair_RCTetaVect",&pair_RCTetaVect,"pair_RCTetaVect[2][10]/I");
  outtree->Branch("pair_RCTphiVect",&pair_RCTphiVect,"pair_RCTphiVect[2][10]/I");
  outtree->Branch("pair_L1isoVect",&pair_L1isoVect,"pair_L1isoVect[2][10]/I");
  outtree->Branch("pair_L1nonisoVect",&pair_L1nonisoVect,"pair_L1nonisoVect[2][10]/I");
  outtree->Branch("pair_L1isoVect_M",&pair_L1isoVect_M,"pair_L1isoVect_M[2][10]/I");
  outtree->Branch("pair_L1nonisoVect_M",&pair_L1nonisoVect_M,"pair_L1nonisoVect_M[2][10]/I");
  outtree->Branch("pair_L1isoVect_To",&pair_L1isoVect_To,"pair_L1isoVect_To[2][10]/I");
  outtree->Branch("pair_L1nonisoVect_To",&pair_L1nonisoVect_To,"pair_L1nonisoVect_To[2][10]/I");

  // USEFUL VARIABLES //
  vector<int> pairIdx;
  int cutEle[2], fidu[2];
  bool cut_HLT_Ele27[2];
  TLorentzVector* cand[2];
  TLorentzVector total;
  bool isGoodRun;
  TString filename;


  // -------------------------------------------------------------------------------
  // LOOP OVER EVENTS
  // -------------------------------------------------------------------------------  
  if(debug) cout << "--- loop over myChain entries (data events)" << endl;

  int numEntries = myChain->GetEntries () ;
  int nProcess = numEntries;
  if(nEntries>=0 && nEntries<numEntries)
    nProcess = nEntries;

  int nCurrentRun = -999;
  outlog << "will process " << nProcess << "/" << numEntries << "entries" << endl;

  for (int iEvent = 0 ; iEvent < nProcess ; iEvent++ )
    { 

      // HLT information
      for(int i=0 ; i<4 ; i++)
	trig_HLT_path[i]=0;

      for(int i=0 ; i<250 ; i++)
	trig_hltInfo[i]=0;

      for(int i=0 ; i<5000 ; i++)
	trig_fired_names[i]=0;

      m_HLT_pathsV.clear();
      m_HLT_triggered.clear();
      m_HLT_pathsV_check.clear();

      // TP Initialization
      trig_tower_N = 0;
      for(int iTow=0 ; iTow<nTow ; iTow++) {
	trig_tower_ieta[iTow] = trig_tower_iphi[iTow]  = -999;
	trig_tower_adc[iTow]  = trig_tower_sFGVB[iTow] = -999;
      }
      trig_tower_N_M = 0;
      for(int iTow=0 ; iTow<nTow ; iTow++) {
	trig_tower_ieta_M[iTow] = trig_tower_iphi_M[iTow]  = -999;
	trig_tower_adc_M[iTow]  = trig_tower_sFGVB_M[iTow] = -999;
      }
      trig_tower_N_E = 0;
      for(int iTow=0 ; iTow<nTow ; iTow++) {
	trig_tower_ieta_E[iTow] = trig_tower_iphi_E[iTow] = -999;
	for(int i=0 ; i<5 ; i++)
	  trig_tower_adc_E[iTow][i] = trig_tower_sFGVB_E[iTow][i] = -999;
      }

      // L1 candidates
      trig_L1emIso_N    = 0; 
      trig_L1emNonIso_N = 0;
      trig_L1emIso_N_M_To  = 0; 
      trig_L1emNonIso_N_M_To = 0;
      trig_L1emIso_N_To  = 0; 
      trig_L1emNonIso_N_To = 0;

      for(int il1=0 ; il1<4 ; il1++) {
	trig_L1emIso_ieta[il1] = 0; 
	trig_L1emIso_iphi[il1] = 0; 
	trig_L1emIso_rank[il1] = 0; 
	trig_L1emNonIso_ieta[il1] = 0; 
	trig_L1emNonIso_iphi[il1] = 0; 
	trig_L1emNonIso_rank[il1] = 0; 
 
	trig_L1emIso_ieta_To[il1] = 0; 
	trig_L1emIso_iphi_To[il1] = 0; 
	trig_L1emIso_rank_To[il1] = 0; 
 	trig_L1emNonIso_ieta_To[il1] = 0; 
	trig_L1emNonIso_iphi_To[il1] = 0; 
	trig_L1emNonIso_rank_To[il1] = 0; 
  
	trig_L1emNonIso_ieta_M_To[il1] = 0; 
	trig_L1emNonIso_iphi_M_To[il1] = 0; 
	trig_L1emNonIso_rank_M_To[il1] = 0; 
 	trig_L1emIso_ieta_M_To[il1] = 0; 
	trig_L1emIso_iphi_M_To[il1] = 0; 
	trig_L1emIso_rank_M_To[il1] = 0; 
      }
	
      for(int i=0 ; i<2 ; i++) {
	pair_L1iso[i]=0;
	pair_L1noniso[i]=0;
	pair_L1iso_M[i]=0;
	pair_L1noniso_M[i]=0;

	for(int j=0 ; j<10 ; j++) {
	  pair_L1isoVect[i][j]=0;
	  pair_L1nonisoVect[i][j]=0;
	  pair_L1isoVect_M[i][j]=0;
	  pair_L1nonisoVect_M[i][j]=0;
	}
      }

      //cout << "GET ENTRY" << endl;
      myChain->GetEntry (iEvent) ;
      //cout << "GOT" << endl;
      
      // show processed file
      if(iEvent==0) {
        filename = myChain->GetFile()->GetName() ;
        outlog << "File : " << filename << endl << endl;
      }
      else if( filename != myChain->GetFile()->GetName() ) {
        filename = myChain->GetFile()->GetName() ;
        outlog << "File : " << myChain->GetFile()->GetName() << endl << endl;
      }
     
      // show current run/iCat processed
      if(iEvent==0) {
	nCurrentRun = nRun ;
	outlog << "nRun=" << nRun << endl;
      }
      else if(nRun!=nCurrentRun) {
	nCurrentRun=nRun ;
	outlog << "nRun=" << nRun << endl;
      }

      // run selection (using both json files)
      //int iJson = detJson(nRun);
      int iJson = 4;
      if(debug) cout << "iJson = " << iJson << endl;
      if( iJson>-1 && iJson<5) {
	isGoodRun = AcceptEventByRunAndLumiSection(nRun, nLumi, jsonMap[iJson]);
	if(!isGoodRun && dojson) {
	  outlog << "failed JSON" << endl;
	  continue;
	}
      }
      else {
	outlog << "failed JSON" << endl;
	continue;
      }

      // at least 2 electrons
      if(ele_N<2) continue;
      else outlog << "ele_N=" << ele_N << endl;


      // Load trigOnly information ////////////////////////////////////////////////////////////////////
      if( GetTrigOnly ) {

	runevtTo = make_pair(nRun,nEvent);
	iterMapTo = mapTo.find( runevtTo );
	//
	if( iterMapTo != mapTo.end() )
	  iEntryTo = mapTo[ runevtTo ];
	else iEntryTo = -1;
	//
	outlog << "trigOnly info : nRun=" << nRun << " | nEvent=" << nEvent << " | iEntryTo=" << iEntryTo << endl;
	//
	if(iEntryTo>=0) {
	  trigOnlyChain->GetEntry( iEntryTo );
	  if(debug) cout << "----> got entry from trigOnlyChain" << endl;
	}
	else {
	  if(debug) cout << "----> entry not mapped" << endl;
	  continue;
	}
      }

      // REDO MATCHING ELE / L1-CANDIDATE //////////////////////////////////////////////////////////////////
      outcheckL1 << "ISO_N_data : " ;
      for(int i=0 ; i<4 ; i++)
	outcheckL1 << "(" << trig_L1emIso_ieta[i] << ";" << trig_L1emIso_iphi[i] << ")=" << trig_L1emIso_rank[i] << " | ";

      outcheckL1 << endl
	   << "NONISO_N_data : ";
      for(int i=0 ; i<4 ; i++)
	outcheckL1 << "(" << trig_L1emNonIso_ieta[i] << ";" << trig_L1emNonIso_iphi[i] << ")=" << trig_L1emNonIso_rank[i] << " | ";

      outcheckL1 << endl
	   << "ISO_N : " ;
      for(int i=0 ; i<4 ; i++)
	outcheckL1 << "(" << trig_L1emIso_ieta_To[i] << ";" << trig_L1emIso_iphi_To[i] << ")=" << trig_L1emIso_rank_To[i] << " | ";

      outcheckL1 << endl
	   << "NONISO_N : ";
      for(int i=0 ; i<4 ; i++)
	outcheckL1 << "(" << trig_L1emNonIso_ieta_To[i] << ";" << trig_L1emNonIso_iphi_To[i] << ")=" << trig_L1emNonIso_rank_To[i] << " | ";

      outcheckL1 << endl
	   << "ISO_M : ";
      for(int i=0 ; i<4 ; i++)
	outcheckL1 << "(" << trig_L1emIso_ieta_M_To[i] << ";" << trig_L1emIso_iphi_M_To[i] << ")=" << trig_L1emIso_rank_M_To[i] << " | ";

      outcheckL1 << endl
	   << "NONISO_M : ";
      for(int i=0 ; i<4 ; i++)
	outcheckL1 << "(" << trig_L1emNonIso_ieta_M_To[i] << ";" << trig_L1emNonIso_iphi_M_To[i] << ")=" << trig_L1emNonIso_rank_M_To[i] << " | ";

      outcheckL1 << endl << endl;

      const int nReg=10;
      for( int iEle=0 ; iEle<ele_N ; iEle++ ) {
	
	outcheckL1 << "BEFORE :  ele_RCTL1iso_M[iEle]=" << ele_RCTL1iso_M[iEle] << endl;

	// Data Iso candidates
	matchL1toEle( trig_L1emIso_N_To, trig_L1emIso_ieta_To, trig_L1emIso_iphi_To, trig_L1emIso_rank_To,
		      ele_RCTeta[iEle],  ele_RCTphi[iEle],  ele_RCTL1iso_To[iEle],
		      ele_RCTetaVect[iEle], ele_RCTphiVect[iEle], ele_RCTL1isoVect_To[iEle], nReg, outcheckL1);

	// Data NonIso candidates
	// cout << "trig_L1emNonIso_N_To=" << trig_L1emNonIso_N_To << endl;
	matchL1toEle( trig_L1emNonIso_N_To, trig_L1emNonIso_ieta_To, trig_L1emNonIso_iphi_To, trig_L1emNonIso_rank_To,
		      ele_RCTeta[iEle],  ele_RCTphi[iEle],  ele_RCTL1noniso_To[iEle],
		      ele_RCTetaVect[iEle], ele_RCTphiVect[iEle], ele_RCTL1nonisoVect_To[iEle], nReg, outcheckL1);
	
	// Emul Iso candidates
	matchL1toEle( trig_L1emIso_N_M_To, trig_L1emIso_ieta_M_To, trig_L1emIso_iphi_M_To, trig_L1emIso_rank_M_To,
		      ele_RCTeta[iEle],  ele_RCTphi[iEle],  ele_RCTL1iso_M[iEle],
		      ele_RCTetaVect[iEle], ele_RCTphiVect[iEle], ele_RCTL1isoVect_M[iEle], nReg, outcheckL1);
	
	// Emul NonIso candidates
	matchL1toEle( trig_L1emNonIso_N_M_To, trig_L1emNonIso_ieta_M_To, trig_L1emNonIso_iphi_M_To, trig_L1emNonIso_rank_M_To,
		      ele_RCTeta[iEle],  ele_RCTphi[iEle],  ele_RCTL1noniso_M[iEle],
		      ele_RCTetaVect[iEle], ele_RCTphiVect[iEle], ele_RCTL1nonisoVect_M[iEle], nReg, outcheckL1);
	
	outcheckL1 << "AFTER :  ele_RCTL1iso_M[iEle]=" << ele_RCTL1iso_M[iEle] << endl;
      }
      
      ///////////////////////////////////////////////////////////////////////////////////////////////////////

      // LOOP OVER ELECTRONS //
      if(debug) cout << "<-- ele_N=" << ele_N << endl
		     << "--- electrons.size=" << electrons->GetSize() << endl;
      for( int iEle1=0 ; iEle1<ele_N ; iEle1++ ) {
	if(debug) cout << "--- get ele #" << iEle1 << endl;
	cand[0] = (TLorentzVector*) (electrons->At (iEle1)) ;
	if(debug) cout << "--- got it" << endl;

	// severity selection
	if( ele_severityLevelSeed[iEle1] >= 3 ) continue;
	
	// check whether electrons of the pair pass HLT_Ele27 Id/Iso cuts
	cut_HLT_Ele27[0] = VBTFcuts( "HLT_Ele27", RunPhase,
				     cand[0]->Pt(), cand[0]->Et(), ele_sclEta[iEle1], cand[0]->Eta(), ele_tkSumPt_dr03[iEle1], ele_ecalRecHitSumEt_dr03[iEle1], 
				     ele_hcalDepth1TowerSumEt_dr03[iEle1], ele_hcalDepth2TowerSumEt_dr03[iEle1], ele_expected_inner_hits[iEle1],
				     ele_deltaphiin[iEle1], ele_deltaetain[iEle1], ele_he[iEle1], ele_sigmaietaieta[iEle1],
				     ele_conv_dist[iEle1], ele_conv_dcot[iEle1], ele_fbrem[iEle1], ele_isConversion[iEle1] ) ;

	// check if ele is a good tag candidate : pass VBTF 95 and has pT>5 GeV
	cutEle[0] = 0;
	cutEle[0] = whichCuts( RunPhase, cand[0]->Pt(), cand[0]->Et(), ele_sclEta[iEle1], cand[0]->Eta(), ele_tkSumPt_dr03[iEle1], ele_ecalRecHitSumEt_dr03[iEle1], 
			       ele_hcalDepth1TowerSumEt_dr03[iEle1], ele_hcalDepth2TowerSumEt_dr03[iEle1], ele_expected_inner_hits[iEle1],
			       ele_deltaphiin[iEle1], ele_deltaetain[iEle1], ele_he[iEle1], ele_sigmaietaieta[iEle1],
			       ele_conv_dist[iEle1], ele_conv_dcot[iEle1], ele_fbrem[iEle1], ele_isConversion[iEle1] ) ;
	fidu[0] = 0;
	if ( fabs(ele_sclEta[iEle1]) < 2.5 && ( fabs(ele_sclEta[iEle1]) > 1.566 || fabs(ele_sclEta[iEle1])<1.4442 ) ) 
	  fidu[0] = 1 ;

	if( cutEle[0]>0 && cand[0]->Et()>=5. ) {

	  // loop to find probe candidates
	  for( int iEle2=0 ; iEle2<ele_N ; iEle2++ ) {
	    cand[1] = (TLorentzVector*) (electrons->At (iEle2)) ;

	    // severity
	    if( ele_severityLevelSeed[iEle2] >= 3 ) continue;

	    // check HLT_Ele27 cuts
	    cut_HLT_Ele27[1] = VBTFcuts( "HLT_Ele27", RunPhase,
					 cand[1]->Pt(), cand[1]->Et(), ele_sclEta[iEle1], cand[1]->Eta(), ele_tkSumPt_dr03[iEle1], ele_ecalRecHitSumEt_dr03[iEle1], 
					 ele_hcalDepth1TowerSumEt_dr03[iEle1], ele_hcalDepth2TowerSumEt_dr03[iEle1], ele_expected_inner_hits[iEle1],
					 ele_deltaphiin[iEle1], ele_deltaetain[iEle1], ele_he[iEle1], ele_sigmaietaieta[iEle1],
					 ele_conv_dist[iEle1], ele_conv_dcot[iEle1], ele_fbrem[iEle1], ele_isConversion[iEle1] ) ;
	    

	    // check cuts passed by probe candidate
	    cutEle[1] = whichCuts( RunPhase, cand[1]->Pt(), cand[0]->Et(), ele_sclEta[iEle2], cand[1]->Eta(), ele_tkSumPt_dr03[iEle2], ele_ecalRecHitSumEt_dr03[iEle2], 
				   ele_hcalDepth1TowerSumEt_dr03[iEle2], ele_hcalDepth2TowerSumEt_dr03[iEle2], ele_expected_inner_hits[iEle2],
				   ele_deltaphiin[iEle2], ele_deltaetain[iEle2], ele_he[iEle2], ele_sigmaietaieta[iEle2],
				   ele_conv_dist[iEle2], ele_conv_dcot[iEle2], ele_fbrem[iEle2], ele_isConversion[iEle2] ) ;
	    fidu[1] = 0;
	    if ( fabs(ele_sclEta[iEle2]) < 2.5 && ( fabs(ele_sclEta[iEle2]) > 1.566 || fabs(ele_sclEta[iEle2])<1.4442 ) ) 
	      fidu[1] = 1 ;

	    if( cutEle[1]>0 && iEle2<=iEle1 ) continue; // prevents to create several times the same pair

	    // get the pair informations
	    total = (*cand[0]) + (*cand[1]) ;

	    // keep only pairs with Mee > 30 GeV
	    if( total.M() < 30. ) continue; 

	    pair_M = total.M() ;

	    pairIdx.clear();
	    pairIdx.push_back(iEle1);
	    pairIdx.push_back(iEle2);

	    for(int iP=0 ; iP<2 ; iP++) {

	      pair_cuts[iP] = cutEle[iP];
	      pair_fidu[iP] = fidu[iP];
	      pair_HLT_Ele27_cut[iP] = cut_HLT_Ele27[iP];
	      //
	      pair_eta[iP] = cand[iP]->Eta();
	      pair_sclEta[iP] = ele_sclEta[pairIdx[iP]];
	      pair_phi[iP] = cand[iP]->Phi();
	      pair_RCTeta[iP] = ele_RCTeta[pairIdx[iP]];
	      pair_RCTphi[iP] = ele_RCTphi[pairIdx[iP]];
	      //
	      pair_charge[iP] = ele_echarge[pairIdx[iP]];
	      pair_pT[iP] = cand[iP]->Pt();
	      pair_eT[iP] = cand[iP]->Et();
	      pair_sclEt[iP] = ele_sclEt[pairIdx[iP]];
	      pair_E[iP] = cand[iP]->E();
	      //
	      pair_L1iso[iP] = ele_RCTL1iso[pairIdx[iP]];
	      pair_L1noniso[iP] = ele_RCTL1noniso[pairIdx[iP]];
	      pair_L1iso_M[iP] = ele_RCTL1iso_M[pairIdx[iP]];
 	      pair_L1noniso_M[iP] = ele_RCTL1noniso_M[pairIdx[iP]];
	      pair_L1iso_To[iP] = ele_RCTL1iso_To[pairIdx[iP]];
 	      pair_L1noniso_To[iP] = ele_RCTL1noniso_To[pairIdx[iP]];
	      //
	      for(int iV=0 ; iV<10 ; iV++) {
		pair_RCTetVect[iP][iV] = ele_RCTetVect[pairIdx[iP]][iV];
		pair_RCTetaVect[iP][iV] = ele_RCTetaVect[pairIdx[iP]][iV];
		pair_RCTphiVect[iP][iV] = ele_RCTphiVect[pairIdx[iP]][iV]; 
		pair_L1isoVect[iP][iV] = ele_RCTL1isoVect[pairIdx[iP]][iV]; 
		pair_L1nonisoVect[iP][iV] = ele_RCTL1nonisoVect[pairIdx[iP]][iV];
 		pair_L1isoVect_M[iP][iV] = ele_RCTL1isoVect_M[pairIdx[iP]][iV];
 		pair_L1nonisoVect_M[iP][iV] = ele_RCTL1nonisoVect_M[pairIdx[iP]][iV];
 		pair_L1isoVect_To[iP][iV] = ele_RCTL1isoVect_To[pairIdx[iP]][iV];
 		pair_L1nonisoVect_To[iP][iV] = ele_RCTL1nonisoVect_To[pairIdx[iP]][iV];
	      } 
	      //
	      for(int iV=0 ; iV<50 ; iV++) {
		pair_TTetaVect[iP][iV] = ele_TTetaVect[pairIdx[iP]][iV];
		pair_TTphiVect[iP][iV] = ele_TTphiVect[pairIdx[iP]][iV];
		pair_TTetVect[iP][iV] = ele_TTetVect[pairIdx[iP]][iV];
	      }
	    }
	    
	    // Fill OUTTREE //
	    if(debug) cout << "--- Fill output tree" << endl;
	    outtree->Fill();

	  } // loop for probe
	} // endif ele1 is good tag candidate
      } // loop over electrons
     
    }//loop over events

  // Record tree
  outlog << "recording tree..." << endl;
  outtree->Write();
  outlog << "recorded !" << endl;
  outfile->Close();
  outlog << "file closed." << endl;

  return 1;

}
/////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////
// LOOP OVER THE TREE LIST //
/////////////////////////////
int looplist(TString nameChain, vector<TString> treeList, vector<TString> listTrees_To, 
	     int n1, int n2, int nEntries, int nEntriesTo,
	     TString dirIn, TString dirOut, ofstream& outcheck, 
	     bool dojson, TString RunPhase, bool GetTrigOnly, bool debug, bool b_mapcheck) {

  if(debug) cout << "- define JSON" << endl;

  // 1- JSON FILE READER //
  // 2010B        :   -> 2010B                           
  // May10        : 160404-163869 -> May10                                                                                      
  // PRV4         : 165071-168437 -> Prompt                                                                                        
  // Aug05 (PRV5) : 170053-172619 -> Aug05                                                                                         
  // Oct03 (PRV6) : 172620-175770 -> Prompt                                                                                       
  // 2011B-PRV1   : 175832-180296 -> Prompt                                                                                        
  //
  string jsonDir = "/data_CMS/cms/ndaci/ndaci_2011A/JSON/" ;
  const int nJson = 5;
  string jsonFile[nJson]; // 2010B, May10, Aug05, Prompt
  MAPJson jsonMap[nJson] ;
  //
  jsonFile[0] = jsonDir + "goodrunlist_json.txt" ;
  jsonFile[1] = jsonDir + "Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.txt" ;
  jsonFile[2] = jsonDir + "Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v3.txt" ;
  jsonFile[3] = jsonDir + "Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt" ;
  jsonFile[4] = jsonDir + "Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt" ;
  //
  for(int i=0 ; i<nJson ; i++)
    jsonMap[i] = readJSONFile(jsonFile[i]);
  //

  // 2- INPUT TrigOnly TREE //
  if(debug) cout << "- get TrigOnly trees --> chain" << endl;

  TChain * trigOnlyChain = new TChain(nameChain);
  MAPTrigOnly mapTo;

  if( GetTrigOnly ) {
    for(int iTo=0 ; iTo<listTrees_To.size() ; iTo++)
      trigOnlyChain->Add(listTrees_To[iTo]);
    
    if(debug) cout << "--> map TrigOnly entries --> MAPTrigOnly" << endl;
    mapTo = MapTrigOnlyEntries(trigOnlyChain, nEntriesTo, outcheck, jsonMap, dojson, debug, b_mapcheck);
  }

  // Loop over data trees
  if(debug) cout << "- loop over trees [" << n1 << ";" << n2 << "]" << endl;

  for(int iFile = n1 ; iFile < n2 ; iFile++) {
    if(debug) cout << "--> tree #" << iFile << endl;
    process(trigOnlyChain, jsonMap, mapTo, nameChain, treeList[iFile], iFile, nEntries, 
	    dirIn, dirOut, dojson, RunPhase, GetTrigOnly, debug);
  }  
  /*
    int process(TChain *trigOnlyChain, MAPJson *jsonMap, MAPTrigOnly mapTo, TString nameChain, TString file, int iFile, int nEntries, 
                TString dirIn, TString dirOut, bool dojson, TString RunPhase, bool debug)
  */

  return 1;
}
//////////////////////////////////////////////////////////////////////////

///////////////
// makePairs //
///////////////
void makePairs(int nEntries=-1, int nEntriesTo=-1, int iHalf=0, int nHalf=1,
	       TString dirOut="",
	       TString dirIn="",
	       TString data="",
	       TString data_To="",
	       TString nameChain="produceNtuple/eIDSimpleTree", 
	       bool GetTrigOnly=false, bool dojson=true, TString RunPhase="2011A", bool debug=false, bool b_mapcheck=false)
{
  // Output Log
  ofstream outcheck(dirOut+"mapTrigOnly.txt",ios::out);

  // Process the tree
  if(debug) cout << "Define tree list" << endl;

  vector<TString> treeList;
  vector<TString> listTrees_To;

  if(data=="SingEle_ReReco08Nov") {
    treeList = list_SingEle_ReReco08Nov();
  }
  else if(data=="SingEle_ReReco19Nov") {
    treeList = list_SingEle_ReReco19Nov();
  }
  else if(data=="SingEle_ReReco08Nov_Aod") {
    treeList = list_SingEle_08Nov_AodHits();
  }
  else if(data=="SingEle_ReReco19Nov_Aod") {
    treeList = list_SingEle_19Nov_AodHits();
  }
  else {
    cout << "asked data tag non recognized ! stopping..." << endl;
    return;
  }

  if( GetTrigOnly ) {
    if(data_To=="TrigOnly_ReReco08Nov") {
      listTrees_To = list_TrigOnly_ReReco08Nov();
    }
    else if(data_To=="TrigOnly_ReReco19Nov") {
      listTrees_To = list_TrigOnly_ReReco19Nov();
    }
    else if(data_To=="TrigOnly_transpaTag_08Nov") {
      listTrees_To = list_TrigOnly_transpaTag_08Nov();
    }
    else if(data_To=="TrigOnly_transpaTag_19Nov") {
      listTrees_To = list_TrigOnly_transpaTag_19Nov();
    }
    else if(data_To=="TrigOnly_transpaTag_08Nov_redo") {
      listTrees_To = list_TrigOnly_transpaTag_08Nov_redo();
    }
    else if(data_To=="TrigOnly_transpaTag_19Nov_redo") {
      listTrees_To = list_TrigOnly_transpaTag_19Nov_redo();
    }
    else if(data_To=="TrigOnly_Newkill_2011A") {
      listTrees_To = list_SingEle_Newkill_2011A();
    }
    else if(data_To=="TrigOnly_Newkill_2011B") {
      listTrees_To = list_SingEle_Newkill_2011B();
    }
  }

  if(debug) cout << "--> cut the list into " << nHalf << " parts" << endl;

  int n,n1,n2;
  n1=n2=n=0;

  n = treeList.size() / nHalf ;
  n1 = iHalf * n ;

  if( iHalf < (nHalf-1) )
    n2 = (iHalf+1) * n ;

  else if( iHalf == nHalf-1 )
    n2 = treeList.size();
 
  if(debug) cout << "Launch looplist" << endl;
  looplist(nameChain, treeList, listTrees_To, n1, n2, nEntries, nEntriesTo, dirIn, dirOut, outcheck, dojson, RunPhase, GetTrigOnly, debug, b_mapcheck);

  if(debug) cout << "DONE" << endl;

}
