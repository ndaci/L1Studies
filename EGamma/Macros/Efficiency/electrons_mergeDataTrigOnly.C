///////////////////////////////////////////////////////////////////////////////
/////////////////////////                                   ///////////////////
/////////////////////////    electrons_mergeDataTrigOnly    ///////////////////
/////////////////////////                                   ///////////////////
///////////////////////////////////////////////////////////////////////////////

// 1.3 : 
// 1.3.1 : makes " *= 2 " optional for 2010 menu
// 1.3.2 : menu = every EG between 0 and 30

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
			       MAPJson *jsonMap, int iJson, bool debug, bool b_mapcheck) {
  
  // Get trigOnlyChain information //
  if(debug) cout << "--- define trigOnlyChain branches" << endl;
  int nEvent_To, nRun_To, nLumi_To;
  bool isGoodRun;
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

    // JSON //
    if( iJson>-1 && iJson<5) {
      isGoodRun = AcceptEventByRunAndLumiSection(nRun_To, nLumi_To, jsonMap[iJson]);
      if(!isGoodRun) {
	outcheck << "trig failed JSON" << endl;
	continue;
      }
    }
    
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
	    TString dirIn, TString dirOut, TString usr_vbtf, double cut_et, TString usr_hlt, int iJson, TString RunPhase, bool GetTrigOnly, bool debug)  {


  // TrigOnly MAP //
  MAPTrigOnly::iterator iterMapTo;
  pair<int,int> runevtTo;
  int iEntryTo=-1;

  // OUTPUT FILE //
  if(debug) cout << "--- define output file : " ;
  ostringstream ossi("");
  ossi << iFile ;
  //
  TFile *outfile = new TFile(dirIn+"electrons_"+ossi.str()+".root","RECREATE");
  ofstream outlog(dirOut+"/logs/log_"+ossi.str()+".txt",ios::out);
  //ofstream outlogEle(dirOut+"/logsEle/log_"+ossi.str()+".txt",ios::out);
  ofstream outcheckL1(dirOut+"/checkL1/checkL1_"+ossi.str()+".txt",ios::out);

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
  int ele_expected_inner_hits_aod[10];
  //int ele_ambiguousGsfTracks[10];
  int ele_isConversion[10];
  int ele_echarge[10];
  //
  int ele_RCTeta[10], ele_RCTphi[10], ele_RCTL1iso[10], ele_RCTL1noniso[10], ele_RCTL1iso_To[10], ele_RCTL1noniso_To[10], ele_RCTL1iso_M[10], ele_RCTL1noniso_M[10];
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
  //myChain->SetBranchAddress ("trig_HLT_triggered", &m_HLT_triggered);
  //myChain->SetBranchAddress ("trig_HLT_pathsV", &m_HLT_pathsV);
  //myChain->SetBranchAddress ("trig_HLT_pathsV_check", &m_HLT_pathsV_check);
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
  myChain->SetBranchAddress("ele_expected_inner_hits",&ele_expected_inner_hits);
  myChain->SetBranchAddress("ele_expected_inner_hits_aod",&ele_expected_inner_hits_aod);
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
  trigOnlyChain->SetBranchAddress("trig_preL1emIso_N",     &trig_preL1emIso_N);
  trigOnlyChain->SetBranchAddress("trig_preL1emIso_ieta",  &trig_preL1emIso_ieta);
  trigOnlyChain->SetBranchAddress("trig_preL1emIso_iphi",  &trig_preL1emIso_iphi);
  trigOnlyChain->SetBranchAddress("trig_preL1emIso_rank",  &trig_preL1emIso_rank);
  //
  trigOnlyChain->SetBranchAddress("trig_preL1emNonIso_N",     &trig_preL1emNonIso_N);
  trigOnlyChain->SetBranchAddress("trig_preL1emNonIso_ieta",  &trig_preL1emNonIso_ieta);
  trigOnlyChain->SetBranchAddress("trig_preL1emNonIso_iphi",  &trig_preL1emNonIso_iphi);
  trigOnlyChain->SetBranchAddress("trig_preL1emNonIso_rank",  &trig_preL1emNonIso_rank);
  //
  trigOnlyChain->SetBranchAddress("trig_postL1emIso_N",     &trig_postL1emIso_N);
  trigOnlyChain->SetBranchAddress("trig_postL1emIso_ieta",  &trig_postL1emIso_ieta);
  trigOnlyChain->SetBranchAddress("trig_postL1emIso_iphi",  &trig_postL1emIso_iphi);
  trigOnlyChain->SetBranchAddress("trig_postL1emIso_rank",  &trig_postL1emIso_rank);
  //
  trigOnlyChain->SetBranchAddress("trig_postL1emNonIso_N",     &trig_postL1emNonIso_N);
  trigOnlyChain->SetBranchAddress("trig_postL1emNonIso_ieta",  &trig_postL1emNonIso_ieta);
  trigOnlyChain->SetBranchAddress("trig_postL1emNonIso_iphi",  &trig_postL1emNonIso_iphi);
  trigOnlyChain->SetBranchAddress("trig_postL1emNonIso_rank",  &trig_postL1emNonIso_rank);

  // Strip masking
  trigOnlyChain->SetBranchAddress("trig_strip_mask_N", &trig_strip_mask_N);
  //trigOnlyChain->SetBranchAddress("trig_strip_mask_status", &trig_strip_mask_status);
  trigOnlyChain->SetBranchAddress("trig_strip_mask_TTieta", &trig_strip_mask_TTieta);
  trigOnlyChain->SetBranchAddress("trig_strip_mask_TTiphi", &trig_strip_mask_TTiphi);
  trigOnlyChain->SetBranchAddress("trig_strip_mask_StripID", &trig_strip_mask_StripID);
  trigOnlyChain->SetBranchAddress("trig_strip_mask_PseudoStripID", &trig_strip_mask_PseudoStripID);
  trigOnlyChain->SetBranchAddress("trig_strip_mask_TccID", &trig_strip_mask_TccID);
  trigOnlyChain->SetBranchAddress("trig_strip_mask_CCU", &trig_strip_mask_CCU);
  trigOnlyChain->SetBranchAddress("trig_strip_mask_xtal_ix", &trig_strip_mask_xtal_ix);
  trigOnlyChain->SetBranchAddress("trig_strip_mask_xtal_iy", &trig_strip_mask_xtal_iy);
  trigOnlyChain->SetBranchAddress("trig_strip_mask_xtal_iz", &trig_strip_mask_xtal_iz);
  //
  // Crystal masking
  trigOnlyChain->SetBranchAddress("trig_xtal_mask_N", &trig_xtal_mask_N);
  trigOnlyChain->SetBranchAddress("trig_xtal_mask_ieta", &trig_xtal_mask_ieta);
  trigOnlyChain->SetBranchAddress("trig_xtal_mask_iphi", &trig_xtal_mask_iphi);
  trigOnlyChain->SetBranchAddress("trig_xtal_mask_TTieta", &trig_xtal_mask_TTieta);
  trigOnlyChain->SetBranchAddress("trig_xtal_mask_TTiphi", &trig_xtal_mask_TTiphi);
  trigOnlyChain->SetBranchAddress("trig_xtal_mask_Rieta", &trig_xtal_mask_Rieta);
  trigOnlyChain->SetBranchAddress("trig_xtal_mask_Riphi", &trig_xtal_mask_Riphi);
  trigOnlyChain->SetBranchAddress("trig_xtal_mask_status", &trig_xtal_mask_status);
  trigOnlyChain->SetBranchAddress("trig_xtal_mask_EBEE", &trig_xtal_mask_EBEE);

  // Masking
  trigOnlyChain->SetBranchAddress("trig_nMaskedRCT",      &trig_nMaskedRCT);      
  trigOnlyChain->SetBranchAddress("trig_iMaskedRCTeta",   &trig_iMaskedRCTeta);                                          
  trigOnlyChain->SetBranchAddress("trig_iMaskedRCTcrate", &trig_iMaskedRCTcrate);
  trigOnlyChain->SetBranchAddress("trig_iMaskedRCTphi",   &trig_iMaskedRCTphi);
  trigOnlyChain->SetBranchAddress("trig_nMaskedCh",       &trig_nMaskedCh);    
  trigOnlyChain->SetBranchAddress("trig_iMaskedTTeta",    &trig_iMaskedTTeta);   
  trigOnlyChain->SetBranchAddress("trig_iMaskedTTphi",    &trig_iMaskedTTphi);      	


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
  //trigOnlyChain->SetBranchAddress("trig_tower_FG",  &trig_tower_FG);
 
  // modified collections
  trigOnlyChain->SetBranchAddress("trig_tower_N_M", &trig_tower_N_M);
  trigOnlyChain->SetBranchAddress("trig_tower_ieta_M",  &trig_tower_ieta_M);
  trigOnlyChain->SetBranchAddress("trig_tower_iphi_M",  &trig_tower_iphi_M);
  trigOnlyChain->SetBranchAddress("trig_tower_adc_M",  &trig_tower_adc_M);
  trigOnlyChain->SetBranchAddress("trig_tower_sFGVB_M",  &trig_tower_sFGVB_M);
  //trigOnlyChain->SetBranchAddress("trig_tower_FG_M",  &trig_tower_FG_M);

  trigOnlyChain->SetBranchAddress("trig_tower_N_E",     &trig_tower_N_E);
  trigOnlyChain->SetBranchAddress("trig_tower_ieta_E",  &trig_tower_ieta_E);
  trigOnlyChain->SetBranchAddress("trig_tower_iphi_E",  &trig_tower_iphi_E);
  trigOnlyChain->SetBranchAddress("trig_tower_adc_E",   &trig_tower_adc_E);
  trigOnlyChain->SetBranchAddress("trig_tower_sFGVB_E", &trig_tower_sFGVB_E);
  //trigOnlyChain->SetBranchAddress("trig_tower_FG_E", &trig_tower_FG_E);
  
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

  const int nColl = 2; // _N / _M
  const int nEG = 30;

  vector<int> menu;
  for(int iEG=0; iEG<nEG ; iEG++) {
    menu.push_back(iEG);
    if(RunPhase=="2010") menu[iEG] *= 2 ; // retablit le bon menu pour 2010 !  
  }

  vector<int> firedEG[2]; // [N/M]

  int ele_trigEG[nEG], ele_trigEG_M[nEG];
  for(int iEG=0;iEG<nEG;iEG++)
    ele_trigEG[iEG] = ele_trigEG_M[iEG] = 0;

  // OutTree branches //
  double ele_out_eta;
  double ele_out_sclEta;
  double ele_out_sclEt;
  double ele_out_phi;
  //double ele_out_charge;
  double ele_out_pT;
  double ele_out_eT;  
  double ele_out_E;
  int ele_out_cut, ele_out_fidu;
  int ele_out_cut_HLT_Ele27;
  //
  int ele_out_severityLevelSeed;
  double ele_out_he, ele_out_sigmaietaieta;
  double ele_out_hcalDepth1TowerSumEt_dr03, ele_out_hcalDepth2TowerSumEt_dr03;
  double ele_out_ecalRecHitSumEt_dr03, ele_out_tkSumPt_dr03;
  double ele_out_deltaphiin, ele_out_deltaetain;
  double ele_out_fbrem;
  int ele_out_expected_inner_hits;
  int ele_out_expected_inner_hits_aod;
  //
  int ele_out_RCTeta, ele_out_RCTphi, ele_out_RCTL1iso, ele_out_RCTL1noniso, ele_out_RCTL1iso_M, ele_out_RCTL1noniso_M, ele_out_RCTL1iso_To, ele_out_RCTL1noniso_To;
  int ele_out_TTetaVect[50], ele_out_TTphiVect[50];
  double ele_out_TTetVect[50];
  int ele_out_RCTetaVect[10], ele_out_RCTphiVect[10], ele_out_RCTL1isoVect[10], 
    ele_out_RCTL1nonisoVect[10],ele_out_RCTL1isoVect_M[10], ele_out_RCTL1nonisoVect_M[10],ele_out_RCTL1isoVect_To[10], ele_out_RCTL1nonisoVect_To[10];
  double ele_out_RCTetVect[10];
  //                  //

  TTree * outtree = new TTree("Electrons","Electrons");

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
//   outtree->Branch ("trig_HLT_triggered", &m_HLT_triggered, 256000,0);
//   outtree->Branch ("trig_HLT_pathsV", &m_HLT_pathsV, 256000,0);
//   outtree->Branch ("trig_HLT_pathsV_check", &m_HLT_pathsV_check, 256000,0);
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
  outtree->Branch("trig_strip_mask_status", &trig_strip_mask_status, "trig_strip_mask_status[trig_strip_mask_N]/I");
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

  outtree->Branch("ele_cut",&ele_out_cut,"ele_cut/I");
  // 0 : noCut | 1 : VBTF 95 | 2 : VBTF 80 | 3 : VBTF 60
  outtree->Branch("ele_cut_HLT_Ele27",&ele_out_cut_HLT_Ele27,"ele_cut_HLT_Ele27/I");
  outtree->Branch("ele_fidu",&ele_out_fidu,"ele_fidu/I");
  //
  outtree->Branch("ele_eta",&ele_out_eta,"ele_eta/D");
  outtree->Branch("ele_sclEta",&ele_out_sclEta,"ele_sclEta/D");
  outtree->Branch("ele_phi",&ele_out_phi,"ele_phi/D");
  outtree->Branch("ele_RCTeta",&ele_out_RCTeta,"ele_RCTeta/I");
  outtree->Branch("ele_RCTphi",&ele_out_RCTphi,"ele_RCTphi/I");
  //
  //outtree->Branch("ele_charge",&ele_out_charge,"ele_charge/I");
  outtree->Branch("ele_pT",&ele_out_pT,"ele_pT/D");
  outtree->Branch("ele_eT",&ele_out_eT,"ele_eT/D");
  outtree->Branch("ele_sclEt",&ele_out_sclEt,"ele_sclEt/D");
  outtree->Branch("ele_E",&ele_out_E,"ele_E/D");

  outtree->Branch("ele_deltaetain",&ele_out_deltaetain,"ele_deltaetain/D");
  outtree->Branch("ele_deltaphiin",&ele_out_deltaphiin,"ele_deltaphiin/D");
  outtree->Branch("ele_tkSumPt_dr03",&ele_out_tkSumPt_dr03,"ele_tkSumPt_dr03/D");
  outtree->Branch("ele_ecalRecHitSumEt_dr03",&ele_out_ecalRecHitSumEt_dr03,"ele_ecalRecHitSumEt_dr03/D");
  outtree->Branch("ele_hcalDepth1TowerSumEt_dr03",&ele_out_hcalDepth1TowerSumEt_dr03,"ele_hcalDepth1TowerSumEt_dr03/D");
  outtree->Branch("ele_hcalDepth2TowerSumEt_dr03",&ele_out_hcalDepth2TowerSumEt_dr03,"ele_hcalDepth2TowerSumEt_dr03/D");
  outtree->Branch("ele_expected_inner_hits",&ele_out_expected_inner_hits,"ele_expected_inner_hits/I");
  outtree->Branch("ele_expected_inner_hits_aod",&ele_out_expected_inner_hits_aod,"ele_expected_inner_hits_aod/I");
  outtree->Branch("ele_sigmaietaieta",&ele_out_sigmaietaieta,"ele_sigmaietaieta/D");
  outtree->Branch("ele_he",&ele_out_he,"ele_he/D");
  outtree->Branch("ele_severityLevelSeed",&ele_out_severityLevelSeed,"ele_severityLevelSeed/I");
  outtree->Branch("ele_fbrem",&ele_out_fbrem,"ele_fbrem/D");

  outtree->Branch("ele_TTetaVect", &ele_out_TTetaVect,"ele_TTetaVect[50]/I");
  outtree->Branch("ele_TTphiVect", &ele_out_TTphiVect,"ele_TTphiVect[50]/I");
  outtree->Branch("ele_TTetVect", &ele_out_TTetVect,"ele_TTetVect[50]/D");

  outtree->Branch("ele_RCTL1iso",&ele_out_RCTL1iso,"ele_RCTL1iso/I");
  outtree->Branch("ele_RCTL1noniso",&ele_out_RCTL1noniso,"ele_RCTL1noniso/I");
  outtree->Branch("ele_RCTL1iso_M",&ele_out_RCTL1iso_M,"ele_RCTL1iso_M/I");
  outtree->Branch("ele_RCTL1noniso_M",&ele_out_RCTL1noniso_M,"ele_RCTL1noniso_M/I");
  outtree->Branch("ele_RCTL1iso_To",&ele_out_RCTL1iso_To,"ele_RCTL1iso_To/I");
  outtree->Branch("ele_RCTL1noniso_To",&ele_out_RCTL1noniso_To,"ele_RCTL1noniso_To/I");
  //
  outtree->Branch("ele_RCTetVect",&ele_out_RCTetVect,"ele_RCTetVect[10]/D");
  outtree->Branch("ele_RCTetaVect",&ele_out_RCTetaVect,"ele_RCTetaVect[10]/I");
  outtree->Branch("ele_RCTphiVect",&ele_out_RCTphiVect,"ele_RCTphiVect[10]/I");
  outtree->Branch("ele_RCTL1isoVect",&ele_out_RCTL1isoVect,"ele_RCTL1isoVect[10]/I");
  outtree->Branch("ele_RCTL1nonisoVect",&ele_out_RCTL1nonisoVect,"ele_RCTL1nonisoVect[10]/I");
  outtree->Branch("ele_RCTL1isoVect_M",&ele_out_RCTL1isoVect_M,"ele_RCTL1isoVect_M[10]/I");
  outtree->Branch("ele_RCTL1nonisoVect_M",&ele_out_RCTL1nonisoVect_M,"ele_RCTL1nonisoVect_M[10]/I");
  outtree->Branch("ele_RCTL1isoVect_To",&ele_out_RCTL1isoVect_To,"ele_RCTL1isoVect_To[10]/I");
  outtree->Branch("ele_RCTL1nonisoVect_To",&ele_out_RCTL1nonisoVect_To,"ele_RCTL1nonisoVect_To[10]/I");

  outtree->Branch("ele_trigEG",&ele_trigEG,"ele_trigEG[30]/I");
  outtree->Branch("ele_trigEG_M",&ele_trigEG_M,"ele_trigEG_M[30]/I");

  // USEFUL VARIABLES //
  vector<int> pairIdx;
  TLorentzVector* cand;
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
	ele_RCTL1iso[i]=0;
	ele_RCTL1noniso[i]=0;
	ele_RCTL1iso_M[i]=0;
	ele_RCTL1noniso_M[i]=0;

	for(int j=0 ; j<10 ; j++) {
	  ele_RCTL1isoVect[i][j]=0;
	  ele_RCTL1nonisoVect[i][j]=0;
	  ele_RCTL1isoVect_M[i][j]=0;
	  ele_RCTL1nonisoVect_M[i][j]=0;
	}
      }

      myChain->GetEntry (iEvent) ;
     
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
      if(debug) cout << "iJson = " << iJson << endl;
      if( iJson>-1 && iJson<5) {
	isGoodRun = AcceptEventByRunAndLumiSection(nRun, nLumi, jsonMap[iJson]);
	if(!isGoodRun) {
	  outlog << "failed JSON" << endl;
	  continue;
	}
      }

      // HLT selection //
      if(usr_hlt=="unbias")
	if(trig_HLT_path[0]==0) continue;


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
      for( int iEle=0 ; iEle<ele_N ; iEle++ ) {
	if(debug) cout << "--- get ele #" << iEle << endl;
	cand = (TLorentzVector*) (electrons->At (iEle)) ;
	if(debug) cout << "--- got it" << endl;

	// remove EE electrons
	//if( fabs(cand->Eta()) >= 1.479 ) continue;
	
	// severity selection
	//if( ele_severityLevelSeed[iEle] >= 3 ) continue;
	
	/*
	double trackIsoRel03 = ele_tkSumPt_dr03[iEle] / (cand->Pt());
	double ecalIsoRel03  = ele_ecalRecHitSumEt_dr03[iEle]/(cand->Pt());
	double hcalIsoRel03  =(ele_hcalDepth1TowerSumEt_dr03[iEle]+ele_hcalDepth2TowerSumEt_dr03[iEle])/(cand->Pt());

	outlogEle << "ele_deltaphiin=" << ele_deltaphiin[iEle]
		  << "   fabs="              << fabs(ele_deltaphiin[iEle])
		  << endl
		  << "ele_deltaetain=" << ele_deltaetain[iEle]
		  << "   fabs="              << fabs(ele_deltaetain[iEle])
		  << endl
		  << "ele_he="         << ele_he[iEle]
		  << endl
		  << "ele_sigmaietaieta=" << ele_sigmaietaieta[iEle]
		  << endl
		  << "trackIsoRel03="   << trackIsoRel03
		  << "   fabs="         << fabs(trackIsoRel03)
		  << endl
		  << "hcalIsoRel03=" << hcalIsoRel03
		  << "   fabs=" << fabs(hcalIsoRel03)
		  << endl
		  << "ecalIsoRel03=" << ecalIsoRel03
		  << "   fabs=" << fabs(ecalIsoRel03)
		  << endl
		  << "ele_expected_inner_hits=" << ele_expected_inner_hits[iEle]
		  << endl;
	*/

	// cuts passed by the electron
	ele_out_cut_HLT_Ele27 = 0;
	ele_out_cut_HLT_Ele27 = (int) VBTFcuts( "HLT_Ele27", RunPhase,
				      cand->Pt(), cand->Et(), ele_sclEta[iEle], cand->Eta(), ele_tkSumPt_dr03[iEle], ele_ecalRecHitSumEt_dr03[iEle], 
				      ele_hcalDepth1TowerSumEt_dr03[iEle], ele_hcalDepth2TowerSumEt_dr03[iEle], ele_expected_inner_hits_aod[iEle],
				      ele_deltaphiin[iEle], ele_deltaetain[iEle], ele_he[iEle], ele_sigmaietaieta[iEle],
				      ele_conv_dist[iEle], ele_conv_dcot[iEle], ele_fbrem[iEle], ele_isConversion[iEle] ) ;

	ele_out_cut = 0;
	ele_out_cut = whichCuts( RunPhase, cand->Pt(), cand->Et(), ele_sclEta[iEle], cand->Eta(), ele_tkSumPt_dr03[iEle], ele_ecalRecHitSumEt_dr03[iEle], 
			     ele_hcalDepth1TowerSumEt_dr03[iEle], ele_hcalDepth2TowerSumEt_dr03[iEle], ele_expected_inner_hits_aod[iEle],
			     ele_deltaphiin[iEle], ele_deltaetain[iEle], ele_he[iEle], ele_sigmaietaieta[iEle],
			     ele_conv_dist[iEle], ele_conv_dcot[iEle], ele_fbrem[iEle], ele_isConversion[iEle] ) ;
	ele_out_fidu = 0;
	if ( fabs(ele_sclEta[iEle]) < 2.5 && ( fabs(ele_sclEta[iEle]) > 1.566 || fabs(ele_sclEta[iEle])<1.4442 ) ) 
	  ele_out_fidu = 1 ;


	// Apply user cuts //
	int cut = 0;
	if(usr_vbtf=="WP95") cut = 1;
	else if(usr_vbtf=="WP80") cut = 2;
	else if(usr_vbtf=="WP60") cut = 3;

	if( !(ele_out_cut>=cut && cand->Et()>=cut_et) )
	  continue;


	// Check L1 triggering //
	for(int iColl=0 ; iColl<2 ; iColl++) {
	  firedEG[iColl].clear();
	  firedEG[iColl].resize(nEG,0);	  
	}
	for(int iEG=0 ; iEG<nEG ; iEG++)
	  ele_trigEG[iEG] = ele_trigEG_M[iEG] = 0;

	for(int iR=0 ; iR < 10 ; iR++) {
	  
	  globalFireL1_Normal( ele_RCTetVect[iEle][iR] , ele_RCTL1nonisoVect[iEle][iR],
			       ele_RCTL1isoVect[iEle][iR], firedEG[0], menu); 
	  
	  globalFireL1_Normal( ele_RCTetVect[iEle][iR] , ele_RCTL1nonisoVect_M[iEle][iR],
			       ele_RCTL1isoVect_M[iEle][iR], firedEG[1], menu); 
	  
	  if(debug) cout << ele_RCTetVect[iEle][iR] << " " << ele_RCTL1nonisoVect[iEle][iR] << " " << ele_RCTL1isoVect[iEle][iR] << endl;
	  
	}
      
	for(int iEG=0 ; iEG<nEG ; iEG++) {
	  ele_trigEG[iEG] = firedEG[0][iEG];
	  ele_trigEG_M[iEG] = firedEG[1][iEG];
	}	  
	
	// Fill variables //
	ele_out_eta = cand->Eta();
	ele_out_phi = cand->Phi();
	ele_out_pT = cand->Pt();
	ele_out_eT = cand->Et();  
	ele_out_E = cand->E();
	//
	ele_out_severityLevelSeed = ele_severityLevelSeed[iEle];
	ele_out_he = ele_he[iEle];
	ele_out_sigmaietaieta = ele_sigmaietaieta[iEle];
	ele_out_hcalDepth1TowerSumEt_dr03 = ele_hcalDepth1TowerSumEt_dr03[iEle];
	ele_out_hcalDepth2TowerSumEt_dr03 = ele_hcalDepth2TowerSumEt_dr03[iEle];
	ele_out_ecalRecHitSumEt_dr03 = ele_ecalRecHitSumEt_dr03[iEle];
	ele_out_tkSumPt_dr03 = ele_tkSumPt_dr03[iEle];
	ele_out_deltaphiin = ele_deltaphiin[iEle]; 
	ele_out_deltaetain = ele_deltaetain[iEle];
	ele_out_fbrem = ele_fbrem[iEle];
	ele_out_expected_inner_hits = ele_expected_inner_hits[iEle];
	ele_out_expected_inner_hits_aod = ele_expected_inner_hits_aod[iEle];
	//
	ele_out_RCTeta = ele_RCTeta[iEle];
	ele_out_RCTphi = ele_RCTphi[iEle]; 
	ele_out_RCTL1iso = ele_RCTL1iso[iEle]; 
	ele_out_RCTL1noniso = ele_RCTL1noniso[iEle]; 
	ele_out_RCTL1iso_M = ele_RCTL1iso_M[iEle]; 
	ele_out_RCTL1noniso_M = ele_RCTL1noniso_M[iEle];
	ele_out_RCTL1iso_To = ele_RCTL1iso_To[iEle]; 
	ele_out_RCTL1noniso_To = ele_RCTL1noniso_To[iEle];
	//
	for(int i=0;i<50;i++) {
	  ele_out_TTetaVect[i] = ele_TTetaVect[iEle][i];
	  ele_out_TTphiVect[i] = ele_TTphiVect[iEle][i];
	  ele_out_TTetVect[i] = ele_TTetVect[iEle][i];
	}
	for(int i=0;i<10;i++) {
	  ele_out_RCTetaVect[i] = ele_RCTetaVect[iEle][i];
	  ele_out_RCTphiVect[i] = ele_RCTphiVect[iEle][i];
	  ele_out_RCTL1isoVect[i] = ele_RCTL1isoVect[iEle][i];
	  ele_out_RCTL1nonisoVect[i] = ele_RCTL1nonisoVect[iEle][i];
	  ele_out_RCTL1isoVect_M[i] = ele_RCTL1isoVect_M[iEle][i];
	  ele_out_RCTL1nonisoVect_M[i] = ele_RCTL1nonisoVect_M[iEle][i];
	  ele_out_RCTL1isoVect_To[i] = ele_RCTL1isoVect_To[iEle][i];
	  ele_out_RCTL1nonisoVect_To[i] = ele_RCTL1nonisoVect_To[iEle][i];
	  ele_out_RCTetVect[i] = ele_RCTetVect[iEle][i];
	}
	
	outtree->Fill();

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
	     TString dirIn, TString dirOut, TString usr_vbtf, double cut_et, TString usr_hlt, ofstream& outcheck, 
	     int iJson, TString RunPhase, bool GetTrigOnly, bool debug, bool b_mapcheck) {

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
    mapTo = MapTrigOnlyEntries(trigOnlyChain, nEntriesTo, outcheck, jsonMap, iJson, debug, b_mapcheck);
  }

  // Loop over data trees
  if(debug) cout << "- loop over trees [" << n1 << ";" << n2 << "]" << endl;

  for(int iFile = n1 ; iFile < n2 ; iFile++) {
    if(debug) cout << "--> tree #" << iFile << endl;
    process(trigOnlyChain, jsonMap, mapTo, nameChain, treeList[iFile], iFile, nEntries, 
	    dirIn, dirOut, usr_vbtf, cut_et, usr_hlt, iJson, RunPhase, GetTrigOnly, debug);
  }  
  /*
    int process(TChain *trigOnlyChain, MAPJson *jsonMap, MAPTrigOnly mapTo, TString nameChain, TString file, int iFile, int nEntries, 
                TString dirIn, TString dirOut, TString usr_vbtf, double cut_et, bool dojson, TString RunPhase, bool debug)
  */

  return 1;
}
//////////////////////////////////////////////////////////////////////////

/////////////////////////////////
// electrons_mergeDataTrigOnly //
/////////////////////////////////
void electrons_mergeDataTrigOnly(int nEntries=-1, int nEntriesTo=-1, int iHalf=0, int nHalf=1,
				 TString dirOut="",
				 TString dirIn="",
				 TString data="",
				 TString data_To="",
				 TString usr_vbtf="WP95", 
				 double cut_et=5, 
				 TString usr_hlt="unbias",
				 TString nameChain="produceNtuple/eIDSimpleTree", 
				 bool GetTrigOnly=false, int iJson=-1, 
				 TString RunPhase="2011A", bool debug=false, bool b_mapcheck=false)
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
  else if(data=="AN2010_Data_EGMonitor") {
    treeList = list_AN2010_Data_EGMonitor();
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
    else if(data_To=="TrigOnly_transpaTag_08Nov_redo") {
      listTrees_To = list_TrigOnly_transpaTag_08Nov_redo();
    }
    else if(data_To=="TrigOnly_transpaTag_19Nov_redo") {
      listTrees_To = list_TrigOnly_transpaTag_19Nov_redo();
    }
    else if(data_To=="AN2010_TrigOnly_EGMonitor_10adc") {
      listTrees_To = list_an2010_Trig_10adc();
    }
    else if(data_To=="AN2010_TrigOnly_EGMonitor_15adc") {
      listTrees_To = list_an2010_Trig_15adc();
    }
    else if(data_To=="AN2010_TrigOnly_EGMonitor_17adc") {
      listTrees_To = list_an2010_Trig_17adc();
    }
    else if(data_To=="AN2010_TrigOnly_EGMonitor_19adc") {
      listTrees_To = list_an2010_Trig_19adc();
    }
    else if(data_To=="AN2010_TrigOnly_EGMonitor_21adc") {
      listTrees_To = list_an2010_Trig_21adc();
    }
    else if(data_To=="AN2010_TrigOnly_EGMonitor_23adc") {
      listTrees_To = list_an2010_Trig_23adc();
    }
    else if(data_To=="AN2010_TrigOnly_EGMonitor_30adc") {
      listTrees_To = list_an2010_Trig_30adc();
    }
    else if(data_To=="AN2010_TrigOnly_EGMonitor_35adc") {
      listTrees_To = list_an2010_Trig_35adc();
    }
    else if(data_To=="AN2010_TrigOnly_EGMonitor_40adc") {
      listTrees_To = list_an2010_Trig_40adc();
    }
    else if(data_To=="AN2010_TrigOnly_EGMonitor_50adc") {
      listTrees_To = list_an2010_Trig_50adc();
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
  looplist(nameChain, treeList, listTrees_To, n1, n2, nEntries, nEntriesTo, dirIn, dirOut, 
	   usr_vbtf, cut_et, usr_hlt, outcheck, iJson, RunPhase, GetTrigOnly, debug, b_mapcheck);

  if(debug) cout << "DONE" << endl;

}
//#include "/home/llr/cms/ndaci/SKWork/macro/skEfficiency/analyzers/COMMON/baseFuncNad.1.5.h"
