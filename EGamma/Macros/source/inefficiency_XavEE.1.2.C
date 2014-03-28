///////////////////////////////////////////////////////////////////////////////
/////////////////////////                                   ///////////////////
/////////////////////////    inefficiency                   ///////////////////
/////////////////////////                                   ///////////////////
///////////////////////////////////////////////////////////////////////////////

#include "/home/llr/cms/ndaci/SKWork/macro/skEfficiency/analyzers/COMMON/baseFuncNad.1.5.h"

// 1.2.7 : all inefficiency cases are implemented
// 1.3 : implements counting
// 1.3.1 :
// 1.3.2 : prioritize counting + change spread definition

typedef vector< vector< int > > VectMask ;

int matchL1toEle(int trig_L1emIso_N, int* trig_L1emIso_ieta, int* trig_L1emIso_iphi, int* trig_L1emIso_rank,
		 int ele_RCTeta, int ele_RCTphi, int &ele_RCTL1iso,
		 int* ele_RCTetaVect, int* ele_RCTphiVect, int* ele_RCTL1isoVect, const int nReg, ofstream& outcheckL1) {

  int nC = trig_L1emIso_N;
  if( nC>4 || nC<0 ) nC=4;

  outcheckL1 << "ELE : (" << ele_RCTeta << ";" << ele_RCTphi << ")=" << ele_RCTL1iso << endl;
  
  for(int iC=0 ; iC<nC ; iC++) {
    if( trig_L1emIso_ieta[iC]== ele_RCTeta && trig_L1emIso_iphi[iC]== ele_RCTphi )
      ele_RCTL1iso = trig_L1emIso_rank[iC];
    else
      ele_RCTL1iso = -999;

    for(int iR=0 ; iR<nReg ; iR++) {
      if( trig_L1emIso_ieta[iC]== ele_RCTetaVect[iR] && trig_L1emIso_iphi[iC]== ele_RCTphiVect[iR] )
	ele_RCTL1isoVect[iR] = trig_L1emIso_rank[iC];
      else
	ele_RCTL1isoVect[iR] = -999;
    }
  }

  return 1;
}

VectMask checkMasking(int trig_nMaskedCh,    int* trig_iMaskedTTeta,      int* trig_iMaskedTTphi,
		      int trig_strip_mask_N, int* trig_strip_mask_TTieta, int* trig_strip_mask_TTiphi,
		      int trig_xtal_mask_N,  int* trig_xtal_mask_TTieta,  int* trig_xtal_mask_TTiphi,
		      int TTieta, int TTiphi)
{

  VectMask answer;
  
  // check tower masking
  vector<int> towmask;
  for(int i=0 ; i<trig_nMaskedCh ; i++) {
    if( trig_iMaskedTTeta[i]==TTieta && trig_iMaskedTTphi[i]==TTiphi )
      towmask.push_back(i);
  }

  // check strip masking
  vector<int> stripmask;
  for(int i=0 ; i<trig_strip_mask_N ; i++) {
    if( trig_strip_mask_TTieta[i]==TTieta && trig_strip_mask_TTiphi[i]==TTiphi )
      stripmask.push_back(i);
  }
  
  // check strip masking
  vector<int> xtalmask;
  for(int i=0 ; i<trig_strip_mask_N ; i++) {
    if( trig_xtal_mask_TTieta[i]==TTieta && trig_xtal_mask_TTiphi[i]==TTiphi )
      xtalmask.push_back(i);
  }
  
  // merge it all
  answer.push_back(towmask);
  answer.push_back(stripmask);
  answer.push_back(xtalmask);

  return answer;

}

int analyzeTree(TChain* myChain, int nEntAsk, double cutet, bool debug, ofstream& outlog)
{

  // VARIABLES //
  int nEvent, nRun, nLumi ; // data
  int trig_HLT_path[4]; // unbias, EG5, EG8, EG12

  // Vertices //
  int _vtx_N;
  double _vtx_x[200], _vtx_y[200], _vtx_z[200];
  double _vtx_normalizedChi2[200], _vtx_ndof[200], _vtx_nTracks[200], _vtx_d0[200];

  // Electrons
  //
  double ele_sclEta, ele_sclEt, ele_eT, ele_pT, ele_E, ele_eta, ele_phi;
  int ele_cut, ele_fidu;
  //
  double ele_he, ele_sigmaietaieta, ele_hcalDepth1TowerSumEt_dr03, ele_hcalDepth2TowerSumEt_dr03,ele_ecalRecHitSumEt_dr03, ele_tkSumPt_dr03,
    ele_deltaphiin, ele_deltaetain, ele_conv_dist, ele_conv_dcot, ele_fbrem;
  int ele_expected_inner_hits, ele_expected_inner_hits_aod,ele_isConversion, ele_severityLevelSeed;
  //
  int ele_RCTL1iso_pre, ele_RCTL1noniso_pre, ele_RCTL1isoVect_pre[10], ele_RCTL1nonisoVect_pre[10];
  int ele_RCTL1iso_post, ele_RCTL1noniso_post, ele_RCTL1isoVect_post[10], ele_RCTL1nonisoVect_post[10];
  //
  int ele_RCTeta, ele_RCTphi, ele_RCTL1iso, ele_RCTL1noniso, ele_RCTL1iso_M, ele_RCTL1noniso_M, ele_RCTL1iso_To, ele_RCTL1noniso_To;
  int ele_TTetaVect[50], ele_TTphiVect[50];
  double ele_TTetVect[50];
  int ele_RCTetaVect[10], ele_RCTphiVect[10], ele_RCTL1isoVect[10], 
    ele_RCTL1nonisoVect[10],ele_RCTL1isoVect_M[10], ele_RCTL1nonisoVect_M[10],ele_RCTL1isoVect_To[10], ele_RCTL1nonisoVect_To[10];
  double ele_RCTetVect[10];

  // Electron triggering
  const int nEG=30;
  //int EG_2011[nEG] = {0,1,2,5,7,8,10,12,13,15,18,20,22,30};
  vector<int> menu;
  for(int i=0; i<nEG ; i++) {
    //EG_2011[i] *= 2 ; // erreur de menu (ele merge 1.3.1 non corrigé)
    menu.push_back(i);
  }

  int ele_trigEG[nEG], ele_trigEG_M[nEG], ele_trigEG_pre[nEG], ele_trigEG_post[nEG];
  for(int iEG=0;iEG<nEG;iEG++)
    ele_trigEG[iEG] = ele_trigEG_M[iEG] = ele_trigEG_pre[iEG] = ele_trigEG_post[iEG] = 0;

  // TP info
  const int nTow = 4032;
  int trig_tower_N,trig_tower_ieta[nTow],trig_tower_iphi[nTow],trig_tower_adc[nTow],trig_tower_sFGVB[nTow],trig_tower_FG[nTow]; 
  int trig_tower_N_M,trig_tower_ieta_M[nTow],trig_tower_iphi_M[nTow],trig_tower_adc_M[nTow],trig_tower_sFGVB_M[nTow],trig_tower_FG_M[nTow]; 
  int trig_tower_N_E,trig_tower_ieta_E[nTow],trig_tower_iphi_E[nTow],trig_tower_adc_E[nTow][5],trig_tower_sFGVB_E[nTow][5],trig_tower_FG_E[nTow][5];

  // HCAL TP
  const int nHtow=5000;
  int trig_tower_hcal_N, trig_tower_hcal_ieta[nHtow], trig_tower_hcal_iphi[nHtow], trig_tower_hcal_FG[nHtow],trig_tower_hcal_et[nHtow];

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
  for(int iTow=0 ; iTow<nHtow ; iTow++) {
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

  ele_RCTL1iso_pre = ele_RCTL1iso_post = -999;
  ele_RCTL1noniso_pre = ele_RCTL1noniso_post = -999;
  for(int i=0;i<10;i++)
    ele_RCTL1isoVect_pre[i] = ele_RCTL1nonisoVect_pre[i] = ele_RCTL1isoVect_post[i] = ele_RCTL1nonisoVect_post[i] = -999;

  // INPUT TREE //

  myChain->SetBranchAddress("nEvent",&nEvent);
  myChain->SetBranchAddress("nRun",&nRun);
  myChain->SetBranchAddress("nLumi",&nLumi);

  myChain->SetBranchAddress("trig_HLT_path",&trig_HLT_path);

  // Electrons
  myChain->SetBranchAddress("ele_fidu",&ele_fidu);
  myChain->SetBranchAddress("ele_cut",&ele_cut);
  // 0 : noCut | 1 : VBTF 95 | 2 : VBTF 80 | 3 : VBTF 60
  //
  myChain->SetBranchAddress("ele_eta",&ele_eta);
  myChain->SetBranchAddress("ele_phi",&ele_phi);
  //
  myChain->SetBranchAddress("ele_sclEta",&ele_sclEta);
  myChain->SetBranchAddress("ele_sclEt",&ele_sclEt);
  //
  myChain->SetBranchAddress("ele_RCTeta",&ele_RCTeta);
  myChain->SetBranchAddress("ele_RCTphi",&ele_RCTphi);
  //
  myChain->SetBranchAddress("ele_pT",&ele_pT);
  myChain->SetBranchAddress("ele_eT",&ele_eT);
  myChain->SetBranchAddress("ele_E",&ele_E);
  //
  myChain->SetBranchAddress("ele_trigEG",&ele_trigEG);
  myChain->SetBranchAddress("ele_trigEG_M",&ele_trigEG_M);
  //
  myChain->SetBranchAddress("ele_severityLevelSeed", &ele_severityLevelSeed);
  myChain->SetBranchAddress("ele_he",&ele_he);
  myChain->SetBranchAddress("ele_sigmaietaieta",&ele_sigmaietaieta);
  myChain->SetBranchAddress("ele_hcalDepth1TowerSumEt_dr03", &ele_hcalDepth1TowerSumEt_dr03);
  myChain->SetBranchAddress("ele_hcalDepth2TowerSumEt_dr03", &ele_hcalDepth2TowerSumEt_dr03);
  myChain->SetBranchAddress("ele_ecalRecHitSumEt_dr03", &ele_ecalRecHitSumEt_dr03);
  myChain->SetBranchAddress("ele_tkSumPt_dr03",&ele_tkSumPt_dr03);
  myChain->SetBranchAddress("ele_expected_inner_hits",&ele_expected_inner_hits);
  myChain->SetBranchAddress("ele_expected_inner_hits_aod",&ele_expected_inner_hits_aod);
  myChain->SetBranchAddress("ele_deltaphiin",&ele_deltaphiin);
  myChain->SetBranchAddress("ele_deltaetain",&ele_deltaetain);
  myChain->SetBranchAddress("ele_fbrem",&ele_fbrem);

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
  myChain->SetBranchAddress("ele_RCTL1iso_To", &ele_RCTL1iso_To);
  myChain->SetBranchAddress("ele_RCTL1noniso_To", &ele_RCTL1noniso_To);

  myChain->SetBranchAddress("ele_RCTetaVect", &ele_RCTetaVect);
  myChain->SetBranchAddress("ele_RCTphiVect", &ele_RCTphiVect);
  myChain->SetBranchAddress("ele_RCTetVect", &ele_RCTetVect);
  myChain->SetBranchAddress("ele_RCTL1isoVect", &ele_RCTL1isoVect);
  myChain->SetBranchAddress("ele_RCTL1nonisoVect", &ele_RCTL1nonisoVect);
  myChain->SetBranchAddress("ele_RCTL1isoVect_M", &ele_RCTL1isoVect_M);
  myChain->SetBranchAddress("ele_RCTL1nonisoVect_M", &ele_RCTL1nonisoVect_M);
  myChain->SetBranchAddress("ele_RCTL1isoVect_To", &ele_RCTL1isoVect_To);
  myChain->SetBranchAddress("ele_RCTL1nonisoVect_To", &ele_RCTL1nonisoVect_To);

  // L1 candidates
  myChain->SetBranchAddress("trig_L1emIso_N", &trig_L1emIso_N);
  myChain->SetBranchAddress("trig_L1emIso_ieta", &trig_L1emIso_ieta);
  myChain->SetBranchAddress("trig_L1emIso_iphi", &trig_L1emIso_iphi);
  myChain->SetBranchAddress("trig_L1emIso_rank", &trig_L1emIso_rank);

  myChain->SetBranchAddress("trig_L1emNonIso_N", &trig_L1emNonIso_N);
  myChain->SetBranchAddress("trig_L1emNonIso_ieta", &trig_L1emNonIso_ieta);
  myChain->SetBranchAddress("trig_L1emNonIso_iphi", &trig_L1emNonIso_iphi); 
  myChain->SetBranchAddress("trig_L1emNonIso_rank", &trig_L1emNonIso_rank);

  myChain->SetBranchAddress("trig_L1emIso_N_To", &trig_L1emIso_N_To);
  myChain->SetBranchAddress("trig_L1emIso_ieta_To", &trig_L1emIso_ieta_To);
  myChain->SetBranchAddress("trig_L1emIso_iphi_To", &trig_L1emIso_iphi_To);
  myChain->SetBranchAddress("trig_L1emIso_rank_To", &trig_L1emIso_rank_To);

  myChain->SetBranchAddress("trig_L1emNonIso_N_To", &trig_L1emNonIso_N_To);
  myChain->SetBranchAddress("trig_L1emNonIso_ieta_To", &trig_L1emNonIso_ieta_To);
  myChain->SetBranchAddress("trig_L1emNonIso_iphi_To", &trig_L1emNonIso_iphi_To); 
  myChain->SetBranchAddress("trig_L1emNonIso_rank_To", &trig_L1emNonIso_rank_To);

  myChain->SetBranchAddress("trig_L1emIso_N_M_To", &trig_L1emIso_N_M_To);
  myChain->SetBranchAddress("trig_L1emIso_ieta_M_To", &trig_L1emIso_ieta_M_To);
  myChain->SetBranchAddress("trig_L1emIso_iphi_M_To", &trig_L1emIso_iphi_M_To);
  myChain->SetBranchAddress("trig_L1emIso_rank_M_To", &trig_L1emIso_rank_M_To);

  myChain->SetBranchAddress("trig_L1emNonIso_N_M_To", &trig_L1emNonIso_N_M_To);
  myChain->SetBranchAddress("trig_L1emNonIso_ieta_M_To", &trig_L1emNonIso_ieta_M_To);
  myChain->SetBranchAddress("trig_L1emNonIso_iphi_M_To", &trig_L1emNonIso_iphi_M_To); 
  myChain->SetBranchAddress("trig_L1emNonIso_rank_M_To", &trig_L1emNonIso_rank_M_To);

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

  // normal collection
  myChain->SetBranchAddress("trig_tower_N", &trig_tower_N);
  myChain->SetBranchAddress("trig_tower_ieta",  &trig_tower_ieta);
  myChain->SetBranchAddress("trig_tower_iphi",  &trig_tower_iphi);
  myChain->SetBranchAddress("trig_tower_adc",  &trig_tower_adc);
  myChain->SetBranchAddress("trig_tower_sFGVB",  &trig_tower_sFGVB);  
  myChain->SetBranchAddress("trig_tower_FG",  &trig_tower_FG);  

  // HCAL Towers
  myChain->SetBranchAddress("trig_tower_hcal_N", &trig_tower_hcal_N);
  myChain->SetBranchAddress("trig_tower_hcal_ieta",  &trig_tower_hcal_ieta);
  myChain->SetBranchAddress("trig_tower_hcal_iphi",  &trig_tower_hcal_iphi);
  myChain->SetBranchAddress("trig_tower_hcal_et",  &trig_tower_hcal_et);
  myChain->SetBranchAddress("trig_tower_hcal_FG",  &trig_tower_hcal_FG);  

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


  ofstream outcheckL1("logs/checkL1_XavEE_short.log",ios::out);

  // Map TT
  MAPTT adcTT;
  MAPTT::iterator iterTT;
  pair<int,int> coords;
  int iTow, iTT, iTTH, adc, het, sFGVB, FG;
  double hOverE;

  // Counters
  const int nECAL=2;
  int nEle[nECAL]={0,0};
  int iECAL;
  double ele_TT_Reta, ele_TT_Rphi;
  VectMask masked;
  bool present;

  // Counters of cases //
  int nEleAbo15GeV[nECAL]={0,0};
  int nEleTrig[nECAL]={0,0};
  int nEleNoTrig[nECAL]={0,0};
  int nEleIneff[nECAL][10]={ {0,0,0,0,0,0,0,0,0,0} , {0,0,0,0,0,0,0,0,0,0} };
  // Pre,Post firing ; FG ; H/E ; Mask:Tow,Strip,Xtal ; Spike ; Etalement

  // Variables
  int iTP_max, iTP_max2, rank_TP_max, rank_TP_max2, iTow_max, iTow_max2, iTow_max3;
  double eT_Tow_max, eT_Tow_max2, eT_Tow_max3, ele_eT_dep_in2TPmax;
  vector<int> firedEG[2]; // [N/M]
  vector<int> firedEG_pre, firedEG_post;

  // Get Issues
  vector< pair<int,int> > stripmaskIssueV, towmaskIssueV, xtalmaskIssueV, spikeIssueV, FG_issueV, HE_issueV;
  bool maskIssue, towmaskIssue, stripmaskIssue, xtalmaskIssue, spikeIssue, FG_issue, HE_issue, prefire, postfire;

  // Number of entries to process
  int nEntDisp = myChain->GetEntries();
  int nProcess = nEntDisp;
  if(nEntAsk>=0 && nEntAsk<nEntDisp) nProcess=nEntAsk;

  int nLooked=0;

  for(int iEntry=0 ; iEntry<nProcess ; iEntry++) {

    if(iEntry%100==0) cout << "process entry : " << iEntry << "/" << nProcess << endl;

    myChain->GetEntry(iEntry);

    // Select ONLY EE+3 (17;9) and EE-4 (4;12) regions to investigate ////////
    int inv_ieta[2]={17,4};
    int inv_iphi[2]={9,12};
    bool skip=true;
    for(int inv=0 ; inv<2 ; inv++)
      if(ele_RCTeta == inv_ieta[inv] && ele_RCTphi == inv_iphi[inv]) skip=false;
    if(skip) continue;	
    //////////////////////////////////////////////////////////////////////////

    if(ele_eta<1.479) iECAL=0;
    else iECAL=1;

    nEle[iECAL]++ ;

    // Check electron eT
    if(ele_eT<15) continue;
    nEleAbo15GeV[iECAL]++ ;

    for(int iColl=0 ; iColl<2 ; iColl++) {
      firedEG[iColl].clear();
      firedEG[iColl].resize(nEG,0);	  
    }
    //
    firedEG_pre.clear();
    firedEG_pre.resize(nEG,0);
    firedEG_post.clear();
    firedEG_post.resize(nEG,0);

    for(int iEG=0 ; iEG<nEG ; iEG++)
      ele_trigEG[iEG] = ele_trigEG_M[iEG] = ele_trigEG_pre[iEG] = ele_trigEG_post[iEG] = 0;
    
    for(int iR=0 ; iR < 10 ; iR++) {
      
      globalFireL1_Normal( ele_RCTetVect[iR] , ele_RCTL1nonisoVect[iR],
			   ele_RCTL1isoVect[iR], firedEG[0], menu); 
      
      globalFireL1_Normal( ele_RCTetVect[iR] , ele_RCTL1nonisoVect_M[iR],
			   ele_RCTL1isoVect_M[iR], firedEG[1], menu); 
      
      if(debug) cout << ele_RCTetVect[iR] << " " << ele_RCTL1nonisoVect[iR] << " " << ele_RCTL1isoVect[iR] << endl;
      
    }
    
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      ele_trigEG[iEG] = firedEG[0][iEG];
      ele_trigEG_M[iEG] = firedEG[1][iEG];
    }	  

    // Check electron triggering
    if(ele_trigEG[15]!=0 && ele_trigEG[15]!=1) continue;
    if(ele_trigEG[15]==1) { 
      nEleTrig[iECAL]++ ;
      continue; 
    }
    nEleNoTrig[iECAL]++ ;

    nLooked++;
    if(nLooked>2) break;

    // Map the trigger tower //
    adcTT.clear();
    for(int t=0 ; t<trig_tower_N ; t++) {
      coords = make_pair( trig_tower_ieta[t] , trig_tower_iphi[t] );
      adcTT[coords].first = t;
    }
    for(int t=0 ; t<trig_tower_hcal_N ; t++) {
      if(t>=nHtow) {
	cout << "FATAL :: PROBLEM WITH trig_tower_hcal_N=" << trig_tower_hcal_N << endl;
	return 3;
      }
      coords = make_pair( trig_tower_hcal_ieta[t] , trig_tower_hcal_iphi[t] );
      iterTT = adcTT.find( coords );

      if( iterTT != adcTT.end() ) {
	adcTT[coords].second = t;
      }
    }

    // General Informations //
    outlog << "Run #" << nRun << "  Lumi #" << nLumi << "  Event #" << nEvent << endl
	   << "Ele : eta=" << ele_eta << "  phi=" << ele_phi << "  eT=" << ele_eT << endl
	   << "RCT : (" << ele_RCTeta << " ; " << ele_RCTphi << ") : iso=" << ele_RCTL1iso << " , noniso=" << ele_RCTL1noniso << endl
	   << "RCT Vect : " ;

    // L1 candidates & RCT regions //
    present=false;
    for(int iR=0 ; iR<10 ; iR++) {
      if( ele_RCTL1isoVect[iR]>0 || ele_RCTL1nonisoVect[iR]>0 ) {
	present=true;
	outlog << "(" << ele_RCTetaVect[iR] << " ; " << ele_RCTphiVect[iR] 
	       << ") : iso=" << ele_RCTL1isoVect[iR] 
	       << " , noniso=" << ele_RCTL1nonisoVect[iR] 
	       << " , eT=" << ele_RCTetVect[iR] << endl;
      }
    }
    if(!present) outlog << endl;


    // Pre/Post Firing //
    matchL1toEle( trig_preL1emIso_N, trig_preL1emIso_ieta, trig_preL1emIso_iphi, trig_preL1emIso_rank,
		  ele_RCTeta,  ele_RCTphi,  ele_RCTL1iso_pre,
		  ele_RCTetaVect, ele_RCTphiVect, ele_RCTL1isoVect_pre, 10, outcheckL1);
    
    matchL1toEle( trig_preL1emNonIso_N, trig_preL1emNonIso_ieta, trig_preL1emNonIso_iphi, trig_preL1emNonIso_rank,
		  ele_RCTeta,  ele_RCTphi,  ele_RCTL1noniso_pre,
		  ele_RCTetaVect, ele_RCTphiVect, ele_RCTL1nonisoVect_pre, 10, outcheckL1);

    matchL1toEle( trig_postL1emIso_N, trig_postL1emIso_ieta, trig_postL1emIso_iphi, trig_postL1emIso_rank,
		  ele_RCTeta,  ele_RCTphi,  ele_RCTL1iso_post,
		  ele_RCTetaVect, ele_RCTphiVect, ele_RCTL1isoVect_post, 10, outcheckL1);

    matchL1toEle( trig_postL1emNonIso_N, trig_postL1emNonIso_ieta, trig_postL1emNonIso_iphi, trig_postL1emNonIso_rank,
		  ele_RCTeta,  ele_RCTphi,  ele_RCTL1noniso_post,
		  ele_RCTetaVect, ele_RCTphiVect, ele_RCTL1nonisoVect_post, 10, outcheckL1);

    for(int iR=0 ; iR < 10 ; iR++) {
      
      globalFireL1_Normal( ele_RCTetVect[iR] , ele_RCTL1nonisoVect_pre[iR],
			   ele_RCTL1isoVect_pre[iR], firedEG_pre, menu); 
      
      globalFireL1_Normal( ele_RCTetVect[iR] , ele_RCTL1nonisoVect_post[iR],
			   ele_RCTL1isoVect_post[iR], firedEG_post, menu); 
    }

    outlog << "Pre-firing : iso=" << ele_RCTL1iso_pre << " noniso=" << ele_RCTL1noniso_pre << endl
	   << "Vect : " ;

    present=false;
    for(int iR=0 ; iR<10 ; iR++) {
      if( ele_RCTL1isoVect_pre[iR]>0 || ele_RCTL1nonisoVect_pre[iR]>0 ) {
	present=true;
	outlog << "(" << ele_RCTetaVect[iR] << " ; " << ele_RCTphiVect[iR] 
	       << ") : iso=" << ele_RCTL1isoVect_pre[iR] 
	       << " , noniso=" << ele_RCTL1nonisoVect_pre[iR] << endl;
      }
    }
    if(!present) outlog << endl;
    
    outlog << "Post-firing : iso=" << ele_RCTL1iso_post << " noniso=" << ele_RCTL1noniso_post << endl
	   << "Vect : " ;

    present=false;
    for(int iR=0 ; iR<10 ; iR++) {
      if( ele_RCTL1isoVect_post[iR]>0 || ele_RCTL1nonisoVect_post[iR]>0 ) {
	present=true;
	outlog << "(" << ele_RCTetaVect[iR] << " ; " << ele_RCTphiVect[iR] 
	       << ") : iso=" << ele_RCTL1isoVect_post[iR] 
	       << " , noniso=" << ele_RCTL1nonisoVect_post[iR] << endl;
      }
    }
    if(!present) outlog << endl;

    // Pre/Post firing //
    prefire=postfire=false;
    if(firedEG_pre[15]==1) {
      outlog << "Pre-Firing" << endl;
      prefire=true;
    }
    if(firedEG_post[15]==1) {
      outlog << "Post-Firing" << endl;
      postfire=true;
    }
    
    // Towers //
    outlog << "TT Vect : " << endl;

    iTP_max=-999; rank_TP_max=0;

    iTow_max=-999; iTow_max2=-999; iTow_max3=-999;
    eT_Tow_max=0;  eT_Tow_max2=0;  eT_Tow_max3=0;

    towmaskIssueV.clear();
    stripmaskIssueV.clear();
    xtalmaskIssueV.clear();
    spikeIssueV.clear();
    FG_issueV.clear();
    HE_issueV.clear();

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Loop over electron's towers ////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //
    int ele_TPrank_dep_all = 0;
    //
    for(int iT=0 ; iT<50 ; iT++) {
      if( ele_TTetVect[iT]>0 ) {

	ele_TT_Reta = getGCTRegionEta(ele_TTetaVect[iT]);
	ele_TT_Rphi = getGCTRegionPhi(ele_TTphiVect[iT]);
	
	masked = checkMasking( trig_nMaskedCh,     trig_iMaskedTTeta,       trig_iMaskedTTphi,
			       trig_strip_mask_N,  trig_strip_mask_TTieta,  trig_strip_mask_TTiphi,
			       trig_xtal_mask_N,   trig_xtal_mask_TTieta,   trig_xtal_mask_TTiphi,
			       ele_TTetaVect[iT],  ele_TTphiVect[iT]);

	// tower status
	iTT=-1; adc=0; sFGVB=0;
	iTTH=-1; het=0; hOverE=0;

	coords = make_pair( ele_TTetaVect[iT] , ele_TTphiVect[iT] );	
	iTT = adcTT[coords].first;
	iTTH = adcTT[coords].second;

	if(iTT>=0) {
	  adc = trig_tower_adc[iTT];
	  ele_TPrank_dep_all += adc;
	  sFGVB = trig_tower_sFGVB[iTT];
	  FG = trig_tower_FG[iTT];
	}
	else {
	  adc= sFGVB = FG = -777;
	}

	if(iTTH>=0 && iTTH<nHtow) {
	  het = trig_tower_hcal_et[iTTH];
	  if(adc!=0) hOverE = ((double)het*4) / ((double)(0.5*adc)) ;
	  else hOverE=-999;
	}
	else {
	  het=hOverE=-777;
	}

	outlog << "(" << ele_TT_Reta << ";" << ele_TT_Rphi << "):"
	       << "(" << ele_TTetaVect[iT] << ";" << ele_TTphiVect[iT] << ")=" << ele_TTetVect[iT] 
	       << " adc=" << adc << " sFGVB=" << sFGVB << " FG=" << FG << " het=" << het << " hOverE=" << hOverE ;

	// Fine Grain veto
	if(FG==1) {
	  outlog << " | FG-veto" ;
	  FG_issueV.push_back( make_pair(iT,iTT) );
	  FG_issue=true;
	}

	// H/E veto
	if(hOverE>0.05) {
	  outlog << " | H/E-veto" ;
	  HE_issueV.push_back( make_pair(iT,iTT) );
	  HE_issue=true;
	}

	// Masking
	maskIssue=false;

	// tower masking
	if( masked[0].size()>0 ) {
	  outlog << " | TowMask ";
	  towmaskIssueV.push_back( make_pair(iT,iTT) );
	  maskIssue=true;
	}

	// strip masking
	if( masked[1].size()>0 ) {
	  outlog << " | StripMask : ";
	  for(int iS=0 ; iS<(int)(masked[1].size()) ; iS++) {
	    outlog << "id=" << trig_strip_mask_StripID[iS] << " ; "
		   << "pseudo=" << trig_strip_mask_PseudoStripID[iS] << " ; "
		   << "TT(" << trig_strip_mask_TTieta[iS] << "," << trig_strip_mask_TTiphi[iS] << ") ; "
		   << "TCC(" << trig_strip_mask_TccID[iS] << ") ; "
		   << "CCU(" << trig_strip_mask_CCU[iS] << ") ; "
		   << "status(" << trig_strip_mask_status[iS] << ") ; ";
	  }		   
	  maskIssue=true;
	  stripmaskIssueV.push_back( make_pair(iT,iTT) );
	}

	// xtal masking
	if( masked[2].size()>0 ) {
	  outlog << " | XtalMask : ";
	  for(int iX=0 ; iX<(int)(masked[2].size()) ; iX++)
	    outlog << "(" <<  trig_xtal_mask_ieta << ";" << trig_xtal_mask_ieta << ") ; ";
	  maskIssue=true;
	  xtalmaskIssueV.push_back( make_pair(iT,iTT) );	  
	}	

	// spike zeroing
	if( sFGVB==0 && ele_TTetVect[iT]>8 ) {
	  outlog << " | Zeroed ";
	  spikeIssueV.push_back( make_pair(iT,iTT) );
	  spikeIssue=true;
	}

	outlog << endl;

	// check maxima //
	if(adc>rank_TP_max) { iTP_max=iTT; rank_TP_max=adc; }
	//
	if(ele_TTetVect[iT]>eT_Tow_max) { iTow_max=iT; eT_Tow_max=ele_TTetVect[iT]; }
	//
	if(ele_TTetVect[iT]<eT_Tow_max && ele_TTetVect[iT]>eT_Tow_max2) { iTow_max2=iT; eT_Tow_max2=ele_TTetVect[iT]; }
	//
	if(ele_TTetVect[iT]<eT_Tow_max && ele_TTetVect[iT]<eT_Tow_max2 && ele_TTetVect[iT]>eT_Tow_max3) { iTow_max3=iT; eT_Tow_max3=ele_TTetVect[iT]; }
	//
      }// endif TT > 0
    } // end loop over TTs
    outlog << endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Find max neighbour of max TP //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    pair<int,int> op[4] = {make_pair(1,0),make_pair(-1,0),make_pair(0,1),make_pair(0,-1)};
    //
    iTT=iTP_max2=-999; rank_TP_max2=0;
    //
    if(iTP_max>=0) {

      for(int iOp=0 ; iOp<4 ; iOp++) {

	coords = make_pair( trig_tower_ieta[iTP_max]+op[iOp].first , trig_tower_iphi[iTP_max]+op[iOp].second );	
	iTT = adcTT[coords].first;

	if(iTT>=0)
	  if(trig_tower_adc[iTT]>rank_TP_max2) { iTP_max2=iTT; rank_TP_max2=trig_tower_adc[iTT];}
	//else outlog << "iTT=" << iTT ;
      }
    }
    if(iTP_max>=0)
      outlog << "MaxTP(" << iTP_max << "):(" 
	     << trig_tower_ieta[iTP_max] << "," << trig_tower_iphi[iTP_max] << ")=" << trig_tower_adc[iTP_max] ;
    else outlog << "NoMaxTPFound" ;
    //
    if(iTP_max2>=0)
      outlog << "   MaxNeighbour("<< iTP_max2 << "):(" 
	     << trig_tower_ieta[iTP_max2] << "," << trig_tower_iphi[iTP_max2] << ")=" << trig_tower_adc[iTP_max2] ;
    else outlog << "   NoNeighbourFound" ;
    outlog << endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Check FG, H/E and Masking : concerns 2 max TP ? ///////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // FG //
    FG_issue = false;
    if(fabs(ele_eta)<1.479) {
      for(u_int i=0 ; i<FG_issueV.size() ; i++) {
      
	iTow  = FG_issueV[i].first ;
	iTT = FG_issueV[i].second ;

	if(iTT>=0 && iTP_max>=0)
	  if(iTP_max == iTT) FG_issue = true;
      }
    }
    //
    // H/E //
    HE_issue = false;
    for(u_int i=0 ; i<HE_issueV.size() ; i++) {
      
      iTow  = HE_issueV[i].first ;
      iTT = HE_issueV[i].second ;

      if(iTT>=0 && iTP_max>=0)
	if(iTP_max == iTT) HE_issue = true;
      
    }
    //
    // Masking : check on electron's max eT towers //
    // Tower
    towmaskIssue = false;
    for(u_int i=0 ; i<towmaskIssueV.size() ; i++) {

      iTow = towmaskIssueV[i].first ;
      iTT = towmaskIssueV[i].second ;

      if(iTow>=0 && iTow_max>=0)
	if(iTow == iTow_max) towmaskIssue = true;
      
      if(iTow>=0 && iTow_max2>=0) {
	if(iTow == iTow_max2) {
	  if( ele_TTetVect[iTow_max]<15 && (ele_TTetVect[iTow_max]+ele_TTetVect[iTow_max2])>15 )
	    towmaskIssue = true;
	  }
      }
    }
    //
    // Strip
    stripmaskIssue = false;
    for(u_int i=0 ; i<stripmaskIssueV.size() ; i++) {

      iTow = stripmaskIssueV[i].first ;
      iTT = stripmaskIssueV[i].second ;
      
      if(iTow>=0 && iTow_max>=0)
	if(iTow == iTow_max) stripmaskIssue = true;
      
      if(iTow>=0 && iTow_max2>=0) {
	if(iTow == iTow_max2) {
	  if( ele_TTetVect[iTow_max]<15 && (ele_TTetVect[iTow_max]+ele_TTetVect[iTow_max2])>15 )
	    stripmaskIssue = true;
	}
      }
    }
    //
    // Xtal
    xtalmaskIssue = false;
    for(u_int i=0 ; i<xtalmaskIssueV.size() ; i++) {

      iTow = xtalmaskIssueV[i].first ;
      iTT = xtalmaskIssueV[i].second ;

      if(iTow>=0 && iTow_max>=0)
	if(iTow == iTow_max) xtalmaskIssue = true;
      
      if(iTow>=0 && iTow_max2>=0) {
	if(iTow == iTow_max2) {
	  if( ele_TTetVect[iTow_max]<15 && (ele_TTetVect[iTow_max]+ele_TTetVect[iTow_max2])>15 )
	    xtalmaskIssue = true;
	  }
      }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    // Spikes //
    //    
    spikeIssue = false;
    if(ele_eta<1.479) {
      for(u_int i=0 ; i<spikeIssueV.size() ; i++) {
	
	iTow = spikeIssueV[i].first ;
	iTT = spikeIssueV[i].second ;
	
	bool quit=false;
	for(u_int j=0 ; j<towmaskIssueV.size() ; j++) {
	  int iT_M = towmaskIssueV[j].first ;
	  int iTT_M = towmaskIssueV[j].second ;
	  
	  if(iT_M>=0)
	    if(iT_M==iTow) quit=true;
	}
	if(quit==true) continue;
	
	if(iTow>=0 && iTow_max>=0)
	  if(iTow == iTow_max) spikeIssue = true;
	
	if(iTow>=0 && iTow_max2>=0) {
	  if(iTow == iTow_max2) {
	    if( ele_TTetVect[iTow_max]<15 && (ele_TTetVect[iTow_max]+ele_TTetVect[iTow_max2])>15 )
	      spikeIssue = true;
	  }
	}
      }
    }
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Max electron's trigger tower deposits /////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(iTow_max>=0)
      outlog << "MaxTow(" << iTow_max << "):(" 
	     << ele_TTetaVect[iTow_max] << "," << ele_TTphiVect[iTow_max] << ")=" << ele_TTetVect[iTow_max] ;
    else outlog << "NoMaxTowFound" ;
    //
    if(iTow_max2>=0)
      outlog << "   MaxTow2(" << iTow_max2 << "):(" 
	     << ele_TTetaVect[iTow_max2] << "," << ele_TTphiVect[iTow_max2] << ")=" << ele_TTetVect[iTow_max2] ;
    else outlog << "   NoMaxTow2Found" ;
    //
    if(iTow_max3>=0)
      outlog << "   MaxTow3(" << iTow_max3 << "):(" 
	     << ele_TTetaVect[iTow_max3] << "," << ele_TTphiVect[iTow_max3] << ")=" << ele_TTetVect[iTow_max3] ;
    else outlog << "   NoMaxTow3Found" ;
    //
    outlog << endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Check Energy deposited by ele /////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // in the 2 max TP //
    ele_eT_dep_in2TPmax = 0;
    for(int iT=0 ; iT<50 ; iT++) {
      if( ele_TTetVect[iT]<=0 ) continue;
	
      if( ele_TTetaVect[iT]==trig_tower_ieta[iTP_max] && ele_TTphiVect[iT]==trig_tower_iphi[iTP_max] ) {
	//outlog << "adding ele_TTetaVect[iT]=" << ele_TTetVect[iT] << endl;
	ele_eT_dep_in2TPmax += ele_TTetVect[iT];
      }
      else if( ele_TTetaVect[iT]==trig_tower_ieta[iTP_max2] && ele_TTphiVect[iT]==trig_tower_iphi[iTP_max2] )
	ele_eT_dep_in2TPmax += ele_TTetVect[iT];
    }
    outlog << "Ele eT deposit in L1 candidate's 2 towers : " << ele_eT_dep_in2TPmax << " GeV" << endl;
    //
    // in the TP corresponding to 3 max Tower eT deposits //
    //
    int ele_TPrank_dep_in3TowMax = 0, ele_TPrank_dep_in2TowMax = 0;
    bool spread=false, spreadAll=false;
    //
    if(iTow_max>=0) {
      coords = make_pair( ele_TTetaVect[iTow_max] , ele_TTphiVect[iTow_max] );
      iTT = adcTT[coords].first;
      if(iTT>=0) {
	ele_TPrank_dep_in3TowMax += trig_tower_adc[iTT];
	ele_TPrank_dep_in2TowMax += trig_tower_adc[iTT];
      }
    }
    //
    if(iTow_max2>=0) {
      coords = make_pair( ele_TTetaVect[iTow_max2] , ele_TTphiVect[iTow_max2] );
      iTT = adcTT[coords].first;
      if(iTT>=0) {
	ele_TPrank_dep_in3TowMax += trig_tower_adc[iTT];
	ele_TPrank_dep_in2TowMax += trig_tower_adc[iTT];
      }
    }
    //
    if(iTow_max3>=0) {
      coords = make_pair( ele_TTetaVect[iTow_max3] , ele_TTphiVect[iTow_max3] );
      iTT = adcTT[coords].first;
      if(iTT>=0) ele_TPrank_dep_in3TowMax += trig_tower_adc[iTT];
    }
    //
    if( 0.5*ele_TPrank_dep_in3TowMax >= 15 ) {
      spread=true;
      outlog << "Spread : " << 0.5*ele_TPrank_dep_in3TowMax << endl;
    }
    if( 0.5*ele_TPrank_dep_in2TowMax < 15 && 0.5*ele_TPrank_dep_all>15 ) {
      spreadAll=true;
      outlog << "SpreadAll : " << 0.5*ele_TPrank_dep_all << endl;
    }
    outlog << endl << endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ///////////////////////////////////////
    // COUNTING ///////////////////////////
    ///////////////////////////////////////
    if(towmaskIssue || stripmaskIssue || xtalmaskIssue) {
      if(towmaskIssue)   nEleIneff[iECAL][4]++ ;  //
      if(stripmaskIssue) nEleIneff[iECAL][5]++ ;  //
      if(xtalmaskIssue)  nEleIneff[iECAL][6]++ ;  //
    }
    else if(spikeIssue)  nEleIneff[iECAL][7]++ ;
    else {
      if(prefire)        nEleIneff[iECAL][0]++ ;        //
      if(postfire)       nEleIneff[iECAL][1]++ ;        //
      if(FG_issue)       nEleIneff[iECAL][2]++ ;        //
      if(HE_issue)       nEleIneff[iECAL][3]++ ;        //
      if(spikeIssue)     nEleIneff[iECAL][7]++ ;  //
      if(spread)         nEleIneff[iECAL][8]++ ;  //
      if(spreadAll)      nEleIneff[iECAL][9]++ ;  //
    }
    ///////////////////////////////////////

  }

  TString ecalname[nECAL]={"BARREL","ENDCAPS"};
  TString countname[10] = {"pre-firing","post-firing", "FG", "H/E", "Mask:Tow", "Mask:Strip","Mask:Xtal", "Spike", "Etalement", "Etalement:Tout"};

  ofstream outreport("logs/log_report_XavEE_short.log",ios::out);
  for(iECAL=0 ; iECAL<2 ; iECAL++) {

    outreport << ecalname[iECAL] << endl
	      << "nEle=" << nEle[iECAL] << "  nEleAbo15GeV=" << nEleAbo15GeV[iECAL] << "  nEleTrig=" << nEleTrig[iECAL] 
	      << "  nEleNoTrig=" << nEleNoTrig[iECAL] << endl;

    for(int i=0 ; i<10 ; i++)
      outreport << countname[i] << "=" << nEleIneff[iECAL][i] << endl;
    outreport << endl << endl;
  }
  
  return 0;
}

int inefficiency(int nEntAsk=-1, TString dir="/data_CMS/cms/ndaci/ndaci_2011A/Inefficiency/Both/", TString nameChain="Electrons", double cutet=15., bool debug=false)
{

  TChain* myChain = new TChain(nameChain);
  myChain->Add(dir+"/*.root");

  ofstream outlog("logs/log_inefficiency_XavEE_short.log",ios::out);
  int ant = analyzeTree( myChain , nEntAsk, cutet, debug, outlog);

  return ant;
}
