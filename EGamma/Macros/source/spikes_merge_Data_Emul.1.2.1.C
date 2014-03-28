#include "/home/llr/cms/ndaci/SKWork/macro/skEfficiency/analyzers/COMMON/baseFuncNad.1.5.h"

#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/CommissioningEmulReco/NewKill/highPU/treelist.h"
#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/HLTEG12/treelist_HighPU.h"

// AN 2010 Repro 2012                                                                                                                
/// data                                                                                                                             
//#include "/data_CMS/cms/ndaci/ndaci_2010B/Repro2012/Data/treelist.h"
#include "/data_CMS/cms/ndaci/ndaci_2010B/EGMonitorModifSpikes2/Skimmed/treelist.h"

/// trigger                                                                                                                          
#include "/data_CMS/cms/ndaci/ndaci_2010B/Repro2012/TrigOnly/TrigOnly_10adc/treelist.h"
#include "/data_CMS/cms/ndaci/ndaci_2010B/Repro2012/TrigOnly/TrigOnly_15adc/treelist.h"
#include "/data_CMS/cms/ndaci/ndaci_2010B/Repro2012/TrigOnly/TrigOnly_17adc/treelist.h"
#include "/data_CMS/cms/ndaci/ndaci_2010B/Repro2012/TrigOnly/TrigOnly_19adc/treelist.h"
#include "/data_CMS/cms/ndaci/ndaci_2010B/Repro2012/TrigOnly/TrigOnly_21adc/treelist.h"
#include "/data_CMS/cms/ndaci/ndaci_2010B/Repro2012/TrigOnly/TrigOnly_23adc/treelist.h"
#include "/data_CMS/cms/ndaci/ndaci_2010B/Repro2012/TrigOnly/TrigOnly_30adc/treelist.h"
#include "/data_CMS/cms/ndaci/ndaci_2010B/Repro2012/TrigOnly/TrigOnly_35adc/treelist.h"
#include "/data_CMS/cms/ndaci/ndaci_2010B/Repro2012/TrigOnly/TrigOnly_40adc/treelist.h"
#include "/data_CMS/cms/ndaci/ndaci_2010B/Repro2012/TrigOnly/TrigOnly_50adc/treelist.h"

int matchL1toSpike(int trig_L1emIso_N, int* trig_L1emIso_ieta, int* trig_L1emIso_iphi, int* trig_L1emIso_rank,
		   int spike_RCTeta, int spike_RCTphi, int &spike_RCTL1iso, int &maxEGtriged, ofstream& outcheckL1) {

  int nC = trig_L1emIso_N;
  if( nC>4 || nC<0 ) nC=4;

  outcheckL1 << "SPIKE : (" << spike_RCTeta << ";" << spike_RCTphi << ")=" << spike_RCTL1iso << endl;

  for(int iC=0 ; iC<nC ; iC++) {
    if( trig_L1emIso_ieta[iC]==spike_RCTeta && trig_L1emIso_iphi[iC]==spike_RCTphi )
      spike_RCTL1iso = trig_L1emIso_rank[iC];
  }

  if(spike_RCTL1iso>maxEGtriged) maxEGtriged = spike_RCTL1iso;

  return 1;
}

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

int process(TChain *trigOnlyChain, MAPJson *jsonMap, MAPTrigOnly mapTo, TString nameChain, TString file, int iFile, int nEntries,
            TString dirIn, TString dirOut, TString usr_hlt, int iJson, 
	    bool GetTrigOnly, bool debug) {

  // TrigOnly MAP //                                                                                          
  MAPTrigOnly::iterator iterMapTo;
  pair<int,int> runevtTo;
  int iEntryTo=-1;

  // OUTPUT FILE //
  if(debug) cout << "--- define output file : " ;
  ostringstream ossi("");
  ossi << iFile ;
  //                                                                                                                                 
  TFile *outfile = new TFile(dirIn+"spikes_"+ossi.str()+".root","RECREATE");
  //ofstream outcheckdata(dirIn+"logcheckdata_"+ossi.str()+".txt",ios::out);
  ofstream outlog(dirOut+"/logs/log_"+ossi.str()+".txt",ios::out);
  ofstream outcheckL1(dirOut+"/checkL1/checkL1_"+ossi.str()+".txt",ios::out);

  //outlog << "TRY TRY TRY" << endl;
  //outcheckL1 << "TRY TRY TRY" << endl;
  //outcheckdata << "TRY TRY TRY" << endl;

  if(debug) cout << "spikes_"+ossi.str()+".root" << endl;
  ossi.str("");

  // INPUT TREE //                                                                                                                   
  if(debug) cout << "--- add data trees to myChain" << endl;
  TChain * myChain = new TChain(nameChain);
  myChain->Add(file);

  // VARIABLES //                                                                                      
  int nEvent, nRun, nLumi ; // data                                                      
  int nEvent_To, nRun_To, nLumi_To; // trigOnly
  // Vertices //
  int vtx_N;
  double _vtx_x[200], _vtx_y[200], _vtx_z[200];
  double _vtx_normalizedChi2[200], _vtx_ndof[200], _vtx_nTracks[200], _vtx_d0[200];

  // Trigger Paths //                                                                           
  int trig_hltInfo[250];
  //int _trig_isEleHLTpath;                                                                                    
  int trig_HLT_path[4]; // unbias, EG5, EG8, EG12
  char trig_fired_names[5000];
  vector<string> m_HLT_pathsV;
  vector<string> m_HLT_triggered;
  vector<int> m_HLT_pathsV_check;

  // Spikes
  int spike_N,spike_TTieta[5000], spike_TTiphi[5000], spike_Rieta[5000], spike_Riphi[5000], spike_severityLevel[5000], 
    spike_outOfTime[5000], spike_TT_sFGVB[5000], spike_TT_adc[5000], spike_TT_sFGVB_E[5000][5], spike_TT_adc_E[5000][5],
    spike_TT_t[5000], spike_TT_tE[5000];
  double  spike_Et[5000], spike_eta[5000], spike_phi[5000], spike_theta[5000];
  int spike_RCTL1iso[5000], spike_RCTL1noniso[5000], spike_RCTL1iso_To[5000], spike_RCTL1noniso_To[5000], 
    spike_RCTL1iso_M_To[5000], spike_RCTL1noniso_M_To[5000],
    spike_maxEGtrig[5000], spike_maxEGtrig_To[5000], spike_maxEGtrig_M_To[5000];
  
  int spike_out_TTieta, spike_out_TTiphi, spike_out_Rieta, spike_out_Riphi, spike_out_severityLevel, 
    spike_out_outOfTime, spike_out_TT_sFGVB, spike_out_TT_adc, spike_out_TT_sFGVB_E[5], spike_out_TT_adc_E[5],
    spike_out_TT_t, spike_out_TT_tE;
  double  spike_out_Et, spike_out_eta, spike_out_phi, spike_out_theta;
  int spike_out_RCTL1iso, spike_out_RCTL1noniso, spike_out_RCTL1iso_To, spike_out_RCTL1noniso_To, 
    spike_out_RCTL1iso_M_To, spike_out_RCTL1noniso_M_To,
    spike_out_maxEGtrig, spike_out_maxEGtrig_To, spike_out_maxEGtrig_M_To;
  
  // TP info
  const int nTow = 4032;
  int trig_tower_N,trig_tower_ieta[nTow],trig_tower_iphi[nTow],trig_tower_adc[nTow],trig_tower_sFGVB[nTow],trig_tower_FG[nTow];
  int trig_tower_N_M,trig_tower_ieta_M[nTow],trig_tower_iphi_M[nTow],trig_tower_adc_M[nTow],trig_tower_sFGVB_M[nTow],
    trig_tower_FG_M[nTow];
  int trig_tower_N_E,trig_tower_ieta_E[nTow],trig_tower_iphi_E[nTow],trig_tower_adc_E[nTow][5],trig_tower_sFGVB_E[nTow][5],
    trig_tower_FG_E[nTow][5];

  // HCAL TP //
  int trig_tower_hcal_N, trig_tower_hcal_ieta[4032], trig_tower_hcal_iphi[4032], trig_tower_hcal_FG[4032],trig_tower_hcal_et[4032];

  // L1 Candidates //                                                                                                             
  int trig_L1emIso_N, trig_L1emNonIso_N; 
  int trig_L1emIso_N_To, trig_L1emNonIso_N_To, trig_L1emIso_N_M_To, trig_L1emNonIso_N_M_To;

  // L1 candidates info //
  int trig_L1emIso_ieta[4], trig_L1emIso_iphi[4], trig_L1emIso_rank[4];
  int trig_L1emNonIso_ieta[4], trig_L1emNonIso_iphi[4], trig_L1emNonIso_rank[4];

  int trig_L1emIso_ieta_To[4], trig_L1emIso_iphi_To[4], trig_L1emIso_rank_To[4];
  int trig_L1emNonIso_ieta_To[4], trig_L1emNonIso_iphi_To[4], trig_L1emNonIso_rank_To[4];
  int trig_L1emIso_ieta_M_To[4], trig_L1emIso_iphi_M_To[4], trig_L1emIso_rank_M_To[4];
  int trig_L1emNonIso_ieta_M_To[4], trig_L1emNonIso_iphi_M_To[4], trig_L1emNonIso_rank_M_To[4];

  // INITIALIZATION //                                                                                                               
  //                                                                                                                                 
  // Global                                                                                                                        
  nEvent = 0;
  nRun = 0;
  nLumi = 0;
  //                                                                                                                          
  // Vertices                                                                                                                   
  vtx_N = 0;
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
  }

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

  // SPIKES //
  spike_N = 0;  
  for(int isp=0 ; isp<5000 ; isp++) {
    spike_outOfTime[isp] = -999;
    spike_severityLevel[isp] = -999 ;
    //spike_SwissCross[isp] = -999 ;
    spike_Et[isp] = -999 ;
    spike_phi[isp] = -999 ;
    spike_eta[isp] = -999 ;
    spike_theta[isp] = -999 ;
    spike_TTiphi[isp] = -999 ;
    spike_TTieta[isp] = -999 ;
    spike_TT_t[isp] = -999 ;
    spike_TT_tE[isp] = -999 ;
    spike_TT_sFGVB[isp] = -999 ;
    spike_TT_adc[isp] = -999 ;
    spike_Riphi[isp] = -999 ;
    spike_Rieta[isp] = -999 ;

    for(int i=0;i<5;i++) {
      spike_TT_sFGVB_E[isp][i] = -999 ;
      spike_TT_adc_E[isp][i]   = -999 ;
    }
  }

  spike_out_outOfTime = -999;
  spike_out_severityLevel = -999 ;
  //spike_out_SwissCross = -999 ;
  spike_out_Et = -999 ;
  spike_out_phi = -999 ;
  spike_out_eta = -999 ;
  spike_out_theta = -999 ;
  spike_out_TTiphi = -999 ;
  spike_out_TTieta = -999 ;
  spike_out_TT_t = -999 ;
  spike_out_TT_tE = -999 ;
  spike_out_TT_sFGVB = -999 ;
  spike_out_TT_adc = -999 ;
  spike_out_Riphi = -999 ;
  spike_out_Rieta = -999 ;
  
  for(int i=0;i<5;i++) {
    spike_out_TT_sFGVB_E[i] = -999 ;
    spike_out_TT_adc_E[i]   = -999 ;
  }

  myChain->SetBranchAddress("nEvent",&nEvent);
  myChain->SetBranchAddress("nRun",&nRun);
  myChain->SetBranchAddress("nLumi",&nLumi);

  // Vertices
  myChain->SetBranchAddress("vtx_N",&vtx_N);

  //myChain->SetBranchAddress("trig_HLT_path",&trig_HLT_path);
  myChain->SetBranchAddress("trig_fired_names",&trig_fired_names);
  myChain->SetBranchAddress("trig_hltInfo",&trig_hltInfo);

  // L1 candidates                             
  myChain->SetBranchAddress("trig_L1emIso_N", &trig_L1emIso_N);
  myChain->SetBranchAddress("trig_L1emIso_ieta", &trig_L1emIso_ieta);
  myChain->SetBranchAddress("trig_L1emIso_iphi", &trig_L1emIso_iphi);
  myChain->SetBranchAddress("trig_L1emIso_rank", &trig_L1emIso_rank);

  myChain->SetBranchAddress("trig_L1emNonIso_N", &trig_L1emNonIso_N);
  myChain->SetBranchAddress("trig_L1emNonIso_ieta", &trig_L1emNonIso_ieta);
  myChain->SetBranchAddress("trig_L1emNonIso_iphi", &trig_L1emNonIso_iphi);
  myChain->SetBranchAddress("trig_L1emNonIso_rank", &trig_L1emNonIso_rank);

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

  // TrigOnly Chain //
  trigOnlyChain->SetBranchAddress("nRun", &nRun_To);
  trigOnlyChain->SetBranchAddress("nLumi", &nLumi_To);
  trigOnlyChain->SetBranchAddress("nEvent", &nEvent_To);

  trigOnlyChain->SetBranchAddress("trig_HLT_path",&trig_HLT_path);

  trigOnlyChain->SetBranchAddress("trig_tower_N", &trig_tower_N);
  trigOnlyChain->SetBranchAddress("trig_tower_ieta",  &trig_tower_ieta);
  trigOnlyChain->SetBranchAddress("trig_tower_iphi",  &trig_tower_iphi);
  trigOnlyChain->SetBranchAddress("trig_tower_adc",  &trig_tower_adc);
  trigOnlyChain->SetBranchAddress("trig_tower_sFGVB",  &trig_tower_sFGVB);
  //trigOnlyChain->SetBranchAddress("trig_tower_FG",  &trig_tower_FG);

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

  TTree * outtree = new TTree("Spikes","Spikes");
  // General informations //
  outtree->Branch("nRun",&nRun,"nRun/I");
  outtree->Branch("nLumi",&nLumi,"nLumi/I");
  outtree->Branch("nEvent",&nEvent,"nEvent/I");

  // Vertices //
  outtree->Branch("vtx_N",&vtx_N,"vtx_N/I");

  outtree->Branch("trig_HLT_path",&trig_HLT_path,"trig_HLT_path[4]/I");
  outtree->Branch("trig_fired_names",&trig_fired_names,"trig_fired_names[5000]/C");
  outtree->Branch("trig_hltInfo",&trig_hltInfo,"trig_hltInfo[250]/I");

  // Trigger towers //
  outtree->Branch("trig_tower_N",&trig_tower_N,"trig_tower_N/I");
  outtree->Branch("trig_tower_ieta",&trig_tower_ieta,"trig_tower_ieta[4032]/I");
  outtree->Branch("trig_tower_iphi",&trig_tower_iphi,"trig_tower_iphi[4032]/I");
  outtree->Branch("trig_tower_adc",&trig_tower_adc,"trig_tower_adc[4032]/I");
  outtree->Branch("trig_tower_sFGVB",&trig_tower_sFGVB,"trig_tower_sFGVB[4032]/I");
  outtree->Branch("trig_tower_FG",&trig_tower_FG,"trig_tower_FG[4032]/I");
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

  // L1 candidates from Data
  outtree->Branch("trig_L1emIso_N",     &trig_L1emIso_N,     "trig_L1emIso_N/I");
  outtree->Branch("trig_L1emIso_ieta",  &trig_L1emIso_ieta,  "trig_L1emIso_ieta[4]/I");
  outtree->Branch("trig_L1emIso_iphi",  &trig_L1emIso_iphi,  "trig_L1emIso_iphi[4]/I");
  outtree->Branch("trig_L1emIso_rank",  &trig_L1emIso_rank,  "trig_L1emIso_rank[4]/I");

  outtree->Branch("trig_L1emNonIso_N",     &trig_L1emNonIso_N,     "trig_L1emNonIso_N/I");
  outtree->Branch("trig_L1emNonIso_ieta",  &trig_L1emNonIso_ieta,  "trig_L1emNonIso_ieta[4]/I");
  outtree->Branch("trig_L1emNonIso_iphi",  &trig_L1emNonIso_iphi,  "trig_L1emNonIso_iphi[4]/I");
  outtree->Branch("trig_L1emNonIso_rank",  &trig_L1emNonIso_rank,  "trig_L1emNonIso_rank[4]/I");
  
  // L1 candidates from TrigOnly 
  outtree->Branch("trig_L1emIso_N_To",     &trig_L1emIso_N_To,     "trig_L1emIso_N_To/I");
  outtree->Branch("trig_L1emIso_ieta_To",  &trig_L1emIso_ieta_To,  "trig_L1emIso_ieta_To[4]/I");
  outtree->Branch("trig_L1emIso_iphi_To",  &trig_L1emIso_iphi_To,  "trig_L1emIso_iphi_To[4]/I");
  outtree->Branch("trig_L1emIso_rank_To",  &trig_L1emIso_rank_To,  "trig_L1emIso_rank_To[4]/I");
  //                                                                                                   
  outtree->Branch("trig_L1emNonIso_N_To", &trig_L1emNonIso_N_To, "trig_L1emNonIso_N_To/I");
  outtree->Branch("trig_L1emNonIso_ieta_To", &trig_L1emNonIso_ieta_To, "trig_L1emNonIso_ieta_To[4]/I");
  outtree->Branch("trig_L1emNonIso_iphi_To", &trig_L1emNonIso_iphi_To, "trig_L1emNonIso_iphi_To[4]/I");
  outtree->Branch("trig_L1emNonIso_rank_To", &trig_L1emNonIso_rank_To, "trig_L1emNonIso_rank_To[4]/I");

  outtree->Branch("trig_L1emIso_N_M_To",     &trig_L1emIso_N_M_To,     "trig_L1emIso_N_M_To/I");
  outtree->Branch("trig_L1emIso_ieta_M_To",  &trig_L1emIso_ieta_M_To,  "trig_L1emIso_ieta_M_To[4]/I");
  outtree->Branch("trig_L1emIso_iphi_M_To",  &trig_L1emIso_iphi_M_To,  "trig_L1emIso_iphi_M_To[4]/I");
  outtree->Branch("trig_L1emIso_rank_M_To",  &trig_L1emIso_rank_M_To,  "trig_L1emIso_rank_M_To[4]/I");
  //                                                                                
  outtree->Branch("trig_L1emNonIso_N_M_To", &trig_L1emNonIso_N_M_To, "trig_L1emNonIso_N_M_To/I");
  outtree->Branch("trig_L1emNonIso_ieta_M_To", &trig_L1emNonIso_ieta_M_To, "trig_L1emNonIso_ieta_M_To[4]/I");
  outtree->Branch("trig_L1emNonIso_iphi_M_To", &trig_L1emNonIso_iphi_M_To, "trig_L1emNonIso_iphi_M_To[4]/I");
  outtree->Branch("trig_L1emNonIso_rank_M_To", &trig_L1emNonIso_rank_M_To, "trig_L1emNonIso_rank_M_To[4]/I");

  outtree->Branch("spike_TTieta",&spike_out_TTieta,"spike_TTieta/I");
  outtree->Branch("spike_TTiphi",&spike_out_TTiphi,"spike_TTiphi/I");
  //
  outtree->Branch("spike_TT_t",&spike_out_TT_t,"spike_TT_t/I");
  outtree->Branch("spike_TT_sFGVB",&spike_out_TT_sFGVB,"spike_TT_sFGVB/I");
  outtree->Branch("spike_TT_adc",&spike_out_TT_adc,"spike_TT_adc/I");
  outtree->Branch("spike_TT_tE",&spike_out_TT_tE,"spike_TT_tE/I");
  outtree->Branch("spike_TT_sFGVB_E",&spike_out_TT_sFGVB_E,"spike_TT_sFGVB[5]/I");
  outtree->Branch("spike_TT_adc_E",&spike_out_TT_adc_E,"spike_TT_adc[5]/I");
  //
  outtree->Branch("spike_Rieta",&spike_out_Rieta,"spike_Rieta/I");
  outtree->Branch("spike_Riphi",&spike_out_Riphi,"spike_Riphi/I");
  outtree->Branch("spike_severityLevel",&spike_out_severityLevel,"spike_severityLevel/I");
  outtree->Branch("spike_outOfTime",&spike_out_outOfTime,"spike_outOfTime/I");
  outtree->Branch("spike_Et",&spike_out_Et,"spike_Et/D");
  outtree->Branch("spike_eta",&spike_out_eta,"spike_eta/D");
  outtree->Branch("spike_phi",&spike_out_phi,"spike_phi/D");
  outtree->Branch("spike_theta",&spike_out_theta,"spike_theta/D");
  //
  outtree->Branch("spike_maxEGtrig",&spike_out_maxEGtrig,"spike_maxEGtrig/I");
  outtree->Branch("spike_maxEGtrig_To",&spike_out_maxEGtrig_To,"spike_maxEGtrig_To/I");
  outtree->Branch("spike_maxEGtrig_M_To",&spike_out_maxEGtrig_M_To,"spike_maxEGtrig_M_To/I");
  //
  outtree->Branch("spike_RCTL1iso",&spike_out_RCTL1iso,"spike_RCTL1iso/I");
  outtree->Branch("spike_RCTL1noniso",&spike_out_RCTL1noniso,"spike_RCTL1noniso/I");
  outtree->Branch("spike_RCTL1iso_To",&spike_out_RCTL1iso_To,"spike_RCTL1iso_To/I");
  outtree->Branch("spike_RCTL1noniso_To",&spike_out_RCTL1noniso_To,"spike_RCTL1noniso_To/I");
  outtree->Branch("spike_RCTL1iso_M_To",&spike_out_RCTL1iso_M_To,"spike_RCTL1iso_M_To/I");
  outtree->Branch("spike_RCTL1noniso_M_To",&spike_out_RCTL1noniso_M_To,"spike_RCTL1noniso_M_To/I");

  // Variables
  bool isGoodRun, evt_trig_EG12;
  TString filename;

  // Map TT
  MAPTT adcTT;
  MAPTT::iterator iterTT;
  pair<int,int> coords;
  int TT_t, TT_tE;

  int numEntries = myChain->GetEntries () ;
  int nProcess = numEntries;
  if(nEntries>=0 && nEntries<numEntries)
    nProcess = nEntries;

  int nCurrentRun = -999;
  outlog << "will process " << nProcess << "/" << numEntries << "entries" << endl;

  // Loop over entries //
  for (int iEvent = 0 ; iEvent < nProcess ; iEvent++ ) {
    
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

    for(int isp=0 ; isp<5000 ; isp++) {
      spike_outOfTime[isp] = -999;
      spike_severityLevel[isp] = -999 ;
      //spike_SwissCross[isp] = -999 ;
      spike_Et[isp] = -999 ;
      spike_phi[isp] = -999 ;
      spike_eta[isp] = -999 ;
      spike_theta[isp] = -999 ;
      spike_TTiphi[isp] = -999 ;
      spike_TTieta[isp] = -999 ;
      spike_TT_t[isp] = -999 ;
      spike_TT_tE[isp] = -999 ;
      spike_TT_sFGVB[isp] = -999 ;
      spike_TT_adc[isp] = -999 ;
      spike_Riphi[isp] = -999 ;
      spike_Rieta[isp] = -999 ;

      for(int i=0;i<5;i++) {
	spike_TT_sFGVB_E[isp][i] = -999 ;
	spike_TT_adc_E[isp][i]   = -999 ;
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

    if(debug) cout << "iJson = " << iJson << endl;
    if( iJson>-1 && iJson<5) {
      isGoodRun = AcceptEventByRunAndLumiSection(nRun, nLumi, jsonMap[iJson]);
      if(!isGoodRun) {
	outlog << "failed JSON" << endl;
	continue;
      }
    }

    // JSON selection ////////////////////////////////////////////
    if(debug) cout << "iJson = " << iJson << endl;
    if( iJson>-1 && iJson<5) {
      isGoodRun = false;
      isGoodRun = AcceptEventByRunAndLumiSection(nRun, nLumi, jsonMap[iJson]);
      if(!isGoodRun) {
	outlog << "failed JSON" << endl;
	continue;
      }
    }
    
    // HLT selection //
    if(usr_hlt=="unbias")
      if(trig_HLT_path[0]==0) continue;
    
    // Map the trigger tower //////////////////////////////////////////////////////////////////////
    adcTT.clear();
    for(int t=0 ; t<trig_tower_N ; t++) {
      coords = make_pair( trig_tower_ieta[t] , trig_tower_iphi[t] );
      adcTT[coords].first = t;
    }
    for(int t=0 ; t<trig_tower_N_E ; t++) {
      coords = make_pair( trig_tower_ieta_E[t] , trig_tower_iphi_E[t] );
      iterTT = adcTT.find( coords );

      if( iterTT != adcTT.end() ) {
	adcTT[coords].second = t;
      }
    }
    
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

    // MATCHING SPIKE / L1-CANDIDATE //////////////////////////////////////////////////////////////////                         
    evt_trig_EG12=false;

    outcheckL1 << "ISO_N_data : " ;
    for(int i=0 ; i<4 ; i++) {
      outcheckL1 << "(" << trig_L1emIso_ieta[i] << ";" << trig_L1emIso_iphi[i] << ")=" << trig_L1emIso_rank[i] << " | ";
      if(trig_L1emIso_rank[i] >= 12) evt_trig_EG12=true;
    }

    outcheckL1 << endl
	       << "NONISO_N_data : ";
    for(int i=0 ; i<4 ; i++) {
      outcheckL1 << "(" << trig_L1emNonIso_ieta[i] << ";" << trig_L1emNonIso_iphi[i] << ")=" << trig_L1emNonIso_rank[i] << " | ";
      if(trig_L1emNonIso_rank[i] >= 12) evt_trig_EG12=true;
    }

    if(!evt_trig_EG12) continue;

    outcheckL1 << endl
	       << "ISO_N : " ;
    for(int i=0 ; i<4 ; i++)
      outcheckL1 << "(" << trig_L1emIso_ieta_To[i] << ";" << trig_L1emIso_iphi_To[i] << ")=" << trig_L1emIso_rank_To[i] << " | ";

    outcheckL1 << endl
	       << "NONISO_N : ";
    for(int i=0 ; i<4 ; i++)
      outcheckL1 << "(" << trig_L1emNonIso_ieta_To[i] << ";" << trig_L1emNonIso_iphi_To[i] << ")=" << trig_L1emNonIso_rank_To[i] 
		 << " | ";

    outcheckL1 << endl
	       << "ISO_M : ";
    for(int i=0 ; i<4 ; i++)
      outcheckL1 << "(" << trig_L1emIso_ieta_M_To[i] << ";" << trig_L1emIso_iphi_M_To[i] << ")=" << trig_L1emIso_rank_M_To[i] 
		 << " | ";

    outcheckL1 << endl
	       << "NONISO_M : ";
    for(int i=0 ; i<4 ; i++)
      outcheckL1 << "(" << trig_L1emNonIso_ieta_M_To[i] << ";" << trig_L1emNonIso_iphi_M_To[i] << ")=" 
		 << trig_L1emNonIso_rank_M_To[i] << " | ";

    outcheckL1 << endl << endl;
  

    ///////////////////////////////////////
    // LOOP OVER SPIKES ///////////////////
    ///////////////////////////////////////

    for(int iS=0 ; iS<spike_N ; iS++) {

      if(iS>=5000) break;
      if(spike_Et[iS]<=0) continue;

      // Find spikes'TT informations
      coords = make_pair( spike_TTieta[iS] , spike_TTiphi[iS] );
      iterTT = adcTT.find( coords );
      if( iterTT != adcTT.end() ) {
	TT_t  = (iterTT->second).first ;
	TT_tE = (iterTT->second).second ;
	
	if(TT_t<0  || TT_t  >= 4032) continue;
	spike_TT_t[iS]  = TT_t ;
	spike_TT_sFGVB[iS] = trig_tower_sFGVB[ TT_t ];
	spike_TT_adc[iS] = trig_tower_adc[ TT_t ];
	
	if(TT_tE<0 || TT_tE >= 4032) continue;
	spike_TT_tE[iS] = TT_tE ;
	for(int iSam=0 ; iSam<5 ; iSam++) {
	  spike_TT_sFGVB_E[iS][iSam] = trig_tower_sFGVB_E[ TT_tE ][iSam];
	  spike_TT_adc_E[iS][iSam] = trig_tower_adc_E[ TT_tE ][iSam];
	}
      }

      // Match spike to L1 candidate
      spike_maxEGtrig[iS] = spike_maxEGtrig_To[iS] = spike_maxEGtrig_M_To[iS] = 0;
      
      matchL1toSpike( trig_L1emIso_N, trig_L1emIso_ieta, trig_L1emIso_iphi, trig_L1emIso_rank,
		      spike_Rieta[iS],  spike_Riphi[iS],  spike_RCTL1iso[iS], spike_maxEGtrig[iS], outcheckL1 );
      
      matchL1toSpike( trig_L1emIso_N, trig_L1emIso_ieta, trig_L1emIso_iphi, trig_L1emIso_rank,
		      spike_Rieta[iS],  spike_Riphi[iS],  spike_RCTL1noniso[iS], spike_maxEGtrig[iS], outcheckL1 );
      
      matchL1toSpike( trig_L1emIso_N_To, trig_L1emIso_ieta_To, trig_L1emIso_iphi_To, trig_L1emIso_rank_To,
		      spike_Rieta[iS],  spike_Riphi[iS],  spike_RCTL1iso_To[iS], spike_maxEGtrig_To[iS], outcheckL1 );
      
      matchL1toSpike( trig_L1emIso_N_To, trig_L1emIso_ieta_To, trig_L1emIso_iphi_To, trig_L1emIso_rank_To,
		      spike_Rieta[iS],  spike_Riphi[iS],  spike_RCTL1noniso_To[iS], spike_maxEGtrig_To[iS], outcheckL1 );
      
      matchL1toSpike( trig_L1emIso_N_M_To, trig_L1emIso_ieta_M_To, trig_L1emIso_iphi_M_To, trig_L1emIso_rank_M_To,
		      spike_Rieta[iS],  spike_Riphi[iS],  spike_RCTL1iso_M_To[iS], spike_maxEGtrig_M_To[iS], outcheckL1 );
      
      matchL1toSpike( trig_L1emIso_N_M_To, trig_L1emIso_ieta_M_To, trig_L1emIso_iphi_M_To, trig_L1emIso_rank_M_To,
		      spike_Rieta[iS],  spike_Riphi[iS],  spike_RCTL1noniso_M_To[iS], spike_maxEGtrig_M_To[iS], outcheckL1 );
      

      // Fill the tree //
      spike_out_outOfTime = spike_outOfTime[iS] ;
      spike_out_severityLevel = spike_severityLevel[iS] ;
      spike_out_Et = spike_Et[iS] ;
      spike_out_phi = spike_phi[iS] ;
      spike_out_eta = spike_eta[iS] ;
      spike_out_theta = spike_theta[iS] ;
      spike_out_TTiphi = spike_TTiphi[iS] ;
      spike_out_TTieta = spike_TTieta[iS] ;
      spike_out_TT_t = spike_TT_t[iS] ;
      spike_out_TT_tE = spike_TT_tE[iS] ;
      spike_out_TT_sFGVB = spike_TT_sFGVB[iS] ;
      spike_out_TT_adc = spike_TT_adc[iS] ;
      spike_out_Riphi = spike_Riphi[iS] ;
      spike_out_Rieta = spike_Rieta[iS] ;

      spike_out_maxEGtrig      = spike_maxEGtrig[iS];
      spike_out_maxEGtrig_To   = spike_maxEGtrig_To[iS];
      spike_out_maxEGtrig_M_To = spike_maxEGtrig_M_To[iS];

      spike_out_RCTL1iso = spike_RCTL1iso[iS];
      spike_out_RCTL1noniso = spike_RCTL1noniso[iS];
      spike_out_RCTL1iso_To = spike_RCTL1iso_To[iS];
      spike_out_RCTL1noniso_To = spike_RCTL1noniso_To[iS];
      spike_out_RCTL1iso_M_To = spike_RCTL1iso_M_To[iS];
      spike_out_RCTL1noniso_M_To = spike_RCTL1noniso_M_To[iS];


      for(int i=0;i<5;i++) {
	spike_out_TT_sFGVB_E[i] = spike_TT_sFGVB_E[iS][i] ;
	spike_out_TT_adc_E[i]   = spike_TT_adc_E[iS][i] ;
      }
      
      outtree->Fill();
    }

  }

    // Record tree
  outlog << "recording tree..." << endl;
  outtree->Write();
  outlog << "recorded !" << endl;
  outfile->Close();
  outlog << "file closed." << endl;

  return 1;

}


/////////////////////////////
// LOOP OVER THE TREE LIST //
/////////////////////////////
int looplist(TString nameChain, TString nameChain_To, vector<TString> treeList, vector<TString> listTrees_To, 
	     int n1, int n2, int nEntries, int nEntriesTo,
	     TString dirIn, TString dirOut, TString usr_hlt, ofstream& outcheck, 
	     int iJson, bool GetTrigOnly, bool debug, bool b_mapcheck) {

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

  TChain * trigOnlyChain = new TChain(nameChain_To);
  MAPTrigOnly mapTo;

  if( GetTrigOnly ) {
    for(u_int iTo=0 ; iTo<listTrees_To.size() ; iTo++)
      trigOnlyChain->Add(listTrees_To[iTo]);
    
    if(debug) cout << "--> map TrigOnly entries --> MAPTrigOnly" << endl;
    mapTo = MapTrigOnlyEntries(trigOnlyChain, nEntriesTo, outcheck, jsonMap, iJson, debug, b_mapcheck);
  }

  // Loop over data trees
  if(debug) cout << "- loop over trees [" << n1 << ";" << n2 << "]" << endl;

  for(int iFile = n1 ; iFile < n2 ; iFile++) {
    if(debug) cout << "--> tree #" << iFile << endl;
    process(trigOnlyChain, jsonMap, mapTo, nameChain, treeList[iFile], iFile, nEntries, 
	    dirIn, dirOut, usr_hlt, iJson, GetTrigOnly, debug);
  }  
  /*
    int process(TChain *trigOnlyChain, MAPJson *jsonMap, MAPTrigOnly mapTo, TString nameChain, TString file, int iFile, int nEntries, 
                TString dirIn, TString dirOut, TString usr_vbtf, double cut_et, bool dojson, TString RunPhase, bool debug)
  */

  return 1;
}
//////////////////////////////////////////////////////////////////////////

void spikes_merge_Data_Emul(int nEntries=-1, int nEntriesTo=-1, int iHalf=0, int nHalf=1,
			    TString dirOut="/home/llr/cms/ndaci/SKWork/macro/skEfficiency/Spikes/Spikes2010/AN2010Repro/test_10adc/",
			    TString dirIn="/data_CMS/cms/ndaci/ndaci_2010B/Repro2012/MergedSpikes/test_10adc/",
			    TString data="AN2010_Data_EGMonitor",
			    TString data_To="AN2010_TrigOnly_EGMonitor_10adc",
			    TString usr_hlt="",
			    TString nameChain="produceNtuple/eIDSimpleTree", 
			    TString nameChain_To="produceNtuple/eIDSimpleTree", 
			    bool GetTrigOnly=true, int iJson=0,
			    bool debug=false, bool b_mapcheck=false)
{
  // Output Log
  ofstream outcheck(dirOut+"mapTrigOnly.txt",ios::out);

  // Process the tree
  if(debug) cout << "Define tree list" << endl;

  vector<TString> treeList;
  vector<TString> listTrees_To;

  if(data=="Spikes_HighPU_HLTEG12") {
    treeList = list_Spikes_HighPU_HLTEG12();
  }
  else if(data=="AN2010_Data_EGMonitor") {
    treeList = list_AN2010_Data_EGMonitor();
  }
  else {
    cout << "asked data tag non recognized ! stopping..." << endl;
    return;
  }

  if( GetTrigOnly ) {
    if(data_To=="Commi_HPU_EmulOnly") {
      listTrees_To = list_Commi_HPU_EmulOnly();
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
    else {
      cout << "asked TrigOnly tag non recognized ! stopping..." << endl;
      return;
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
 
  if(debug) cout << "Launch looplist : n=" << n << " treeList.size()=" << treeList.size() << endl;
  looplist(nameChain, nameChain_To, treeList, listTrees_To, n1, n2, nEntries, nEntriesTo, dirIn, dirOut, 
	   usr_hlt, outcheck, iJson, GetTrigOnly, debug, b_mapcheck);

  if(debug) cout << "DONE" << endl;

}
