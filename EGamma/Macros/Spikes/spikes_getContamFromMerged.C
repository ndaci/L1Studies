#include "../Common/baseFuncNad.1.5.h"

int spikes_getContamFromMerged(int nEntries=-1,
			       TString dirOut="/home/llr/cms/ndaci/SKWork/macro/skEfficiency/Spikes/HighPuContamination/try2/",
			       TString dirIn="/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/CommissioningEmulReco/NewKill/highPU/MergedToData/",
			       TString file="all_spikes.root",
			       TString usr_hlt="noHLT",
			       TString nameChain="Spikes", 
			       float spikeCut=8.,
			       bool debug=false)
{

  ofstream outlog(dirOut+"/log_getContam.txt",ios::out);

  // INPUT TREE //                                                                                                                   
  if(debug) cout << "--- add data trees to myChain" << endl;
  TChain * myChain = new TChain(nameChain);
  myChain->Add(dirIn+file);

  // VARIABLES //                                                                                      
  int nEvent, nRun, nLumi ; // data                                                      
  //int nEvent_To, nRun_To, nLumi_To; // trigOnly
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
    spike_outOfTime[5000];
  double  spike_Et[5000], spike_eta[5000], spike_phi[5000], spike_theta[5000];
  int spike_RCTL1iso[5000], spike_RCTL1noniso[5000], spike_RCTL1iso_To[5000], spike_RCTL1noniso_To[5000], 
    spike_RCTL1iso_M_To[5000], spike_RCTL1noniso_M_To[5000],
    spike_maxEGtrig[5000], spike_maxEGtrig_To[5000], spike_maxEGtrig_M_To[5000];
  
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
    spike_Riphi[isp] = -999 ;
    spike_Rieta[isp] = -999 ;
    //
    spike_RCTL1iso[isp] = -999 ;
    spike_RCTL1noniso[isp] = -999 ;
    spike_RCTL1iso_To[isp] = -999 ;
    spike_RCTL1noniso_To[isp] = -999 ;
    spike_RCTL1iso_M_To[isp] = -999 ;
    spike_RCTL1noniso_M_To[isp] = -999 ;
    //spike_time[isp] = -999 ;
  }


  myChain->SetBranchAddress("nEvent",&nEvent);
  myChain->SetBranchAddress("nRun",&nRun);
  myChain->SetBranchAddress("nLumi",&nLumi);

  myChain->SetBranchAddress("vtx_N",&vtx_N);

  myChain->SetBranchAddress("trig_HLT_path",&trig_HLT_path);
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
//
  myChain->SetBranchAddress("spike_maxEGtrig",&spike_maxEGtrig);
  myChain->SetBranchAddress("spike_maxEGtrig_To",&spike_maxEGtrig_To);
  myChain->SetBranchAddress("spike_maxEGtrig_M_To",&spike_maxEGtrig_M_To);

  myChain->SetBranchAddress("trig_tower_N", &trig_tower_N);
  myChain->SetBranchAddress("trig_tower_ieta",  &trig_tower_ieta);
  myChain->SetBranchAddress("trig_tower_iphi",  &trig_tower_iphi);
  myChain->SetBranchAddress("trig_tower_adc",  &trig_tower_adc);
  myChain->SetBranchAddress("trig_tower_sFGVB",  &trig_tower_sFGVB);
  myChain->SetBranchAddress("trig_tower_FG",  &trig_tower_FG);

  myChain->SetBranchAddress("trig_tower_N_E",     &trig_tower_N_E);
  myChain->SetBranchAddress("trig_tower_ieta_E",  &trig_tower_ieta_E);
  myChain->SetBranchAddress("trig_tower_iphi_E",  &trig_tower_iphi_E);
  myChain->SetBranchAddress("trig_tower_adc_E",   &trig_tower_adc_E);
  myChain->SetBranchAddress("trig_tower_sFGVB_E", &trig_tower_sFGVB_E);
  myChain->SetBranchAddress("trig_tower_FG_E", &trig_tower_FG_E);

  // HCAL TP                                                        
  myChain->SetBranchAddress("trig_tower_hcal_N", &trig_tower_hcal_N);
  myChain->SetBranchAddress("trig_tower_hcal_ieta",  &trig_tower_hcal_ieta);
  myChain->SetBranchAddress("trig_tower_hcal_iphi",  &trig_tower_hcal_iphi);
  myChain->SetBranchAddress("trig_tower_hcal_et",  &trig_tower_hcal_et);
  myChain->SetBranchAddress("trig_tower_hcal_FG",  &trig_tower_hcal_FG);

  myChain->SetBranchStatus("trig_tower*",0);
  
  // L1 candidates collections                        
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

  int numEntries = myChain->GetEntries () ;
  int nProcess = numEntries;
  if(nEntries>=0 && nEntries<numEntries)
    nProcess = nEntries;

  TString filename;
  bool evt_trig_EG12, evt_trig, evt_trig_M;
  int spike_trig, spike_trig_M;

  // Numbers //
  const int nM=2;
  const int nEG=8;
  int thresh[nEG] = {2,5,8,10,12,15,20,30} ;
  TString trigname[nEG] = {"2","5","8","10","12","15","20","30"};
  TString collname[nM]={"N","M"};

  // HISTOGRAMS //
  TH1F * h_spikes;
  TH1F * h_triggering_spikes[nEG][nM];
  TH1F * h_evts_trigBy_any[nEG][nM];
  TH1F * h_evts_trigBy_spikes[nEG][nM];

  const int nVtx = 40;
  TString name_histo;

  h_spikes = new TH1F("h_spikes","h_spikes",nVtx,0,nVtx); // abscisse = vtx_N

  for(int iEG=0 ; iEG<nEG ; iEG++ ) {
    for(int iM=0 ; iM<nM ; iM++ ) {
      
      name_histo = "h_triggering_spikes_EG"+trigname[iEG]+"_"+collname[iM];
      h_triggering_spikes[iEG][iM]=new TH1F(name_histo,name_histo,nVtx,0,nVtx);
      
      name_histo = "h_evts_trigBy_any_EG"+trigname[iEG]+"_"+collname[iM];
      h_evts_trigBy_any[iEG][iM]=new TH1F(name_histo,name_histo,nVtx,0,nVtx);
      
      name_histo = "h_evts_trigBy_spikes_EG"+trigname[iEG]+"_"+collname[iM];
      h_evts_trigBy_spikes[iEG][iM]=new TH1F(name_histo,name_histo,nVtx,0,nVtx);     

    }
  }
  
  // Loop over entries //
  int nCurrentRun = -999;
  outlog << "will process " << nProcess << "/" << numEntries << "entries" << endl;
  //
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

    // HLT selection //
    if(usr_hlt=="unbias")
      if(trig_HLT_path[0]==0) continue;
    
    // CHECK EG12 FIRED //
    evt_trig_EG12=false;
    for(int i=0 ; i<4 ; i++) {
      if(trig_L1emIso_rank[i] >= 12) evt_trig_EG12=true;
      if(trig_L1emNonIso_rank[i] >= 12) evt_trig_EG12=true;
    }
    if(!evt_trig_EG12) continue;

    // EXTRACT CONTAMINATION //
    //
    spike_trig=spike_trig_M=0;
    // Loop over Spikes
    for(int iS=0 ; iS<spike_N ; iS++) {

      h_spikes->Fill(vtx_N);

      if( spike_Et[iS] < spikeCut ) continue;
      
      if( spike_maxEGtrig[iS] != spike_maxEGtrig_To[iS] )
	if(debug) cout << "Spike #" << iS << " discripancy between OnlineFromData and OnlineFromEmul !" << endl;
    
      if( spike_maxEGtrig[iS]>spike_trig ) spike_trig = spike_maxEGtrig[iS];
      if( spike_maxEGtrig_M_To[iS]>spike_trig_M ) spike_trig_M = spike_maxEGtrig_M_To[iS];

      for(int iEG=0 ; iEG<nEG ; iEG++) {
	if( spike_maxEGtrig[iS]>thresh[iEG] ) h_triggering_spikes[iEG][0]->Fill(vtx_N);
	if( spike_maxEGtrig_M_To[iS]>thresh[iEG] ) h_triggering_spikes[iEG][1]->Fill(vtx_N);
      }
    }
    // 
    // Object triggering event
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      
      // EG trig by any object
      evt_trig=evt_trig_M=false;
      for(int i=0 ; i<4 ; i++) {
	if(trig_L1emIso_rank_To[i] >= thresh[iEG]) evt_trig=true;
	if(trig_L1emNonIso_rank_To[i] >= thresh[iEG]) evt_trig=true;
	//
	if(trig_L1emIso_rank_M_To[i] >= thresh[iEG]) evt_trig_M=true;
	if(trig_L1emNonIso_rank_M_To[i] >= thresh[iEG]) evt_trig_M=true;	
      }
      
      if(evt_trig)
	h_evts_trigBy_any[iEG][0]->Fill(vtx_N);
      
      if(evt_trig_M)
	h_evts_trigBy_any[iEG][1]->Fill(vtx_N);
      
      // EG trig by spike
      if(spike_trig>=thresh[iEG])
	h_evts_trigBy_spikes[iEG][0]->Fill(vtx_N);

      if(spike_trig_M>=thresh[iEG])
	h_evts_trigBy_spikes[iEG][1]->Fill(vtx_N);
    }
     
  }
    
  TFile *outplot = new TFile(dirOut+"/spike_plots.root","RECREATE");
  h_spikes->Write();

  for(int iEG=0 ; iEG<nEG ; iEG++) {
    for(int iM=0; iM<nM ; iM++) {
      h_triggering_spikes[iEG][iM]->Write();
      h_evts_trigBy_any[iEG][iM]->Write();
      h_evts_trigBy_spikes[iEG][iM]->Write();
    }
  }
  outplot->Close();

  cout << "DONE, BABY" << endl;

  return 0;
}
