////////////////////////////////////////                                                                                            
/// HTauTau TriggerStudies : MuMu T&P //                                                                                            
////////////////////////////////////////                                                                   

// v1.0 : Unifying muTau and muMu macros

//#include "LumiReWeighting.cc"
#include "selection.h"

// Treelists                                                                                                                     
// Include here your treelists

int process(TString triggersource, TString tagfilter, TString tagsanity, TString probefilter, bool probeL1, 
	    TString probeL2, TString probeL25, TString obj_tag, TString obj_probe, 
	    TString call_tag, TString call_probe, float massCutLow, float massCutHigh,
	    float dRcut, float dR_trigmatch, float pTcut_tag, float pTcut_probe, 
	    float pTcut_probeL1jet, float pTcut_probeL1tau,
	    float pTcut_loose, float mTcut, float etacut_tag, float etacut_probe,
	    TString dirStore, TString dirRes, TString file, 
	    int iFile, int nEntries, bool debug,
	    int& nEvents, int& nEventsJSON, int& nEventsMuTau, int& nPairsChargeOpp, 
	    int& nPairsMT, int& nPairsMass, int& nPairsTrigProbe, int& nPairsTrigTag)
{

  // INPUT TREE //            
  TChain * myChain = new TChain("produceNtuple/eIDSimpleTree");
  myChain->Add(file);
  
  // Variables
  int nCurrentRun, nCurrentLumi;
  nCurrentRun = nCurrentLumi = -1;

  // Declaration of leaf types
  Int_t           nEvent;
  Int_t           nRun;
  Int_t           nLumi;
  Int_t           PU_N;
  Double_t        PU_rhoCorr;
  Double_t        PU_sigmaCorr;
  Int_t           vtx_N;
  Double_t        vtx_normalizedChi2[25];   //[vtx_N]
  Double_t        vtx_ndof[25];   //[vtx_N]
  Double_t        vtx_nTracks[25];   //[vtx_N]
  Double_t        vtx_d0[25];   //[vtx_N]
  Double_t        vtx_x[25];   //[vtx_N]
  Double_t        vtx_y[25];   //[vtx_N]
  Double_t        vtx_z[25];   //[vtx_N]
  Char_t          trig_fired_names[50000];

  std::vector<std::string> * trig_HLT_algoStudied = 0 ;
  std::vector<std::string> * trig_HLT_filterStudied = 0 ;

  // L1
  int _trig_L1tau_N,_trig_L1jet_N; 
  double _trig_L1tau_eta[4], _trig_L1tau_phi[4],_trig_L1tau_energy[4],_trig_L1tau_et[4]; 
  double _trig_L1jet_eta[4], _trig_L1jet_phi[4],_trig_L1jet_energy[4],_trig_L1jet_et[4]; 

  Int_t           trig_HLT_N;
  Double_t        trig_HLT_eta[50];   //[trig_HLT_N]
  Double_t        trig_HLT_phi[50];   //[trig_HLT_N]
  Double_t        trig_HLT_energy[50];   //[trig_HLT_N]
  Double_t        trig_HLT_pt[50];   //[trig_HLT_N]
  Int_t           trig_HLT_name[50];   //[trig_HLT_N]
  Int_t           trig_HLT_id[50];   //[trig_HLT_N]
  Int_t           trig_HLT_vids[50];   //[trig_HLT_N]
  int             trig_HLT_NPath ;
  int             trig_HLT_pathname[128] ;   //[trig_HLT_NPath]

  // Electrons
  Int_t           ele_N;
  Int_t           ele_echarge[10];   //[ele_N]
  Double_t        ele_mva[10];   //[ele_N]
  Int_t           ele_hasMatchedConversion[10];   //[ele_N]
  Int_t           ele_missing_hits[10];   //[ele_N]
  Int_t           ele_expected_inner_hits[10];   //[ele_N]
  Double_t        ele_dxy[10];   //[ele_N]
  Double_t        ele_dz[10];   //[ele_N]
  Double_t        ele_dzPV[10];   //[ele_N]
  Double_t        ele_dxyPV[10];   //[ele_N]
  Int_t           ele_isMCEle[10];   //[ele_N]
  Int_t           ele_isMCPhoton[10];   //[ele_N]
  Int_t           ele_isMCHadron[10];   //[ele_N]
  Int_t           ele_isSIM[10];   //[ele_N]
  Int_t           ele_isSIMEle[10];   //[ele_N]
  Int_t           ele_idPDGMatch[10];   //[ele_N]
  Int_t           ele_idPDGmother_MCEle[10];   //[ele_N]
  Int_t           ele_idPDGMatchSim[10];   //[ele_N]
  Int_t           ele_mvaId_noTrg[10];
  Int_t           ele_mvaId_trg[10];
  double          ele_pfIsolCh[10];   //[ele_N]
  double          ele_pfIsolNe[10];   //[ele_N]
  double          ele_pfIsolPh[10];   //[ele_N]
  double          ele_pfIsolChPU[10];   //[ele_N]
  TClonesArray    *electrons;

  // Muons
  Int_t _muons_N;
  //TClonesArray * m_muons;
  Int_t _muons_charge[20];
  Int_t _muons_istracker[20];
  Int_t _muons_isglobal[20];
  Int_t _muons_ispflow[20];
  Double_t _muons_dxy[20];
  Double_t _muons_dz[20];
  Double_t _muons_dxyPV[20];
  Double_t _muons_dzPV[20];
  Double_t _muons_normalizedChi2[20];
  Int_t  _muons_NtrackerHits[20];
  Int_t _muons_NpixelHits[20];
  Int_t _muons_NmuonHits[20];
  Int_t _muons_Nmatches[20];
  Double_t _muons_pfIsolCh[20];
  Double_t _muons_pfIsolNe[20];
  Double_t _muons_pfIsolPh[20];
  Double_t _muons_pfIsolChPU[20];

  // Vector for muons
  TClonesArray * m_muons ;	
  TLorentzVector myvector ; 

  // Taus
  Int_t _nhpsTau;
  float _hpsTauEta[5];
  float _hpsTauPhi[5];
  float _hpsTauPt[5];
  float _hpsTauJetPt[5];
  float _hpsTauLeadPionPt[5];
  float _hpsTauLeadTrackPt[5];
  Int_t _hpsTauCharge[5];
  float _hpsTauChargedIso[5];
  float _hpsTauPhotonsIso[5];
  float _hpsTauDiscrByDecMode[5];
  float _hpsTauDiscrByLooseIso[5];
  float _hpsTauDiscrByMediumIso[5];
  float _hpsTauDiscrByLooseIsoMVA[5];
  float _hpsTauDiscrByMediumIsoMVA[5];
  float _hpsTauDiscrAgainstMuonLoose[5];
  float _hpsTauDiscrAgainstMuonTight[5];
  float _hpsTauDiscrAgainstElecLoose[5];
  float _hpsTauDiscrAgainstElecMedium[5];
  float _hpsTauDiscrAgainstElecTight[5];
  float _hpsTauDiscrAgainstElecMVA[5];

  // MET
  Double_t _met_pf_et;
  Double_t _met_pf_px;
  Double_t _met_pf_py;
  Double_t _met_pf_phi;
  Double_t _met_pf_set;
  Double_t _met_pf_sig;

  myChain->SetBranchAddress("nEvent", &nEvent);
  myChain->SetBranchAddress("nRun", &nRun);
  myChain->SetBranchAddress("nLumi", &nLumi);
  myChain->SetBranchAddress("PU_N", &PU_N);
  myChain->SetBranchAddress("PU_rhoCorr", &PU_rhoCorr);
  myChain->SetBranchAddress("PU_sigmaCorr", &PU_sigmaCorr);
  myChain->SetBranchAddress("vtx_N", &vtx_N);
  myChain->SetBranchAddress("vtx_normalizedChi2", vtx_normalizedChi2);
  myChain->SetBranchAddress("vtx_ndof", vtx_ndof);
  myChain->SetBranchAddress("vtx_nTracks", vtx_nTracks);
  myChain->SetBranchAddress("vtx_d0", vtx_d0);
  myChain->SetBranchAddress("vtx_x", vtx_x);
  myChain->SetBranchAddress("vtx_y", vtx_y);
  myChain->SetBranchAddress("vtx_z", vtx_z);

  myChain->SetBranchAddress("trig_fired_names", trig_fired_names);
  myChain->SetBranchAddress("trig_HLT_algoStudied", &trig_HLT_algoStudied);
  myChain->SetBranchAddress("trig_HLT_filterStudied", &trig_HLT_filterStudied);

  myChain->SetBranchAddress("trig_HLT_N", &trig_HLT_N);
  myChain->SetBranchAddress("trig_HLT_eta", trig_HLT_eta);
  myChain->SetBranchAddress("trig_HLT_phi", trig_HLT_phi);
  myChain->SetBranchAddress("trig_HLT_energy", trig_HLT_energy);
  myChain->SetBranchAddress("trig_HLT_pt", trig_HLT_pt);
  myChain->SetBranchAddress("trig_HLT_name", trig_HLT_name);
  myChain->SetBranchAddress("trig_HLT_id", trig_HLT_id);
  myChain->SetBranchAddress("trig_HLT_vids", trig_HLT_vids);
  myChain->SetBranchAddress("trig_HLT_NPath", &trig_HLT_NPath);
  myChain->SetBranchAddress("trig_HLT_pathname", trig_HLT_pathname);

  // Electrons
  if(obj_tag=="ele" || obj_probe=="ele") {
    electrons = new TClonesArray ("TLorentzVector");
    myChain->SetBranchAddress("electrons", &electrons);
    myChain->SetBranchAddress("ele_N", &ele_N);
    myChain->SetBranchAddress("ele_echarge", ele_echarge);
    myChain->SetBranchAddress("ele_mva", ele_mva);
    myChain->SetBranchAddress("ele_hasMatchedConversion", ele_hasMatchedConversion);
    myChain->SetBranchAddress("ele_missing_hits", ele_missing_hits);
    myChain->SetBranchAddress("ele_expected_inner_hits", ele_expected_inner_hits);
    myChain->SetBranchAddress("ele_dxy", ele_dxy);
    myChain->SetBranchAddress("ele_dz", ele_dz);
    myChain->SetBranchAddress("ele_dxyPV", ele_dxyPV);
    myChain->SetBranchAddress("ele_dzPV", ele_dzPV);
    myChain->SetBranchAddress("ele_isMCEle", ele_isMCEle);
    myChain->SetBranchAddress("ele_isMCPhoton", ele_isMCPhoton);
    myChain->SetBranchAddress("ele_isMCHadron", ele_isMCHadron);
    myChain->SetBranchAddress("ele_isSIM", ele_isSIM);
    myChain->SetBranchAddress("ele_isSIMEle", ele_isSIMEle);
    myChain->SetBranchAddress("ele_idPDGMatch", ele_idPDGMatch);
    myChain->SetBranchAddress("ele_idPDGmother_MCEle", ele_idPDGmother_MCEle);
    myChain->SetBranchAddress("ele_idPDGMatchSim", ele_idPDGMatchSim);
    myChain->SetBranchAddress("ele_pfIsolCh", ele_pfIsolCh);
    myChain->SetBranchAddress("ele_pfIsolNe", ele_pfIsolNe);
    myChain->SetBranchAddress("ele_pfIsolPh", ele_pfIsolPh);
    myChain->SetBranchAddress("ele_pfIsolChPU", ele_pfIsolChPU);
    myChain->SetBranchAddress("ele_mvaId_noTrg", ele_mvaId_noTrg);
    myChain->SetBranchAddress("ele_mvaId_trg", ele_mvaId_trg);
  }

  // Muons
  if(obj_tag=="mu" || obj_probe=="mu") {
    m_muons = new TClonesArray ("TLorentzVector");
    myChain->SetBranchAddress("muons", &m_muons);
    myChain->SetBranchAddress("muons_N",&_muons_N);
    myChain->SetBranchAddress("muons_charge",&_muons_charge);
    myChain->SetBranchAddress("muons_istracker",&_muons_istracker);
    myChain->SetBranchAddress("muons_isglobal",&_muons_isglobal);
    myChain->SetBranchAddress("muons_ispflow",&_muons_ispflow);
    myChain->SetBranchAddress("muons_dxy",&_muons_dxy);
    myChain->SetBranchAddress("muons_dz",&_muons_dz);
    myChain->SetBranchAddress("muons_dxyPV",&_muons_dxyPV);
    myChain->SetBranchAddress("muons_dzPV",&_muons_dzPV);
    myChain->SetBranchAddress("muons_normalizedChi2",&_muons_normalizedChi2);
    myChain->SetBranchAddress("muons_NtrackerHits",&_muons_NtrackerHits);
    myChain->SetBranchAddress("muons_NpixelHits",&_muons_NpixelHits);
    myChain->SetBranchAddress("muons_NmuonHits",&_muons_NmuonHits);
    myChain->SetBranchAddress("muons_Nmatches",&_muons_Nmatches);
    myChain->SetBranchAddress("muons_pfIsolCh",&_muons_pfIsolCh);
    myChain->SetBranchAddress("muons_pfIsolNe",&_muons_pfIsolNe);
    myChain->SetBranchAddress("muons_pfIsolPh",&_muons_pfIsolPh);
    myChain->SetBranchAddress("muons_pfIsolChPU",&_muons_pfIsolChPU);
  }
  //else myChain->SetBranchStatus("muons*",0);

  // Taus
  if(obj_tag=="tau" || obj_probe=="tau") {
    myChain->SetBranchAddress("hpsTau_N", &_nhpsTau);
    myChain->SetBranchAddress("hpsTau_eta", &_hpsTauEta);
    myChain->SetBranchAddress("hpsTau_phi", &_hpsTauPhi);
    myChain->SetBranchAddress("hpsTau_pt", &_hpsTauPt);
    myChain->SetBranchAddress("hpsTau_jet_pt", &_hpsTauJetPt);
    myChain->SetBranchAddress("hpsTau_leadPion_pt", &_hpsTauLeadPionPt);
    myChain->SetBranchAddress("hpsTau_leadTrack_pt", &_hpsTauLeadTrackPt);
    myChain->SetBranchAddress("hpsTau_charge", &_hpsTauCharge);
    myChain->SetBranchAddress("hpsTau_chIso", &_hpsTauChargedIso);
    myChain->SetBranchAddress("hpsTau_phIso", &_hpsTauPhotonsIso);
    myChain->SetBranchAddress("hpsTau_decayMode", &_hpsTauDiscrByDecMode);
    myChain->SetBranchAddress("hpsTau_isoL", &_hpsTauDiscrByLooseIso);
    myChain->SetBranchAddress("hpsTau_isoM", &_hpsTauDiscrByMediumIso);
    myChain->SetBranchAddress("hpsTau_isoMVAL", &_hpsTauDiscrByLooseIsoMVA);
    myChain->SetBranchAddress("hpsTau_isoMVAM", &_hpsTauDiscrByMediumIsoMVA);
    myChain->SetBranchAddress("hpsTau_antiMuL", &_hpsTauDiscrAgainstMuonLoose);
    myChain->SetBranchAddress("hpsTau_antiMuT", &_hpsTauDiscrAgainstMuonTight);
    myChain->SetBranchAddress("hpsTau_antiElL", &_hpsTauDiscrAgainstElecLoose);
    myChain->SetBranchAddress("hpsTau_anitElM", &_hpsTauDiscrAgainstElecMedium);
    myChain->SetBranchAddress("hpsTau_antiElT", &_hpsTauDiscrAgainstElecTight);
    myChain->SetBranchAddress("hpsTau_antiElMVA", &_hpsTauDiscrAgainstElecMVA);

    // MET
    myChain->SetBranchAddress("met_pf_et",&_met_pf_et);
    myChain->SetBranchAddress("met_pf_px",&_met_pf_px);
    myChain->SetBranchAddress("met_pf_py",&_met_pf_py);
    myChain->SetBranchAddress("met_pf_phi",&_met_pf_phi);
    myChain->SetBranchAddress("met_pf_set",&_met_pf_set);
    myChain->SetBranchAddress("met_pf_sig",&_met_pf_sig);

    // L1
    myChain->SetBranchAddress("trig_L1tau_N",     &_trig_L1tau_N);
    myChain->SetBranchAddress("trig_L1tau_eta",   &_trig_L1tau_eta);
    myChain->SetBranchAddress("trig_L1tau_phi",   &_trig_L1tau_phi);
    myChain->SetBranchAddress("trig_L1tau_energy",&_trig_L1tau_energy);
    myChain->SetBranchAddress("trig_L1tau_et",    &_trig_L1tau_et);
    myChain->SetBranchAddress("trig_L1jet_N",     &_trig_L1jet_N);
    myChain->SetBranchAddress("trig_L1jet_eta",   &_trig_L1jet_eta);
    myChain->SetBranchAddress("trig_L1jet_phi",   &_trig_L1jet_phi);
    myChain->SetBranchAddress("trig_L1jet_energy",&_trig_L1jet_energy);
    myChain->SetBranchAddress("trig_L1jet_et",    &_trig_L1jet_et);
    
  }

  // OUTPUT FILE // 
  ostringstream ossi("");
  ossi << iFile ;
  //
  TString name=(TString)(dirStore+"tree_MuMu_"+ossi.str()+".root");
  TFile *outfile = new TFile(name,"RECREATE");
  // 
  name=(TString)(dirRes+"/json/MuMu_"+ossi.str()+".json");
  ofstream out_json(name,ios::out);
  out_json << "{" << endl;
  // 
  ossi.str("");

  /////////////////
  // OUTPUT TREE //
  ////////////////

  TTree* treeTnP = new TTree("treeTnP", "treeTnP");

  float t_mass, mass;
  float t_etraw,t_et, t_pt, t_eta,t_et_tag, t_pt_tag, t_eta_tag , t_weight;
  float t_pt_HLT_tag, t_pt_HLT_tag_sanity, t_pt_L3, t_pt_L25, t_pt_L2, t_et_L1_tau, t_et_L1_jet;
  int t_match, t_L1match, t_L1L2match, t_L1L2L25match, t_L1L2L25L3match, 
    t_Run_match, t_Run_nomatch, t_trig_HLT_N , t_NRun, t_NEvent; //, t_LooseTau, t_MediumTau, t_TightMu;

  treeTnP->Branch("mass",&t_mass, "mass/F");
  treeTnP->Branch("etraw",&t_etraw, "etraw/F");
  treeTnP->Branch("et",&t_et, "et/F");
  treeTnP->Branch("pt",&t_pt, "pt/F");

  treeTnP->Branch("pt_HLT_tag",&t_pt_HLT_tag, "pt_HLT_tag/F");
  treeTnP->Branch("pt_HLT_tag_sanity",&t_pt_HLT_tag_sanity, "pt_HLT_tag_sanity/F");
  treeTnP->Branch("pt_L3",&t_pt_L3, "pt_L3/F");
  treeTnP->Branch("pt_L25",&t_pt_L25, "pt_L25/F");
  treeTnP->Branch("pt_L2",&t_pt_L2, "pt_L2/F");
  treeTnP->Branch("et_L1_tau",&t_et_L1_tau, "et_L1_tau/F");
  treeTnP->Branch("et_L1_jet",&t_et_L1_jet, "et_L1_jet/F");

  treeTnP->Branch("eta" ,&t_eta,  "eta/F");
  treeTnP->Branch("et_tag",&t_et_tag, "et_tag/F");
  treeTnP->Branch("pt_tag",&t_pt_tag, "pt_tag/F");
  treeTnP->Branch("eta_tag" ,&t_eta_tag,  "eta_tag/F");

//   treeTnP->Branch("loosetau" ,&t_LooseTau,  "loosetau/I");
//   treeTnP->Branch("mediumtau" ,&t_MediumTau,  "mediumtau/I");
//   treeTnP->Branch("tightmu" ,&t_TightMu,  "tightmu/I");

  treeTnP->Branch("match" ,&t_match,  "match/I");
  treeTnP->Branch("L1match" ,&t_L1match,  "L1match/I");
  //treeTnP->Branch("L1match" ,&t_L1match,  "L1match/I");
  treeTnP->Branch("L1L2match" ,&t_L1L2match,  "L1L2match/I");
  treeTnP->Branch("L1L2L25match" ,&t_L1L2L25match,  "L1L2L25match/I");
  treeTnP->Branch("L1L2L25L3match" ,&t_L1L2L25L3match,  "L1L2L25L3match/I");

  treeTnP->Branch("weight" ,&t_weight,  "weight/F");
  treeTnP->Branch("Run_match" ,&t_Run_match,  "Run_match/I");
  treeTnP->Branch("Run_nomatch" ,&t_Run_nomatch,  "Run_nomatch/I");
  treeTnP->Branch("NRun" ,&t_NRun,  "NRun/I");
  treeTnP->Branch("NEvent" ,&t_NEvent,  "NEvent/I");
  treeTnP->Branch("trig_HLT_N", &t_trig_HLT_N, "trig_HLT_N/I");

  // -------------------------------------------------------------------------------   
  // JSON FILE READER                                                                   
  // -------------------------------------------------------------------------------          
  if(debug) cout << "Getting JSON file" << endl;
  std::string jsonFile = "/data_CMS/cms/htautau/JSON/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt" ;
  std::map<int, std::vector<std::pair<int, int> > > jsonMap = readJSONFile(jsonFile);
  if(debug) cout << "Got JSON file" << endl;

  int numEntries = myChain->GetEntries () ;
  int nProcess = numEntries;
  if(nEntries>=0 && nEntries<numEntries)
    nProcess = nEntries;
  if(debug) cout << "will process " << nProcess << "/" << numEntries << " entries" << endl;
  
  // Variables //
  TLorentzVector* particle;
  vector< vector<float> > tagprobesColl[2]; // 0:probe 1:tag
  vector<float> tagprobe;
  float eta[2], phi[2], et[2], pt[2], mT[2], idx[2], charge[2];
  //int triggers[2];
  bool sourcetriggered, isGoodRun,selected;
  int nLooseGood;
  int idtriggersource, idtagfilter, idtagsanity, idprobefilter, idprobeL2, idprobeL25;
  idtriggersource=idtagfilter=idtagsanity=idprobefilter=idprobeL2=idprobeL25=-1;

  for (int iEvent = 0 ; iEvent < nProcess ; ++iEvent ) {

    // Initialize //
    for(int iTP=0;iTP<2;iTP++)
      tagprobesColl[iTP].clear();

    tagprobe.clear();

    if (iEvent%1000 == 0) cout<<iEvent<<" processed"<<endl ;
    myChain->GetEntry (iEvent) ;

    nEvents++ ;

    if(iEvent==0) cout << "iFile=" << iFile << endl;

    // HLT FILTERS ////////////////////////////////////////////////////////////////

    if (iEvent == 0 && iFile==0) {
      cout <<"==============================" << endl ;
      cout << "HLT paths : " << endl ;
      for (u_int itr=0; itr<trig_HLT_algoStudied->size() ; ++itr) {
	cout << (*trig_HLT_algoStudied)[itr] << " , index = " << itr << endl;
      }
      cout << "==============================" << endl;
      cout << "HLT filters : "<<endl ;

      for (u_int itr=0; itr<trig_HLT_filterStudied->size() ; ++itr)
	cout << (*trig_HLT_filterStudied)[itr] << " , index = " << itr << endl ;
    }

    idtriggersource = idtagfilter = idprobefilter = idprobeL2 = idprobeL25 = -1;
    
    for (u_int itr=0; itr<trig_HLT_algoStudied->size(); itr++)
      if( (*trig_HLT_algoStudied)[itr]== triggersource ) idtriggersource=itr ; 

    for (u_int itr=0; itr<trig_HLT_filterStudied->size() ; ++itr) {
      //cout << (*trig_HLT_filterStudied)[itr] << " , index = " << itr << endl ;
      if((*trig_HLT_filterStudied)[itr]==tagfilter) idtagfilter=itr;
      if((*trig_HLT_filterStudied)[itr]==tagsanity) idtagsanity=itr;
      if((*trig_HLT_filterStudied)[itr]==probefilter) idprobefilter=itr;
      if((*trig_HLT_filterStudied)[itr]==probeL2) idprobeL2=itr;
      if((*trig_HLT_filterStudied)[itr]==probeL25) idprobeL25=itr;
    }
    
    if(probeL2=="noL2")           idprobeL2=999;
    if(probeL25=="noL25")         idprobeL25=999;
    if(probefilter=="noProbeHLT") idprobefilter=999;
    if(tagsanity=="nosanity")     idtagsanity=999;    
    if(tagfilter=="noTagHLT")     idtagfilter=999;

    if(idtriggersource==-1 || idtagfilter==-1 || idtagsanity==-1 || idprobefilter==-1 || idprobeL2==-1 || idprobeL25==-1) {
      cout << "Bad choice of filters : " << endl;
      return 2;
    }
    ///////////////////////////////////////////////////////////////////////////////
    
    if (iEvent == 0 && iFile==0) {
      
      if(debug) cout << "idtagfilter=" << idtagfilter
		     << "  idtagsanity=" << idtagsanity
		     << "  idprobefilter=" << idprobefilter << endl;
      
      cout <<"=============================="<<endl ;
      
      cout << endl ;
      cout << "==============================" << endl;
      cout << "HLT path selection : " << triggersource << endl;
      cout << "HLT filter tag : " << tagfilter << " index =" << idtagfilter << endl;
      cout << "HLT filter tag sanity : " << tagsanity << " index =" << idtagsanity << endl;
      cout << "HLT filter probe : " << probefilter << " index =" << idprobefilter << endl;
      cout <<"=============================="<<endl ;
      cout <<endl ;
    }

    ////////////RUN SELECTION/////////////////////////////////                
    isGoodRun = AcceptEventByRunAndLumiSection(nRun, nLumi, jsonMap);
    if(!isGoodRun) {
      if(debug) cout << "Fails JSON selection" << endl
		     << nRun << " " << nLumi << " " << nEvent << endl;
      //continue;
    }
    else       nEventsJSON++ ;

    // Generates custom JSON file                 
    if( nRun!=nCurrentRun ) {
      if(nCurrentRun!=-1) {
	out_json << "]," << endl;
      }
      nCurrentLumi=-1;
      out_json << '"' << nRun << '"' << ": [" ;
      nCurrentRun = nRun;
    }
    if( nLumi!=nCurrentLumi ) {
      if( nCurrentLumi!=-1 ) {
	out_json << "," ;
      }
      out_json << "[" << nLumi << "," << nLumi << "]" ;
      nCurrentLumi=nLumi;
    }
    /////                 

    // HLT trigger source must be responsible for acquisition //
    sourcetriggered=false;
    for (u_int itr=0; itr<trig_HLT_NPath; itr++) {
      if( trig_HLT_pathname[itr]==idtriggersource ) sourcetriggered = true ; 
    }
    if (!sourcetriggered) continue ;
    //nEventsPath++ ;
    //////////////////////////////////////////////////////////

    // Check objects presence //
    if(obj_tag=="mu" || obj_probe=="mu") { if( _muons_N < 1 ) continue; }
    if(obj_tag=="ele" || obj_probe=="ele") { if( ele_N < 1 ) continue; }
    if(obj_tag=="tau" || obj_probe=="tau") { if(_nhpsTau < 1) continue; }
    ///////////////////////////
    nEventsMuTau++ ;

    /////////////////////
    // LOOP OVER MUONS //
    /////////////////////

    if(obj_tag=="ele" || obj_probe=="ele") {

      if(debug) cout << "look ele" << endl;

      selected=false;

      for (int iEle = 0 ; iEle < ele_N ; ++iEle) {
	//
	if(debug) cout << "Ele #" << iEle << endl;

	selected=false;
	tagprobe.clear();
	mT[0]=0;
	//
	particle=(TLorentzVector*)(electrons->At (iEle));
	mT[0] = sqrt( 2*particle->Pt()*_met_pf_et * ( 1-cos(fabs(particle->Phi()-_met_pf_phi)) ) );

	tagprobe.push_back(particle->Eta());// eta                                                                               
	tagprobe.push_back(particle->Phi());// phi                                                                               
	tagprobe.push_back(particle->Et()); // eT                                                                                
	tagprobe.push_back(particle->Pt()); // pT                                                                                
	tagprobe.push_back(mT[0]); // mT                                                                                         
	tagprobe.push_back(ele_echarge[iEle]); // charge                                                                        
	tagprobe.push_back(iEle); // index                    

	selected=false;
	if(obj_probe=="ele")
	  selected = selectEle("probe", pTcut_probe, particle->Pt(), 
			       ele_pfIsolNe[iEle], ele_pfIsolPh[iEle], ele_pfIsolChPU[iEle], ele_pfIsolCh[iEle],
			       ele_dxyPV[iEle], ele_dzPV[iEle], particle->Eta(), 
			       ele_expected_inner_hits[iEle], ele_mvaId_noTrg[iEle], ele_hasMatchedConversion[iEle]);
	if(selected) {
	  if(debug) cout << "selected ele probe" << endl;
	  tagprobesColl[0].push_back(tagprobe);
   	}

	selected=false;
	if(obj_tag=="ele")
	  selected = selectEle("tag", pTcut_tag, particle->Pt(), 
			       ele_pfIsolNe[iEle], ele_pfIsolPh[iEle], ele_pfIsolChPU[iEle], ele_pfIsolCh[iEle],
			       ele_dxyPV[iEle], ele_dzPV[iEle], particle->Eta(), 
			       ele_expected_inner_hits[iEle], ele_mvaId_noTrg[iEle], ele_hasMatchedConversion[iEle]);
	
	if(selected) {
	  if(debug) cout << "tagprobe[0].size()=" << tagprobe.size() << endl;
	  tagprobesColl[1].push_back(tagprobe);	  
	  if(debug) cout << "selected electrons tag" << endl;
	}
	
      } // end of loop over electrons
    } // close electrons case

    /////////////////////
    // LOOP OVER MUONS //
    /////////////////////

    if(obj_tag=="mu" || obj_probe=="mu") {

      if(debug) cout << "look muons" << endl;

      nLooseGood=0;
      selected=false;

      for (int iMu = 0 ; iMu < _muons_N ; ++iMu) {
	//
	if(debug) cout << "Mu #" << iMu << endl;

	selected=false;
	tagprobe.clear();
	mT[0]=0;
	//
	particle=(TLorentzVector*)(m_muons->At (iMu));
	mT[0] = sqrt( 2*particle->Pt()*_met_pf_et * ( 1-cos(fabs(particle->Phi()-_met_pf_phi)) ) );

	tagprobe.push_back(particle->Eta());// eta                                                                               
	tagprobe.push_back(particle->Phi());// phi                                                                               
	tagprobe.push_back(particle->Et()); // eT                                                                                
	tagprobe.push_back(particle->Pt()); // pT                                                                                
	tagprobe.push_back(mT[0]); // mT                                                                                         
	tagprobe.push_back(_muons_charge[iMu]); // charge                                                                        
	tagprobe.push_back(iMu); // index                    

	if( selectMu("loose_good", pTcut_loose, _muons_isglobal[iMu], _muons_ispflow[iMu], 
		     _muons_NmuonHits[iMu], _muons_NpixelHits[iMu], 
		     _muons_Nmatches[iMu], _muons_NtrackerHits[iMu], 
		     particle->Pt(), _muons_pfIsolNe[iMu], _muons_pfIsolPh[iMu],
		     _muons_pfIsolChPU[iMu], _muons_pfIsolCh[iMu],
		     _muons_normalizedChi2[iMu], _muons_dxyPV[iMu],
		     _muons_dzPV[iMu], particle->Eta()) )
	  nLooseGood++ ;

	selected=false;
	if(obj_probe=="mu")
	  selected = selectMu(call_probe, pTcut_probe, _muons_isglobal[iMu], _muons_ispflow[iMu], 
			      _muons_NmuonHits[iMu], _muons_NpixelHits[iMu], 
			      _muons_Nmatches[iMu], _muons_NtrackerHits[iMu], 
			      particle->Pt(), _muons_pfIsolNe[iMu], _muons_pfIsolPh[iMu],
			      _muons_pfIsolChPU[iMu], _muons_pfIsolCh[iMu],
			      _muons_normalizedChi2[iMu], _muons_dxyPV[iMu],
			      _muons_dzPV[iMu], particle->Eta());
	if(selected) {
	  if(debug) cout << "selected mu probe" << endl;
	  tagprobesColl[0].push_back(tagprobe);
   	}

	selected=false;
	if(obj_tag=="mu")
	  selected = selectMu(call_tag, pTcut_tag, _muons_isglobal[iMu], _muons_ispflow[iMu], 
			      _muons_NmuonHits[iMu], _muons_NpixelHits[iMu], 
			      _muons_Nmatches[iMu], _muons_NtrackerHits[iMu], 
			      particle->Pt(), _muons_pfIsolNe[iMu], _muons_pfIsolPh[iMu],
			      _muons_pfIsolChPU[iMu], _muons_pfIsolCh[iMu],
			      _muons_normalizedChi2[iMu], _muons_dxyPV[iMu],
			      _muons_dzPV[iMu], particle->Eta());
	
	if(selected || true) {
	  if(debug) cout << "tagprobe[0].size()=" << tagprobe.size() << endl;
	  tagprobesColl[1].push_back(tagprobe);	  
	  if(debug) cout << "selected mu tag" << endl;
	}
	
      } // end of loop over muons
      if(obj_probe=="tau" && nLooseGood>1) continue;
    } // close muon case

    /////////////////////
    // LOOP OVER TAU   //
    /////////////////////

    if(obj_tag=="tau" || obj_probe=="tau") {

      if(debug) cout << "Look Taus" << endl;
      
      for (int iTau = 0 ; iTau < _nhpsTau ; iTau++) {
	
	// Initialize //
	tagprobe.clear();
	mT[0]=0;
	mT[0] = sqrt( 2 * _hpsTauPt[iTau] * _met_pf_et * ( 1-cos(fabs(_hpsTauPhi[iTau]-_met_pf_phi)) ) );

	tagprobe.push_back(_hpsTauEta[iTau]);// eta                                                                              
	tagprobe.push_back(_hpsTauPhi[iTau]);// phi                                                                              
	tagprobe.push_back(_hpsTauPt[iTau]); // eT                                                                               
	tagprobe.push_back(_hpsTauPt[iTau]); // pT                                                                               
	tagprobe.push_back(mT[0]); // mT                                                                                         
	tagprobe.push_back(_hpsTauCharge[iTau]); // charge                                                                       
	tagprobe.push_back(iTau); // index                 

	selected=false;
	if(obj_probe=="tau")
	  selected = selectTau(call_probe, pTcut_probe, _hpsTauPt[iTau], _hpsTauDiscrByDecMode[iTau], 
			       _hpsTauDiscrByLooseIsoMVA[iTau], _hpsTauDiscrByMediumIsoMVA[iTau],
			       _hpsTauDiscrAgainstElecLoose[iTau], _hpsTauDiscrAgainstElecMedium[iTau], _hpsTauDiscrAgainstElecTight[iTau],
			       _hpsTauDiscrAgainstElecMVA[iTau], _hpsTauDiscrAgainstMuonLoose[iTau], _hpsTauDiscrAgainstMuonTight[iTau],
			       _hpsTauEta[iTau], etacut_probe);
	if(selected) {
	  if(debug) cout << "selected tau probe" << endl;
	  tagprobesColl[0].push_back(tagprobe);
	  if(debug) cout << "has pushed back" << endl;
	}

	selected=false;
	if(obj_tag=="tau")
	  selected = selectTau(call_probe, pTcut_probe, _hpsTauPt[iTau], _hpsTauDiscrByDecMode[iTau], 
			       _hpsTauDiscrByLooseIsoMVA[iTau], _hpsTauDiscrByMediumIsoMVA[iTau],
			       _hpsTauDiscrAgainstElecLoose[iTau], _hpsTauDiscrAgainstElecMedium[iTau], _hpsTauDiscrAgainstElecTight[iTau],
			       _hpsTauDiscrAgainstElecMVA[iTau], _hpsTauDiscrAgainstMuonLoose[iTau], _hpsTauDiscrAgainstMuonTight[iTau],
			       _hpsTauEta[iTau],etacut_tag);
	
	if(selected) {
	  tagprobesColl[1].push_back(tagprobe);
	  if(debug) cout << "selected tau tag" << endl;
	}

      } // end of loop over taus
    } // close taus case


    ///////////
    // PAIRS //
    ///////////

    int idxTP[2]={0,0};
    mass=0;
    TLorentzVector object[2];

    if(debug) cout << "Explore pairs : "
		   << "  tagprobesColl[0].size()=" << tagprobesColl[0].size()
		   << "  tagprobesColl[1].size()=" << tagprobesColl[1].size()
		   << endl;

    // Get one candidate pair //
    for(u_int iTag=0 ; iTag<tagprobesColl[1].size() ; iTag++) {
      for(u_int iProbe=0 ; iProbe<tagprobesColl[0].size() ; iProbe++) {
	
	if(debug) cout << "iTag=" << iTag << "  iProbe=" << iProbe << endl;

	idxTP[0]=iProbe;
	idxTP[1]=iTag;
	mass=0;

	if(debug) cout << "Size probe #" << iProbe << " : " << tagprobesColl[0][iProbe].size() << endl;
	//if(debug) cout << tagprobesColl[1][0].size() << endl;

	for(int iTP=0 ; iTP<2 ; iTP++) {
	  eta[iTP] = tagprobesColl[iTP][idxTP[iTP]][0];
	  phi[iTP] = tagprobesColl[iTP][idxTP[iTP]][1];
	  et[iTP] = tagprobesColl[iTP][idxTP[iTP]][2];
	  pt[iTP] = tagprobesColl[iTP][idxTP[iTP]][3];
	  mT[iTP] = tagprobesColl[iTP][idxTP[iTP]][4];
	  charge[iTP] = tagprobesColl[iTP][idxTP[iTP]][5];
	  idx[iTP] = tagprobesColl[iTP][idxTP[iTP]][6];
	}

	if(debug) cout << "Got variables" << endl;

	// Check objects are not the same
	if(obj_tag==obj_probe)
	  if(idx[0]==idx[1]) continue;

	// Check charge opposition
	if(charge[0]*charge[1] != -1) continue;
	nPairsChargeOpp++ ;

	// Check transverse mass if needed
	if(obj_tag=="mu" && obj_probe=="tau")
	  if(mT[1]>mTcut) continue;
	nPairsMT++ ;

	if(debug) cout << "mT checked" << endl;

	// Check invariant mass of the pair
	if(obj_tag=="tau") { 
	  object[1].SetPtEtaPhiM(pt[1], eta[1], phi[1], 1.777);
	}
	if(obj_probe=="tau") { 
	  object[0].SetPtEtaPhiM(pt[0], eta[0], phi[0], 1.777);
	}

	if(obj_tag=="mu")
	  object[1]= *((TLorentzVector*)( m_muons->At((int)idx[1]) ));
	if(obj_probe=="mu")
	  object[0]= *((TLorentzVector*)( m_muons->At((int)idx[0]) ));

	if(obj_tag=="ele")
	  object[1]= *((TLorentzVector*)( electrons->At((int)idx[1]) ));
	if(obj_probe=="ele")
	  object[0]= *((TLorentzVector*)( electrons->At((int)idx[0]) ));


	if(debug) cout << pt[0] << " " << eta[0] << " " << phi[0] << endl;
	//
	mass = (object[0]+object[1]).M();
	if(mass>massCutHigh || mass<massCutLow) continue;
	
	if(debug) cout << "mass checked m=" << mass << endl;
	nPairsMass++ ;

	// Check dR between objects
	if( deltaR(eta[0],phi[0],eta[1],phi[1]) < dRcut) continue;
	
	// Check triggering properties ///////////////////////////////////////
	bool trigL1, trigL2, trigL25, trigL3, trig_tag, trig_tag_sanity;
	int trig_probe=0;
	trigL1 = trigL2 = trigL25 = trigL3 = trig_tag = trig_tag_sanity = false;
	t_pt_HLT_tag = t_pt_HLT_tag_sanity = t_pt_L3 = t_pt_L25 = t_pt_L2 = t_et_L1_tau = t_et_L1_jet = 0;

	if( probefilter=="noProbeHLT" ) trigL3=true;
	if( probeL25=="noL25")          trigL25=true;
	if( probeL2=="noL2")            trigL2=true;
	if( tagsanity=="nosanity")      trig_tag_sanity=true;
	if (tagfilter=="noTagHLT")              trig_tag = true;

	// HLT filters //
	for (int itr=0; itr<trig_HLT_N; itr++) {

	  //if(debug) cout << "MATCHED HLT OBJECT : triggers filter #" << trig_HLT_name[itr] << endl;

	    // TAG //
	    if( deltaR(trig_HLT_eta[itr], trig_HLT_phi[itr], eta[1], phi[1] ) < dR_trigmatch ) {

	      if(debug) cout << "trig_HLT_name[itr]=" << trig_HLT_name[itr] << " | idtagfilter=" << idtagfilter << endl;

	      if( trig_HLT_name[itr]==idtagfilter) {
		trig_tag=true;
		t_pt_HLT_tag = trig_HLT_pt[itr];
		if(debug) cout << "matched" << endl;
	      }

	      if( trig_HLT_name[itr]==idtagsanity) {
		trig_tag_sanity=true;
		t_pt_HLT_tag_sanity = trig_HLT_pt[itr];
	      }

	    }

	    // Probe //
	    if( probefilter!="no" ) {

	      if(debug) cout << "trig_HLT_name[itr]=" << trig_HLT_name[itr] << " | idprobefilter=" << idprobefilter << endl;

	      if( deltaR(trig_HLT_eta[itr], trig_HLT_phi[itr], eta[0], phi[0] ) < dR_trigmatch ) {
		if( trig_HLT_name[itr]==idprobefilter) {
		  trigL3=true;
		  t_pt_L3 = trig_HLT_pt[itr];
		  if(debug) cout << "matched" << endl;
		}
	      }
	    }

	    if( probeL25!="no") {
	      if( deltaR(trig_HLT_eta[itr], trig_HLT_phi[itr], eta[0], phi[0] ) < dR_trigmatch )
		if( trig_HLT_name[itr]==idprobeL25 )
		  trigL25=true;
		  t_pt_L25 = trig_HLT_pt[itr];
	    }

	    if( probeL2!="no") {
	      if( deltaR(trig_HLT_eta[itr], trig_HLT_phi[itr], eta[0], phi[0] ) < dR_trigmatch )
		if( trig_HLT_name[itr]==idprobeL2 )
		  trigL2=true;
		  t_pt_L2 = trig_HLT_pt[itr];
	    }
	}
	
	if(!trig_tag) continue;
	if(obj_probe=="tau" && !trig_tag_sanity) continue;
	//if(debug) cout << "Tag candidate triggers tag filter : is really a tag" << endl;
	nPairsTrigTag++ ;

	if( !probeL1 ) trigL1=true;
	else {
	  // L1 tau candidates
	  for(int itr=0; itr<_trig_L1tau_N; itr++)
	    if( deltaR(_trig_L1tau_eta[itr], _trig_L1tau_phi[itr], eta[0], phi[0] ) < dR_trigmatch )
	      if( _trig_L1tau_et[itr] > pTcut_probeL1tau ) {
		trigL1=true;
		t_et_L1_tau = _trig_L1tau_et[itr];
	      }

	  // L1 jet candidates
	  for(int itr=0; itr<_trig_L1jet_N; itr++)
	    if( deltaR(_trig_L1jet_eta[itr], _trig_L1jet_phi[itr], eta[0], phi[0] ) < dR_trigmatch )
	      if( _trig_L1jet_et[itr] > pTcut_probeL1jet ) {
		trigL1=true;
		t_et_L1_jet = _trig_L1jet_et[itr];
	      }

	}
	//////////////////////////////////////////////////////////////////////
	if(debug) cout << "Pair selected" << endl;
	//if(trigL1 && trigL2 && trigL25 && trigL3) { trig_probe=1; nPairsTrigProbe++; }

	t_L1match = t_L1L2match = t_L1L2L25match = t_L1L2L25L3match = t_match = 0;
	
	
	if(trigL3)
	  t_match=1;

	if(trigL1) {
	  t_L1match=1;
	  if(trigL2) {
	    t_L1L2match=1;
	    if(trigL25) {
	      t_L1L2L25match=1;
	      if(trigL3) {
		t_L1L2L25L3match=1;
		nPairsTrigProbe++;
	      }
	    }
	  }
	}

	if(debug) cout << "trig_probe=" << trig_probe << endl;

	t_Run_match = -1 ;
	t_Run_nomatch = -1 ;
	
	if(trig_tag==1) t_Run_match = nRun ;
	else t_Run_nomatch = nRun ;
	
	t_mass = mass ;
	t_etraw = -1;
	t_pt_tag = pt[1]  ;
	t_eta_tag = eta[1] ;
	t_et = et[0];
	t_pt = pt[0]  ;
	t_eta = eta[0] ;
	t_trig_HLT_N = trig_HLT_N ;
	t_NRun = nRun;
	t_NEvent = nEvent;
	t_weight = 1. ;
	treeTnP->Fill() ;
	
      }
    }
    
  } // loop entries
  
  if(debug) cout << "treeTnP->Write()" << endl;
  treeTnP->Write() ;
  
  return 0;
  
}

int tagAndProbe(int nEntries=-1, int iHalf=0, int nHalf=1,
		TString triggersource="", TString tagfilter="", TString tagsanity="", TString probefilter="no",
		bool probeL1=true, TString probeL2="no", TString probeL25="no",
		TString obj_tag="mu", TString obj_probe="tau",
		TString call_tag="tag", TString call_probe="mutau-probe", 
		float massCutLow=45, float massCutHigh=70, float dRcut=0.52, float dR_trigmatch=0.3,
		float pTcut_tag=24, float pTcut_probe=10, float pTcut_loose=15,
		float pTcut_probeL1jet=24, float pTcut_probeL1tau=24,
		float mTcut=40, float etacut_tag=2.4, float etacut_probe=2.4,
		TString dirRes="/home/llr/cms/azabi/CMSSW/CMSSW_535_TAUS/src/Htautau/TriggerStudies/test/tagAndProbe/Macros/",
		TString dirStore="/data_CMS/cms/azabi/",
		TString data="HTT_Alex",
		bool debug=false)
{

  int nEvents=0, nEventsJSON=0, nEventsPath=0, nEventsMuTau=0, nPairsChargeOpp=0, nPairsMT=0, nPairsMass=0, nPairsTrigTag=0, nPairsTrigProbe=0;

  vector<TString> treeList;

  // Insert here your functions list_... from your treelists .h
  /*
  if(data=="HTT_MuTau_2012A_PRV1") {
    treeList = list_HTT_MuTau_2012A_PRV1();
  }
  else if(data=="HTT_MuTau_2012B_PRV1") {
    treeList = list_HTT_MuTau_2012B_PRV1();
  }
  else {
    cout << "asked data tag non recognized ! stopping..." << endl;
    return -1;
  }
  */

  cout << "USE " << data << " TAG" << endl;

  int n,n1,n2;
  n=n1=n2=0;

  n = treeList.size() / nHalf ;
  n1 = iHalf * n ;

  if( iHalf < (nHalf-1) )
    n2 = (iHalf+1) * n ;

  else if( iHalf == nHalf-1 )
    n2 = treeList.size();

  for(int i = n1 ; i < n2 ; i++) {

    if(debug) cout << "process file #" << i << " (n1=" << n1 << " ; n2=" << n2 << ")" << endl;

    /*
    */

    process(triggersource,  tagfilter,  tagsanity, probefilter, 
	    probeL1, probeL2, probeL25,
	    obj_tag,  obj_probe,  call_tag,  call_probe, 
	    massCutLow, massCutHigh, 
	    dRcut, dR_trigmatch,
	    pTcut_tag, pTcut_probe, pTcut_probeL1jet, pTcut_probeL1tau,
	    pTcut_loose, mTcut, etacut_tag, etacut_probe,
	    dirStore,  dirRes, treeList[i],
	    i,  nEntries,  debug,
	    nEvents, nEventsJSON, nEventsMuTau, nPairsChargeOpp, 
	    nPairsMT, nPairsMass, nPairsTrigProbe, nPairsTrigTag);
    
    if(debug) cout << "ended file #" << i << endl;
  }

  cout << "========================================="
       << "COUNTERS" << "  nEvents=" << nEvents << "  nEventsJSON=" << nEventsJSON 
       << "  nEventsPath=" << nEventsPath << "  nEventsMuTau=" << nEventsMuTau
       << "  nPairsChargeOpp=" << nPairsChargeOpp << "  nPairsMT=" << nPairsMT
       << "  nPairsMass=" << nPairsMass << "  nPairsTrigTag=" << nPairsTrigTag
       << "  nPairsTrigProbe=" << nPairsTrigProbe
       << endl;



  return 0;
}
