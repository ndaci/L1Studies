// VERSION DESIGNED TO RUN ON EMUL-RECO DATA (2 Collections + TPs)

// WARNING !!!!!!!!! INVERSION !!!! //
// in source trees : _modif <=> NoKill  |  _normal <=> Kill
// For L1 candidates :
// ==> source tree variable _modif --> in macro variable _normal  (and vice-versa)
// ==> to keep same comparison variables than 2010B ones
//
// For trigger primitives :
// inversion done directly for comparison variables
//
// adding index [iVtx] and changing a bit output
//
// instead of counting wellId (and others) within if(output) [which corresponds to cases the spike is matched to a L1 candidate]
// --> counting wellId whatever matched to a L1 candidate the spike is or not...

#include "/home/llr/cms/ndaci/SKWork/macro/skEfficiency/analyzers/COMMON/baseFuncNad.1.4.h"

// Tree lists
//#include "/data_CMS/cms/ndaci/ndaci_2010B/EGMonitorModifSpikes2/treelist3.h"
//#include "/data_CMS/cms/ndaci/ndaci_2010B/EGMonitorModif8GeV/treelist.h"
//#include "/data_CMS/cms/ndaci/ndaci_2010B/EGMonitorModif12GeV3/treelistTotal.h"
//#include "/data_CMS/cms/ndaci/ndaci_2010B/FullEmul/test8G17a/new/treelist1.h"
//#include "/data_CMS/cms/ndaci/ndaci_2010B/FullEmul/test8G35a_2/treelist1.h"
//#include "/data_CMS/cms/ndaci/ndaci_2010B/FullEmul/test8G1735_2/treelist1.h"
/*
#include "/data_CMS/cms/ndaci/ndaci_2010B/FullEmul/test8G17a/new/treelist4.h"
#include "/data_CMS/cms/ndaci/ndaci_2010B/FullEmul/test8G1735_2/treelist4.h"
#include "/data_CMS/cms/ndaci/ndaci_2010B/FullEmul/test8G35a_2/treelist2.h"
#include "/data_CMS/cms/ndaci/ndaci_2010B/FullEmul/test1517/treelist3.h"
#include "/data_CMS/cms/ndaci/ndaci_2010B/FullEmul/test2117/treelist1.h"
#include "/data_CMS/cms/ndaci/ndaci_2010B/FullEmul/test5017/treelist2.h"
*/

// Run2011A Commissioning EmulReco data
#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/CommissioningEmulReco/REPRO/treelist_EmulRecoNoKill_Run2011A.h"
#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/CommissioningEmulReco/REPRO2/treelist.h"

// 2011A ReallyUnbias
#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/CommissioningEmulReco/ReallyUnbias/treelist.h"

// 2011A ZeroBias
#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/CommissioningEmulReco/ZeroBias/treelist.h"

// 2011A Newkill
#include "/data_CMS/cms/ndaci/ndaci_2011A/Commissioning/CommissioningEmulReco/NewKill/Run2011A/treelist.h"

void TreeToHistos(TChain* myChain, TString dirOut, ofstream& outlog, int nEntries, TString halftag, float spikeCut, TString hlt_sel, int sev_user, bool debug=false)  {
  
  // get the tree
  int nEvent, nRun, nLumi, trig_isUnbiased, trig_isL1SingleEG2, trig_isL1SingleEG5, trig_isL1SingleEG8, vtx_N;
  int trig_HLT_path[4];

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
  //myChain->SetBranchAddress("trig_HLT_path",&trig_HLT_path);

  myChain->SetBranchAddress("vtx_N",&vtx_N);

  myChain->SetBranchAddress("trig_isUnbiased",&trig_isUnbiased);
  myChain->SetBranchAddress("trig_isL1SingleEG2",&trig_isL1SingleEG2);
  myChain->SetBranchAddress("trig_isL1SingleEG5",&trig_isL1SingleEG5);
  myChain->SetBranchAddress("trig_isL1SingleEG8",&trig_isL1SingleEG8);
  

  // L1 candidates
  myChain->SetBranchAddress("trig_L1emIso_N", &trig_L1emIso_N_modif);
  myChain->SetBranchAddress("trig_L1emIso_ieta", &trig_L1emIso_ieta_modif);
  myChain->SetBranchAddress("trig_L1emIso_iphi", &trig_L1emIso_iphi_modif);
  myChain->SetBranchAddress("trig_L1emIso_rank", &trig_L1emIso_rank_modif);

  myChain->SetBranchAddress("trig_L1emNonIso_N", &trig_L1emNonIso_N_modif);
  myChain->SetBranchAddress("trig_L1emNonIso_ieta", &trig_L1emNonIso_ieta_modif);
  myChain->SetBranchAddress("trig_L1emNonIso_iphi", &trig_L1emNonIso_iphi_modif); 
  myChain->SetBranchAddress("trig_L1emNonIso_rank", &trig_L1emNonIso_rank_modif);

  myChain->SetBranchAddress("trig_L1emIso_N_M", &trig_L1emIso_N);
  myChain->SetBranchAddress("trig_L1emIso_ieta_M", &trig_L1emIso_ieta);
  myChain->SetBranchAddress("trig_L1emIso_iphi_M", &trig_L1emIso_iphi);
  myChain->SetBranchAddress("trig_L1emIso_rank_M", &trig_L1emIso_rank);

  myChain->SetBranchAddress("trig_L1emNonIso_N_M", &trig_L1emNonIso_N);
  myChain->SetBranchAddress("trig_L1emNonIso_ieta_M", &trig_L1emNonIso_ieta);
  myChain->SetBranchAddress("trig_L1emNonIso_iphi_M", &trig_L1emNonIso_iphi);
  myChain->SetBranchAddress("trig_L1emNonIso_rank_M", &trig_L1emNonIso_rank);

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

  MAPTTS adcTT;
  MAPTTS::iterator iterTT;
  pair<int,int> coords;

  vector< pair<int,int> > maskedTT;

  const int nEG = 8;

  int trigthresh[nEG] = {2,5,8,10,12,15,20} ;
  bool triggEG[nEG], triggEGM[nEG];
  for(int i=0 ; i<7 ; i++) {
    triggEG[i] = false ;
    triggEGM[i] = false ;
  }
  vector<int> firedEG_N; // EG 2/5/8/10/12/15/20/30
  vector<int> firedEG_M;

  vector<int> matchIso, matchNonIso, matchIsoM, matchNonIsoM;

  const int nVtx=40;
  int iVtx=0;
  vector<int> IdxVtx;

  int nSpikyL1[nEG][nVtx], nLostSpikyL1[nEG][nVtx], nSpikes[nVtx], nSpikyL1sat[nVtx], nSpikesWellId[nVtx], 
    nDiffSFGVB[nVtx], nMatchedTT[nVtx], nLostTT[nVtx];
  for(iVtx=0 ; iVtx<nVtx ; iVtx++) {
    nSpikes[iVtx] = nSpikyL1sat[iVtx] = nSpikesWellId[iVtx] = nDiffSFGVB[iVtx] = nMatchedTT[iVtx] = nLostTT[iVtx] = 0;
    for(int itg=0 ; itg<nEG ; itg++)
      nSpikyL1[itg][iVtx] = nLostSpikyL1[itg][iVtx] = 0 ;
  }
  bool isGoodRun,output,saturates,wellId,missed,eliminated,zeroedTT,trigEG8,trigEG8M;
  string flag;

  // -------------------------------------------------------------------------------
  // JSON FILE READER
  // -------------------------------------------------------------------------------
  // define map of run/LS
  //string jsonFile = "/data_CMS/cms/ndaci/ndaci_2011A/JSON/Cert_160404-177515_7TeV_PromptReco_Collisions11_JSON.txt";
  string jsonFile = "/data_CMS/cms/ndaci/ndaci_2011A/JSON/Cert_160404-178677_7TeV_PromptReco_Collisions11_JSON.txt";
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
      IdxVtx.clear();
      IdxVtx.push_back(0);
      if(vtx_N>0 && vtx_N<40) IdxVtx.push_back(vtx_N);
      
      //if(iEvent%500 != 0) continue;
      myChain->GetEntry (iEvent) ;

      // run selection
      isGoodRun = AcceptEventByRunAndLumiSection(nRun, nLumi, jsonMap);
      if(!isGoodRun) continue;
      
      // HLT selection
      if( hlt_sel=="unbias" ) {
	if( trig_isUnbiased==0 ) {
	  continue;
	}
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

      // show which run/category is being processed
      if(iEvent==0) {
	nCurrentRun = nRun ;
	outlog << "nRun=" << nRun << endl;
      }
      else if(nRun!=nCurrentRun) {
	nCurrentRun=nRun ;
	outlog << "nRun=" << nRun << endl;
      }
      //cout << nRun << " " << nEvent << endl;

      // map the towers
      adcTT.clear();
      for(int t=0 ; t<trig_tower_N ; t++) {
	coords = make_pair( trig_tower_ieta[t] , trig_tower_iphi[t] );
	adcTT[coords].first.first = trig_tower_adc[t];
	adcTT[coords].second.first = trig_tower_sFGVB[t];
	  //if(trig_tower_ieta[t]>-17 && trig_tower_ieta[t]<18)
	  //for(int i=0 ; i<IdxCat.size() ; i++)
	  //h_adc[0][IdxCat[i]]->Fill(0.25*trig_tower_adc[t]);
      }
      for(int t=0 ; t<trig_tower_N_emul ; t++) {
	coords = make_pair( trig_tower_ieta_emul[t] , trig_tower_iphi_emul[t] );
	iterTT = adcTT.find( coords );
	if( iterTT != adcTT.end() ) {
	  adcTT[coords].first.second = trig_tower_adc_emul[t][2];
	  adcTT[coords].second.second = trig_tower_sFGVB_emul[t][2];
	    //if(trig_tower_ieta_emul[t]>-17 && trig_tower_ieta_emul[t]<18) 
	    //for(int i=0 ; i<IdxCat.size() ; i++)
	    //h_adc[3][iCat]->Fill(0.25*trig_tower_adc_emul[t][2]);
	}
	else {
	  outlog << "mapping problem" << endl;
	  adcTT[coords].first.first = -777;
	  adcTT[coords].first.second = trig_tower_adc_emul[t][2];
	  adcTT[coords].second.first = -777;
	  adcTT[coords].second.second = trig_tower_sFGVB_emul[t][2];
	  
	}
      }

      // get masked towers
      maskedTT.clear();
      for(int iTT=0 ; iTT<trig_nMaskedCh ; iTT++) {
	maskedTT.push_back( make_pair( trig_iMaskedTTeta[iTT] , trig_iMaskedTTphi[iTT] ) );
      }

      // loop over evil rechits     
      if(debug) cout << "loop over evil rechits" << endl;
      saturates = false;
      for(int irh=0 ; irh<spike_N ; irh++) {
	if( spike_severityLevel[irh]>=sev_user && spike_severityLevel[irh]!=5 ) { // (kTime) + kWeird but not kProblematic

	  // selection on the energy
	  if( spike_Et[irh]<spikeCut ) continue;
	  
	  //if(vtx_N>0 && vtx_N<40) IdxVtx.push_back(vtx_N);

	  for(iVtx=0 ; iVtx<IdxVtx.size() ; iVtx++)
	    nSpikes[IdxVtx[iVtx]]++ ;
	  
	  matchIso.clear();
	  matchNonIso.clear();
	  matchIsoM.clear();
	  matchNonIsoM.clear();
	  
	  if(spike_Et[irh]>0) {

	    output=false;
	    saturates = false;
	    for(int itg=0 ; itg<nEG ; itg++) triggEG[itg] = triggEGM[itg] = false;
	    
	    for(int ic=0 ; ic<4 ; ic++) {
	      if(trig_L1emIso_ieta[ic]==spike_Rieta[irh] && trig_L1emIso_iphi[ic]==spike_Riphi[irh]) {
		output = true;
		matchIso.push_back(ic);
		if(trig_L1emIso_rank[ic]>=63) saturates = true;
		for(int itg=0 ; itg<nEG ; itg++)
		  if(trig_L1emIso_rank[ic]>=trigthresh[itg]) triggEG[itg] = true;		
	      }
	      if(trig_L1emNonIso_ieta[ic]==spike_Rieta[irh] && trig_L1emNonIso_iphi[ic]==spike_Riphi[irh]) {
		output = true;
		matchNonIso.push_back(ic);
		if(trig_L1emNonIso_rank[ic]>=63) saturates = true;
		for(int itg=0 ; itg<nEG ; itg++)
		  if(trig_L1emNonIso_rank[ic]>=trigthresh[itg]) triggEG[itg] = true;		
	      }
	      if(trig_L1emIso_ieta_modif[ic]==spike_Rieta[irh] && trig_L1emIso_iphi_modif[ic]==spike_Riphi[irh]) {
		output = true;
		matchIsoM.push_back(ic);
		if(trig_L1emIso_rank_modif[ic]>=63) saturates = true;
		for(int itg=0 ; itg<nEG ; itg++)
		  if(trig_L1emIso_rank_modif[ic]>=trigthresh[itg]) triggEGM[itg] = true;		
	      }
	      if(trig_L1emNonIso_ieta_modif[ic]==spike_Rieta[irh] && trig_L1emNonIso_iphi_modif[ic]==spike_Riphi[irh]) {
		output = true;
		matchNonIsoM.push_back(ic);
		if(trig_L1emNonIso_rank_modif[ic]>=63) saturates = true;
		for(int itg=0 ; itg<nEG ; itg++)
		  if(trig_L1emNonIso_rank_modif[ic]>=trigthresh[itg]) triggEGM[itg] = true;		
	      }
	    }

	    coords = make_pair(spike_TTieta[irh] , spike_TTiphi[irh]);

	    bool b_mask = false;
	    for(int iTT=0 ; iTT<maskedTT.size() ; iTT++) {
	      if( maskedTT[iTT] == coords )
		b_mask = true;
	    }

	    //if(output) {
	    if(true) {
	      outlog << "Spike   nRun=" << nRun  << "   nEvent=" << nEvent << "("<<iEvent<<")" << endl
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
		  outlog << "IsoM rank=" << trig_L1emIso_rank_modif[matchIsoM[ic]] << "      ";
	      if(matchNonIsoM.size()>0)
		for(int ic=0 ; ic<matchNonIsoM.size() ; ic++)
		  outlog << "NonIsoM rank=" << trig_L1emNonIso_rank_modif[matchNonIsoM[ic]] ;
	      
	      outlog<<endl;
	      
	      outlog << "Trigger Tower : adc=" << adcTT[coords].first.first
		     << "   adc_emul="         << adcTT[coords].first.second
		     << "   sFGVB="            << adcTT[coords].second.first
		     << "   sFGVB_emul="       << adcTT[coords].second.second ;
	      if(adcTT[coords].second.first != adcTT[coords].second.second) {
		for(iVtx=0 ; iVtx<IdxVtx.size() ; iVtx++)
		  nDiffSFGVB[IdxVtx[iVtx]]++ ;
		outlog << "   !! sFGVB diff !!" ;
	      }
	      outlog << endl;
	      
	      wellId = zeroedTT = eliminated = missed = false;

	      // non L1-matching-dependant properties
	      if(adcTT[coords].second.first==0) {
		wellId = true;
		for(iVtx=0 ; iVtx<IdxVtx.size() ; iVtx++)
		  nSpikesWellId[IdxVtx[iVtx]]++ ;
	      }
	      if(adcTT[coords].first.second>0) {
		for(iVtx=0 ; iVtx<IdxVtx.size() ; iVtx++)
		  nMatchedTT[IdxVtx[iVtx]]++;
		if(adcTT[coords].first.first==0) {
		  for(iVtx=0 ; iVtx<IdxVtx.size() ; iVtx++)
		    nLostTT[IdxVtx[iVtx]]++;
		  zeroedTT = true;
		}
	      }		  
	      // non L1-matching
		
	      // L1-matching-dependant properties
	      if( matchIso.size()>0 || matchNonIso.size()>0 ) {
		if(saturates) {
		  //saturates = true;
		  for(iVtx=0 ; iVtx<IdxVtx.size() ; iVtx++)
		    nSpikyL1sat[IdxVtx[iVtx]]++ ;
		}
	      } // L1-matching
		
		// EGX
	      for(int itg=0 ; itg<nEG ; itg++) {
		if( triggEG[itg] ) {
		  for(iVtx=0 ; iVtx<IdxVtx.size() ; iVtx++)
		    nSpikyL1[itg][IdxVtx[iVtx]]++ ;
		    
		  if( !triggEGM[itg] ) {
		    for(iVtx=0 ; iVtx<IdxVtx.size() ; iVtx++)
		      nLostSpikyL1[itg][IdxVtx[iVtx]]++ ;
		    eliminated = true;
		  }
		} // EGX-L1-matching
	      }
	      // Modified L1-matching
		
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
	      if( b_mask ) {
		outlog << "Trigger towers of the region" << endl;
		vector<int> ietaTT = getEtaTT( spike_Rieta[irh] );
		vector<int> iphiTT = getPhiTT( spike_Riphi[irh] );
		pair<int,int> etaphi;
		for(int i=0 ; i<ietaTT.size() ; i++) {
		  for(int j=0 ; j<iphiTT.size() ; j++) {
		    etaphi = make_pair( ietaTT[i] , iphiTT[j] );
		    if( adcTT[etaphi].first.second != 0 ) {
		      outlog << "ietaTT="    << ietaTT[i]
			     << "   iphiTT=" << iphiTT[j]
			     << "   adc="    << adcTT[etaphi].first.first
			     << "   adc_emul="   << adcTT[etaphi].first.second
			     << "   sFGVB="    << adcTT[etaphi].second.first
			     << "   sFGVB_emul="   << adcTT[etaphi].second.second
			     << endl;
		    }
		  }
		}
	      }
	      outlog << endl;

	    } // endif output
	  } //endif rechit et > 8 GeV
	}//endif
      } // loop over rechits 
    
    }//loop over events  


  // show counters
  outlog << "COUNTING RESULTS" << endl;

  for(iVtx=0 ; iVtx<nVtx ; iVtx++) {
    outlog << "nSpikes="      << nSpikes[iVtx]
	   << " | nSpikesWellId="<< nSpikesWellId[iVtx]
      //<< " | nSpikyL1EG8="  << nSpikyL1[2][iVtx]
	   << " | nSpikyL1sat="  << nSpikyL1sat[iVtx]
	   << endl
      //<< " | nLostSpikyL1EG8=" << nLostSpikyL1[2][iVtx]
	   << " | nMatchedTT="      << nMatchedTT[iVtx]
	   << " | nLostTT="         << nLostTT[iVtx]
	   << " | nDiffSFGVB="      << nDiffSFGVB[iVtx]
	   << endl
	   << "nSpikyL1EG=[";
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outlog << nSpikyL1[iEG][iVtx] ;
      if(iEG<nEG-1) outlog << " ; " ;
    }
    outlog << endl
	   << "nLostSpikyL1EG=[";
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outlog << nLostSpikyL1[iEG][iVtx] ;
      if(iEG<nEG-1) outlog << " ; " ;
    }
    outlog << endl;
  }

  ofstream outcount(dirOut+"counters"+halftag+".csv",ios::out);
  for(iVtx=0 ; iVtx<nVtx ; iVtx++) {
    outcount << nSpikes[iVtx] 
	     << " | " << nSpikesWellId[iVtx]
	     << " | " << nMatchedTT[iVtx]
	     << " | " << nLostTT[iVtx]
	     << " | " << nDiffSFGVB[iVtx]
	     << " | " << nSpikyL1sat[iVtx]
	     << endl;
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outcount << nSpikyL1[iEG][iVtx] ;
      if(iEG<nEG-1) outcount << " | ";
    }
    outcount << endl;
    for(int iEG=0 ; iEG<nEG ; iEG++) {
      outcount << nLostSpikyL1[iEG][iVtx] ;
      if(iEG<nEG-1) outcount << " | ";
    }
    outcount << endl;
  }
  
  // plots
  if(debug) cout << "is gonna plot in file" << endl;
//   TFile *outplot = new TFile(dirOut+"spike_plots"+halftag+".root","RECREATE");
//   for(iCat=0 ; iCat<nCat ; iCat++) {
//     for(int i=0 ; i<nHistos ; i++) {
//       h_adc[i][iCat]->Write();
//       h_sp[i][iCat]->Write();     
//       h_sptp[i][iCat]->Write();
//     }
//   }
//   outplot->Close();

  if(debug) cout << "plotted in file" << endl;
  gStyle->SetPalette(1);
  if(debug) cout << "is gonna plot" << endl;
 
}
  

int spikes(int iHalf=0, int nHalf=1, TString data_tag="Run2011A_EmulReco", int nEntries=-1, 
	   TString dirOut="/home/llr/cms/ndaci/SKWork/macro/skEfficiency/tagAndProbe/Spike2011A/CommiEmulReco_Run2011A/",
	   TString hlt_sel="", int sev_user=4, float spikeCut=8, bool debug=false)
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

  if( data_tag=="Run2011A_EmulReco_Newkill") {
    listTrees = list_Spikes_Run2011A_Newkill();
  }
  else if( data_tag=="Run2011A_EmulReco_ReallyUnbias") {
    listTrees = list_Spikes_EmulRecoNoKill_2011A_ReallyUnbias();
  }
  else if( data_tag=="Run2011A_EmulReco_ZeroBias" ) {
    listTrees = list_Spikes_EmulRecoNoKill_Run2011A_ZeroBias();
  }
  else if( data_tag=="Run2011A_EmulReco" ) {
    listTrees = list_Spikes_EmulRecoNoKill_Run2011A();
  }
  else if( data_tag=="Run2011A_EmulReco_REPRO2" ) {
    listTrees = list_Spikes_Run2011A_EmulReco_REPRO2();
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
  TreeToHistos(myChain, dirOut, outlog, nEntries, halftag, spikeCut, hlt_sel, sev_user, debug);

  return 1;
}
