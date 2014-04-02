#include "../Common/baseFuncNad.h"

bool selectMu(TString kindsel, double pTcut, int global, int pflow, int valMuHits, int valPixHits, int matchStations, int tkLayHits, 
	      double pT, double neutriso, double phiso, double dB, double chiso, double nX2, double dxy, double dz, double eta)
{

  // same selection for Tag & Probe Muons (TightMuon)
  if(kindsel=="tag" || kindsel=="probe" || kindsel=="loose_good") { 
    // ISO
    double isovar, niso;
    isovar=niso=0;
    //
    niso = neutriso + phiso - 0.5*dB ; // neutral iso + photon iso - 0.5*dB
    if(niso < 0) niso=0; // taking max(0,niso)
    //
    if(pT!=0) isovar = (chiso + niso)/pT ; // chiso=charged iso
    else return false;
    //
    if( isovar >= 0.3 ) return false; // 
    //
    // pT cut
    if( pT < pTcut ) return false; // pT
    //
    // ID
    if( global != 1 ) return false; // is Global Mu
    if( fabs(eta) > 2.4 )   return false; // eta (covering all det)
    //
    if(kindsel=="tag" || kindsel=="probe") {
      if( isovar >= 0.1 ) return false; // 
      if( pflow  != 1 ) return false; // is PF Mu
      if( valMuHits  <= 0 ) return false; // # valid mu hits
      if( valPixHits <= 0 ) return false; // # valid pixel hits
      if( matchStations <= 1 ) return false; // # matched stations 
      if( tkLayHits     <= 5 ) return false; // # tracker layer with hits
      if( nX2       >= 10. ) return false; // norm Chi2
      if( fabs(dxy) >= 0.2 ) return false; // |dxy| 2 mm
      if( fabs(dz)  >= 0.2 ) return false; // |dz| 5 mm
      if( fabs(eta) >= 2.1 ) return false; // |eta| 2p1
    }
    return true;
  }
  return false;

}

bool selectEle(TString kindsel, double pTcut, double pT, double neutriso, 
	       double phiso, double dB, double chiso, 
	       double dxy, double dz, double eta,
	       int miss_hits, int mvaID_noTrg, int isConv)
{

  if(kindsel=="tag" || kindsel=="probe") {
    // ISO
    double isovar, niso;
    isovar=niso=0;
    //
    niso = neutriso + phiso - 0.5*dB ; // neutral iso + photon iso - 0.5*dB
    if(niso < 0) niso=0; // taking max(0,niso)
    if(pT!=0) isovar = (chiso + niso)/pT ; // chiso=charged iso
    else return false;
    if( isovar >= 0.1 ) return false; // 
    // pT cut
    if( pT < pTcut ) return false; // pT
    // ID
    if( fabs(dxy) >= 0.045 ) return false; // |dxy| 2 mm
    if( fabs(dz)  >= 0.2 ) return false; // |dz| 5 mm
    // ID MVA : Tight ID (eTau)
    if(pT<pTcut) return false; // pT cut should be 20 GeV
    //
    double mvacut=0;
    if( fabs(eta) < 0.8 ) mvacut=0.925;
    if( fabs(eta) >= 0.8 && fabs(eta) < 1.479 ) mvacut=0.975;
    if( fabs(eta) >= 1.479 ) mvacut=0.985;
    if( mvaID_noTrg <= mvacut ) return false;
    //
    if( miss_hits>0 ) return false;
    if( isConv==1 ) return false;
    //
    return true;
  }
  return false;
}

bool selectTau(TString kindsel, double pTcut, double pT, float decayMode, float isoMVAloose, float isoMVAmedium,
	       float discrEleLoose, float discrEleMedium, float discrEleTight, float discrEleMVA,
	       float discrMuLoose, float discrMuTight, float eta, float etacut)
{

  if(kindsel=="etau-probe") {

    return true;
  }
  if(kindsel=="mutau-probe") {
    if(decayMode!=1) return false;
    if(isoMVAloose!=1) return false;
    if(discrMuTight!=1) return false;
    if(discrEleLoose!=1) return false;
    return true;
  }
  if(kindsel=="tautau-probe" || kindsel=="tautau-tag") {
    if(pT<pTcut) return false;
    if(fabs(eta)>etacut) return false;
    if(decayMode!=1) return false;
    if(isoMVAmedium!=1) return false;
    if(discrMuLoose!=1) return false;

    if(kindsel=="tautau-probe")
      if(discrEleMVA!=1) return false;

    if(kindsel=="tautau-tag")
      if(discrEleLoose!=1) return false;

    return true;
  }
  cout << "selectTau() used with bad kindsel !!" << endl;
  return false;
}
