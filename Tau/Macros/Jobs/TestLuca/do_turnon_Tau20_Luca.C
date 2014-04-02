void do_turnon_Tau20_Luca(int threshold, int iECAL1, int iColl1, int iECAL2, int iColl2)
{
  
  gROOT->ProcessLine(".L ../../FitEfficiency/FuncCB_cxx.so");
  
  gROOT->ProcessLine(".L ../../FitEfficiency/fitEfficiency_L1Tau_C.so");
  
  //[nECAL] = {"Barrel","Endcaps"};
  //[nColl] = {"Current Algo","Current Algo With Iso","New Algo","New Algo With Iso"};
  
  fastEfficiencyNadir(threshold , iECAL1, iColl1, iECAL2, iColl2,"../","18 fb",4, kBlack, kFullCircle, kRed, kOpenSquare, "WP80", "WP80");
  
  gSystem->Exit(0);
  
  return;
}
