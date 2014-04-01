{

  gROOT->ProcessLine(".L ../../../FitEfficiency/FuncCB_cxx.so");

  gROOT->ProcessLine(".L ../../../FitEfficiency/fitEfficiency_C.so");

  fastEfficiencyNadir( 20 , 0,0,1,0,"../","495 pb",4, kBlack, kFullCircle, kRed, kOpenSquare, "WP80", "WP80", "tree_effi_TagProbe.root");

  gSystem->Exit(0);

}
