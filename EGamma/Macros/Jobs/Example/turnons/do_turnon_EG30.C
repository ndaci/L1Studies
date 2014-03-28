{

  gROOT->ProcessLine(".L /home/llr/cms/ndaci/SKWork/macro/skEfficiency/analyzers/tagAndProbe/FastEfficiency/FuncCB_428_cxx.so");

  gROOT->ProcessLine(".L /home/llr/cms/ndaci/SKWork/macro/skEfficiency/analyzers/tagAndProbe/FastEfficiency/fastEfficiency_2012_style.1.0_C.so");

  fastEfficiencyNadir( 30 , 0,0,1,0,"/home/llr/cms/ndaci/SKWork/macro/skEfficiency/tagAndProbe/EfficiencyStudy/Study2012/Run2012C_PRV2/","6.242 fb",4,
		       kBlack, kFullCircle, kRed, kOpenSquare, "WP80", "WP80");

  gSystem->Exit(0);

}
