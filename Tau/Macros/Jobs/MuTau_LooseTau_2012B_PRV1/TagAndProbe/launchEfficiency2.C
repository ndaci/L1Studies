{

  gROOT->ProcessLine(".L /home/llr/cms/ndaci/SKWork/macro/HTauTau/analyzers/tagAndProbe/Efficiency/FuncCB_428_cxx.so");

  gROOT->ProcessLine(".L /home/llr/cms/ndaci/SKWork/macro/HTauTau/analyzers/tagAndProbe/Efficiency/fitEfficiency.1.1.1_C.so");

  float thres=20, iECAL1=0, iColl1=0, iECAL2=0, iColl2=0;
  int nCPU=4;
  TString lumi="200 pb" ;
  TString cutdef="";
  TString cut_choice="HLT";
  float mLow=45, mHigh=70;

  //TString dirIn = "/data_CMS/cms/ndaci/ndaci_2012/HTauTau/TriggerStudy/SingleMu/MuMu/Run2012A_PRV1/Pairs/" ;
  TString dirIn = "/data_CMS/cms/ndaci/ndaci_2012/HTauTau/TriggerStudy/SingleMu/MuTau/Run2012B_PRV1_2/PairsLooseTau/Update/";
  TString dirResults = "/home/llr/cms/ndaci/SKWork/macro/HTauTau/results/TriggerStudies/MuTau/LooseTau/Run2012B_PRV1_2/Test2/" ;
  
  fitEfficiency( thres,  iECAL1,  iColl1,  iECAL2,  iColl2,
		 mLow, mHigh,
		 cutdef, cut_choice,
		 dirIn, dirResults, lumi,  nCPU, 
		 kBlack, kFullCircle, kRed, kOpenSquare,
		 "",  "", "/*.root", "eff_HLT_MuTau_2012B_PRV1_LooseTau");

  //gROOT->ProcessLine(".q");
}
