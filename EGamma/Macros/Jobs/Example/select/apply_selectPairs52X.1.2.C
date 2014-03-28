{

  gSystem->Load("/home/llr/cms/ndaci/SKWork/macro/skEfficiency/analyzers/tagAndProbe/SelectPairs/selectPairs52X.1.2_C.so");

  selectPairs(-1,0,0,"L1_EG_MENU_2011",
	      "/home/llr/cms/ndaci/SKWork/macro/skEfficiency/tagAndProbe/EfficiencyStudy/Study2012/Run2012D_PRV1_up/",
	      "/data_CMS/cms/ndaci/ndaci_2012/SingleElectron/Update/Run2012D_PRV1_up/Elepairs/",
	      "WP80", 20, "WP80", 5,
	      60, 120, "noHLT", "ElePairs",
	      true, true, false) ;
  
  gSystem->Exit(0);

  /*
void selectPairs(int nEntries=-1, int iEff=0, int iEff_M=0, TString eg_menu="L1_EG_MENU_2011",
		 TString dirOut="",
		 TString dirIn="",
		 TString tagCuts="WP80", double tagEt=20, TString probeCuts="WP95", double probeEt=5, 
		 double massCut1=60, double massCut2=120, TString tag_hlt="noHLT", TString nameChain="ElePairs",
		 bool checkCharge=true, bool data2011=true, bool debug=false)

  */

}
