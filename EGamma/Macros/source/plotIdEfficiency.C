void loadPresentationStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(0.01);
  gStyle->SetLabelOffset(0.005, "XYZ");
  gStyle->SetTitleSize(0.07, "XYZ");
  gStyle->SetTitleFont(22,"X");
  gStyle->SetTitleFont(22,"Y");
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetHistLineWidth(2);

  gROOT->ProcessLine(".L /home/llr/cms/ndaci/SKWork/FIT/turnonfit/tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
}

void plotIdEfficiency(float size=1)
{

  loadPresentationStyle();
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetTitleYOffset(1.2);

  /*
  const int n=10;
  float adc[n] = {10,15,17,19,21,23,30,35,40,50};
  float eff[n]={76.5 , 50.6 , 96.7 , 97.0 , 96.0 , 98.1 , 
	    98.6,98.2,99.2,98.6};
  */

  const int n=9;
  float adc[n] = {10,17,19,21,23,30,35,40,50};
  float eff[n]={76.5 , 96.7 , 97.0 , 96.0 , 98.1 ,
		98.6,98.2,99.2,98.6};

  for(int i=0 ; i<n ; i++)
    eff[i] /= 100 ;

  TGraph gr(n,adc,eff);
  
  gr.SetMinimum(0.75);
  gr.SetMaximum(1.0);
  gr.GetXaxis()->Set(11,5,55);

  gr.GetXaxis()->SetTitle("sFGVB threshold (adc)");
  gr.GetYaxis()->SetTitle("Spike id efficiency");

  gr.SetMarkerStyle(kFullTriangleUp);
  gr.SetMarkerSize(size);
  gr.SetMarkerColor(kBlue);
  gr.SetFillColor(kWhite);

  TCanvas c("c","c",0,0,800,600);
  gr.Draw("AP");
  
  c.Print("plotIdEfficiency.eps");
  c.Print("plotIdEfficiency.ps");
  c.Print("plotIdEfficiency.gif");
  c.Print("plotIdEfficiency.png");
  c.Print("plotIdEfficiency.pdf");

}
