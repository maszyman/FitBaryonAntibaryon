void preparehist(TH1D *hist, int ht)
{
  int cols[3] = { kRed, kGreen+2, kBlue };
  int mars[3] = { 24, 20, 25 };

  hist->SetLineColor(cols[ht]);
  hist->SetMarkerColor(cols[ht]);
  hist->SetMarkerStyle(mars[ht]);
  hist->SetMarkerSize(1.1);
  //  hist->GetXaxis()->SetRange(1,149);


  hist->GetXaxis()->SetTitleSize(0.07);
  hist->GetXaxis()->SetLabelSize(0.07);

  hist->GetYaxis()->SetTitleSize(0.07);
  hist->GetYaxis()->SetLabelSize(0.07);

  hist->GetXaxis()->CenterTitle(kTRUE);
  hist->GetYaxis()->CenterTitle(kTRUE);

  hist->GetYaxis()->SetTitleOffset(0.7);

  hist->GetXaxis()->SetNdivisions(4,5,0);
  hist->GetYaxis()->SetNdivisions(4,5,0);

  //  hist->SetTitle(";q_{LCMS} [GeV/c];C(q)");

}

void preparepad()
{
  gPad->SetFillColor(0);
  gPad->SetFillStyle(4000);
  gPad->SetTopMargin(0.005);
  gPad->SetRightMargin(0.005);
  gPad->SetBottomMargin(0.150);
  gPad->SetLeftMargin(0.13);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  
}

void apapalfig2()
{
  gStyle->SetOptStat(1000000000);
  gStyle->SetStatBorderSize(0);
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42,"X");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(0);
  gStyle->SetFrameLineWidth(1);

  TFile *infs[3];

  TFile *infs[0] = new TFile("divp4cfPAP_11hcen0.root");
  TFile *infs[1] = new TFile("divp4cfPAP_11hcen2.root");
  TFile *infs[2] = new TFile("divp4cfPAP_11hcen4.root");
  
  TH1D *hpaps[3];
  
  const char *cents[3] = { "0", "2", "4" };

  for (int iter=0; iter<3; iter++) {
    hpaps[iter] = (TH1D *) infs[iter]->Get(Form("divp4NumOutPckstarPAPtpcM%sPsi3", cents[iter]));
    preparehist(hpaps[iter], iter);
  }

  TH1D *hramka = new TH1D("hramka",";k* (GeV/#it{c});C(k*)",100,0.0,0.5);
  preparehist(hramka, 0);
  hramka->SetMinimum(0.86);
  hramka->SetMaximum(1.26);

  TCanvas *canpapcfsdiv = new TCanvas ("canpapcfsdiv","canpapcfsdiv",1200,800);
  canpapcfsdiv->cd();
  preparepad();

  hramka->Draw();
  hpaps[0]->Draw("SAME");
  hpaps[1]->Draw("SAME");
  hpaps[2]->Draw("SAME");

  TLatex lat;
  lat.SetTextSize(0.065);
  lat.SetTextFont(42);
  lat.DrawLatex(0.05, 1.22, "ALICE Pb-Pb #sqrt{#it{s}_{NN}}=2.76");
  lat.DrawLatex(0.05, 1.18, "p#bar{p} correlation function");
  lat.DrawLatex(0.05, 1.14, "Statistical uncertainty");

  TLegend *leg = new TLegend(0.7,0.6,0.95,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(4000);
  leg->AddEntry(hpaps[0], "0-10\%", "p");
  leg->AddEntry(hpaps[1], "10-30\%", "p");
  leg->AddEntry(hpaps[2], "30-50\%", "p");

  leg->Draw();

}
