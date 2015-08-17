void preparehist(TH1D *hist, int ht)
{
  int cols[6] = { kRed, kGreen+2, kBlue, kCyan+1, kViolet+1, kOrange+7 };
  int mars[6] = { 24, 20, 25, 21, 29, 27 };

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

void apapalfig3()
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

  TFile *infs[0] = new TFile("divgcfaplam.root");
  // TFile *infs[1] = new TFile("divp4cfPAP_11hcen2.root");
  // TFile *infs[2] = new TFile("divp4cfPAP_11hcen4.root");
  
  TH1D *hapls[3];
  TH1D *hpals[3];
  
  const char *cents[3] = { "0", "2", "4" };

  for (int iter=0; iter<3; iter++) {
    hpals[iter] = (TH1D *) infs[0]->Get(Form("divgNumOutPcnonidV0PALtpcM%s", cents[iter]));
    preparehist(hpals[iter], iter);
    hapls[iter] = (TH1D *) infs[0]->Get(Form("divgNumOutPcnonidV0APLtpcM%s", cents[iter]));
    preparehist(hapls[iter], iter+3);
  }

  TH1D *hramka = new TH1D("hramka",";k* (GeV/#it{c});C(k*)",100,0.0,0.5);
  preparehist(hramka, 0);
  hramka->SetMinimum(0.66);
  hramka->SetMaximum(1.03);

  TCanvas *canpalcfsdiv = new TCanvas ("canpalcfsdiv","canpalcfsdiv",1200,800);
  canpalcfsdiv->cd();
  preparepad();

  hramka->Draw();
  hpals[0]->Draw("SAME");
  hpals[1]->Draw("SAME");
  hpals[2]->Draw("SAME");

  hapls[0]->Draw("SAME");
  hapls[1]->Draw("SAME");
  hapls[2]->Draw("SAME");

  TLatex lat;
  lat.SetTextSize(0.065);
  lat.SetTextFont(42);
  lat.DrawLatex(0.15, 0.76, "ALICE Pb-Pb #sqrt{#it{s}_{NN}}=2.76 TeV/#it{c}");
  lat.DrawLatex(0.15, 0.72, "p#bar{#Lambda} and #bar{p}#Lambda CF");
  lat.DrawLatex(0.15, 0.68, "Statistical uncertainty");

  TLegend *leg = new TLegend(0.5,0.5,0.73,0.8,"     p#bar{#Lambda}");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(4000);
  
  leg->AddEntry(hpals[0], "0-10\%", "p");
  leg->AddEntry(hpals[1], "10-30\%", "p");
  leg->AddEntry(hpals[2], "30-50\%", "p");

  leg->Draw();

  TLegend *leg2 = new TLegend(0.75,0.5,0.98,0.8,"     #bar{p}#Lambda");
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(4000);
  
  leg2->AddEntry(hapls[0], "0-10\%", "p");
  leg2->AddEntry(hapls[1], "10-30\%", "p");
  leg2->AddEntry(hapls[2], "30-50\%", "p");

  leg2->Draw();

}
