void plotpaptherm(TString file1) {

  TFile* ifile1_ = new TFile(Form("%s",file1.Data()), "read");

  TH1D* hist1;

  TString histname1="cnum1da";
  
  hist1 = (TH1D*)ifile1_->Get(Form("%s",histname1.Data()));

  TCanvas *myCan = new TCanvas("asd","asd");
  myCan->Draw();
  myCan->cd();
  hist1->SetMarkerColor(kBlack);
  hist1->SetMarkerStyle(20);
  hist1->SetMarkerSize(1.6);
  hist1->SetTitle("");
  hist1->GetXaxis()->SetTitle("k* (GeV/c)");
  hist1->GetYaxis()->SetTitle("C(k*)");
  hist1->GetXaxis()->SetRangeUser(0.,0.44);
  gStyle->SetOptStat(0);
   hist1->SetMaximum(1.15);
  // hist1->SetMinimum(0.88);
  hist1->Draw("lp");
  TLatex Tl;
  Tl.SetTextAlign(23);
  Tl.SetTextSize(0.08);
  Tl.SetTextColor(kBlack);
  Tl.DrawLatex(0.2,1.1,"Therminator, b=3fm");

  myCan->SaveAs(Form("paptherm.eps"));
}
