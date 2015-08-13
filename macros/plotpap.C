void plotpap(TString file1) {

  TFile* ifile1_ = new TFile(Form("%s",file1.Data()), "read");

  TH1D* hist1;
  TH1D* hist2;

  TString histname1="cpapR3.0";
  TString histname2="cpapR4.0";
  
  hist1 = (TH1D*)ifile1_->Get(Form("%s",histname1.Data()));
  hist2 = (TH1D*)ifile1_->Get(Form("%s",histname2.Data()));

  TCanvas *myCan = new TCanvas("asd","asd");
  myCan->Draw();
  myCan->cd();
  hist1->SetMarkerColor(kBlack);
  hist2->SetMarkerColor(kRed);
  hist1->SetMarkerStyle(20);
  hist2->SetMarkerStyle(20);
  hist1->SetMarkerSize(1.6);
  hist2->SetMarkerSize(1.6);
  hist1->SetTitle("");
  hist1->GetXaxis()->SetTitle("k* (GeV/c)");
  hist1->GetYaxis()->SetTitle("C(k*)");
  hist1->GetXaxis()->SetRangeUser(0.,0.44);
  hist2->GetXaxis()->SetRangeUser(0.,0.44);
  gStyle->SetOptStat(0);
   hist1->SetMaximum(1.15);
  // hist1->SetMinimum(0.88);
  hist1->Draw("lp");
  hist2->Draw("lpsame");
  TLatex Tl;
  Tl.SetTextAlign(23);
  Tl.SetTextSize(0.08);
  Tl.SetTextColor(kBlack);
  Tl.DrawLatex(0.2,1.1,"R=3.0 fm");
  Tl.SetTextColor(kRed);
  Tl.DrawLatex(0.2,1.05,"R=4.0 fm");

  myCan->SaveAs(Form("pap.eps"));
}
