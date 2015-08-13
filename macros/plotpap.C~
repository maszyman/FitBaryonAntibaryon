void comparecf(TString file1, TString file2, TString histname1, TString histname2,TString sys) {

  TFile* ifile1_ = new TFile(Form("%s",file1.Data()), "read");
  TFile* ifile2_ = new TFile(Form("%s",file2.Data()), "read");

  TH1D* hist1;
  TH1D* hist2;

  // if ( (strcmp(sys,"V0LL") == 0) || (strcmp(sys,"V0ALAL") == 0) || (strcmp(sys,"V0LAL") == 0)  ) {
  //   hist1 = (TH1D*)ifile1_->Get(Form("Numcqinv%stpcM%d%s",sys,mult,kT));
  //   hist2 = (TH1D*)ifile2_->Get(Form("Numcqinv%stpcM%d%s",sys,mult,kT));
  // }
  // else {
  hist1 = (TH1D*)ifile1_->Get(Form("%s",histname1.Data()));
  hist2 = (TH1D*)ifile2_->Get(Form("%s",histname2.Data()));
  // }

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
  hist1->GetXaxis()->SetRangeUser(0.,0.6);
  hist2->GetXaxis()->SetRangeUser(0.,0.6);
  gStyle->SetOptStat(0);
  hist1->SetMaximum(1.04);
  hist1->SetMinimum(0.88);
  hist1->Draw("p");
  hist2->Draw("psame");
  // TLatex Tl;
  // Tl.SetTextAlign(23);
  // Tl.SetTextSize(0.08);
  // Tl.SetTextColor(kBlack);
  // Tl.DrawLatex(0.2,1.06,"field ++");
  // Tl.SetTextColor(kRed);
  // Tl.DrawLatex(0.19,1.03,"field --");

  myCan->SaveAs(Form("comp%s.eps",sys.Data()));
}
