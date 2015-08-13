void comparecfback(TString file1, TString file2, TString histname1, TString histname2, TString histname3,TString sys) {

  TFile* ifile1_ = new TFile(Form("%s",file1.Data()), "read");
  TFile* ifile2_ = new TFile(Form("%s",file2.Data()), "read");

  TH1D* hist1;
  TH1D* hist2;
  TH1D* hist3;

  // if ( (strcmp(sys,"V0LL") == 0) || (strcmp(sys,"V0ALAL") == 0) || (strcmp(sys,"V0LAL") == 0)  ) {
  //   hist1 = (TH1D*)ifile1_->Get(Form("Numcqinv%stpcM%d%s",sys,mult,kT));
  //   hist2 = (TH1D*)ifile2_->Get(Form("Numcqinv%stpcM%d%s",sys,mult,kT));
  // }
  // else {
  hist1 = (TH1D*)ifile1_->Get(Form("%s",histname1.Data()));
  hist2 = (TH1D*)ifile2_->Get(Form("%s",histname2.Data()));
  hist3 = (TH1D*)ifile2_->Get(Form("%s",histname3.Data()));
  // }

  TCanvas *myCan = new TCanvas("asd","asd");
  myCan->Draw();
  myCan->cd();
  hist1->SetMarkerColor(kBlack);
  hist2->SetMarkerColor(kRed);
  hist3->SetMarkerColor(kBlue);
  hist1->SetMarkerStyle(20);
  hist2->SetMarkerStyle(20);
  hist3->SetMarkerStyle(20);
  hist1->SetMarkerSize(1.6);
  hist2->SetMarkerSize(1.6);
  hist3->SetMarkerSize(1.6);
  hist1->SetTitle("");
  hist1->GetXaxis()->SetRangeUser(0.,1.0);
  hist2->GetXaxis()->SetRangeUser(0.,1.0);
  hist3->GetXaxis()->SetRangeUser(0.,1.0);
  gStyle->SetOptStat(0);

  //hist1->Scale(1./1.026);//1.026, 1.05

    hist1->SetMaximum(1.015);
  hist1->SetMinimum(0.9425);
  hist1->GetXaxis()->SetTitleSize(0.05);
  hist1->GetYaxis()->SetTitleSize(0.05);
  hist1->Draw("p");
  hist2->Draw("psame");
  hist3->Draw("psame");
  TLatex Tl;
  Tl.SetTextAlign(23);
  Tl.SetTextSize(0.08);
  Tl.SetTextColor(kBlack);
  Tl.DrawLatex(0.5,0.995,"p#bar{p}");
  Tl.SetTextColor(kRed);
  Tl.DrawLatex(0.5,0.980,"p#bar{#Lambda}");
  Tl.SetTextColor(kBlue);
  Tl.DrawLatex(0.5,0.965,"#bar{p}#Lambda");

  myCan->SaveAs(Form("comp%s.eps",sys.Data()));
}
