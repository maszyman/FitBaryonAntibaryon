void fitpapvar(const char* infilename, const char* histname) {
  TFile* infile = new TFile(Form("%s",infilename),"read");
  TH1D* hcf = (TH1D*)infile->Get(histname);
  hcf->GetXaxis()->SetRangeUser(0,2);

   TF1 *fc2 = new TF1("fc2","[3]+[2]*TMath::Gaus(x,[0],[1])",0.4,1);
   fc2->SetParameters(1.5,0.3,-0.2,1);

  // TF1 *fc2 = new TF1("fc2","[1]+[0]*x*x",0,1);
  // fc2->SetParameters(-0.01,1.0);

  // TF1 *fc2 = new TF1("fc2","[2]+[1]*x*x*x*x*x",0,1);
  // fc2->SetParameters(-0.01,0.01,1.0);

  // TF1 *fc2 = new TF1("fc2","[3]+[2]*x*x+[1]*x*x*x+[0]*x*x*x*x",0.,1);
  // fc2->SetParameters(-0.01,0.01,1.0,1.0);

  hcf->Fit("fc2","r","",0.3,1);
  TH1D* hnew = new TH1D("hnew","hnew",hcf->GetNbinsX(),0,1);
  for(int i=1;i<=hcf->GetNbinsX();++i){
    // cout << i << endl;
    // cout << hcf->GetBinContent(i)/fc2->Eval(2.*i/hcf->GetNbinsX()) << endl;
    // cout << hcf->GetBinContent(i) << endl;
    // cout << fc2->Eval(2.*i/hcf->GetNbinsX()) << endl << endl;
    hnew->SetBinContent(i, hcf->GetBinContent(i)/fc2->Eval(1.*i/hcf->GetNbinsX()));
    hnew->SetBinError(i, hcf->GetBinError(i));

  }
  hnew->Draw("same");
  fc2->Draw("same");
  hnew->SetName(Form("divp4%s",histname));
  TFile* ofile = new TFile(Form("divp4%s",infilename),"update");
  hnew->Write();
  hcf->Write();
  ofile->Close();
}
