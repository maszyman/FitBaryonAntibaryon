#include <TH1D.h>
#include <TVirtualFitter.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH2D.h>
#include <TVectorD.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

double f0re = -1.0/0.197327;
double f0im = 2.2/0.197327;
double d0re = 1.0/0.197327;
double f0ia = 1.0/0.197327;

double radius0 = 3.2;
double radius1 = 3.0;
double radius2 = 2.8;
double radius3 = 2.5;
double radius4 = 2.5;
double radius5 = 2.5;

double lampl = 0.70;
double lamps = 0.30;
double lampp = 0.00;
double lamll = 0.00;
double lamls = 0.00;
double lamsl = 0.00;
double lamss = 0.00;
double lampx = 0.00;
double lamlx = 0.00;
double lamsx = 0.00;

double lamapl = 0.70;
double lamaps = 0.30;
double lamapp = 0.00;
double lamall = 0.00;
double lamals = 0.00;
double lamasl = 0.00;
double lamass = 0.00;
double lamapx = 0.00;
double lamalx = 0.00;
double lamasx = 0.00;

// double radiuspap = 3.5;
// double lambdapap = 0.7;
// double lambdapal = 0.3;

double radiuspap[3] = {3.5,3.0,2.5};
double lambdapap = 0.52;
double lambdapal = 0.3;

double kstarmin = 0.015;
double kstarmax = 0.495;
double kstarstep = 0.01;
double stepR = 0.1;
double minR = 1.0;
double maxR = 6.0;

double bkgmean = 1.0;
double bkgwidth = 0.4;
double bkgscale = -0.001;
double bkgoffset = 1.0;

double norm0;
double norm1;
double norm2;
double norm3;

double radius = 3.5;
double radiuspp = 0.0;

int binbeg = 1;
int binend = 35;

int NumOfSystems; // 1 - 6

const int npar = 43;
const int nopts = 6;

TFile *files[7];
TH1D *prims[6];
TH1D *cfppdata[3];
TGraph *fitpap[3];
TGraph *fitpapmr[3];
TF1* fc2;
TGraph* fpap;
TGraph* fpal;
TGraph* CFapLakpapVal;
TGraph* CFapLakaplVal;

int dochimap = 0;
int dofscale = 0;
int domomres = 0;
int dodivg = 0;

TH2D *htr;
TH2D *htrpp;
TH2D *htrll;
TH2D *htrsl;
TH2D *htrls;
TH2D *htrss;
TH2D *htrpx;
TH2D *htrlx;
TH2D *htrsx;

TVirtualFitter *fitter;

Int_t dofix[npar];

Double_t parsg[npar];
Double_t parmin[npar];
Double_t parmax[npar];
Double_t opts[nopts];

TGraph *grreout;
TGraph *grimout;
TGraph *gralout;

TGraph *graf1;
TGraph *graf2;
TGraph *graf3;
TGraph *graf4;
TGraph *graf5;
TGraph *graf6;
TGraph *graf1mr;
TGraph *graf2mr;
TGraph *graf3mr;
TGraph *graf4mr;
TGraph *graf5mr;
TGraph *graf6mr;

TGraph *grpf1;
TGraph *grpf2;
TGraph *grpf3;
TGraph *grpf4;

TGraph *grrf1;
TGraph *grrf2;
TGraph *grrf3;
TGraph *grrf4;

TGraph *gruf1;
TGraph *gruf2;
TGraph *gruf3;
TGraph *gruf4;

TGraph *grlf1;
TGraph *grlf2;
TGraph *grlf3;
TGraph *grlf4;

TGraph *grsf1;
TGraph *grsf2;
TGraph *grsf3;
TGraph *grsf4;

TGraph *gropl;
TGraph *gropp;
TGraph *grops;
TGraph *groll;
TGraph *grols;
TGraph *grosl;
TGraph *gross;
TGraph *gropx;
TGraph *grolx;
TGraph *grosx;


TFile *infpp;
TFile *fileInCFppkpp;
TFile *fileInCFplakpp;

double intexp(double z)
{
  double val = 0;
  double vx;
  for (int ix=0; ix<1000; ix++) {
    vx = (ix+0.5)*z/1000.0;

    val += TMath::Exp(+vx*vx);

  }
  val /= 1000;
  val *= z;

  return val;
}

TH1D* CFpapkpap (Double_t R)
{
  TH1D* cfpapkpap = (TH1D*)fileInCFppkpp->Get(Form("cpapR%.1f",R));
  return cfpapkpap;
}

TGraphErrors* CFppkpp (Double_t R)
{
  TGraphErrors* cfppkpp = (TGraphErrors*)fileInCFppkpp->Get(Form("cppR%.1f",R));
  return cfppkpp;
}

TGraphErrors* CFplakpp (Double_t R)
{
  TGraphErrors* cfplakpp = (TGraphErrors*)fileInCFplakpp->Get(Form("grlamR%.1f",R));
  return cfplakpp;
}

Double_t getCFppVal (Double_t R, Double_t kstar)
{
  Double_t Ra = TMath::Floor(R*10.0)/10.0-stepR;
  Double_t Rb = TMath::Floor(R*10.0)/10.0;
  Double_t Rc = TMath::Floor(R*10.0)/10.0+stepR;
  Double_t Rd = TMath::Floor(R*10.0)/10.0+2*stepR;

    Double_t eps = 0.000001;
    if(TMath::Abs((Ra-minR)/stepR - 1.0*TMath::Nint((Ra-minR)/stepR)) > eps) {
      Ra += 0.1;
      Rb = Ra + stepR;
      Rc = Ra + 2*stepR;
      Rd = Ra + 3*stepR;
    }

    //if (Ra == Rc) Rc += 2*stepR;

    if (Ra < minR) {
        Ra = minR;
        Rb = minR + stepR;
        Rc = minR + 2*stepR;
	Rd = minR + 3*stepR;

    }
    else if (Rd > maxR) {
        Ra = maxR - 3*stepR;
        Rb = maxR - 2*stepR;
        Rc = maxR - stepR;
        Rd = maxR;
    }

    // return (CFppkpp(Ra)->Eval(kstar)) * ( (R-Rb)*(R-Rc) ) / ( (Ra-Rb)*(Ra-Rc) )
    //   +  (CFppkpp(Rb)->Eval(kstar)) * ( (R-Ra)*(R-Rc) ) / ( (Rb-Ra)*(Rb-Rc) )
    //   +  (CFppkpp(Rc)->Eval(kstar)) * ( (R-Ra)*(R-Rb) ) / ( (Rc-Ra)*(Rc-Rb) );
    // return 0.5*(CFppkpp(Ra)->Eval(kstar)+CFppkpp(Ra)->Eval(-1.0*kstar)) * ( (R-Rb)*(R-Rc) ) / ( (Ra-Rb)*(Ra-Rc) )
    //   +  0.5*(CFppkpp(Rb)->Eval(kstar)+CFppkpp(Rb)->Eval(-1.0*kstar)) * ( (R-Ra)*(R-Rc) ) / ( (Rb-Ra)*(Rb-Rc) )
    //   +  0.5*(CFppkpp(Rc)->Eval(kstar)+CFppkpp(Rc)->Eval(-1.0*kstar)) * ( (R-Ra)*(R-Rb) ) / ( (Rc-Ra)*(Rc-Rb) );

    // new pap functions
    TH1D* cfra = CFpapkpap(Ra);
    TH1D* cfrb = CFpapkpap(Rb);
    TH1D* cfrc = CFpapkpap(Rc);
    
    // return (cfra->GetBinContent(cfra->FindFixBin(kstar))) * ( (R-Rb)*(R-Rc) ) / ( (Ra-Rb)*(Ra-Rc) )
    //   +  (cfrb->GetBinContent(cfrb->FindFixBin(kstar))) * ( (R-Ra)*(R-Rc) ) / ( (Rb-Ra)*(Rb-Rc) )
    //   +  (cfrc->GetBinContent(cfrc->FindFixBin(kstar))) * ( (R-Ra)*(R-Rb) ) / ( (Rc-Ra)*(Rc-Rb) );

    // 3rd order interpolation
    TH1D* cfrd = CFpapkpap(Rd);

    // if (kstar < 0.02) {
    // cout << R << " " << Ra << " " << Rb << " " << Rc << " " <<Rd <<endl;
    // cout << (cfra->GetBinContent(cfra->FindFixBin(kstar))) * ( (R-Rb)*(R-Rc)*(R-Rd) ) / ( (Ra-Rb)*(Ra-Rc)*(Ra-Rd) )
    //   +  (cfrb->GetBinContent(cfrb->FindFixBin(kstar))) * ( (R-Ra)*(R-Rc)*(R-Rd) ) / ( (Rb-Ra)*(Rb-Rc)*(Rb-Rd) )
    //   +  (cfrc->GetBinContent(cfrc->FindFixBin(kstar))) * ( (R-Ra)*(R-Rb)*(R-Rd) ) / ( (Rc-Ra)*(Rc-Rb)*(Rc-Rd) )
    //   +  (cfrd->GetBinContent(cfrd->FindFixBin(kstar))) * ( (R-Ra)*(R-Rb)*(R-Rc) ) / ( (Rd-Ra)*(Rd-Rb)*(Rd-Rc) ) << " " <<  cfra->GetBinContent(cfra->FindFixBin(kstar)) << " " << cfrb->GetBinContent(cfrb->FindFixBin(kstar)) << " " << cfrc->GetBinContent(cfrc->FindFixBin(kstar)) << " " << cfrd->GetBinContent(cfrd->FindFixBin(kstar)) << " " << endl;
    // }
    return (cfra->GetBinContent(cfra->FindFixBin(kstar))) * ( (R-Rb)*(R-Rc)*(R-Rd) ) / ( (Ra-Rb)*(Ra-Rc)*(Ra-Rd) )
      +  (cfrb->GetBinContent(cfrb->FindFixBin(kstar))) * ( (R-Ra)*(R-Rc)*(R-Rd) ) / ( (Rb-Ra)*(Rb-Rc)*(Rb-Rd) )
      +  (cfrc->GetBinContent(cfrc->FindFixBin(kstar))) * ( (R-Ra)*(R-Rb)*(R-Rd) ) / ( (Rc-Ra)*(Rc-Rb)*(Rc-Rd) )
      +  (cfrd->GetBinContent(cfrd->FindFixBin(kstar))) * ( (R-Ra)*(R-Rb)*(R-Rc) ) / ( (Rd-Ra)*(Rd-Rb)*(Rd-Rc) );

}

Double_t getCFpLaVal (Double_t R, Double_t kstar)
{
    Double_t Ra = TMath::Floor(R*10.0)/10.0;
    Double_t Rb = TMath::Floor(R*10.0)/10.0+stepR;
    Double_t Rc = TMath::Floor(R*10.0)/10.0+2*stepR;

    if (Ra == Rc) Rc += 2*stepR;

    if (Ra < minR) {
        Ra = minR;
        Rb = minR + stepR;
        Rc = minR + 2*stepR;
    }
    else if (Rc > maxR) {
        Ra = maxR - 2*stepR;
        Rb = maxR - stepR;
        Rc = maxR;
    }

    return CFplakpp(Ra)->Eval(kstar) * ( (R-Rb)*(R-Rc) ) / ( (Ra-Rb)*(Ra-Rc) )
        + CFplakpp(Rb)->Eval(kstar) * ( (R-Ra)*(R-Rc) ) / ( (Rb-Ra)*(Rb-Rc) )
        + CFplakpp(Rc)->Eval(kstar) * ( (R-Ra)*(R-Rb) ) / ( (Rc-Ra)*(Rc-Rb) );

}
TGraph *getgrk(double R=radius)
{
  double lam = 1.0;

  TF1 *funf2 = new TF1("funf2","(1-exp(-[0]*[0]*x*x))/([0]*x)");
  funf2->SetParameter(0,2*radius/0.197327);
  funf2->SetRange(0.0,0.5);

  const int npoints = 100;
  double xy[npoints];
  double yk[npoints];

  double ykre[npoints];
  double ykim[npoints];
  double ykal[npoints];

  double rad = radius/0.197327;

  double fzm, fzr, fzi, d0v;
  double fsd, fsr, fsi, fsm;

  for (int ix=0; ix<npoints; ix++) {
    xy[ix] = (ix+0.5)*0.005;
    yk[ix] = intexp(xy[ix])*exp(-xy[ix]*xy[ix])/xy[ix];

    d0v = 0.5*d0re*xy[ix]*xy[ix];

    fzm = (f0re*f0re + f0im*f0im);
    fzr = f0re/fzm;
    fzi = -f0im/fzm;

    //    cout << fzr << " " << d0v << endl;

    fzr += d0v;
    fzi -= xy[ix];

    fsd = (fzr*fzr + fzi*fzi);
    fsr = fzr/fsd;
    fsi = -fzi/fsd;
    fsm = (fsr*fsr + fsi*fsi);

    //    cout << xy[ix] << " " << fsm << " " << fsr << " " << fsi << endl;

    ykre[ix] = 1.0 +
      1.0 * (0.5 * (fsm/(rad*rad))*(1-d0re/(2*TMath::Sqrt(TMath::Pi())*rad)) +
             (2*(fsr)/(TMath::Sqrt(TMath::Pi())*rad)*(intexp(2*xy[ix]*rad)*exp(-2*xy[ix]*rad*2*xy[ix]*rad)/(2*xy[ix]*rad))));
    ykim[ix] = 1.0 +
      1.0 * (-(fsi)*funf2->Eval(xy[ix])/rad);
    ykal[ix] = 1.0 +
      1.0 * (0.5 * (fsm/(rad*rad))*(1-d0re/(2*TMath::Sqrt(TMath::Pi())*rad)) +
             (2*(fsr)/(TMath::Sqrt(TMath::Pi())*rad)*(intexp(2*xy[ix]*rad)*exp(-2*xy[ix]*rad*2*xy[ix]*rad)/(2*xy[ix]*rad))) -
             (fsi)*funf2->Eval(xy[ix])/rad);

    //    cout << ykre[ix] << " " << ykim[ix] << " " << ykal[ix] << endl;

  }

//   double rad = radius/0.197327;
//   double den, f0rk, f0ik, f0md;

//   for (int ix=0; ix<100; ix++) {
//     xy[ix] = (ix+0.5)*0.005;
//     yk[ix] = intexp(xy[ix])*exp(-xy[ix]*xy[ix])/xy[ix];

//     den = (1+f0im*xy[ix])*(1+f0im*xy[ix]) + f0re*f0re*xy[ix]*xy[ix];
//     f0rk = f0re;
//     f0ik = f0im + xy[ix]*(f0re*f0re + f0im*f0im);
//     f0md = (f0rk*f0rk + f0ik*f0ik)/(den*den);

//     //    cout << xy[ix] << " " << f0md << " " << f0rk/den << " " << f0ik/den << endl;

//     ykre[ix] = 1.0 +
//       1.0 * (0.5 * (f0md/(rad*rad))*(1-d0re/(2*TMath::Sqrt(TMath::Pi())*rad)) +
// 	     (2*(f0rk/den)/(TMath::Sqrt(TMath::Pi())*rad)*(intexp(2*xy[ix]*rad)*exp(-2*xy[ix]*rad*2*xy[ix]*rad)/(2*xy[ix]*rad))));
//     ykim[ix] = 1.0 +
//       1.0 * (-(f0ik/den)*funf2->Eval(xy[ix])/rad);
//     ykal[ix] = 1.0 +
//       1.0 * (0.5 * (f0md/(rad*rad))*(1-d0re/(2*TMath::Sqrt(TMath::Pi())*rad)) +
// 	     (2*(f0rk/den)/(TMath::Sqrt(TMath::Pi())*rad)*(intexp(2*xy[ix]*rad)*exp(-2*xy[ix]*rad*2*xy[ix]*rad)/(2*xy[ix]*rad))) -
// 	     (f0ik/den)*funf2->Eval(xy[ix])/rad);

// //     den = (1+xy[ix]*f0im)*(1+xy[ix]*f0im) + xy[ix]*xy[ix]*f0re*f0re;
// //     f0rk = f0re;
// //     f0ik = f0im + xy[ix]*(f0re*f0re+f0im*f0im);
// //     f0md = (f0rk*f0rk + f0ik*f0ik)/(den*den);

// //     ykre[ix] = 1.0 + (0.5 * f0md / (rad*rad) +
// // 		      2*(f0re/den)*intexp(2*xy[ix]*rad)*exp(-2*xy[ix]*rad*2*xy[ix]*rad)/(2*xy[ix]*rad)/(TMath::Sqrt(TMath::Pi())*rad));
// //     ykim[ix] = 1.0 -(f0ik/den)*funf2->Eval(xy[ix])/(rad);
// //     ykal[ix] = 1 + lam*(ykim[ix] + ykre[ix] - 2.0);
//   }

  TGraph *grre = new TGraph(npoints, xy, ykre);
  TGraph *grim = new TGraph(npoints, xy, ykim);
  TGraph *gral = new TGraph(npoints, xy, ykal);

  grreout = grre;
  grimout = grim;
  gralout = gral;

//   grre->SetMaximum(1.2);
//   grre->SetMinimum(0.0);

//   grre->Draw("P");
//   grim->SetLineColor(2);
//   grim->SetMarkerColor(2);
//   grim->Draw("P");
//   gral->SetLineColor(4);
//   gral->SetMarkerColor(4);
//   gral->Draw("P");

  return gral;
}

TGraph* getCFapLakpapVal (Double_t radius)
{

  const int npoints = 100;

  TGraph *cfth = getgrk(radius);
  double xy[npoints];
  double yk[npoints];

  int bcy, bcx;
  double val, valsum, wgtsum;
  double wgt;

  for (int ibpl=0; ibpl<npoints; ibpl++) {
    xy[ibpl] = cfth->GetX()[ibpl];

    bcx = (ibpl+1);

    wgtsum = 0.0;
    valsum = 0.0;

    for (int ibpp = 0; ibpp<npoints; ibpp++) {

      bcy = ibpp+1;

      val = cfth->GetY()[ibpp];
      wgt = htrpp->GetBinContent(bcx, bcy);

      wgtsum += wgt;
      valsum += val*wgt;
    }
    yk[ibpl] = valsum/wgtsum;
  }

  TGraph * grres = new TGraph(npoints, xy, yk);

  return grres;

}

void GetPar(ifstream *inf, Double_t *parval, Int_t *isfixed, Double_t *parmin, Double_t *parmax)
{
  char linebuf[1000];
  char buf[100];
  memset(buf, 0, 100);

  inf->getline(linebuf, 1000);
  istringstream *istr = new istringstream(linebuf);
  (*istr) >> (*parval);
  (*istr) >> buf;

  cout << "Read value |" << *parval << "| |" << buf << "|";
  if (strstr(buf, "F")) {
    *isfixed = 1;
    cout << " Fixed" << endl;
  }
  else if (strstr(buf, "L")) {
    *isfixed = 2;
    cout << " Limits ";
    (*istr) >> (*parmin) >> (*parmax);
    cout << "|" << *parmin << "| |" << (*parmax) << "|" << endl;
  }
  else {
    *isfixed = 0;
    cout << endl;
  }

  delete istr;
}

void GetOpts(ifstream *inf, Double_t *parval)
{
  char linebuf[1000];
  char buf[100];
  memset(buf, 0, 100);

  inf->getline(linebuf, 1000,  '#');
  cout << "linebuf = " << linebuf << endl;
  istringstream *istr = new istringstream(linebuf);
  (*istr) >> (*parval);
  (*istr) >> buf;
  cout << "parval = " << *parval << endl;

  delete istr;
}

void rescaleGraph(TGraph *grp, double pur)
{
  for (int ix=0; ix<grp->GetN(); ix++) {
    grp->GetY()[ix] /= pur;
  }
}



TGraph *getfsq()
{
  double lam = 1.0;

//   TF1 *funf2 = new TF1("funf2","(1-exp(-[0]*[0]*x*x))/([0]*x)");
//   funf2->SetParameter(0,2*radius/0.197327);
//   funf2->SetRange(0.0,0.5);

  double xy[100];
  double yk[100];

//   double ykre[100];
//   double ykim[100];
//   double ykal[100];

  double rad = radius/0.197327;

  double fzm, fzr, fzi, d0v;
  double fsd, fsr, fsi, fsm;

  for (int ix=0; ix<100; ix++) {
    xy[ix] = (ix+0.5)*0.005;

    d0v = 0.5*d0re*xy[ix]*xy[ix];

    fzm = (f0re*f0re + f0im*f0im);
    fzr = f0re/fzm;
    fzi = -f0im/fzm;

    //    cout << fzr << " " << d0v << endl;

    fzr += d0v;
    fzi -= xy[ix];

    fsd = (fzr*fzr + fzi*fzi);
    fsr = fzr/fsd;
    fsi = -fzi/fsd;
    fsm = (fsr*fsr + fsi*fsi);


    yk[ix] = 4*TMath::Pi()*fsm;


  }

  TGraph *grre = new TGraph(100, xy, yk);

  return grre;
}

TGraph* getresidualkpLLL()
{

//   TFile *infpp = new TFile("CFplam.root");
//   TGraphErrors *cfth = infpp->Get("cplR3.5");

//  f0im = f0ia*0.48;
  TGraph *cfth = getgrk();
  //  f0im = f0ia;

  double xy[100];
  double yk[100];

  int bcy, bcx;
  double val, valsum, wgtsum;
  double wgt;

  for (int ibpl=0; ibpl<100; ibpl++) {
    xy[ibpl] = cfth->GetX()[ibpl];

    bcx = (ibpl+1);

    wgtsum = 0.0;
    valsum = 0.0;

    for (int ibpp = 0; ibpp<100; ibpp++) {

      bcy = ibpp+1;

      val = cfth->GetY()[ibpp];
      wgt = htrll->GetBinContent(bcx, bcy);
      //      cout << ibpl << " " << ibpp << "    " <<val << " " << wgt << endl;

      wgtsum += wgt;
      valsum += val*wgt;
    }

    yk[ibpl] = valsum/wgtsum;

  }

  TGraph * grres = new TGraph(100, xy, yk);

  return grres;
}

TGraph* getresidualkpLLS()
{
  //  f0im = f0ia*1.5;
  TGraph *cfth = getgrk();
  //  f0im = f0ia;

  double xy[100];
  double yk[100];

  int bcy, bcx;
  double val, valsum, wgtsum;
  double wgt;

  for (int ibpl=0; ibpl<100; ibpl++) {
    xy[ibpl] = cfth->GetX()[ibpl];

    bcx = (ibpl+1);

    wgtsum = 0.0;
    valsum = 0.0;

    for (int ibpp = 0; ibpp<100; ibpp++) {

      bcy = ibpp+1;

      val = cfth->GetY()[ibpp];
      wgt = htrls->GetBinContent(bcx, bcy);
      //      cout << ibpl << " " << ibpp << "    " <<val << " " << wgt << endl;

      wgtsum += wgt;
      valsum += val*wgt;
    }

    yk[ibpl] = valsum/wgtsum;

  }

  TGraph * grres = new TGraph(100, xy, yk);

  return grres;
}

TGraph* getresidualkpLSL()
{
  //  f0im = f0ia;
  TGraph *cfth = getgrk();

  double xy[100];
  double yk[100];

  int bcy, bcx;
  double val, valsum, wgtsum;
  double wgt;

  for (int ibpl=0; ibpl<100; ibpl++) {
    xy[ibpl] = cfth->GetX()[ibpl];

    bcx = (ibpl+1);

    wgtsum = 0.0;
    valsum = 0.0;

    for (int ibpp = 0; ibpp<100; ibpp++) {

      bcy = ibpp+1;

      val = cfth->GetY()[ibpp];
      wgt = htrsl->GetBinContent(bcx, bcy);
      //      cout << ibpl << " " << ibpp << "    " <<val << " " << wgt << endl;

      wgtsum += wgt;
      valsum += val*wgt;
    }

    yk[ibpl] = valsum/wgtsum;

  }

  TGraph * grres = new TGraph(100, xy, yk);

  return grres;
}

TGraph* getresidualkpLpX()
{
  //  f0im = f0ia*0.5;
  TGraph *cfth = getgrk();
  //  f0im = f0ia;

  double xy[100];
  double yk[100];

  int bcy, bcx;
  double val, valsum, wgtsum;
  double wgt;

  for (int ibpl=0; ibpl<100; ibpl++) {
    xy[ibpl] = cfth->GetX()[ibpl];

    bcx = (ibpl+1);

    wgtsum = 0.0;
    valsum = 0.0;

    for (int ibpp = 0; ibpp<100; ibpp++) {

      bcy = ibpp+1;

      val = cfth->GetY()[ibpp];
      wgt = htrpx->GetBinContent(bcx, bcy);
      //      cout << ibpl << " " << ibpp << "    " <<val << " " << wgt << endl;

      wgtsum += wgt;
      valsum += val*wgt;
    }

    yk[ibpl] = valsum/wgtsum;

  }

  TGraph * grres = new TGraph(100, xy, yk);

  return grres;
}

TGraph* getresidualkpLLX()
{
  //  f0im = f0ia;
  TGraph *cfth = getgrk();

  double xy[100];
  double yk[100];

  int bcy, bcx;
  double val, valsum, wgtsum;
  double wgt;

  for (int ibpl=0; ibpl<100; ibpl++) {
    xy[ibpl] = cfth->GetX()[ibpl];

    bcx = (ibpl+1);

    wgtsum = 0.0;
    valsum = 0.0;

    for (int ibpp = 0; ibpp<100; ibpp++) {

      bcy = ibpp+1;

      val = cfth->GetY()[ibpp];
      wgt = htrlx->GetBinContent(bcx, bcy);
      //      cout << ibpl << " " << ibpp << "    " <<val << " " << wgt << endl;

      wgtsum += wgt;
      valsum += val*wgt;
    }

    yk[ibpl] = valsum/wgtsum;

  }

  TGraph * grres = new TGraph(100, xy, yk);

  return grres;
}

TGraph* getresidualkpLSX()
{
  //  f0im = f0ia;
  TGraph *cfth = getgrk();

  double xy[100];
  double yk[100];

  int bcy, bcx;
  double val, valsum, wgtsum;
  double wgt;

  for (int ibpl=0; ibpl<100; ibpl++) {
    xy[ibpl] = cfth->GetX()[ibpl];

    bcx = (ibpl+1);

    wgtsum = 0.0;
    valsum = 0.0;

    for (int ibpp = 0; ibpp<100; ibpp++) {

      bcy = ibpp+1;

      val = cfth->GetY()[ibpp];
      wgt = htrsx->GetBinContent(bcx, bcy);
      //      cout << ibpl << " " << ibpp << "    " <<val << " " << wgt << endl;

      wgtsum += wgt;
      valsum += val*wgt;
    }

    yk[ibpl] = valsum/wgtsum;

  }

  TGraph * grres = new TGraph(100, xy, yk);

  return grres;
}

TGraph* getresidualkpLSS()
{
  //  f0im = f0ia;
  TGraph *cfth = getgrk();

  double xy[100];
  double yk[100];

  int bcy, bcx;
  double val, valsum, wgtsum;
  double wgt;

  for (int ibpl=0; ibpl<100; ibpl++) {
    xy[ibpl] = cfth->GetX()[ibpl];

    bcx = (ibpl+1);

    wgtsum = 0.0;
    valsum = 0.0;

    for (int ibpp = 0; ibpp<100; ibpp++) {

      bcy = ibpp+1;

      val = cfth->GetY()[ibpp];
      wgt = htrss->GetBinContent(bcx, bcy);
      //      cout << ibpl << " " << ibpp << "    " <<val << " " << wgt << endl;

      wgtsum += wgt;
      valsum += val*wgt;
    }

    yk[ibpl] = valsum/wgtsum;

  }

  TGraph * grres = new TGraph(100, xy, yk);

  return grres;
}


TGraph* getresidualkpLpS()
{

//   TFile *infpp = new TFile("CFplam.root");
//   TGraphErrors *cfth = infpp->Get("cplR3.5");

//  f0im = f0ia;
  TGraph *cfth = getgrk();

  double xy[100];
  double yk[100];

  int bcy, bcx;
  double val, valsum, wgtsum;
  double wgt;

  for (int ibpl=0; ibpl<100; ibpl++) {
    xy[ibpl] = cfth->GetX()[ibpl];

    bcy = (ibpl+1);

    wgtsum = 0.0;
    valsum = 0.0;

    for (int ibpp = 0; ibpp<100; ibpp++) {

      bcx = ibpp+1;

      val = cfth->GetY()[ibpp];
      wgt = htr->GetBinContent(bcx, bcy);
      //      cout << ibpl << " " << ibpp << "    " <<val << " " << wgt << endl;

      wgtsum += wgt;
      valsum += val*wgt;

//       if (ibpl == 6) {
// 	cout << bcx << " " << bcy << " " << val << " " << wgt << " " << valsum/wgtsum << endl;
//       }
    }

    yk[ibpl] = valsum/wgtsum;

  }

  TGraph * grres = new TGraph(100, xy, yk);

  return grres;
//   grres->Draw("ACP");
//   cfth->Draw("CP");

//   TFile *fth = new TFile("CFplamres.root","RECREATE");
//   fth->cd();
//   grres->SetName("CFplamthres");
//   grres->SetTitle("CFplamthres");

//   grres->Write();
//   fth->Close();
}

TGraph* getresidualppgraph()
{
  // Old pp plots
  if (0) {
//     TGraphErrors* cfth = (TGraphErrors *) infpp->Get(Form("cppR%.1f", radiuspp));

//     double xy[60];
//     double yk[60];

//     int bcy, bcx;
//     double val, valsum, wgtsum;
//     double wgt;

//     for (int ibpl=0; ibpl<30; ibpl++) {
//       if (cfth) {

// 	xy[ibpl] = cfth->GetX()[ibpl];

// 	bcy = (ibpl+1);

// 	wgtsum = 0.0;
// 	valsum = 0.0;

// 	for (int ibpp = 0; ibpp<30; ibpp++) {

// 	  bcx = ibpp+1;

// 	  val = (cfth->GetY()[31+ibpp] + cfth->GetY()[32-ibpp])/2.0;
// 	  wgt = (htrpp->GetBinContent(bcx*2-1, bcy*2-1) +
// 		 htrpp->GetBinContent(bcx*2-1, bcy*2) +
// 		 htrpp->GetBinContent(bcx*2, bcy*2-1) +
// 		 htrpp->GetBinContent(bcx*2, bcy*2));

// 	  wgtsum += wgt;
// 	  valsum += val*wgt;
// 	}

// 	yk[ibpl] = valsum/wgtsum;
//       }
//       else {
// 	cout << "No cfth !!! " << radiuspp << endl;
//       }

//     TGraph * grres = new TGraph(30, xy, yk);

//     return grres;
  }

  TGraphErrors* cfth = (TGraphErrors *) infpp->Get(Form("cppR%.1f", radiuspp));

  double xy[100];
  double yk[100];

  int bcy, bcx;
  double val, valsum, wgtsum;
  double wgt;

  for (int ibpl=0; ibpl<100; ibpl++) {
    if (cfth) {

      // xy[ibpl] = cfth->GetX()[ibpl];
      xy[ibpl] = cfth->GetX()[TMath::FloorNint(cfth->GetN()/2.0)+ibpl];

      bcy = (ibpl+1);

      wgtsum = 0.0;
      valsum = 0.0;

      for (int ibpp = 0; ibpp<100; ibpp++) {

        bcx = ibpp+1;

        // val = cfth->GetY()[ibpp];
	if (ibpp > 97)
	  val=1;
	else
	  val = 0.5*(cfth->GetY()[TMath::FloorNint(cfth->GetN()/2.0)+ibpp]+cfth->GetY()[TMath::FloorNint(cfth->GetN()/2.0)-ibpp]);
        wgt = htrpp->GetBinContent(bcx, bcy);

        wgtsum += wgt;
        valsum += val*wgt;
      }

      yk[ibpl] = valsum/wgtsum;
    }
    else {
      cout << "No cfth !!! " << radiuspp << endl;
    }
  }

  TGraph *grres = new TGraph(100, xy, yk);

  return grres;
}

TGraph* getresidualpp()
{
  //  TFile *infpp = new TFile("CFppkpp4.root");
  //  TGraphErrors *cfth = infpp->Get("cppR3.0");

  //  TFile *inftr = new TFile("kplakpp_newbins.root");
  //  TH2D *htr = inftr->Get("kplakpp");

//   TGraphErrors* cfth = (TGraphErrors *) infpp->Get(Form("cppR%.1f", radius));
//   //  cout << "Got cfth " << Form("cppR%.1f", radius) << endl;

  double xy[100];
  double yk[100];

//   int bcy, bcx;
//   double val, valsum, wgtsum;
//   double wgt;

//   for (int ibpl=0; ibpl<30; ibpl++) {
//     xy[ibpl] = cfth->GetX()[ibpl+31];

//     bcy = (ibpl+1);

//     wgtsum = 0.0;
//     valsum = 0.0;

//     for (int ibpp = 0; ibpp<30; ibpp++) {

//       bcx = ibpp+1;

//       val = (cfth->GetY()[31+ibpp] + cfth->GetY()[32-ibpp])/2.0;
//       wgt = (htrpp->GetBinContent(bcx, bcy) +
// 	     htrpp->GetBinContent(bcx, bcy+1) +
// 	     htrpp->GetBinContent(bcx+1, bcy) +
// 	     htrpp->GetBinContent(bcx+1, bcy+1));

//       wgtsum += wgt;
//       valsum += val*wgt;
//     }

//     yk[ibpl] = valsum/wgtsum;

//   }

  if (0) {
    double xv, x1, x2, x3;
    xv = radius*TMath::Sqrt(2.0);
    x1 = floor(xv*10)/10 - 0.1;
    x2 = floor(xv*10)/10;
    x3 = floor(xv*10)/10 + 0.1;

    if (x1 < 0.99) { x1 = 1.0; x2 = 1.2; x3 = 1.4; }
//   if (x2 < 1.01) x2 = 1.0;
//   if (x3 < 1.01) x3 = 1.0;

    radiuspp = x1;
    TGraph *grp1 = getresidualppgraph();
    radiuspp = x2;
    TGraph *grp2 = getresidualppgraph();
    radiuspp = x3;
    TGraph *grp3 = getresidualppgraph();

    if (!grp1) cout << "No grp1" << endl;
    if (!grp2) cout << "No grp2" << endl;
    if (!grp3) cout << "No grp3" << endl;

    for (int ix=0; ix<30; ix++) {
      xy[ix] = grp1->GetX()[ix];
      yk[ix] = (grp1->GetY()[ix]*(x2-xv)*(x3-xv)/((x2-x1)*(x3-x1)) +
                grp2->GetY()[ix]*(x1-xv)*(x3-xv)/((x1-x2)*(x3-x2)) +
                grp3->GetY()[ix]*(x1-xv)*(x2-xv)/((x1-x3)*(x2-x3)));

      //    cout << grp1->GetY()[ix] << " " << grp2->GetY()[ix] << " " << grp3->GetY()[ix] << "   " << yk[ix] << endl;
    }

    TGraph * grres = new TGraph(30, xy, yk);

    return grres;
  }

  double xv, x1, x2, x3;
  xv = radius;
  x1 = 1.0 + floor((xv - 1.2)/0.2)*0.2;
  x2 = 1.0 + floor((xv - 1.0)/0.2)*0.2;
  x3 = 1.0 + floor((xv - 0.8)/0.2)*0.2;

  if (x1 < 1.01) x1 = 1.0;
  if (x2 < 1.01) x2 = 1.0;
  if (x3 < 1.01) x3 = 1.0;

  //  cout << "Interpolating " << x1 << " " << x2 << " "<< x3 << endl;

  radiuspp = x1;
  TGraph *grp1 = getresidualppgraph();
  radiuspp = x2;
  TGraph *grp2 = getresidualppgraph();
  radiuspp = x3;
  TGraph *grp3 = getresidualppgraph();

  if (!grp1) cout << "No grp1" << endl;
  if (!grp2) cout << "No grp2" << endl;
  if (!grp3) cout << "No grp3" << endl;

  for (int ix=0; ix<100; ix++) {
    xy[ix] = grp1->GetX()[ix];
    yk[ix] = (grp1->GetY()[ix]*(x2-xv)*(x3-xv)/((x2-x1)*(x3-x1)) +
              grp2->GetY()[ix]*(x1-xv)*(x3-xv)/((x1-x2)*(x3-x2)) +
              grp3->GetY()[ix]*(x1-xv)*(x2-xv)/((x1-x3)*(x2-x3)));

    //    cout << grp1->GetY()[ix] << " " << grp2->GetY()[ix] << " " << grp3->GetY()[ix] << "   " << yk[ix] << endl;
  }

  TGraph * grres = new TGraph(100, xy, yk);

  return grres;
}

void fitk()
{
//   TFile *fil1 = new TFile("shmout.recalc.APL.M0.cdep.root");
//   TH1D *prim1 = (TH1D *) fil1->Get("CfnReYlm00cylmV0APLtpcM0");
//   Double_t scale = prim1->Integral(71,90)/20.0;
//   prim1->Scale(1.0/scale);
//   prim1->Draw();

  TFile *fil1 = new TFile("NewestResults.root");
  TH1D *prim1 = (TH1D *) fil1->Get("hpal");
//   Double_t scale = prim1->Integral(71,90)/20.0;
//   prim1->Scale(1.0/scale);
  prim1->Draw();

  //  f0im = f0ia;
  TGraph *grpl = getgrk();
  TGraph *grps = getresidualkpLpS();
  //  f0im = f0ia;
  TGraph *gral = getgrk();
  TGraph *grpp = getresidualpp();

  Double_t valpp;

  for (int iter=0; iter<100; iter++) {
    grpl->GetY()[iter] = (grpl->GetY()[iter]-1.0)*lampl + 1.0;
    grps->GetY()[iter] = (grps->GetY()[iter]-1.0)*lamps + 1.0;
    if (iter<60) {
      if (iter == 0)
        valpp = grpp->GetY()[0];
      else if (iter%2)
        valpp = 0.25*grpp->GetY()[iter/2] + 0.75*grpp->GetY()[(iter+1)/2];
      else
        valpp = 0.75*grpp->GetY()[iter/2] + 0.25*grpp->GetY()[(iter+1)/2];
    }
    else
      valpp = 1.0 + (grpp->GetY()[29] - 1.0)*(1.0-(iter-60.0)/50);
    valpp = (valpp-1.0)*lampp + 1.0;

    gral->GetY()[iter] = 1.0 + (grpl->GetY()[iter]-1.0) + (grps->GetY()[iter]-1.0);
  }

  grpl->SetMarkerColor(2);
  grps->SetMarkerColor(kGreen+2);
  gral->SetMarkerColor(kBlue);
  grpl->Draw("CP");
  grps->Draw("CP");
  gral->Draw("CP");

  double chi2 = 0;
  double val, err;
  for (int ix=1; ix<49; ix++) {
    val = gral->GetY()[ix] - prim1->GetBinContent(ix+1);
    err = prim1->GetBinError(ix+1);
    chi2 += val*val/(err*err);
  }

  cout << "Chi2 is " << chi2 << "   " << endl;
}

double fitkcent()
{
  //  prims[0]->Draw();

  // Function 1

  TGraph *grpl;
  TGraph *grps;
  TGraph *gral;
  TGraph *grppf;
  TGraph *grpp;
  TGraph *grll;
  TGraph *grsl;
  TGraph *grls;
  TGraph *grss;
  TGraph *grpx;
  TGraph *grlx;
  TGraph *grsx;

//   TCanvas *canpl = new TCanvas("canpl","canpl",1300,1000);
//   canpl->Divide(2,2,0.0001,0.0001);

  double chi2 = 0;
  double norm;

  // background
  //TF1 *fc2 = new TF1("fc2","[3]+[2]*TMath::Gaus(x,[0],[1])",0.4,1);
  // // fc2->SetParameters(1.5,0.3,-0.2,1);
  // fc2->SetParameters(bkgmean,bkgwidth,bkgscale,bkgoffset);
  // fc2 = new TF1("fc2","[1]+[0]*x",0.5,1);
  // // fc2->SetParameters(1.5,0.3,-0.2,1);
  // fc2->SetParameters(bkgmean,bkgwidth,bkgscale,bkgoffset);

  for (int ifil=0; ifil<NumOfSystems; ifil++) {
    if (ifil == 0) radius = radius0;
    if (ifil == 1) radius = radius1;
    if (ifil == 2) radius = radius2;
    if (ifil == 3) radius = radius3;
    if (ifil == 4) radius = radius4;
    if (ifil == 5) radius = radius5;

    if (ifil == 0) norm = norm0;
    if (ifil == 1) norm = norm1;
    if (ifil == 2) norm = norm2;
    if (ifil == 3) norm = norm3;

    //    double f0imsav = f0im;

//     canpl->cd(ifil+1);
    //    cout << "Drawing " << ifil << endl;
    prims[ifil]->SetMaximum(1.04);
    prims[ifil]->SetMinimum(0.84);
    //    prims[ifil]->Draw();

    gral = getgrk();
    grppf = getresidualpp();
    grpp = getgrk();

    if (dofscale == 1) {
      double f0resav = f0re;
      double f0imsav = f0im;

      f0re = f0resav;
      f0im = f0imsav;
      grpl = getgrk();

      f0re = f0resav*0.20;
      f0im = f0imsav*0.48;
      grll = getresidualkpLLL();

      f0re = f0resav*0.20;
      f0im = f0imsav*0.42;
      grsl = getresidualkpLSL();

      f0re = f0resav*0.23;
      f0im = f0imsav*0.60;
      grps = getresidualkpLpS();

      f0re = f0resav*0.20;
      f0im = f0imsav*0.42;
      grls = getresidualkpLLS();

      f0re = f0resav*0.17;
      f0im = f0imsav*0.37;
      grss = getresidualkpLSS();

      f0re = f0resav*0.20;
      f0im = f0imsav*0.45;
      grpx = getresidualkpLpX();

      f0re = f0resav*0.14;
      f0im = f0imsav*0.33;
      grlx = getresidualkpLLX();

      f0re = f0resav*0.17;
      f0im = f0imsav*0.32;
      grsx = getresidualkpLSX();

      f0re = f0resav;
      f0im = f0imsav;
    }
    else if (dofscale == 2) {
      double f0resav = f0re;
      double f0imsav = f0im;

      f0re = f0resav;
      f0im = f0imsav;
      grpl = getgrk();

      f0re = f0resav*1.5;
      f0im = f0imsav*1.5;
      grll = getresidualkpLLL();

      f0re = f0resav*1.0;
      f0im = f0imsav*1.0;
      grsl = getresidualkpLSL();

      f0re = f0resav*1.0;
      f0im = f0imsav*1.0 ;
      grps = getresidualkpLpS();

      f0re = f0resav*1.5 ;
      f0im = f0imsav*1.5 ;
      grls = getresidualkpLLS();

      f0re = f0resav*1.0 ;
      f0im = f0imsav*1.0 ;
      grss = getresidualkpLSS();

      f0re = f0resav*0.5 ;
      f0im = f0imsav*0.5 ;
      grpx = getresidualkpLpX();

      f0re = f0resav*1.0 ;
      f0im = f0imsav*1.0 ;
      grlx = getresidualkpLLX();

      f0re = f0resav*1.0 ;
      f0im = f0imsav*1.0 ;
      grsx = getresidualkpLSX();

      f0re = f0resav;
      f0im = f0imsav;
    }
    else if (dofscale == 3) { // AQM scaling
      double f0resav = f0re;
      double f0imsav = f0im;

      f0re = f0resav;
      f0im = f0imsav;
      grpl = getgrk();

      f0re = f0resav*TMath::Sqrt(0.865);
      f0im = f0imsav*TMath::Sqrt(0.865);
      grll = getresidualkpLLL();

      f0re = f0resav*TMath::Sqrt(0.865);
      f0im = f0imsav*TMath::Sqrt(0.865);
      grsl = getresidualkpLSL();

      f0re = f0resav*1.0;
      f0im = f0imsav*1.0 ;
      grps = getresidualkpLpS();

      f0re = f0resav*TMath::Sqrt(0.865) ;
      f0im = f0imsav*TMath::Sqrt(0.865) ;
      grls = getresidualkpLLS();

      f0re = f0resav*TMath::Sqrt(0.865) ;
      f0im = f0imsav*TMath::Sqrt(0.865) ;
      grss = getresidualkpLSS();

      f0re = f0resav*TMath::Sqrt(0.844) ;
      f0im = f0imsav*TMath::Sqrt(0.844) ;
      grpx = getresidualkpLpX();

      f0re = f0resav*TMath::Sqrt(0.732) ;
      f0im = f0imsav*TMath::Sqrt(0.732) ;
      grlx = getresidualkpLLX();

      f0re = f0resav*TMath::Sqrt(0.732) ;
      f0im = f0imsav*TMath::Sqrt(0.732) ;
      grsx = getresidualkpLSX();

      f0re = f0resav;
      f0im = f0imsav;
    }
    else {
      grpl = getgrk();
      grll = getresidualkpLLL();
      grsl = getresidualkpLSL();
      grps = getresidualkpLpS();
      grls = getresidualkpLLS();
      grss = getresidualkpLSS();
      grpx = getresidualkpLpX();
      grlx = getresidualkpLLX();
      grsx = getresidualkpLSX();
    }


    Double_t valpp = 0.0;
    Double_t thscale = 0.0;
    Double_t wgtth = 0.0;
    Double_t wgt;

    for (int iter=0; iter<100; iter++) {
      if (ifil < 3) {
	grpl->GetY()[iter] = (grpl->GetY()[iter]-1.0)*lampl + 1.0;
	grps->GetY()[iter] = (grps->GetY()[iter]-1.0)*lamps + 1.0;
	grll->GetY()[iter] = (grll->GetY()[iter]-1.0)*lamll + 1.0;
	grls->GetY()[iter] = (grls->GetY()[iter]-1.0)*lamls + 1.0;
	grsl->GetY()[iter] = (grsl->GetY()[iter]-1.0)*lamsl + 1.0;
	grss->GetY()[iter] = (grss->GetY()[iter]-1.0)*lamss + 1.0;
	grpx->GetY()[iter] = (grpx->GetY()[iter]-1.0)*lampx + 1.0;
	grlx->GetY()[iter] = (grlx->GetY()[iter]-1.0)*lamlx + 1.0;
	grsx->GetY()[iter] = (grsx->GetY()[iter]-1.0)*lamsx + 1.0;
      }
      else {
	grpl->GetY()[iter] = (grpl->GetY()[iter]-1.0)*lamapl + 1.0;
	grps->GetY()[iter] = (grps->GetY()[iter]-1.0)*lamaps + 1.0;
	grll->GetY()[iter] = (grll->GetY()[iter]-1.0)*lamall + 1.0;
	grls->GetY()[iter] = (grls->GetY()[iter]-1.0)*lamals + 1.0;
	grsl->GetY()[iter] = (grsl->GetY()[iter]-1.0)*lamasl + 1.0;
	grss->GetY()[iter] = (grss->GetY()[iter]-1.0)*lamass + 1.0;
	grpx->GetY()[iter] = (grpx->GetY()[iter]-1.0)*lamapx + 1.0;
	grlx->GetY()[iter] = (grlx->GetY()[iter]-1.0)*lamalx + 1.0;
	grsx->GetY()[iter] = (grsx->GetY()[iter]-1.0)*lamasx + 1.0;
      }
      
//       if (iter<58) {
//  	if (iter == 0)
//  	  valpp = grppf->GetY()[0];
//  	else if (iter%2)
//  	  valpp = 0.25*grppf->GetY()[iter/2] + 0.75*grppf->GetY()[(iter+1)/2];
//  	else
//  	  valpp = 0.75*grppf->GetY()[iter/2] + 0.25*grppf->GetY()[(iter+1)/2];
//       }
//       else
//  	valpp = 1.0 + (grppf->GetY()[29] - 1.0)*(1.0-(iter-60.0)/50);
//       valpp = (valpp-1.0)*lampp + 1.0;
//      grpp->GetY()[iter] = valpp;
      grpp->GetY()[iter] = (grppf->GetY()[iter]-1.0)*lampp + 1.0;

      gral->GetY()[iter] = 1.0
        + (grpl->GetY()[iter]-1.0)
        + (grps->GetY()[iter]-1.0)
        + (grll->GetY()[iter]-1.0)
        + (grls->GetY()[iter]-1.0)
        + (grsl->GetY()[iter]-1.0)
        + (grss->GetY()[iter]-1.0)
        + (grpx->GetY()[iter]-1.0)
        + (grlx->GetY()[iter]-1.0)
        + (grsx->GetY()[iter]-1.0)
        + (grpp->GetY()[iter]-1.0);
      //	+ (valpp-1.0);

      if ((iter >= 72) && (iter <= 97)) {
        wgt = iter*iter;
        thscale += wgt*gral->GetY()[iter];
        wgtth += wgt;
// 	//	cout << "    Th scale " << thscale << endl;
      }

      //      cout << iter << " " << (grpl->GetY()[iter]-1.0) << " " << (grps->GetY()[iter]-1.0) << " " << (grll->GetY()[iter]-1.0) << " " << (valpp-1.0) << endl;
    }

    //    thscale /= 29.0;
    thscale /= wgtth;

    //    cout << "  Th scale " << thscale << endl;

    if (ifil == 0) {
      //       if (grpf1) delete grpf1;
      //       if (grrf1) delete grrf1;
      if (graf1mr) delete graf1mr;
      if (graf1) delete graf1;
      //       if (gruf1) delete gruf1;
      //       if (grlf1) delete grlf1;
      //       if (grsf1) delete grsf1;
      if (gropl) delete gropl;
      if (gropp) delete gropp;
      if (grops) delete grops;
      if (groll) delete groll;
      if (grols) delete grols;
      if (grosl) delete grosl;
      if (gross) delete gross;
      if (gropx) delete gropx;
      if (grolx) delete grolx;
      if (grosx) delete grosx;
    }
    if (ifil == 1) {
      if (grpf2) delete grpf2;
      if (grrf2) delete grrf2;
      if (graf2) delete graf2;
      if (graf2mr) delete graf2mr;
      if (gruf2) delete gruf2;
    }
    if (ifil == 2) {
      if (grpf3) delete grpf3;
      if (grrf3) delete grrf3;
      if (graf3) delete graf3;
      if (graf3mr) delete graf3mr;
      if (gruf3) delete gruf3;
    }
    if (ifil == 3) {
      if (grpf4) delete grpf4;
      if (grrf4) delete grrf4;
      if (graf4) delete graf4;
      if (graf4mr) delete graf4mr;
      if (gruf4) delete gruf4;
    }
    if (ifil == 4) {
      if (graf5) delete graf5;
      if (graf5mr) delete graf5mr;
    }
    if (ifil == 5) {
      if (graf6) delete graf6;
      if (graf6mr) delete graf6mr;
    }

    // momentum resolution
    Double_t xymr[100];
    Double_t ykmr[100];
    Double_t momspread = 0.01;
    Double_t valmr, wgtmr;
    Double_t wgtsummr, valsummr;

    for (int ix=0; ix<100; ix++) {
      xymr[ix] = ix*0.005 + 0.0025;
    
      valsummr = 0.0;
      wgtsummr = 0.0;
      for (int iy=-6; iy<7; iy++) {
	// Weight from momentum resolution spread
	wgtmr = exp(-(iy*0.005)*(iy*0.005)/(2*momspread*momspread)); 
	// Weight from q^2 phase space
	wgtmr *= TMath::Power(TMath::Abs(ix-iy)*0.005+0.0025,2);
	valmr = gral->GetY()[TMath::Abs(ix-iy)];
	//if (valmr > 1.1) valmr=1;
	wgtsummr += wgtmr;
	valsummr += valmr*wgtmr;
      }
      if (ix > 80)
      	ykmr[ix] = gral->GetY()[TMath::Abs(ix)];
      else
      	ykmr[ix] = valsummr/wgtsummr;
      if (ykmr[ix] > 100)
	cout << "mom res" <<ykmr[ix] << " " <<gral->GetY()[TMath::Abs(ix)] << " " << valsummr << " " << valmr << " " << wgtmr << " " <<ix << " " << endl;
    }
    
    TGraph *grplamom = new TGraph(100, xymr, ykmr);
    grplamom->SetName("grplamomres");
    grplamom->SetTitle("plamomres");
    
    if (ifil == 0) {  //grpf1 = grpl; grrf1 = grps; graf1 = gral; gruf1 = grpp; grlf1 = grll; grsf1 = grsl;  }
      gropl = grpl;
      gropp = grpp;
      grops = grps;
      groll = grll;
      grols = grls;
      grosl = grsl;
      gross = grss;
      gropx = grpx;
      grolx = grlx;
      grosx = grsx;
      graf1 = gral;
      graf1mr = grplamom;

    }
    if (ifil == 1) {  grpf2 = grpl; grrf2 = grps; graf2 = gral; graf2mr = grplamom; gruf2 = grpp;  }
    if (ifil == 2) {  grpf3 = grpl; grrf3 = grps; graf3 = gral; graf3mr = grplamom; gruf3 = grpp;  }
    if (ifil == 3) {  grpf4 = grpl; grrf4 = grps; graf4 = gral; graf4mr = grplamom; gruf4 = grpp;  }
    if (ifil == 4) {  graf5 = gral; graf5mr = grplamom;  }
    if (ifil == 5) {  graf6 = gral; graf6mr = grplamom;  }

    grpl->SetMarkerColor(2);
    grps->SetMarkerColor(kGreen+2);
    gral->SetMarkerColor(kBlue);
//     grpl->Draw("CP");
//     grps->Draw("CP");
//     gral->Draw("CP");


    double val, err;
    for (int ix=binbeg; ix<binend; ix++) {
//       val = (1.0/thscale)*(gral->GetY()[ix*2] + gral->GetY()[ix*2+1])/2.0 - norm*(prims[ifil]->GetBinContent(ix+1));
      if (prims[ifil]->GetBinError(ix+1) < 0.000000000001) continue;
      if (domomres)
	val = (grplamom->GetY()[ix*2] + grplamom->GetY()[ix*2+1])/2.0 - norm*(prims[ifil]->GetBinContent(ix+1));
      else
	val = (gral->GetY()[ix*2] + gral->GetY()[ix*2+1])/2.0 - norm*(prims[ifil]->GetBinContent(ix+1));
      err = prims[ifil]->GetBinError(ix+1)*norm;
      chi2 += val*val/(err*err);
    }

  }
  norm=norm0;
  //proton-antiproton
  for (int irad=0; irad < 3; ++irad) {
    
    CFapLakpapVal = getCFapLakpapVal(radiuspap[irad]);
    CFapLakaplVal = getgrk(radiuspap[irad]);
    int n=0;

    for (Double_t kstar = kstarmin; kstar < kstarmax; kstar += kstarstep) {
      Double_t valData = norm*cfppdata[irad]->GetBinContent(cfppdata[irad]->GetXaxis()->FindFixBin(kstar));
      Double_t valFunNoMomRes = (1 + lambdapap * (getCFppVal(radiuspap[irad], kstar)-1)  + lambdapal*(CFapLakpapVal->Eval(kstar)-1));
      Double_t valFunMomRes = 0.;
      Double_t valFun = 0.;
      //cout << endl << kstar << endl << "before momres = " << valFun << endl;
      if (domomres) {
	Double_t momspread = 0.01;
	Double_t valmr=0., wgtmr=0.;
	Double_t wgtsummr=0., valsummr=0.;

	// for (int iy=-6; iy<7; iy++) {
	for (Double_t iy=-3*momspread; iy<3*momspread+kstarstep; iy+=kstarstep) {
	  // Weight from momentum resolution spread
	  wgtmr = exp(-(iy)*(iy)/(2*momspread*momspread)); 
	  // Weight from q^2 phase space
	  wgtmr *= TMath::Power(TMath::Abs(kstar-iy),2);
	  valmr = (1 + lambdapap * (getCFppVal(radiuspap[irad], TMath::Abs(kstar-iy))-1)  + lambdapal*(CFapLakpapVal->Eval(TMath::Abs(kstar-iy))-1));
	  wgtsummr += wgtmr;
	  valsummr += valmr*wgtmr;

	  //cout << kstar << " " << iy << " " << valmr << " " << wgtmr << endl;
	}
	valFunMomRes = valsummr/wgtsummr;
	valFun = valFunMomRes;
      }
      else {
	valFun = valFunNoMomRes;
      }
      
      //cout << "after momres = " << valFun << endl;
      
      //Double_t valFun = (1 + lambdapap * (getCFppVal(radiuspap, kstar)-1)  + lambdapal*(getCFpLaVal(radiuspap, TMath::Abs(kstar))-1)) ;
      chi2 += (valData-valFun)*(valData-valFun) / (cfppdata[irad]->GetBinError(cfppdata[irad]->GetXaxis()->FindFixBin(kstar)) * cfppdata[irad]->GetBinError(cfppdata[irad]->GetXaxis()->FindFixBin(kstar))*norm*norm);

      // if (TMath::Abs(getCFpLaVal(radiuspap, TMath::Abs(kstar)) - getCFapLaVal(radiuspap, TMath::Abs(kstar))) > 0.05)
      //cout << kstar << " " << getCFpLaVal(radiuspap, TMath::Abs(kstar)) << " " << CFapLakpapVal->Eval(TMath::Abs(kstar)) << " " << CFapLakaplVal->Eval(kstar) <<endl;
    
      fitpap[irad]->SetPoint(n,kstar,valFunNoMomRes);
      fitpapmr[irad]->SetPoint(n,kstar,valFunMomRes);
      fpap->SetPoint(n,kstar,1 + lambdapap * (getCFppVal(radiuspap[irad], kstar)-1));
      fpal->SetPoint(n++,kstar,1 + lambdapal * (getCFpLaVal(radiuspap[irad], TMath::Abs(kstar))-1));
 
      // cout << kstar << " " << valData << " " << valFun << " " << radiuspap << endl;
    }
    delete CFapLakaplVal;
    delete CFapLakpapVal;
  }
  cout << "Chi2 is " << chi2 << "   " << radius << endl;

  return chi2;
}

void myfuncf(Int_t& i, Double_t *x, Double_t &f, Double_t *par, Int_t iflag)
{
  norm0 = par[0];
  norm1 = par[1];
  norm2 = par[2];
  norm3 = par[3];

  lampl = par[4];
  lamps = par[5];

  radius0 = par[6];
  radius1 = par[7];
  radius2 = par[8];
  radius3 = par[9];
  radius4 = par[10];
  radius5 = par[11];

  f0re = par[12]/0.197327;
  f0im = par[13]/0.197327;
  d0re = par[15]/0.197327;
  f0ia = par[18]/0.197327;

  lampp = par[14];
  lamll = par[16];
  lamls = par[17];
  lamsl = par[19];
  lamss = par[20];
  lampx = par[21];
  lamlx = par[22];
  lamsx = par[23];

  lamapl = par[24];
  lamaps = par[25];

  lamapp = par[26];
  lamall = par[27];
  lamals = par[28];
  lamasl = par[29];
  lamass = par[30];
  lamapx = par[31];
  lamalx = par[32];
  lamasx = par[33];

  radiuspap[0] = par[34];
  lambdapap = par[35];
  lambdapal = par[36];

  bkgmean = par[37];
  bkgwidth = par[38];
  bkgscale = par[39];
  bkgoffset = par[40];

  radiuspap[1] = par[41];
  radiuspap[2] = par[42];

  f = fitkcent();
}


int main(int argc, char **argv)
{
//   files[0] = new TFile("shmout.recalc.APL.M0.cdep.root");
//   files[1] = new TFile("shmout.recalc.APL.M1.cdep.root");
//   files[2] = new TFile("shmout.recalc.APL.M2.cdep.root");
//   files[3] = new TFile("shmout.recalc.APL.M3.cdep.root");

//   prims[0] = (TH1D *) files[0]->Get("CfnReYlm00cylmV0APLtpcM0");
//   prims[1] = (TH1D *) files[1]->Get("CfnReYlm00cylmV0APLpcM1");
//   prims[2] = (TH1D *) files[2]->Get("CfnReYlm00cylmV0APLtpcM2");
//   prims[3] = (TH1D *) files[3]->Get("CfnReYlm00cylmV0APLtpcM3");


  ifstream *inf;

  if (argc < 2)
    inf = new ifstream("fitpar.in");
  else
    inf = new ifstream(argv[1]);

  for (int iter=0; iter<npar; iter++) {
    GetPar(inf, &parsg[iter], &dofix[iter], &parmin[iter], &parmax[iter]);
    cout << "Got fix " << iter << " " << dofix[iter] << endl;
  }

  // for (int iter=0; iter<nopts; iter++) {
  //   GetOpts(inf, &opts[iter]);
  // }
  // binbeg = opts[0];
  // binend = opts[1];
  // dochimap = opts[2];
  // dofscale = opts[3];
  // domomres = opts[4];
  // dodivg = opts[5];

  char linebuf[100];
  inf->getline(linebuf, 100);
  binbeg = atoi( linebuf );
  inf->getline(linebuf, 100);
  binend = atoi( linebuf );
  inf->getline(linebuf, 100);
  dochimap = atoi( linebuf );
  inf->getline(linebuf, 100);
  dofscale = atoi( linebuf );
  inf->getline(linebuf, 100);
  domomres = atoi( linebuf );
  inf->getline(linebuf, 100);
  dodivg = atoi( linebuf );

  cout <<  "binbeg = " << binbeg;
  cout <<  binend;
  cout <<  dochimap;
  cout <<  dofscale;
  cout <<  domomres;
  cout <<  dodivg;

  // (*inf) >> binbeg;
  // (*inf) >> binend;
  // (*inf) >> dochimap;
  // (*inf) >> dofscale;
  // (*inf) >> domomres;
  // (*inf) >> dodivg;

  // files[0] = new TFile("cfssdc.root");
  // files[1] = new TFile("cfssdc.root");
  // files[2] = new TFile("cfssdc.root");
  // files[3] = new TFile("cfssdc.root");

  // files[0] = new TFile("divp4divgcfaplam.root");
  // files[1] = new TFile("divp4divgcfaplam.root");
  // files[2] = new TFile("divp4divgcfaplam.root");
  // files[3] = new TFile("divp4divgcfaplam.root");
  
  files[0] = new TFile("divgcfaplam.root");
  files[1] = new TFile("divgcfaplam.root");
  files[2] = new TFile("divgcfaplam.root");
  files[3] = new TFile("divgcfaplam.root");
  
  files[4] = new TFile("divp4cfPAP_11hcen0.root");
  files[5] = new TFile("divp4cfPAP_11hcen2.root");
  files[6] = new TFile("divp4cfPAP_11hcen4.root");
  // cfppdata[0] = (TH1D *) files[4]->Get("NumOutPckstarPAPtpcM0Psi3");
  // cfppdata[1] = (TH1D *) files[5]->Get("NumOutPckstarPAPtpcM2Psi3");
  // cfppdata[2] = (TH1D *) files[6]->Get("NumOutPckstarPAPtpcM4Psi3");
  cfppdata[0] = (TH1D *) files[4]->Get("divp4NumOutPckstarPAPtpcM0Psi3");
  cfppdata[1] = (TH1D *) files[5]->Get("divp4NumOutPckstarPAPtpcM2Psi3");
  cfppdata[2] = (TH1D *) files[6]->Get("divp4NumOutPckstarPAPtpcM4Psi3");

  // files[4] = new TFile("PosM0M1_0727_corr.root");
  // cfppdata[0] = (TH1D *) files[4]->Get("OutPAPM0M1");
  
  // files[4] = new TFile("fitresults_pap_cen0_3.root");
  // cfppdata = (TH1D *) files[4]->Get("cfdata");
  
  fileInCFplakpp = new TFile("CFaplakpap_out.root");
  fileInCFppkpp = new TFile("CFpap.root");
  for (int i=0;i<3;++i){
    fitpap[i] = new TGraph();
    fitpap[i]->SetName(Form("fitpap%d",i));
    fitpap[i]->SetMarkerStyle(20);
    fitpapmr[i] = new TGraph();
    fitpapmr[i]->SetName(Form("fitpapmr%d",i));
    fitpapmr[i]->SetMarkerStyle(20);
  }
  fpap = new TGraph();
  fpap->SetName("fpap");
  fpap->SetMarkerStyle(21);
  fpal = new TGraph();
  fpal->SetName("fpal");
  fpal->SetMarkerStyle(22);
  
  if (0 == strcmp(argv[2],"global")) {
    NumOfSystems = 0;
    if (dodivg) {
      prims[0] = (TH1D *) files[0]->Get("divgNumOutPcnonidV0PALtpcM0");
      prims[1] = (TH1D *) files[1]->Get("divgNumOutPcnonidV0PALtpcM2");
      prims[2] = (TH1D *) files[2]->Get("divgNumOutPcnonidV0PALtpcM4");
      prims[3] = (TH1D *) files[3]->Get("divgNumOutPcnonidV0APLtpcM0");
      prims[4] = (TH1D *) files[3]->Get("divgNumOutPcnonidV0APLtpcM2");
      prims[5] = (TH1D *) files[3]->Get("divgNumOutPcnonidV0APLtpcM4");
    }
    else {
      prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0PALtpcM0");
      prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0PALtpcM2");
      prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0PALtpcM4");
      prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
      prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
      prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");      
    }
  }
  else if (0 == strcmp(argv[2],"globalAQM")) {
    NumOfSystems = 6;
    prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0PALtpcM0");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0PALtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0PALtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }
  else   if (0 == strcmp(argv[2],"PALfreeRf")) {
    NumOfSystems = 3;
    prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0PALtpcM0");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0PALtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0PALtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }
  else   if (0 == strcmp(argv[2],"APLfreeRf")) {
    NumOfSystems = 3;
    prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0APLtpcM0");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0APLtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0APLtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }
  else   if (0 == strcmp(argv[2],"APLfreeRfM0")) {
    NumOfSystems = 1;
    prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0APLtpcM0");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0APLtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0APLtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }
  else   if (0 == strcmp(argv[2],"APLfreeRfM2")) {
    NumOfSystems = 1;
    prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0APLtpcM2");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0APLtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0APLtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }
  else   if (0 == strcmp(argv[2],"APLfreeRfM4")) {
    NumOfSystems = 1;
    prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0APLtpcM4");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0APLtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0APLtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }
  else   if (0 == strcmp(argv[2],"PALfreeRfM0")) {
    NumOfSystems = 1;
    prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0PALtpcM0");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0APLtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0APLtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }
  else   if (0 == strcmp(argv[2],"PALfreeRfM2")) {
    NumOfSystems = 1;
    prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0PALtpcM2");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0APLtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0APLtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }
  else   if (0 == strcmp(argv[2],"PALfreeRfM4")) {
    NumOfSystems = 1;
    prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0PALtpcM4");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0APLtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0APLtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }
  else   if (0 == strcmp(argv[2],"PALavgRf")) {
    NumOfSystems = 3;
    prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0PALtpcM0");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0PALtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0PALtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }
  else   if (0 == strcmp(argv[2],"APLavgRf")) {
    NumOfSystems = 3;
    prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0APLtpcM0");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0APLtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0APLtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }
  else   if (0 == strcmp(argv[2],"PALfreeRavgf")) {
    NumOfSystems = 3;
    prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0PALtpcM0");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0PALtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0PALtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }
  else   if (0 == strcmp(argv[2],"APLfreeRavgf")) {
    NumOfSystems = 3;
    prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0APLtpcM0");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0APLtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0APLtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }
  else   if (0 == strcmp(argv[2],"globalfreeRavgf")) {
    NumOfSystems = 6;
    prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0PALtpcM0");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0PALtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0PALtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }
  else   if (0 == strcmp(argv[2],"globalavgRf")) {
    NumOfSystems = 6;
    prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0PALtpcM0");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0PALtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0PALtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }
  else   if (0 == strcmp(argv[2],"PALfreeRfixstarf")) {
    NumOfSystems = 3;
    prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0PALtpcM0");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0PALtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0PALtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }
  else   if (0 == strcmp(argv[2],"APLfreeRfixstarf")) {
    NumOfSystems = 3;
    prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0APLtpcM0");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0APLtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0APLtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }
  else   if (0 == strcmp(argv[2],"divg")) {
    NumOfSystems = 2;
    prims[0] = (TH1D *) files[0]->Get("divgNumOutPcnonidV0PALtpcM0");
    prims[1] = (TH1D *) files[1]->Get("divgNumOutPcnonidV0APLtpcM0");
    prims[2] = (TH1D *) files[2]->Get("divgNumOutPcnonidV0PALtpcM0");
    prims[3] = (TH1D *) files[3]->Get("divgNumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("divgNumOutPcnonidV0APLtpcM0");
    prims[5] = (TH1D *) files[3]->Get("divgNumOutPcnonidV0APLtpcM0");
  }
  else   if (0 == strcmp(argv[2],"globaldivg")) {
    NumOfSystems = 6;
    prims[0] = (TH1D *) files[0]->Get("divgNumOutPcnonidV0PALtpcM0");
    prims[1] = (TH1D *) files[1]->Get("divgNumOutPcnonidV0PALtpcM2");
    prims[2] = (TH1D *) files[2]->Get("divgNumOutPcnonidV0PALtpcM4");
    prims[3] = (TH1D *) files[3]->Get("divgNumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("divgNumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("divgNumOutPcnonidV0APLtpcM4");
  }
  else  {
    NumOfSystems = 1;
    prims[0] = (TH1D *) files[0]->Get("cnum1da");
    prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0APLtpcM2");
    prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0APLtpcM4");
    prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
    prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
    prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");
  }

  // prims[0] = (TH1D *) files[0]->Get("NumOutPcnonidV0PALtpcM0");
  // prims[1] = (TH1D *) files[1]->Get("NumOutPcnonidV0PALtpcM2");
  // prims[2] = (TH1D *) files[2]->Get("NumOutPcnonidV0PALtpcM4");
  // prims[3] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM0");
  // prims[4] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM2");
  // prims[5] = (TH1D *) files[3]->Get("NumOutPcnonidV0APLtpcM4");

  if (dofscale)
    cout << "Doing the f scaling with sqrt{s}" << endl;

//   files[0] = new TFile("outalicepal.root");
//   files[1] = new TFile("outalicepal.root");
//   files[2] = new TFile("outalicepal.root");
//   files[3] = new TFile("outalicepal.root");

//   prims[0] = (TH1D *) files[0]->Get("CfnReYlm00cylmV0PALtpcM0");
//   prims[1] = (TH1D *) files[1]->Get("CfnReYlm00cylmV0PALtpcM0");
//   prims[2] = (TH1D *) files[2]->Get("CfnReYlm00cylmV0PALtpcM0");
//   prims[3] = (TH1D *) files[3]->Get("CfnReYlm00cylmV0PALtpcM0");

  // // Uncorrect for purity
  // for (int ib = 1; ib<=prims[0]->GetNbinsX(); ib++) {
  //   prims[0]->SetBinContent(ib, 1.0 + (prims[0]->GetBinContent(ib)-1)*0.15);
  //   prims[0]->SetBinError(ib, prims[0]->GetBinError(ib)*0.15);
  // }

  // for (int ib = 1; ib<=prims[1]->GetNbinsX(); ib++) {
  //   prims[1]->SetBinContent(ib, 1.0 + (prims[1]->GetBinContent(ib)-1)*0.15);
  //   prims[1]->SetBinError(ib, prims[1]->GetBinError(ib)*0.15);
  // }

  // for (int ib = 1; ib<=prims[2]->GetNbinsX(); ib++) {
  //   prims[2]->SetBinContent(ib, 1.0 + (prims[2]->GetBinContent(ib)-1)*0.15);
  //   prims[2]->SetBinError(ib, prims[2]->GetBinError(ib)*0.15);
  //   //    cout << "Is " << prims[2]->GetBinContent(ib) << endl;
  // }


  Double_t scale;
//    scale = prims[0]->Integral(45,49)/5.0;
//    prims[0]->Scale(1./scale);
  //   scale = prims[1]->Integral(81,100)/20.0;
//   prims[1]->Scale(1./scale);
//   scale = prims[2]->Integral(81,100)/20.0;
//   prims[2]->Scale(1./scale);
//   scale = prims[3]->Integral(81,100)/20.0;
//   prims[3]->Scale(1./scale);

  TFile *inftr = new TFile("output1new.root");
  htr = (TH2D *) inftr->Get("h");

  TFile *inftrll = new TFile("output2new.root");
  htrll = (TH2D *) inftrll->Get("h");

  TFile *inftrsl = new TFile("output3new.root");
  htrsl = (TH2D *) inftrsl->Get("h");

  TFile *inftrls = new TFile("output4new.root");
  htrls = (TH2D *) inftrls->Get("h");

  TFile *inftrss = new TFile("output5new.root");
  htrss = (TH2D *) inftrss->Get("h");

  TFile *inftrpx = new TFile("output6new.root");
  htrpx = (TH2D *) inftrpx->Get("h");

  TFile *inftrlx = new TFile("output7new.root");
  htrlx = (TH2D *) inftrlx->Get("h");

  TFile *inftrsx = new TFile("output8new.root");
  htrsx = (TH2D *) inftrsx->Get("h");

  //  infpp = new TFile("CFpapkpap2.root");
  // infpp = new TFile("CFpapraw.root");
  infpp = new TFile("CFpapkpap.root");

//   TFile *inftrpp = new TFile("kplakpp_newbins.root");
//   htrpp = (TH2D *) inftrpp->Get("kplakpp");

  TFile *inftrpp = new TFile("output1pppla.root");
  htrpp = (TH2D *) inftrpp->Get("h");

  graf1mr = 0;
  graf2mr = 0;
  graf3mr = 0;
  graf4mr = 0;
  graf5mr = 0;
  graf6mr = 0;

  graf1 = 0;
  graf2 = 0;
  graf3 = 0;
  graf4 = 0;
  graf5 = 0;
  graf6 = 0;

  grpf1 = 0;
  grpf2 = 0;
  grpf3 = 0;
  grpf4 = 0;

  grrf1 = 0;
  grrf2 = 0;
  grrf3 = 0;
  grrf4 = 0;

  gruf1 = 0;
  gruf2 = 0;
  gruf3 = 0;
  gruf4 = 0;

//   files[0] = new TFile("shmout.recalc.APL.M0.cdep.root");
//   files[1] = new TFile("shmout.recalc.APL.M1.cdep.root");
//   files[2] = new TFile("shmout.recalc.APL.M2.cdep.root");
//   files[3] = new TFile("shmout.recalc.APL.M3.cdep.root");

//   prims[0] = (TH1D *) files[0]->Get("CfnReYlm00cylmV0APLtpcM0");
//   prims[1] = (TH1D *) files[1]->Get("CfnReYlm00cylmV0APLtpcM1");
//   prims[2] = (TH1D *) files[2]->Get("CfnReYlm00cylmV0APLtpcM2");
//   prims[3] = (TH1D *) files[3]->Get("CfnReYlm00cylmV0APLtpcM3");

//  Double_t scale;
//   scale = prims[0]->Integral(75,94)/20.0;
//   prims[0]->Scale(1.0/scale);
//   scale = prims[1]->Integral(75,94)/20.0;
//   prims[1]->Scale(1.0/scale);
//   scale = prims[2]->Integral(75,94)/20.0;
//   prims[2]->Scale(1.0/scale);
//   scale = prims[3]->Integral(75,94)/20.0;
//   prims[3]->Scale(1.0/scale);

  if (0) {
    parsg[0] = 1.0;     parmin[0] = 0.95;     parmax[0] = 1.05;    dofix[0] = 1;
    parsg[1] = 1.0;     parmin[1] = 0.95;     parmax[1] = 1.05;    dofix[1] = 1;
    parsg[2] = 1.0;     parmin[2] = 0.95;     parmax[2] = 1.05;    dofix[2] = 1;
    parsg[3] = 1.0;     parmin[3] = 0.95;     parmax[3] = 1.05;    dofix[3] = 1;
    parsg[4] = 0.15 ;     parmin[4] = 0.10;     parmax[4] = 1.60;    dofix[4] = 1;
    parsg[5] = 0.10 ;     parmin[5] = 0.00;     parmax[5] = 0.90;    dofix[5] = 1;
    parsg[6] = 2.9;     parmin[6] = 1.0;      parmax[6] = 6.0;     dofix[6] = 0;
    parsg[7] = 3.7;     parmin[7] = 1.0;      parmax[7] = 6.0;     dofix[7] = 1;
    parsg[8] = 3.5;     parmin[8] = 1.0;      parmax[8] = 6.0;     dofix[8] = 1;
    parsg[9] = 3.1;     parmin[9] = 1.0;      parmax[9] = 6.0;     dofix[9] = 1;
    parsg[10] = -0.10;   parmin[10] = -10.0;   parmax[10] = 10.0;   dofix[10] = 0;
    parsg[11] = 0.53;    parmin[11] = -10.0;   parmax[11] = 10.0;   dofix[11] = 0;
//   parsg[10] = -0.294;   parmin[10] = -10.0;   parmax[10] = 10.0;   dofix[10] = 0;
//   parsg[11] = 0.595;    parmin[11] = -10.0;   parmax[11] = 10.0;   dofix[11] = 0;
    parsg[12] = 0.07;   parmin[12] = 0.0;     parmax[12] = 0.30;   dofix[12] = 1;
    parsg[13] = 3.35;    parmin[13] = -10.0;   parmax[13] = 10.0;   dofix[13] = 1;
    parsg[14] = 0.11;   parmin[14] = 0.0;     parmax[14] = 0.30;   dofix[14] = 1;
    parsg[15] = 0.07;   parmin[15] = 0.0;     parmax[15] = 0.30;   dofix[15] = 1;
  }

  fitter=TVirtualFitter::Fitter(0, npar);
  fitter->SetFCN(myfuncf);

  fitter->SetParameter(0, "Norm1"    ,parsg[0],  0.001, parmin[0],  parmax[0]);
  fitter->SetParameter(1, "Norm2"    ,parsg[1],  0.001, parmin[1],  parmax[1]);
  fitter->SetParameter(2, "Norm3"    ,parsg[2],  0.001, parmin[2],  parmax[2]);
  fitter->SetParameter(3, "Norm4"    ,parsg[3],  0.001, parmin[3],  parmax[3]);
  fitter->SetParameter(4, "Lampl"    ,parsg[4],  0.001, parmin[4],  parmax[4]);
  fitter->SetParameter(5, "LampS"    ,parsg[5],  0.001, parmin[5],  parmax[5]);
  fitter->SetParameter(6, "Rad1"     ,parsg[6],  0.001, parmin[6],  parmax[6]);
  fitter->SetParameter(7, "Rad2"     ,parsg[7],  0.001, parmin[7],  parmax[7]);
  fitter->SetParameter(8, "Rad3"     ,parsg[8],  0.001, parmin[8],  parmax[8]);
  fitter->SetParameter(9, "Rad4"     ,parsg[9],  0.001, parmin[9],  parmax[9]);
  fitter->SetParameter(10, "Rad5"     ,parsg[10],  0.001, parmin[10],  parmax[10]);
  fitter->SetParameter(11, "Rad6"     ,parsg[11],  0.001, parmin[11],  parmax[11]);
  fitter->SetParameter(12,"f0re"     ,parsg[12], 0.001, parmin[12], parmax[12]);
  fitter->SetParameter(13,"f0im"     ,parsg[13], 0.001, parmin[13], parmax[13]);
  fitter->SetParameter(14,"Lampp"    ,parsg[14], 0.001, parmin[14], parmax[14]);
  fitter->SetParameter(15,"d0re"     ,parsg[15], 0.001, parmin[15], parmax[15]);
  fitter->SetParameter(16,"Lamll"    ,parsg[16], 0.001, parmin[16], parmax[16]);
  fitter->SetParameter(17,"Lamls"    ,parsg[17], 0.001, parmin[17], parmax[17]);
  fitter->SetParameter(18,"f0ia"     ,parsg[18], 0.001, parmin[18], parmax[18]);
  fitter->SetParameter(19,"Lamsl"    ,parsg[19], 0.001, parmin[19], parmax[19]);
  fitter->SetParameter(20,"Lamss"    ,parsg[20], 0.001, parmin[20], parmax[20]);
  fitter->SetParameter(21,"Lampx"    ,parsg[21], 0.001, parmin[21], parmax[21]);
  fitter->SetParameter(22,"Lamlx"    ,parsg[22], 0.001, parmin[22], parmax[22]);
  fitter->SetParameter(23,"Lamsx"    ,parsg[23], 0.001, parmin[23], parmax[23]);
  fitter->SetParameter(24, "Lamapl"    ,parsg[24],  0.001, parmin[24],  parmax[24]);
  fitter->SetParameter(25, "LamapS"    ,parsg[25],  0.001, parmin[25],  parmax[25]);
  fitter->SetParameter(26,"Lamapp"    ,parsg[26], 0.001, parmin[26], parmax[26]);
  fitter->SetParameter(27,"Lamall"    ,parsg[27], 0.001, parmin[27], parmax[27]);
  fitter->SetParameter(28,"Lamals"    ,parsg[28], 0.001, parmin[28], parmax[28]);
  fitter->SetParameter(29,"Lamasl"    ,parsg[29], 0.001, parmin[29], parmax[29]);
  fitter->SetParameter(30,"Lamass"    ,parsg[30], 0.001, parmin[30], parmax[30]);
  fitter->SetParameter(31,"Lamapx"    ,parsg[31], 0.001, parmin[31], parmax[31]);
  fitter->SetParameter(32,"Lamalx"    ,parsg[32], 0.001, parmin[32], parmax[32]);
  fitter->SetParameter(33,"Lamasx"    ,parsg[33], 0.001, parmin[33], parmax[33]);

  fitter->SetParameter(34,"radiuspap0"    ,parsg[34], 0.001, parmin[34], parmax[34]);
  fitter->SetParameter(35,"lambdapap"    ,parsg[35], 0.001, parmin[35], parmax[35]);
  fitter->SetParameter(36,"lambdapal"    ,parsg[36], 0.001, parmin[36], parmax[36]);

  fitter->SetParameter(37,"bkgmean"    ,parsg[37], 0.001, parmin[37], parmax[37]);
  fitter->SetParameter(38,"bkgwidth"    ,parsg[38], 0.001, parmin[38], parmax[38]);
  fitter->SetParameter(39,"bkgscale"    ,parsg[39], 0.001, parmin[39], parmax[39]);
  fitter->SetParameter(40,"bkgoffset"    ,parsg[40], 0.001, parmin[40], parmax[40]);

  fitter->SetParameter(41,"radiuspap2"    ,parsg[41], 0.001, parmin[41], parmax[41]);
  fitter->SetParameter(42,"radiuspap4"    ,parsg[42], 0.001, parmin[42], parmax[42]);

//   fitter->SetParameter(12,"RSideS"   ,parsg[12], 0.001, parmin[12], parmax[12]);
//   fitter->SetParameter(13,"RlongS"   ,parsg[13], 0.001, parmin[13], parmax[13]);
//   fitter->SetParameter(14,"lambdaS"  ,parsg[14], 0.001, parmin[14], parmax[14]);
//   fitter->SetParameter(15,"c20lam"   ,parsg[15], 0.001, parmin[15], parmax[15]);
//   fitter->SetParameter(16,"c20rad"   ,parsg[16], 0.001, parmin[16], parmax[16]);
//   fitter->SetParameter(17,"c20gval"  ,parsg[17], 0.001, parmin[17], parmax[17]);
//   fitter->SetParameter(18,"c20mval"  ,parsg[18], 0.001, parmin[18], parmax[18]);
//   fitter->SetParameter(19,"c20wval"  ,parsg[19], 0.001, parmin[19], parmax[19]);

  for (int ipar=0; ipar<npar; ipar++) {
    if (dofix[ipar] == 1)
      fitter->FixParameter(ipar);
  }

  Double_t arglist[100];
  arglist[0] = 1;
  fitter->ExecuteCommand("CALL FCN", arglist, 1);
  //    fitter->FixParameter(0);
  //    fitter->FixParameter(1);
  //    fitter->FixParameter(2);
  //   fitter->FixParameter(3);
//    fitter->FixParameter(4);
//    fitter->FixParameter(5);
  //   fitter->FixParameter(6);
  //     fitter->FixParameter(11);
//   fitter->FixParameter(14);

//   fitter->FixParameter(3);
//   fitter->FixParameter(4);
//   fitter->FixParameter(15);

//   fitter->FixParameter(7);
//   fitter->FixParameter(8);
//   fitter->FixParameter(9);
//    fitter->FixParameter(11);
//    fitter->FixParameter(12);
//    fitter->FixParameter(13);
//    fitter->FixParameter(14);


  //fitter->FixParameter(15);
  arglist[0] = 0;
  fitter->ExecuteCommand("SET PRINT", arglist, 1);
  fitter->ExecuteCommand("MIGRAD", arglist, 0);
//   fitter->ExecuteCommand("MIGRAD", arglist, 0);
//   fitter->ExecuteCommand("MIGRAD", arglist, 0);

//   f1 ->SetParameters(fitter->GetParameter(0), fitter->GetParameter(1), fitter->GetParameter(2), 1.0, n0);
//   f20->SetParameters(Ro, Rs, Rl, 1.0, n0*n1);
//   f22->SetParameters(Ro, Rs, Rl, 1.0, n0*n2);

  double errs[npar];
  double pars[npar];

  for (int ipar=0; ipar<npar; ipar++) {
    pars[ipar] = fitter->GetParameter(ipar);    errs[ipar] = fitter->GetParError(ipar);
  }

  // pars[0] = fitter->GetParameter(0);    errs[0] = fitter->GetParError(0);
  // pars[1] = fitter->GetParameter(1);  	errs[1] = fitter->GetParError(1);
  // pars[2] = fitter->GetParameter(2);  	errs[2] = fitter->GetParError(2);
  // pars[3] = fitter->GetParameter(3);  	errs[3] = fitter->GetParError(3);
  // pars[4] = fitter->GetParameter(4);  	errs[4] = fitter->GetParError(4);
  // pars[5] = fitter->GetParameter(5);  	errs[5] = fitter->GetParError(5);
  // pars[6] = fitter->GetParameter(6);  	errs[6] = fitter->GetParError(6);
  // pars[7] = fitter->GetParameter(7);  	errs[7] = fitter->GetParError(7);
  // pars[8] = fitter->GetParameter(8);  	errs[8] = fitter->GetParError(8);
  // pars[9] = fitter->GetParameter(9);  	errs[9] = fitter->GetParError(9);
  // pars[10] = fitter->GetParameter(10);	errs[10] = fitter->GetParError(10);
  // pars[11] = fitter->GetParameter(11);	errs[11] = fitter->GetParError(11);
  // pars[12] = fitter->GetParameter(12);	errs[12] = fitter->GetParError(12);
  // pars[13] = fitter->GetParameter(13);	errs[13] = fitter->GetParError(13);
  // pars[14] = fitter->GetParameter(14);	errs[14] = fitter->GetParError(14);
  // pars[15] = fitter->GetParameter(15);	errs[15] = fitter->GetParError(15);
  // pars[16] = fitter->GetParameter(16);	errs[16] = fitter->GetParError(16);
  // pars[17] = fitter->GetParameter(17);	errs[17] = fitter->GetParError(17);
  // pars[18] = fitter->GetParameter(18);	errs[18] = fitter->GetParError(18);
  // pars[19] = fitter->GetParameter(19);	errs[19] = fitter->GetParError(19);
  // pars[20] = fitter->GetParameter(20);	errs[20] = fitter->GetParError(20);
  // pars[21] = fitter->GetParameter(21);	errs[21] = fitter->GetParError(21);
  // pars[22] = fitter->GetParameter(22);	errs[22] = fitter->GetParError(22);
  // pars[23] = fitter->GetParameter(23);	errs[23] = fitter->GetParError(23);

  TVectorD VecRes;
  VecRes.ResizeTo(npar);
  VecRes.SetElements(pars);

  TVectorD VecErr;
  VecErr.ResizeTo(npar);
  VecErr.SetElements(errs);

  double  chimult;
  int iter;
  double xyv[2];
  myfuncf(iter, xyv, chimult, pars, 1);
  cout << "chi2 is " << chimult;
  chimult /= (1*(binend-binbeg+1)-3);
  chimult = TMath::Sqrt(chimult);
  cout << "   " << chimult << endl;

//   cout << "Norm   " << fitter->GetParameter(6) << " +/- " << fitter->GetParError(6)*chimult << endl;
//   cout << "Lambda " << fitter->GetParameter(3) << " +/- " << fitter->GetParError(3)*chimult << endl;
//   cout << "Rout   " << fitter->GetParameter(0)*0.197327 << " +/- " << fitter->GetParError(0)*0.197327*chimult << endl;
//   cout << "Rside  " << fitter->GetParameter(1)*0.197327 << " +/- " << fitter->GetParError(1)*0.197327*chimult << endl;
//   cout << "Rlong  " << fitter->GetParameter(2)*0.197327 << " +/- " << fitter->GetParError(2)*0.197327*chimult << endl;

//   cout << "Dev20  " << fitter->GetParameter(7) << " +/- " << fitter->GetParError(7)*chimult << endl;
//   cout << "Dev22  " << fitter->GetParameter(8) << " +/- " << fitter->GetParError(8)*chimult << endl;

//   cout << "QBeg   " << fitter->GetParameter(9) << " +/- " << fitter->GetParError(9)*chimult << endl;
//   cout << "QSlope " << fitter->GetParameter(10) << " +/- " << fitter->GetParError(10)*chimult << endl;

//   cout << "LambdaS " << fitter->GetParameter(14) << " +/- " << fitter->GetParError(14)*chimult << endl;
//   cout << "RoutS  " << fitter->GetParameter(11)*0.197327 << " +/- " << fitter->GetParError(11)*0.197327*chimult << endl;
//   cout << "RsideS " << fitter->GetParameter(12)*0.197327 << " +/- " << fitter->GetParError(12)*0.197327*chimult << endl;
//   cout << "RlongS " << fitter->GetParameter(13)*0.197327 << " +/- " << fitter->GetParError(13)*0.197327*chimult << endl;

//   cout << "c20lam  " << fitter->GetParameter(15) << " +/- " << fitter->GetParError(15)*chimult << endl;
//   cout << "c20rad  " << fitter->GetParameter(16) << " +/- " << fitter->GetParError(16)*chimult << endl;
//   cout << "c20hval " << fitter->GetParameter(17) << " +/- " << fitter->GetParError(17)*chimult << endl;
//   cout << "c20mval " << fitter->GetParameter(18) << " +/- " << fitter->GetParError(18)*chimult << endl;
//   cout << "c20wval " << fitter->GetParameter(19) << " +/- " << fitter->GetParError(19)*chimult << endl;

  cout << "Radius  " << fitter->GetParameter(6) << " +/- " << fitter->GetParError(6)*chimult << endl;
  cout << "f0 real " << fitter->GetParameter(10) << " +/- " << fitter->GetParError(10)*chimult << endl;
  cout << "f0 imag " << fitter->GetParameter(11) << " +/- " << fitter->GetParError(11)*chimult << endl;

  cout << endl;
  cout << "Covariance matrix " << endl;

  for (int iter=0; iter<3; iter++) cout << fitter->GetCovarianceMatrixElement(0,iter) << " " ;
  cout << endl;

  for (int iter=0; iter<3; iter++) cout << fitter->GetCovarianceMatrixElement(1,iter) << " " ;
  cout << endl;

  for (int iter=0; iter<3; iter++) cout << fitter->GetCovarianceMatrixElement(2,iter) << " " ;
  cout << endl;

  cout << endl;
  cout << "Correlation matrix " << endl;

  for (int iter=0; iter<3; iter++) cout << fitter->GetCovarianceMatrixElement(0,iter)/TMath::Sqrt(fitter->GetCovarianceMatrixElement(0,0)*fitter->GetCovarianceMatrixElement(iter,iter)) << " " ;
  cout << endl;

  for (int iter=0; iter<3; iter++) cout << fitter->GetCovarianceMatrixElement(1,iter)/TMath::Sqrt(fitter->GetCovarianceMatrixElement(1,1)*fitter->GetCovarianceMatrixElement(iter,iter)) << " " ;
  cout << endl;

  for (int iter=0; iter<3; iter++) cout << fitter->GetCovarianceMatrixElement(2,iter)/TMath::Sqrt(fitter->GetCovarianceMatrixElement(2,2)*fitter->GetCovarianceMatrixElement(iter,iter)) << " " ;
  cout << endl;

  double radv = fitter->GetParameter(6);
  double lamv = fitter->GetParameter(4);
  double f0rv = fitter->GetParameter(10);
  double f0iv = fitter->GetParameter(11);

  TH2D *radlamcorr = new TH2D("radlamcorr",";radius [fm];\\lambda",21,radv*0.79,radv*1.21,21,lamv*0.79,lamv*1.21);
  TH2D *f0ref0imcorr = new TH2D("f0ref0im",";Re(f_{0}) (fm);Im(f_{0}) (fm)",51,f0rv*0.49,f0rv*1.51,31,f0iv*0.69,f0iv*1.31);
  if (dochimap==10) {
    for (int ibx=1; ibx<= radlamcorr->GetNbinsX(); ibx++) {
      pars[6] = radlamcorr->GetXaxis()->GetBinCenter(ibx);
      for (int iby=1; iby <= radlamcorr->GetNbinsY(); iby++) {
        pars[4] = radlamcorr->GetYaxis()->GetBinCenter(iby);
        myfuncf(iter, xyv, chimult, pars, 1);
        radlamcorr->SetBinContent(ibx,iby,chimult);
      }
    }
  }
  else if (dochimap==1) {
    for (int ibx=1; ibx<= f0ref0imcorr->GetNbinsX(); ibx++) {
      pars[10] = f0ref0imcorr->GetXaxis()->GetBinCenter(ibx);
      for (int iby=1; iby <= f0ref0imcorr->GetNbinsY(); iby++) {
        pars[11] = f0ref0imcorr->GetYaxis()->GetBinCenter(iby);
        myfuncf(iter, xyv, chimult, pars, 1);
        f0ref0imcorr->SetBinContent(ibx,iby,chimult);
      }
    }
  }


  TFile *ofile = new TFile(Form("palfitoutstarda%s%d%d.root",argv[2],domomres,dodivg),"RECREATE");
  ofile->cd();

  prims[0]->Write();
  prims[1]->Write();
  prims[2]->Write();
  prims[3]->Write();
  prims[4]->Write();
  prims[5]->Write();
  for (int i=0;i<3;++i){
    cfppdata[i]->Write();
    fitpap[i]->Write();
    fitpapmr[i]->Write();
  }
  fpap->Write();
  fpal->Write();
  //fc2->Write();
  // CFapLakpapVal->SetName("CFapLakpapVal");
  // CFapLakaplVal->SetName("CFapLakaplVal");
  // CFapLakpapVal->Write();
  // CFapLakaplVal->Write();

  TGraph *grcs = getfsq();
  grcs->SetName("grtotcs");
  grcs->Write();

  if (graf1) {
    graf1mr->SetName("grpalaf1mr");
    graf1->SetName("grpalaf1");
//     grpf1->SetName("grpalpf1");
//     grrf1->SetName("grpalrf1");
//     gruf1->SetName("grpaluf1");
//     grlf1->SetName("grpallf1");
//     grsf1->SetName("grpalsf1");
    gropl->SetName("gropl");
    gropp->SetName("gropp");
    grops->SetName("grops");
    groll->SetName("groll");
    grols->SetName("grols");
    grosl->SetName("grosl");
    gross->SetName("gross");
    gropx->SetName("gropx");
    grolx->SetName("grolx");
    grosx->SetName("grosx");

    rescaleGraph(graf1mr, pars[0]);
    rescaleGraph(graf1, pars[0]);
//     rescaleGraph(grpf1, pars[0]);
//     rescaleGraph(grrf1, pars[0]);
//     rescaleGraph(gruf1, pars[0]);
//     rescaleGraph(grlf1, pars[0]);
//     rescaleGraph(grsf1, pars[0]);

    graf1mr->Write();
    graf1->Write();
    gropl->Write();
    gropp->Write();
    grops->Write();
    groll->Write();
    grols->Write();
    grosl->Write();
    gross->Write();
    gropx->Write();
    grolx->Write();
    grosx->Write();

//     grpf1->Write();
//     grrf1->Write();
//     gruf1->Write();
//     grlf1->Write();
//     grsf1->Write();
  }
  if (graf2) {
    graf2->SetName("grpalaf2");
    graf2mr->SetName("grpalaf2mr");
    grpf2->SetName("grpalpf2");
    grrf2->SetName("grpalrf2");
    gruf2->SetName("grpaluf2");

    rescaleGraph(graf2, pars[1]);
    rescaleGraph(graf2mr, pars[1]);
    rescaleGraph(grpf2, pars[1]);
    rescaleGraph(grrf2, pars[1]);
    rescaleGraph(gruf2, pars[1]);

    graf2->Write();
    graf2mr->Write();
    grpf2->Write();
    grrf2->Write();
    gruf2->Write();
  }
  if (graf3) {
    graf3->SetName("grpalaf3");
    graf3mr->SetName("grpalaf3mr");
    grpf3->SetName("grpalpf3");
    grrf3->SetName("grpalrf3");
    gruf3->SetName("grpaluf3");

    rescaleGraph(graf3, pars[2]);
    rescaleGraph(graf3mr, pars[2]);
    rescaleGraph(grpf3, pars[2]);
    rescaleGraph(grrf3, pars[2]);
    rescaleGraph(gruf3, pars[2]);

    graf3->Write();
    graf3mr->Write();
    grpf3->Write();
    grrf3->Write();
    gruf3->Write();
  }
  if (graf4) {
    graf4->SetName("grpalaf4");
    graf4mr->SetName("grpalaf4mr");
    grpf4->SetName("grpalpf4");
    grrf4->SetName("grpalrf4");
    gruf4->SetName("grpaluf4");

    rescaleGraph(graf4, pars[3]);
    rescaleGraph(graf4mr, pars[3]);
    rescaleGraph(grpf4, pars[3]);
    rescaleGraph(grrf4, pars[3]);
    rescaleGraph(gruf4, pars[3]);

    graf4->Write();
    graf4mr->Write();
    grpf4->Write();
    grrf4->Write();
    gruf4->Write();
  }
  if (graf5) {
    graf5->SetName("grpalaf5");
    rescaleGraph(graf5, pars[3]);
    graf5->Write();
    graf5mr->SetName("grpalaf5mr");
    rescaleGraph(graf5mr, pars[3]);
    graf5mr->Write();
  }
  if (graf6) {
    graf6->SetName("grpalaf6");
    rescaleGraph(graf6, pars[3]);
    graf6->Write();
    graf6mr->SetName("grpalaf6mr");
    rescaleGraph(graf6mr, pars[3]);
    graf6mr->Write();
  }

  VecRes.Write();
  VecErr.Write();

  radlamcorr->Write();
  f0ref0imcorr->Write();

  ofile->Close();



  return 0;
}
