#define _USE_MATH_DEFINES
//#include "myTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>
#include <TPaveText.h>
#include "RooCategory.h"
#include "RooNumIntConfig.h"
#include "RooPlotable.h"
#include <TUnfold.h>

void plotting()
{
  Double_t ptbins []={6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 24, 27, 30, 35, 40, 45, 50};
  Int_t npb = sizeof(ptbins)/sizeof(double);

  TFile *f1 = TFile::Open("prdataplots.root");
  TFile *f2 = TFile::Open("prMcplots.root");
  TFile *f3 = TFile::Open("nprMcplots.root");
  TFile *f4 = TFile::Open("prEff.root");
  TFile *f5 = TFile::Open("nprEff.root");

  TH1F *heff = new TH1F ("heff","; p_{t}^{#mu#mu}(GeV/c); Eff #times Acc", npb-1, ptbins);
  heff->GetYaxis()->SetLimits(0,1.1);
  TH1F *hres = new TH1F ("hres", "; z(J/#psi); dN/N",10, 0, 1);
  hres->GetYaxis()->SetLimits(0,0.35);
  TCanvas* c = new TCanvas ("c", "", 1000, 800);

  TLine * lpty1 = new TLine (6.5,1,50,1);
  lpty1->SetLineColor(kRed);
  lpty1->SetLineStyle(2);
  lpty1->SetLineWidth(2);
  
  ////////eff plots/////////
  TEfficiency *hpreff = (TEfficiency*) f4->Get("pteff");
  TEfficiency *hnpreff = (TEfficiency*) f5->Get("pteff");

  hpreff->SetMarkerColor(kRed);
  hpreff->SetMarkerStyle(33);
  hnpreff->SetMarkerColor(kBlue);
  hnpreff->SetMarkerStyle(33);

  TLegend* l = new TLegend (0.5,0.7,0.7,0.8);
  l->AddEntry(hpreff, "Prompt J/#psi", "lp");
  l->AddEntry(hnpreff, "Nonprompt J/#psi", "lp");
  l->SetBorderSize(0);

  TPaveText* tbox0 = new TPaveText(0.15,0.6,0.3,0.8, "BRNDC");
  tbox0->AddText("p_{t}(J/#psi) > 6.5 GeV");
  tbox0->AddText("p_{t}(jets) > 20 GeV");
  tbox0->AddText("|y| < 2.4");
  tbox0->SetBorderSize(0);
  tbox0->SetFillColor(0);

  gStyle->SetOptStat(0);

  heff->Draw();
  lpty1->Draw("same");
  hpreff->Draw("same");
  hnpreff->Draw("same");
  tbox0->Draw("same");
  l->Draw("same");
  c->SaveAs("prNprPtEff.png");

  ////////results///////////
  TH1F *prd = (TH1F*) f1->Get("prz");
  TH1F *nprd = (TH1F*) f1->Get("nprz");
  TH1F *prm = (TH1F*) f2->Get("nz");
  TH1F *nprm = (TH1F*) f3->Get("nz");

  prd->SetMarkerColor(kMagenta);
  prd->SetMarkerStyle(33);
  prd->SetLineColor(kMagenta+2);
  prm->SetLineColor(kBlue);
  prm->SetLineWidth(2);

  TLegend* l1 = new TLegend (0.75,0.7,0.9,0.8);
  l1->AddEntry(prd, "Data", "lp");
  l1->AddEntry(prm, "MC", "lp");
  l1->SetBorderSize(0);
  l1->SetFillStyle(0);

  TPaveText* tbox1 = new TPaveText(0.15,0.6,0.3,0.8, "BRNDC");
  tbox1->AddText("Prompt J/#psi");
  tbox1->AddText("p_{t}(J/#psi) > 6.5 GeV");
  tbox1->AddText("p_{t}(jets) > 20 GeV");
  tbox1->AddText("|y| < 2.4");
  tbox1->SetBorderSize(0);
  tbox1->SetFillColor(0);

  TPaveText* tbox2 = new TPaveText(0.15,0.6,0.3,0.8, "BRNDC");
  tbox2->AddText("Nonprompt J/#psi");
  tbox2->AddText("p_{t}(J/#psi) > 6.5 GeV");
  tbox2->AddText("p_{t}(jets) > 20 GeV");
  tbox2->AddText("|y| < 2.4");
  tbox2->SetBorderSize(0);
  tbox2->SetFillColor(0);

  gStyle->SetOptStat(0);

  //hres->Draw();
  //hres->GetYaxis()->SetLimits(0,0.35);
  prm->GetYaxis()->SetTitle("dN/N");
  prm->GetXaxis()->SetTitle("z(J/#psi)");
  prm->SetLineColor(0);
  prm->Draw();
  prd->Draw("same");
  tbox1->Draw("same");
  c->SaveAs("resPrData.png");
  prm->SetLineColor(kBlue);
  prm->Draw("same");
  l1->Draw("same");
  c->SaveAs("resPrDataMc.png");

  nprd->SetMarkerColor(kMagenta);
  nprd->SetMarkerStyle(33);
  nprd->SetLineColor(kMagenta+2);
  nprm->SetLineColor(kBlue);
  nprm->SetLineWidth(2);

  //hres->Draw();
  //hres->GetYaxis()->SetLimits(0,0.35);
  nprd->GetYaxis()->SetTitle("dN/N");
  nprd->GetXaxis()->SetTitle("z(J/#psi)");
  nprd->Draw();
  tbox2->Draw("same");
  c->SaveAs("resNPrData.png");
  nprm->Draw("same");
  l1->Draw("same");
  c->SaveAs("resNPrDataMc.png");
}
