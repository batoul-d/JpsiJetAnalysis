#define myTree_cxx
#define _USE_MATH_DEFINES
#include "myTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>
#include <TPaveText.h>
#include "RooCategory.h"
#include "RooNumIntConfig.h"
#include "RooPlotable.h"

using namespace std;
using namespace  RooFit;
//using namespace RooPlot;

TLorentzVector* matchReco;
Float_t jpsi_m;
Float_t jpsi_pt;
Float_t jpsi_eta;
Float_t jpsi_phi;
Float_t jpsi_rap;
Float_t dr;
Float_t dphi;
Float_t dphimin;
Float_t deta;
Float_t drmin; 
Float_t z=100;
int triggerIndex_PP =0;
int k;
Float_t weight;
Double_t ptbins []={6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 24, 27, 30, 35, 40, 45, 50};
Double_t etabins []={-2.4, -1.6, -0.8,  0, 0.8, 1.6,  2.4};

RooRealVar* mass = new RooRealVar("Mass","#mu#mu mass", 2.6, 3.4, "GeV/c^{2}");
RooRealVar* genmass = new RooRealVar("genMass","#mu#mu mass", 2.6, 3.4, "GeV/c^{2}");
RooRealVar* ctau = new RooRealVar("ctau","c_{#tau}", -2.0, 5.0, "mm");
RooRealVar* genctau = new RooRealVar("genctau","gen c_{#tau}", -2.0, 5.0, "mm");
RooRealVar* zed = new RooRealVar("z", "z_{J/#psi}", 0, 1);
RooRealVar* r = new RooRealVar("r", "#DeltaR", 0, 0.5);
RooRealVar* genzed = new RooRealVar("genz", "genz_{J/#psi}", 0, 1);
RooRealVar* w = new RooRealVar("w","MC weight", 0.0, 10000000.0);
RooRealVar* genr = new RooRealVar("genr", "#genDeltaR", 0, 0.5);
RooRealVar* pt = new RooRealVar("pt", "pt", 6.5, 25.5);
RooRealVar* rap = new RooRealVar("rap", "rapidity", -2.4, 2.4);
RooRealVar* gen_pt = new RooRealVar("gen_pt", "pt for gen j/psi", 6.5, 25.5);
RooRealVar* gen_rap = new RooRealVar("gen_rap", "rapidity for gen j/psi", -2.4, 2.4);

RooArgSet* recoset = new RooArgSet(*mass, *ctau, *zed, *r, *pt, *rap, *w);
RooArgSet*  unwset = new RooArgSet(*mass, *ctau, *zed, *r, *pt, *rap);
RooArgSet*  genset = new RooArgSet (*genmass, *genctau, *genzed, *genr, *gen_pt, *gen_rap);
RooArgSet* prset = new RooArgSet(*mass, *ctau, *zed, *r, *pt, *rap, *w);
RooArgSet* nprset = new RooArgSet(*mass, *ctau, *zed, *r, *pt, *rap, *w);

RooDataSet* data = new RooDataSet ("data", "data for reconstructed J/#psi", *recoset, WeightVar(*w));
RooDataSet* unwdata = new RooDataSet ("unwdata", "unweighted data", *unwset);
RooDataSet* gendata = new RooDataSet ("gendata", "data for generated J/#psi", *genset);
RooDataSet* prdata = new RooDataSet ("prdata", "data for prompt J/#psi", *prset, WeightVar(*w));
RooDataSet* nprdata = new RooDataSet ("nprdata", "data for non prompt J/#psi", *nprset, WeightVar(*w));

void myTree::EffCalc()
{
  if (isMc)
    {
      TH1F* ptg= new TH1F ("ptg", "N_{gen} vs p_{T}; p_{T}; N_{total}", 25, ptbins);
      TH1F* ptr= new TH1F ("ptr", "N_{reco} vs p_{T}; p_{T}; N_{reco}", 25, ptbins);
      TH1F* rapr= new TH1F ("rapr","N_{reco} vs #eta; #eta; N_{reco}", 6, etabins);
      TH1F* rapg= new TH1F ("rapg","N_{gen} vs #eta; #eta; N_{total}", 6, etabins);
      TH2F* ptrapg= new TH2F ("ptrapg", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", 6, etabins, 25, ptbins);
      TH2F* ptrapr= new TH2F ("ptrapr", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", 6, etabins, 25, ptbins);

      Long64_t nentries =fChain->GetEntries();
      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;
	  for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	    {
	      TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	      jpsi_m=GenQQ4mom->M();
	      jpsi_pt = GenQQ4mom->Pt();
	      jpsi_rap = GenQQ4mom->Rapidity();
	      if (jpsi_pt>6.5 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5)
		{
		  ptg->Fill(jpsi_pt);
		  rapg->Fill(jpsi_rap);
		  ptrapg->Fill(jpsi_rap, jpsi_pt);
		}
	    }

	  if ( HLT_HIL1DoubleMu0ForPPRef_v1 && pPAprimaryVertexFilter)
	    {
	      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++)
		{
		  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
		  jpsi_pt = RecoQQ4mom->Pt();
		  jpsi_rap = RecoQQ4mom->Rapidity();
		  jpsi_m=RecoQQ4mom->M();
		  if (
		      jpsi_pt > 6.5  &&
		      (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		      (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		      (isTriggerMatch(iQQ, triggerIndex_PP)) //&& if it matches the trigger 
		      //(isMatchedRecoDiMuon(iQQ))
		      )
		    {
		      if (Reco_QQ_sign[iQQ]==0 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5) 
			{			    
			  ptr->Fill(jpsi_pt);
			  rapr->Fill(jpsi_rap);
			  ptrapr->Fill(jpsi_rap, jpsi_pt);
			}
		    }
		}
	    }
	}
      TEfficiency* gptef = new TEfficiency("gptef", "reconstruction efficiency fct of pt", 25, ptbins);
      if(TEfficiency::CheckConsistency(*ptr,*ptg))
      gptef = new TEfficiency (*ptr,*ptg);

      TEfficiency* grapef = new TEfficiency("grapef", "reconstruction efficiency fct of rapidity", 6, etabins);
      if(TEfficiency::CheckConsistency(*rapr,*rapg))
      grapef = new TEfficiency (*rapr,*rapg);

      TEfficiency* gptrapef = new TEfficiency("gptetaef", "reconstruction efficiency fct of pt and rapidity; y; pt; eff", 6, etabins, 25, ptbins);
      if(TEfficiency::CheckConsistency(*ptrapr, *ptrapg))
      gptrapef = new TEfficiency (*ptrapr, *ptrapg);

      TFile* ef (0x0);
      if (isPr)
	ef = new TFile ("prEff.root", "RECREATE");
      else
	ef = new TFile ("nprEff.root", "RECREATE");
      gptef->Write("pteff");
      gptrapef->Write("ptrap");
      grapef->Write("rapeff");
      ef->Close();
    }
  else 
    cout<< "this is data and not MC"<<endl;
}

void myTree::ClosureTest()
{
  TH1F* rpt = new TH1F ("rpt", "pt distribution at reco level", 10, 6.5, 26.5);
  TH1F* gpt = new TH1F ("gpt", "pt distribution at gen level", 10, 6.5, 26.5);
  TFile* f (0x0);
  if (isPr)
    f= TFile::Open("prEff.root");
  else
    f= TFile::Open ("nprEff.root");
  TEfficiency* eff = (TEfficiency*) f->Get("ptrap");

  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  if (isMc)
    {
      for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;


	  for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	    {
	      TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	      jpsi_m = GenQQ4mom->M();
	      jpsi_pt = GenQQ4mom->Pt();
	      jpsi_rap = GenQQ4mom->Rapidity();

	      if (jpsi_pt>6.5 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5)
		{
		  gpt->Fill(jpsi_pt);
		}
	    }
	  if ( HLT_HIL1DoubleMu0ForPPRef_v1 && pPAprimaryVertexFilter)
	    {
	      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) 
		{
		  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
		  jpsi_pt = RecoQQ4mom->Pt();
		  jpsi_rap = RecoQQ4mom->Rapidity();
		  jpsi_m = RecoQQ4mom->M();

		  if (
		      jpsi_pt > 6.5  &&
		      (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		      (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		      (isTriggerMatch(iQQ, triggerIndex_PP)) //&&// if it matches the trigger 
		      )
		    {
		      if (Reco_QQ_sign[iQQ]==0 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5) 
			{
			  weight=1.0/(eff->GetEfficiency(eff->FindFixBin(jpsi_rap, jpsi_pt)));
			  rpt->Fill(jpsi_pt, weight);
			}
		    }
		}
	    }
	}
      TFile* testfile (0x0);
      if (isPr)
	testfile = new TFile ("prClosureTest.root","RECREATE");
      else
	testfile = new TFile ("nprClosureTest.root","RECREATE");
      rpt->Write("recopt");
      gpt->Write("genpt");
      testfile->Close();
    }
  else
    cout<< "this is data and not MC"<<endl;
}

void myTree::Loop()
{
  TF1* lcut = new TF1 ("lcut", "0.013+0.25/x", 6.5, 50);
  if (fChain == 0) return;
  //genzed->setBins(10);
  //mass->setBins(40);
  //zed->setBins(10);
  //rap->setBins(12);
  //gen_rap->setBins(12);
  Bool_t ismatch;
  TFile* f (0x0);
  TFile* prf = TFile::Open("prEff.root");
  TFile* nprf = TFile::Open ("nprEff.root");

  TEfficiency* preff = (TEfficiency*) prf->Get("ptrap");
  TEfficiency* npreff = (TEfficiency*) prf->Get("ptrap");
  TEfficiency* eff(0x0);
  if (isPr)
    eff=preff;
  else
    eff=npreff;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  if (isMc)
    {
      for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
	  ismatch=false;
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;


	  for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	    {
	      TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	      jpsi_m=GenQQ4mom->M();
	      jpsi_pt = GenQQ4mom->Pt();
	      jpsi_rap = GenQQ4mom->Rapidity();
	      drmin = 1000;

	      if (jpsi_pt>6.5 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5)
		{
  
		  for (Long64_t ijet=0; ijet<ngen; ijet++)
		    {
		      if (abs(geny[ijet])<2.4 && genpt[ijet]>20)
			{
			  TLorentzVector v_jet;
			  v_jet.SetPtEtaPhiM(genpt[ijet], geneta[ijet], genphi[ijet], genm[ijet]);
			  dphi= GenQQ4mom->DeltaPhi(v_jet);
			  dr = GenQQ4mom->DeltaR (v_jet);
			  if (dr<=drmin)
			    {
			      drmin=dr;
			      dphimin=dphi;
			      deta=(jpsi_eta-geneta[ijet]);
			      z= jpsi_pt/genpt[ijet];
			    }
			}
		    }
		  if (drmin<0.5)
		    {
		      genmass->setVal(jpsi_m);
		      gen_pt->setVal(jpsi_pt);
		      gen_rap->setVal(jpsi_rap);  
		      genzed->setVal(z);
		      genr->setVal(drmin);
		      genctau->setVal(Gen_QQ_ctau[iQQ]);
		      gendata->add(*genset);
		      ismatch=true;
		    }
		}
	    }
	  if (!ismatch) continue;
	  if ( HLT_HIL1DoubleMu0ForPPRef_v1 && pPAprimaryVertexFilter)
	    {
	      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) 
		{
		  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
		  jpsi_pt = RecoQQ4mom->Pt();
		  jpsi_eta = RecoQQ4mom->Eta();
		  jpsi_rap = RecoQQ4mom->Rapidity();
		  jpsi_m = RecoQQ4mom->M();
		  drmin = 1000;

		  if (
		      jpsi_pt > 6.5  &&
		      (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		      (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		      (isTriggerMatch(iQQ, triggerIndex_PP)) //&&// if it matches the trigger 
		      )
		    {
		      if (Reco_QQ_sign[iQQ]==0 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5) 
			{
			  for (Long64_t ijet=0; ijet<nref; ijet++)
			    {
	  
			      if (abs(jty[ijet])<2.4 && jtpt[ijet]>20)
				{
				  TLorentzVector v_jet;
				  v_jet.SetPtEtaPhiM(jtpt[ijet], jteta[ijet], jtphi[ijet], jtm[ijet]);
				  dphi= RecoQQ4mom->DeltaPhi(v_jet);
				  dr = RecoQQ4mom->DeltaR (v_jet);
				  if (dr<=drmin)
				    {
				      drmin=dr;
				      dphimin=dphi;
				      deta=(jpsi_eta-jteta[ijet]);
				      z= jpsi_pt/jtpt[ijet];
				    }
				}
			    }
			  if (drmin<0.5)
			    {
			      mass->setVal(jpsi_m);
			      pt->setVal(jpsi_pt);
			      rap->setVal(jpsi_rap);
			      zed->setVal(z);
			      r->setVal(drmin);
			      ctau->setVal(Reco_QQ_ctau[iQQ]);
			      weight=(eff->GetEfficiency(eff->FindFixBin(jpsi_rap, jpsi_pt)));
			      if (weight==0)
				cout <<"rap= "<<jpsi_rap<<"  pt= "<<jpsi_pt<< "  e= "<<weight<<endl;
			      w->setVal(1.0/weight);
			      data->add(*recoset, w->getVal());
			      unwdata->add(*unwset);
			    }
			}
		    }
		}
	    }
	}
    }
  else
    {
      for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;

	  if ( HLT_HIL1DoubleMu0_v1 && pPAprimaryVertexFilter)
	    {
	      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++)
		{
		  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
		  jpsi_pt = RecoQQ4mom->Pt();
		  jpsi_eta = RecoQQ4mom->Eta();
		  jpsi_rap = RecoQQ4mom->Rapidity();
		  jpsi_m = RecoQQ4mom->M();
		  drmin = 1000;

		  if (
		      jpsi_pt > 6.5  &&
		      (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		      (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		      (isTriggerMatch(iQQ, triggerIndex_PP)) // if it matches the trigger 
		      )
		    {
		      if (Reco_QQ_sign[iQQ]==0 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5) 
			{
			  for (Long64_t ijet=0; ijet<nref; ijet++)
			    {
	  
			      if (abs(jty[ijet])<2.4 && jtpt[ijet]>20)
				{
				  TLorentzVector v_jet;
				  v_jet.SetPtEtaPhiM(jtpt[ijet], jteta[ijet], jtphi[ijet], jtm[ijet]);
				  dphi= RecoQQ4mom->DeltaPhi(v_jet);
				  dr = RecoQQ4mom->DeltaR (v_jet);
				  if (dr<=drmin)
				    {
				      drmin=dr;
				      dphimin=dphi;
				      deta=(jpsi_eta-jteta[ijet]);
				      z= jpsi_pt/jtpt[ijet];
				    }
				}
			    }
			  if (drmin<0.5)
			    {
			      mass->setVal(jpsi_m);
			      pt->setVal(jpsi_pt);
			      rap->setVal(jpsi_rap);
			      zed->setVal(z);
			      r->setVal(drmin);
			      weight=(eff->GetEfficiency(eff->FindFixBin(jpsi_rap, jpsi_pt)));
			      if (weight==0)
				{
				  cout <<"rap= "<<jpsi_rap<<"  pt= "<<jpsi_pt<< "  e= "<<weight<<endl;
				  weight=1;}
			      w->setVal(1.0/weight);
			      data->add(*recoset, w->getVal());
			      unwdata->add(*unwset);
			      if(Reco_QQ_ctau[iQQ] < lcut->Eval(jpsi_pt))
				{
				  weight=preff->GetEfficiency(preff->FindFixBin(jpsi_rap, jpsi_pt));
				  if (weight==0)
				    weight=1;
				  w->setVal(1.0/weight);
				  prdata->add(*prset, w->getVal());
				}
			      else
				{
				  weight=npreff->GetEfficiency(npreff->FindFixBin(jpsi_rap, jpsi_pt));
				  if (weight==0)
				    weight=1;
				  w->setVal(1.0/weight);
				  nprdata->add(*nprset, w->getVal());
				}
			    }
			}
		    }
		}
	    }
	}
    }	    
  TFile *datafile (0x0);
  if (isMc)
    {
      if (isPr)
	datafile = new TFile ("prMcDatasets","RECREATE");
      else 
	datafile = new TFile ("nprMcDatasets","RECREATE");
      datafile->cd();
      gendata->Write("genData");
      data->Write("recoData");
      unwdata->Write("unwData");
    }
  else
    {
      if (isPr)
	datafile = new TFile ("datasets_pr","RECREATE");
      else
	datafile = new TFile ("datasets_npr","RECREATE");
      datafile->cd();
      data->Write("data");
      unwdata->Write("unwData");
      prdata->Write("prData");
      nprdata->Write("nprData");
    }
      datafile->Write();
      datafile->Close();
      delete datafile;
}

void myTree::Plot()
{
  Float_t nv[12];
  Float_t av[12];
  TFile* mcfile (0x0);
  TFile* datafile (0x0);
  TH1F* gn = new TH1F ("gn", "values of n from mc", 12, 0, 1.1);
  TH1F* ga = new TH1F ("ga", "values of alpha from mc", 12, 0, 1.1);
  TH1F* nz = new TH1F ("nz", "N_{J/#psi} in z bins", 10, 0, 1);
  TH1F* np = new TH1F ("np", "N_{J/#psi} in z bins", 10, 0, 1);
  TAxis* xaxis(0x0);
  if (isPr)
    {
      mcfile = TFile::Open("prMcDatasets");
      datafile = TFile::Open("datasets_pr");
    }
  else
    {
      mcfile = TFile::Open("nprMcDatasets");
      datafile = TFile::Open("datasets_npr");
    }
  gendata= (RooDataSet*) mcfile->Get("genData");
  prdata= (RooDataSet*) datafile->Get("prData");
  nprdata= (RooDataSet*) datafile->Get("nprData");
  if (isMc)
    {
      data= (RooDataSet*) mcfile->Get("recoData");
      unwdata= (RooDataSet*) mcfile->Get("unwData");
    }
  else
    {  
      data= (RooDataSet*) datafile->Get("data");
      unwdata= (RooDataSet*) datafile->Get("unwData");
    }
  genzed->setBins(10);
  mass->setBins(40);
  zed->setBins(10);
  rap->setBins(12);
  gen_rap->setBins(12);
  pt->setBins(19);
  gen_pt->setBins(19);	
  RooDataSet* d = (RooDataSet*) data->reduce("Mass>2.6 && Mass<3.5 && r<0.5");
  RooDataSet* unwd = (RooDataSet*) unwdata->reduce("Mass>2.6 && Mass<3.5 && r<0.5");
  RooDataSet* gend = (RooDataSet*) gendata->reduce("genMass>2.6 && genMass<3.5 && genr<0.5");

  RooPlot* mframe = mass->frame();
  RooPlot* comp = mass->frame();
  RooPlot* zframe = zed->frame();
  RooPlot* rframe = r->frame();
  RooPlot* ptframe = pt->frame();
  RooPlot* rapframe = rap->frame();

  RooRealVar* mean = new RooRealVar ("mean", "mean", 3.096, 3, 3.2);
  RooRealVar* sigma = new RooRealVar ("sigma", "sigma", 0.04, 0, 0.1);
  RooRealVar* sigma1 = new RooRealVar("sigma1", "sigma1", 0.02, 0, 0.1);
  RooRealVar* alpha (0x0);
  RooRealVar* n (0x0);


  if (isMc) 
    {
      alpha = new RooRealVar ("alpha", "alpha", 1.5, 0, 4);
      n = new RooRealVar("n", "n", 1.5, 0.2, 4);
    }
  else
    {
      alpha = new RooRealVar ("alpha", "alpha", 1.94);
      n = new RooRealVar("n", "n", 1.39);
    }

  RooCBShape cb ("cb", "cb", *mass, *mean, *sigma, *alpha, *n);
  RooCBShape cb1 ("cb1", "cb1", *mass, *mean, *sigma1, *alpha, *n);

  RooRealVar c0("c0","coefficient #0", -0.5,-1,1);
  RooRealVar c1("c1","coefficient #1", 0.06,-1,1);
  RooChebychev bkg("bkg","background p.d.f.", *mass, RooArgList(c0,c1)) ;

  RooRealVar fb ("fb", "background nb", 10, 0, data->numEntries());
  RooRealVar fs ("fs", "signal nb", 100, 0, data->numEntries());
  RooRealVar fcb ("fcb", "cb", 0.1, 0, 1);

  mass->setRange("window", 2.9, 3.2) ;

  RooExtendPdf ecb ("ecb", "ecb", cb, fcb);

  RooAddPdf masspeak ("masspeak", "fit of mass peak", RooArgList(cb, cb1), fcb);
  RooExtendPdf es ("es", "es", masspeak, fs, "window");
  RooExtendPdf eb ("eb", "eb", bkg, fb);
  RooAddPdf* model;
  if (isMc)
    model = new RooAddPdf ("model", "model", RooArgList(cb, cb1), fcb);
  else 
    model = new RooAddPdf ("model", "model", RooArgList(es, eb));

  d->plotOn(mframe,DataError(RooAbsData::SumW2));
  d->plotOn(rframe);
  d->plotOn(ptframe, DataError(RooAbsData::SumW2));
  d->plotOn(rapframe, DataError(RooAbsData::SumW2));
  //if (isMc)
  //d->plotOn(zframe, Name("recod"), DataError(RooAbsData::SumW2));

  model->fitTo(*d, Extended(!isMc), SumW2Error(kTRUE));
  if (isMc)
    {
      xaxis=gn->GetXaxis();
      gn->SetBinContent(xaxis->FindBin(1), n->getValV());
      xaxis=ga->GetXaxis();
      ga->SetBinContent(xaxis->FindBin(1), alpha->getValV());
    }
  model->plotOn(mframe, Name("background"), Components(RooArgSet(bkg)),DrawOption("F"), FillColor(kGray));
  model->plotOn(mframe, Name("signal"), Components(RooArgSet(masspeak)), LineStyle(kDashed));
  model->plotOn(mframe, Name("total"), LineColor(kRed));
  model->paramOn(mframe);


  TFile* fsave (0x0);
  if (isMc)
    {
      if (isPr)
	fsave = new TFile ("prMcplots.root","RECREATE");
      else
	fsave = new TFile ("nprMcplots.root", "RECREATE");
    }
  else
    {
      if (isPr)
	fsave = new TFile("prdataplots.root", "RECREATE");
      else 
	fsave = new TFile("nprdataplots.root", "RECREATE");
    }

  RooWorkspace* w (0x0);
  if (isMc)
    {
      for (int i =0; i<10; i++)
	{
	  RooPlot* mzframe = mass->frame();
	  w = new RooWorkspace(Form("w%d",i));
	  RooDataSet* dz = (RooDataSet*) d->reduce(Form("z >= %f && z<%f", i*0.1, (i+1)*0.1));
	  w->import(*dz);
	  w->import(*model);
	  w->data("data")->plotOn(mzframe, DataError(RooAbsData::SumW2));
	  w->pdf("model")->fitTo(*w->data("data"),  Extended(!isMc), SumW2Error(kTRUE));
	  w->pdf("model")->plotOn(mzframe);
	  if (isMc)
	    {
	      xaxis=gn->GetXaxis();
	      gn->SetBinContent(xaxis->FindBin(i*0.1), w->var("n")->getValV());
	      xaxis=ga->GetXaxis();
	      ga->SetBinContent(xaxis->FindBin(i*0.1), w->var("alpha")->getValV());
	      xaxis=nz->GetXaxis();
	      nz->SetBinContent(xaxis->FindBin(i*0.1), w->data("data")->sumEntries());
	      nz->SetBinError(xaxis->FindBin(i*0.1), sqrt(w->data("data")->sumEntries())); 
	    }
	  TPaveText* tbox = new TPaveText(0.15,0.6,0.5,0.7, "BRNDC");
	  tbox-> AddText (Form("%.1f #leq z < %.1f", i*0.1, (i+1)*0.1));
	  mzframe->addObject(tbox);
	  tbox->SetBorderSize(0);
	  tbox->SetFillColor(0);
	  w->pdf("model")->paramOn(mzframe);
	  mzframe->Write(Form("massframe%d%d",i,(i+1)));
	}
    }

  else
    {
      for (int i =0; i<10; i++)
	{
	  w = new RooWorkspace(Form("prw%d",i));
	  RooPlot* mzframe = mass->frame();
	  RooDataSet* prdz = (RooDataSet*) prdata->reduce(Form("z >= %f && z<%f", i*0.1, (i+1)*0.1));
	  w->import(*prdz);
	  w->import(*model);
	  w->data("prdata")->plotOn(mzframe, DataError(RooAbsData::SumW2));
	  w->pdf("model")->fitTo(*w->data("prdata"),  Extended(!isMc), SumW2Error(kTRUE));
	  w->pdf("model")->plotOn(mzframe);

	  xaxis=nz->GetXaxis();
	  nz->SetBinContent(xaxis->FindBin(i*0.1), w->var("fs")->getValV());
	  nz->SetBinError(xaxis->FindBin(i*0.1), w->var("fs")->getError());

	  TPaveText* tbox = new TPaveText(0.15,0.6,0.5,0.7, "BRNDC");
	  tbox-> AddText (Form("%.1f #leq z < %.1f", i*0.1, (i+1)*0.1));
	  mzframe->addObject(tbox);
	  tbox->SetBorderSize(0);
	  tbox->SetFillColor(0);
	  w->pdf("model")->paramOn(mzframe);
	  mzframe->Write(Form("prmfr%d%d",i,(i+1)));
	}

      for (int i =0; i<10; i++)
	{
	  w = new RooWorkspace(Form("nprw%d",i));
	  RooPlot* mzframe = mass->frame();
	  RooDataSet* nprdz = (RooDataSet*) nprdata->reduce(Form("z >= %f && z<%f", i*0.1, (i+1)*0.1));
	  w->import(*nprdz);
	  w->import(*model);
	  w->data("nprdata")->plotOn(mzframe, DataError(RooAbsData::SumW2));
	  w->pdf("model")->fitTo(*w->data("nprdata"),  Extended(!isMc), SumW2Error(kTRUE));
	  w->pdf("model")->plotOn(mzframe);
	  xaxis=np->GetXaxis();
	  np->SetBinContent(xaxis->FindBin(i*0.1), w->var("fs")->getValV());
	  np->SetBinError(xaxis->FindBin(i*0.1), w->var("fs")->getError ());
	  TPaveText* tbox = new TPaveText(0.15,0.6,0.5,0.7, "BRNDC");
	  tbox-> AddText (Form("%.1f #leq z < %.1f", i*0.1, (i+1)*0.1));
	  mzframe->addObject(tbox);
	  tbox->SetBorderSize(0);
	  tbox->SetFillColor(0);
	  w->pdf("model")->paramOn(mzframe);
	  mzframe->Write(Form("nprmfr%d%d",i,(i+1)));
	}
    }



  mframe->Write("massframe");
  rframe->Write("rframe");
  ptframe->Write("recopt");
  rapframe->Write("recorap");
 if (isMc)
   {
     RooPlot* genzframe= genzed->frame();
     gend->plotOn(genzframe, MarkerColor(kBlue), LineColor(kOrange+2), LineStyle(1));
     //TLegend* l = new TLegend (0.2,0.6,0.4,0.8);
     //l->AddEntry(zframe->RooPlot::findObject("recod"), "N_{reco J/#psi} with corrections", "L");
     //l->AddEntry(genzframe->RooPlot::findObject("gend"), "N_{reco J/#psi} with corrections", "L");
     //l->SetBorderSize(0);
     //genzframe->addObject(l);
     //genzframe->Write("genz");
     gn->Write("nv");
     ga->Write("av");
   }
 else
   {
     RooPlot* prmass=mass->frame();
     RooPlot* nprmass=mass->frame();
     prdata->plotOn(prmass, DataError(RooAbsData::SumW2));
     nprdata->plotOn(nprmass, DataError(RooAbsData::SumW2));
     prmass->Write("prmassframe");
     nprmass->Write("nprmassframe");
   }

 nz->Scale(1.0/(nz->Integral()));
 RooDataHist* zd = new RooDataHist ("zd", "N_{J/#psi} in z bins", *zed, nz);
 if (isMc)
   {
   zd->plotOn(zframe, Name("zmc"), MarkerStyle(kFullSquare), MarkerColor(kBlue));
   zframe->Write("zframe");
   }
 else
   {
     RooPlot* nprzfr=zed->frame();
     np->Scale(1.0/(np->Integral()));
     zd->plotOn(zframe, Name("prz"), MarkerColor(kRed));
     zd= new RooDataHist ("nzd", "N_{J/#psi} in z bins", *zed, np);
     zframe->Write("prz");
     zd->plotOn(nprzfr, Name("nprz"), MarkerColor(kRed));
     nprzfr->Write("nprz");
   }
 fsave->Close();
}



Bool_t myTree::isTriggerMatch (Int_t iRecoQQ, Int_t TriggerBit) 
  {
    Bool_t cond = true;
    cond = cond && ( (HLTriggers&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) ); 
    cond = cond && ( (Reco_QQ_trig[iRecoQQ]&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) );
    return cond;
  };

Bool_t myTree::isGlobalMuonInAccept2015 (TLorentzVector* Muon) 
  {
  return (fabs(Muon->Eta()) < 2.4 &&
          ((fabs(Muon->Eta()) < 1.2 && Muon->Pt() >= 3.5) ||
           (1.2 <= fabs(Muon->Eta()) && fabs(Muon->Eta()) < 2.1 && Muon->Pt() >= 5.77-1.89*fabs(Muon->Eta())) ||
           (2.1 <= fabs(Muon->Eta()) && Muon->Pt() >= 1.8)));
  };

Bool_t myTree::areMuonsInAcceptance2015 (Int_t iRecoQQ)
  {
    TLorentzVector *RecoQQmupl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iRecoQQ);
    TLorentzVector *RecoQQmumi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iRecoQQ);
    return ( isGlobalMuonInAccept2015(RecoQQmupl) && isGlobalMuonInAccept2015(RecoQQmumi) );
  };

Bool_t myTree::areGenMuonsInAcceptance2015 (Int_t iGenQQ)
  {
    TLorentzVector *GenQQmupl = (TLorentzVector*) Gen_QQ_mupl_4mom->At(iGenQQ);
    TLorentzVector *GenQQmumi = (TLorentzVector*) Gen_QQ_mumi_4mom->At(iGenQQ);
    return (isGlobalMuonInAccept2015(GenQQmupl) && isGlobalMuonInAccept2015(GenQQmumi));
  };  
  
Bool_t myTree::passQualityCuts2015 (Int_t iRecoQQ) 
  {
    Bool_t cond = true;
    
    // cond = cond && (Reco_QQ_mumi_highPurity[iRecoQQ]);
    cond = cond && (Reco_QQ_mumi_isGoodMuon[iRecoQQ]==1);
    cond = cond && (Reco_QQ_mumi_nTrkWMea[iRecoQQ] > 5);
    cond = cond && (Reco_QQ_mumi_nPixWMea[iRecoQQ] > 0);
    cond = cond && (fabs(Reco_QQ_mumi_dxy[iRecoQQ]) < 0.3);
    cond = cond && (fabs(Reco_QQ_mumi_dz[iRecoQQ]) < 20.);
    
    // cond = cond && (Reco_QQ_mupl_highPurity[iRecoQQ]);
    cond = cond && (Reco_QQ_mupl_isGoodMuon[iRecoQQ]==1);
    cond = cond && (Reco_QQ_mupl_nTrkWMea[iRecoQQ] > 5);
    cond = cond && (Reco_QQ_mupl_nPixWMea[iRecoQQ] > 0);
    cond = cond && (fabs(Reco_QQ_mupl_dxy[iRecoQQ]) < 0.3);
    cond = cond && (fabs(Reco_QQ_mupl_dz[iRecoQQ]) < 20.);    
    cond = cond && (Reco_QQ_VtxProb[iRecoQQ] > 0.01);
    
    return cond;
  };


Double_t myTree::deltaR(TLorentzVector* GenMuon, TLorentzVector* RecoMuon)
{
  double dEta = RecoMuon->Eta() - GenMuon->Eta();
  double dPhi = TVector2::Phi_mpi_pi(RecoMuon->Phi() - GenMuon->Phi());
  return ((double) TMath::Sqrt( (dEta*dEta) + (dPhi*dPhi) ) );
};

Bool_t myTree::isMatchedRecoDiMuon(int iRecoDiMuon, double maxDeltaR =0.03)
{
  TLorentzVector* RecoMuonpl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iRecoDiMuon);
  TLorentzVector* RecoMuonmi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iRecoDiMuon);
  
  bool isMatched(false);
  int iGenMuon(0);
  while ( !isMatched && (iGenMuon < Gen_QQ_size) )
  {
    TLorentzVector *GenMuonpl = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGenMuon);
    TLorentzVector *GenMuonmi = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGenMuon);
    double dRpl = deltaR(GenMuonpl,RecoMuonpl);
    double dRmi = deltaR(GenMuonmi,RecoMuonmi);
    if ( (dRpl < maxDeltaR) && (dRmi < maxDeltaR)  ) isMatched = true;
    iGenMuon++;
  }
  
  return isMatched;
};



Bool_t myTree::isMatchedGenDiMuon(int iGenDiMuon, double maxDeltaR =0.03)
{
  TLorentzVector* GenMuonpl = (TLorentzVector*) Gen_QQ_mupl_4mom->At(iGenDiMuon);
  TLorentzVector* GenMuonmi = (TLorentzVector*) Gen_QQ_mumi_4mom->At(iGenDiMuon);
  TLorentzVector* GenMuon = (TLorentzVector*) Gen_QQ_4mom->At(iGenDiMuon);

  bool isMatched(false);
  int iRecoMuon(0);
  while ( !isMatched && (iRecoMuon < Reco_QQ_size) )
  {
    TLorentzVector *RecoMuonpl = (TLorentzVector*)Reco_QQ_mupl_4mom->At(iRecoMuon);
    TLorentzVector *RecoMuonmi = (TLorentzVector*)Reco_QQ_mumi_4mom->At(iRecoMuon);
    double dRpl = deltaR(GenMuonpl,RecoMuonpl);
    double dRmi = deltaR(GenMuonmi,RecoMuonmi);
    if ( (dRpl < maxDeltaR) && (dRmi < maxDeltaR) ) 
      {
	isMatched = true;
	matchReco = (TLorentzVector*) Reco_QQ_4mom->At(iRecoMuon);
      }
    iRecoMuon++;
  }
  
  return isMatched;
};
