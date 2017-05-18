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
Double_t ptbins []={6.5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 21, 24, 27, 30,40, 50};
Double_t etabins []={-2.4, -2, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2, 2.4};

RooRealVar* mass = new RooRealVar("Mass","#mu#mu mass", 2.6, 3.4, "GeV/c^{2}");
RooRealVar* genmass = new RooRealVar("genMass","#mu#mu mass", 2.6, 3.4, "GeV/c^{2}");
RooRealVar* ctau = new RooRealVar("ctau","c_{#tau}", -100000.0, 100000.0, "mm");
RooRealVar* zed = new RooRealVar("z", "z_{J/#psi}", 0, 1);
RooRealVar* r = new RooRealVar("r", "#DeltaR", 0, 5);
RooRealVar* genzed = new RooRealVar("genz", "genz_{J/#psi}", 0, 1);
RooRealVar* w = new RooRealVar("w","MC weight", 0.0, 10000000.0);
RooRealVar* genr = new RooRealVar("genr", "#genDeltaR", 0, 5);
RooRealVar* pt = new RooRealVar("pt", "pt", 6.5, 25.5);
RooRealVar* rap = new RooRealVar("rap", "rapidity", -2.4, 2.4);
RooRealVar* gen_pt = new RooRealVar("gen_pt", "pt for gen j/psi", 6.5, 25.5);
RooRealVar* gen_rap = new RooRealVar("gen_rap", "rapidity for gen j/psi", -2.4, 2.4);

RooArgSet* recoset = new RooArgSet(*mass, *ctau, *zed, *r, *pt, *rap, *w);
RooArgSet*  unwset = new RooArgSet(*mass, *ctau, *zed, *r, *pt, *rap);
RooArgSet*  genset = new RooArgSet (*genmass, *genzed, *genr, *gen_pt, *gen_rap);
RooDataSet* data = new RooDataSet ("data", "data for reconstructed J/#psi", *recoset, WeightVar(*w));
RooDataSet* unwdata = new RooDataSet ("unwdata", "unweighted data", *unwset);
RooDataSet* gendata = new RooDataSet ("gendata", "data for generated J/#psi", *genset);

void myTree::EffCalc()
{
  if (isMc)
    {
      TH1F* ptg= new TH1F ("ptg", "N_{gen} vs p_{T}; p_{T}; N_{total}", 16, ptbins);
      TH1F* ptr= new TH1F ("ptr", "N_{reco} vs p_{T}; p_{T}; N_{reco}", 16, ptbins);
      TH1F* rapr= new TH1F ("rapr","N_{reco} vs #eta; #eta; N_{reco}", 12, etabins);
      TH1F* rapg= new TH1F ("rapg","N_{gen} vs #eta; #eta; N_{total}", 12, etabins);
      TH2F* ptrapg= new TH2F ("ptrapg", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", 12, etabins, 16, ptbins);
      TH2F* ptrapr= new TH2F ("ptrapr", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", 12, etabins, 16, ptbins);

      Long64_t nentries = fChain->GetEntries();
      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;

	  for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	    {
	      TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	      TLorentzVector *GenQQmupl4mom = (TLorentzVector*) Gen_QQ_mupl_4mom->At(iQQ);
	      TLorentzVector *GenQQmumi4mom = (TLorentzVector*) Gen_QQ_mumi_4mom->At(iQQ);
	      jpsi_pt = GenQQ4mom->Pt();
	      jpsi_rap = GenQQ4mom->Rapidity();

	      if (jpsi_pt>6.5 && abs(jpsi_rap)<2.4)
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
		  TLorentzVector *RecoQQmupl4mom = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iQQ);
		  TLorentzVector *RecoQQmumi4mom = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iQQ);
		  jpsi_pt = RecoQQ4mom->Pt();
		  jpsi_rap = RecoQQ4mom->Rapidity();
		  if (
		      jpsi_pt > 6.5  &&
		      (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		      (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		      (isTriggerMatch(iQQ, triggerIndex_PP)) //&&// if it matches the trigger 
		      //(isMatchedRecoDiMuon(iQQ))
		      )
		    {
		      if (Reco_QQ_sign[iQQ]==0 && abs(jpsi_rap)<2.4) 
			{
			  ptr->Fill(jpsi_pt);
			  rapr->Fill(jpsi_rap);
			  ptrapr->Fill(jpsi_rap, jpsi_pt);

			}
		    }
		}
	    }
	}
      TEfficiency* gptef = new TEfficiency("gptef", "reconstruction efficiency fct of pt", 16, ptbins);
      if(TEfficiency::CheckConsistency(*ptr,*ptg))
	gptef = new TEfficiency (*ptr,*ptg);

      TEfficiency* grapef = new TEfficiency("grapef", "reconstruction efficiency fct of rapidity", 12, etabins);
      if(TEfficiency::CheckConsistency(*rapr,*rapg))
	grapef = new TEfficiency (*rapr,*rapg);

      TEfficiency* gptrapef = new TEfficiency("gptetaef", "reconstruction efficiency fct of pt and rapidity; y; pt; eff", 12, etabins, 16, ptbins);
      if(TEfficiency::CheckConsistency(*ptrapr, *ptrapg))
	gptrapef = new TEfficiency (*ptrapr, *ptrapg);
      TFile* ef (0x0);
      if (isPr)
	ef = new TFile ("prEff.root", "RECREATE");
      else
	ef = new TFile ("nprEff.root", "RECREATE");

      gptef->Write("pt");
      gptrapef->Write("ptrap");
      grapef->Write("rap");
      ef->Close();
    }
  else 
    cout<< "this is data and not MC"<<endl;
}



void myTree::Loop()
{
  if (fChain == 0) return;
  genzed->setBins(10);
  mass->setBins(40);
  zed->setBins(10);
  rap->setBins(12);
  gen_rap->setBins(12);
  TFile* f (0x0);
  if (isPr)
    f=new TFile ("prEff.root");
  else
    f= new TFile ("nprEff.root");
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
	      TLorentzVector *GenQQmupl4mom = (TLorentzVector*) Gen_QQ_mupl_4mom->At(iQQ);
	      TLorentzVector *GenQQmumi4mom = (TLorentzVector*) Gen_QQ_mumi_4mom->At(iQQ);
	      jpsi_pt = GenQQ4mom->Pt();
	      jpsi_rap = GenQQ4mom->Rapidity();

	      if (jpsi_pt>6.5 && abs(jpsi_rap)<2.4)
		{
		  genmass->setVal(GenQQ4mom->M());
		  gen_pt->setVal(jpsi_pt);
		  gen_rap->setVal(jpsi_rap);    
		  drmin = 10000;
		  for (Long64_t ijet=0; ijet<ngen; ijet++)
		    {
	  
		      if (geneta[ijet]>(-2.4) && geneta[ijet]<(2.4) && genpt[ijet]>=20)
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
			      //if (drmin < 0.5)
			      z= jpsi_pt/genpt[ijet];
			    }
			}
		    }
		  genzed->setVal(z);
		  genr->setVal(drmin);
		  gendata->add(*genset);
		}
	    }
	  if ( HLT_HIL1DoubleMu0ForPPRef_v1 && pPAprimaryVertexFilter)
	    {
	      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) 
		{
		  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
		  TLorentzVector *RecoQQmupl4mom = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iQQ);
		  TLorentzVector *RecoQQmumi4mom = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iQQ);
		  jpsi_pt = RecoQQ4mom->Pt();
		  jpsi_eta = RecoQQ4mom->Eta();
		  jpsi_rap = RecoQQ4mom->Rapidity();
		  if (
		      jpsi_pt > 6.5  &&
		      (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		      (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		      (isTriggerMatch(iQQ, triggerIndex_PP)) //&&// if it matches the trigger 
		      )
		    {
		      if (Reco_QQ_sign[iQQ]==0 && abs(jpsi_rap)<2.4) 
			{
			  drmin = 1000;
			  for (Long64_t ijet=0; ijet<nref; ijet++)
			    {
	  
			      if (jteta[ijet]>(-2.4) && jteta[ijet]<(2.4) && jtpt[ijet]>=20)
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
				      //if (drmin < 0.5)
				      z= jpsi_pt/jtpt[ijet];
				    }
				}
			    }
			  mass->setVal(RecoQQ4mom->M());
			  pt->setVal(jpsi_pt);
			  rap->setVal(jpsi_rap);
			  zed->setVal(z);
			  r->setVal(drmin);
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
		  TLorentzVector *RecoQQmupl4mom = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iQQ);
		  TLorentzVector *RecoQQmumi4mom = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iQQ);
		  jpsi_pt = RecoQQ4mom->Pt();
		  jpsi_eta = RecoQQ4mom->Eta();
		  jpsi_rap = RecoQQ4mom->Rapidity();
		  if (
		      jpsi_pt > 6.5  &&
		      (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		      (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		      (isTriggerMatch(iQQ, triggerIndex_PP)) // if it matches the trigger 
		      )
		    {
		      if (Reco_QQ_sign[iQQ]==0 && abs(jpsi_rap)<2.4) 
			{
			  drmin = 1000;
			  for (Long64_t ijet=0; ijet<nref; ijet++)
			    {
	  
			      if (jteta[ijet]>(-2.4) && jteta[ijet]<(2.4) && jtpt[ijet]>=20)
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
				      //if (drmin < 0.5)
				      z= jpsi_pt/jtpt[ijet];
				    }
				}
			    }
			  mass->setVal(RecoQQ4mom->M());
			  pt->setVal(jpsi_pt);
			  rap->setVal(jpsi_rap);
			  zed->setVal(z);
			  r->setVal(drmin);
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
  TH1F* gn = new TH1F ("gn", "values of n from mc", 12, 0, 1.2);
  TH1F* ga = new TH1F ("ga", "values of alpha from mc", 12, 0, 1.2);
  TH1F* nz = new TH1F ("nz", "N_{J/#psi} in z bins", 10, 0, 1);
  if (isPr)
    {
      mcfile = new TFile ("prMcDatasets");
      datafile = new TFile ("datasets_pr");
    }
  else
    {
      mcfile = new TFile ("nprMcDatasets");
      datafile = new TFile ("datasets_npr");
    }
  gendata= (RooDataSet*) mcfile->Get("genData");

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
  RooDataSet* d = (RooDataSet*) data->reduce("Mass>2.6 && Mass<3.4 && r<0.5");
  RooDataSet* unwd = (RooDataSet*) unwdata->reduce("Mass>2.6 && Mass<3.4 && r<0.5");
  RooDataSet* gend = (RooDataSet*) gendata->reduce("genMass>2.6 && genMass<3.4 && genr<0.5");

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
      alpha = new RooRealVar ("alpha", "alpha", 1, 1, 3);
      n = new RooRealVar("n", "n", 1, 1, 3);
    }
  else
    {
      alpha = new RooRealVar ("alpha", "alpha", 1.8);
      n = new RooRealVar("n", "n", 1.7);
    }

  RooCBShape cb ("cb", "cb", *mass, *mean, *sigma, *alpha, *n);
  RooCBShape cb1 ("cb1", "cb1", *mass, *mean, *sigma1, *alpha, *n);

  RooRealVar c0("c0","coefficient #0", -0.5,-1,1);
  RooRealVar c1("c1","coefficient #1", 0.06,-1,1);
  RooRealVar c2("c2","coefficient #2",-0.04,-10,10);
  RooChebychev bkg("bkg","background p.d.f.", *mass, RooArgList(c0,c1,c2)) ;
  RooRealVar fb ("fb", "background nb", 100000, 0, data->numEntries());
  RooRealVar fs ("fs", "signal nb", 5000, 0, data->numEntries());
  RooRealVar fcb ("fcb", "cb", 0.2, 0, 1);
  //RooRealVar fcb1 ("fcb1", "cb1", 0.8, 0, 1);

  mass->setRange("window",2.9, 3.2) ;

  RooExtendPdf ecb ("ecb", "ecb", cb, fcb);
  //RooExtendPdf ecb1 ("ecb1", "ecb1", cb1, fcb1);

  RooAddPdf masspeak ("masspeak", "fit of mass peak", RooArgList(cb, cb1), fcb);
  RooExtendPdf es ("es", "es", masspeak, fs, "window");
  RooExtendPdf eb ("eb", "eb", bkg, fb);
  RooAddPdf* model;
  if (isMc)
    model = new RooAddPdf ("model", "model", RooArgList(cb, cb1), fcb);
  else 
    model = new RooAddPdf ("model", "model", RooArgList(es, eb));


  d->plotOn(mframe,DataError(RooAbsData::SumW2));
  data->plotOn(zframe);
  data->plotOn(rframe);
  data->plotOn(ptframe, Name("recod"));
  data->plotOn(rapframe);
 
  model->fitTo(*d, Extended(!isMc), SumW2Error(kTRUE));
  if (isMc)
    {
      gn->Fill(1, n->getValV());
      ga->Fill(1, alpha->getValV());
    }
  model->plotOn(mframe, Name("background"), Components(RooArgSet(bkg)),DrawOption("F"), FillColor(kGray));
  model->plotOn(mframe, Name("signal"), Components(RooArgSet(masspeak)), LineStyle(kDashed));
  model->plotOn(mframe, Name("total"), LineColor(kRed));
  model->paramOn(mframe);

  RooWorkspace* unw = new RooWorkspace ("unw");
  unw->import(*unwd);
  unw->import(*model);
  unw->data("unwdata")->plotOn(comp);
  unw->pdf("model")->fitTo(*unw->data("unwdata"), Extended(!isMc));
  unw->pdf("model")->plotOn(comp, Name("background"), Components(RooArgSet(bkg)),DrawOption("F"), FillColor(kGray));
  unw->pdf("model")->plotOn(comp, LineColor(kBlue), LineStyle(kDashed));
  unw->pdf("model")->paramOn(comp);
  if (isMc)
    {
      gn->Fill(1.1, n->getValV());
      ga->Fill(1.1, alpha->getValV());
    }


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
  mframe->Write("massframe");
  rframe->Write("rframe");

  for (int i =0; i<10; i++)
    {
      RooPlot* mzframe = mass->frame();
      RooWorkspace* w = new RooWorkspace(Form("w%d",i));
      RooDataSet* dz = (RooDataSet*) d->reduce(Form("z >= %f && z<%f", i*0.1, (i+1)*0.1));
      w->import(*dz);
      w->import(*model);
      w->data("data")->plotOn(mzframe, DataError(RooAbsData::SumW2));
      w->pdf("model")->fitTo(*w->data("data"),  Extended(!isMc), SumW2Error(kTRUE));
      w->pdf("model")->plotOn(mzframe);
      if (isMc)
	{
	  gn->Fill(i*0.1, w->var("n")->getValV());
	  ga->Fill(i*0.1, w->var("alpha")->getValV()); 
	}
      else 
	nz->Fill(i*0.1, w->var("fs")->getValV());
      TPaveText* tbox = new TPaveText(0.15,0.6,0.5,0.7, "BRNDC");
      tbox-> AddText (Form("%.1f #leq z < %.1f", i*0.1, (i+1)*0.1));
      mzframe->addObject(tbox);
      tbox->SetBorderSize(0);
      tbox->SetFillColor(0);
      w->pdf("model")->paramOn(mzframe);
      mzframe->Write(Form("massframe%d%d",i,(i+1)));
    }
  zframe->Write("zframe");
  comp->Write("unweighted");
  ptframe->Write("recopt");
  rapframe->Write("recorap");
 if (isMc)
   {
     RooPlot* genzframe= genzed->frame();
     RooPlot* genptfr = gen_pt->frame();
     RooPlot* genrapfr = gen_rap->frame();
     gendata->plotOn(genzframe,  MarkerColor(kBlue));
     gendata->plotOn(genptfr, Name("gend"), MarkerColor(kBlue), LineColor(kOrange+2), LineStyle(1), Precision(1e-4));
     gendata->plotOn(genrapfr, MarkerColor(kBlue), LineColor(kOrange+2), LineStyle(1), Precision(1e-4));
     TLegend* l = new TLegend (0.2,0.6,0.4,0.8);
     l->AddEntry(ptframe->RooPlot::findObject("recod"), "N_{reco J/#psi}", "LPE");
     l->AddEntry(genptfr->RooPlot::findObject("gend"), "N_{gen J/#psi}", "LPE");
     l->SetBorderSize(0);
     genptfr->addObject(l);
     genzframe->Write("genz");
     genptfr->Write("genpt");
     genrapfr->Write("genrap");
     gn->Write("nv");
     ga->Write("av");
   }
 else
   nz->Write("nvsz");
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
