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

void myTree::Loop(int q=0)
{
  if (fChain == 0) return;
  Float_t  jpsi_pt;
  Float_t  jpsi_eta;
  Float_t  jpsi_phi;
  Float_t dr;
  Float_t dphi;
  Float_t dphimin;
  Float_t deta;
  Float_t drmin; 
  Float_t z=100;
  int triggerIndex_PP =0;
  Float_t jz[10];
  Float_t zb[10];

  RooDataSet* data = NULL;
  RooRealVar* mass = new RooRealVar("invMass","#mu#mu mass", 2.6, 3.4, "GeV/c^{2}");
  RooRealVar* ctau = new RooRealVar("ctau","c_{#tau}", -100000.0, 100000.0, "mm");
  RooRealVar* zed = new RooRealVar("z", "z_{J/#psi}", 0, 1);
  RooRealVar* r = new RooRealVar("r", "#DeltaR", 0, 5);
  RooArgSet* set = NULL;
  set = new RooArgSet(*mass, *ctau, *zed, *r);
  data = new RooDataSet ("data", "data", *set);
  Long64_t nentries = 100000; //fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  if (q==0)
    {
      for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;

	  //if (Cut(ientry) < 0) continue;
      
	  if (HLT_HIL1DoubleMu0_v1 && pPAprimaryVertexFilter)
	    {
	      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) {
		TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
		TLorentzVector *RecoQQmupl4mom = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iQQ);
		TLorentzVector *RecoQQmumi4mom = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iQQ);
		mass->setVal(RecoQQ4mom->M());
		jpsi_pt = RecoQQ4mom->Pt();
		if (
		    jpsi_pt > 6  &&
		    (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		    (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		    (isTriggerMatch(iQQ, triggerIndex_PP)) // if it matches the trigger 
		    )
		  {
		    if (Reco_QQ_sign[iQQ]==0 && jpsi_eta>(-2.4) && jpsi_eta < 2.4) 
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
				    if (drmin < 0.6 )
				      z= jpsi_pt/jtpt[ijet];
				  }
			      }
			  }
			zed->setVal(z);
			r->setVal(drmin);
		      }
		  }
	      }
	    }
	  data->add(*set);
	}

      TFile *datafile = TFile::Open("datas","RECREATE");
      datafile->cd();
      data->Write("data");
      datafile->Write();
      datafile->Close();
      delete datafile;
    }

  else 
    {
      TFile datafile ("datas");
      data= (RooDataSet*) datafile.Get("data");

      RooDataSet* d = (RooDataSet*) data->reduce("invMass>2.6 && invMass<3.4");
      RooPlot* mframe = mass->frame();
      RooPlot* zframe =zed->frame();
      RooPlot* rframe =r->frame();
      RooRealVar mean ("mean", "mean", 3.096, 3, 3.2);
      RooRealVar sigma ("sigma", "sigma", 0.03, 0, 0.1);
      RooRealVar alpha ("alpha", "alpha", 1.94);
      RooRealVar n ("n", "n", 1.64);
      RooCBShape masspeak ("masspeak", "fit of invariant mass", *mass, mean, sigma, alpha, n);
      RooRealVar c0("c0","coefficient #0", 1.0,-1,1) ;
      RooRealVar c1("c1","coefficient #1", 0.1,-1,1);
      RooRealVar c2("c2","coefficient #2",-0.1,-100,100);
      RooChebychev bkg("bkg","background p.d.f.",*mass,RooArgList(c0,c1,c2)) ;
      RooRealVar fb ("fb", "background nb", 5000, 0, data->numEntries());
      RooRealVar fs ("fs", "signal nb", 5000, 0,  data->numEntries());
      RooExtendPdf es ("es", "es", masspeak, fs);
      RooExtendPdf eb ("eb", "eb", bkg, fb);
      RooAddPdf model ("model", "model", RooArgList(es, eb));
      d->plotOn(mframe);
      data->plotOn(zframe);
      data->plotOn(rframe);
      model.fitTo(*d);
      model.plotOn(mframe, Name("background"), Components(RooArgSet(bkg)),DrawOption("F"), FillColor(kGray));
      model.plotOn(mframe, Name("signal"), Components(RooArgSet(masspeak)), LineStyle(kDashed));
      model.plotOn(mframe, Name("total"), LineColor(kRed));
      TLegend *leg = new TLegend(0.2,0.6,0.4,0.8);
      leg->AddEntry(mframe->RooPlot::findObject("signal"), "J/#psi distribution", "l");
      leg->AddEntry(mframe->RooPlot::findObject("background"), "background", "f");
      leg->AddEntry(mframe->RooPlot::findObject("total"), "total distribution", "l");
      leg->SetBorderSize(0);
      mframe->addObject(leg);
      masspeak.paramOn(mframe);
      TPaveText* tb = new TPaveText(0.2,0.35,0.4,0.45, "BRNDC");
      tb-> AddText (Form("#splitline{N_{J/#psi}=%.2f}{N_{bkg}=%.2f}", fs.getValV(),fb.getValV()));
      mframe->addObject(tb);
      tb->SetBorderSize(0);
      tb->SetFillColor(0);
      TFile fsave("f.root","RECREATE");
      mframe->Write();
      rframe->Write();


      for (int i =0; i<10; i++)
      {
	RooPlot* mzframe = mass->frame();
	RooWorkspace* w = new RooWorkspace(Form("w%d",i));
	RooDataSet* dz = (RooDataSet*) data->reduce(Form("z >= %f && z<%f && invMass>2.6 && invMass<3.4", i*0.1, (i+1)*0.1));
	w->import(*dz);
	w->import(model);
	(w->data("data"))->plotOn(mzframe);
	w->pdf("model")->fitTo(*w->data("data"));
	w->pdf("model")->plotOn(mzframe);
	TPaveText* tbox = new TPaveText(0.15,0.6,0.5,0.7, "BRNDC");
	tbox-> AddText (Form("%.1f #leq z < %.1f", i*0.1, (i+1)*0.1));
	mzframe->addObject(tbox);
	tbox->SetBorderSize(0);
	tbox->SetFillColor(0);
	model.paramOn(mzframe);
	mzframe->Write();
	zb[i]=i*0.1+0.05;
	jz[i]=(w->var("fs"))->getValV();
      }
      zframe->Write();
      fsave.Close();

      //TGraph* gr = new TGraph(10, zb, jz);
      TH1F *gr =new TH1F ("gr", "gr", 10, 0, 1);
      for (int i =0; i<10; i++)
	{
	  cout <<jz[i]<<"  ";
	  gr->Fill(i*0.1, jz[i]);
	}
      gr->Draw();
    }
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
