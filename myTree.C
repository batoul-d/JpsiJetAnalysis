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
  int k;
  Double_t ptbins []={6, 9, 12, 15, 20, 30, 50};
  Double_t etabins []={-2.4, -2, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2, 2.4};
  TLorentzVector matchjet;

  TH1F* zz= new TH1F ("z","z gen after", 15, 0, 1.5);
  TH1F* zmr= new TH1F ("zmr", "matched reco after", 15, 0, 1.5);
  TH1F* zmg= new TH1F ("zmg", "matched gen after", 15, 0, 1.5);
  TH1F* ptg= new TH1F ("ptg", "N_{gen} vs p_{T}; p_{T}; N_{total}", 6, ptbins);
  TH1F* ptgm= new TH1F ("ptgm", "N_{matched gen} vs p_{T}; p_{T}; N_{matched}", 6, ptbins);
  TH1F* etagm= new TH1F ("etagm","N_{matched gen} vs #eta; #eta; N_{matched}",6, 0, 2.4);
  TH1F* etag= new TH1F ("etag","N_{gen} vs #eta; #eta; N_{total}",6, 0, 2.4);
  TH2F* ptetag= new TH2F ("ptetag", "N_{matched gen} vs p_{T} and #eta; #eta; p_{T}; N_{matched}", 12, etabins, 6, ptbins);
  TH2F* ptetamg= new TH2F ("ptetamg", "N_{gen} vs p_{T} and #eta; #eta; p_{T}; N_{total}", 12, etabins, 6, ptbins);
  TH1F* opan= new TH1F ("opan", "opening angle between J/#psi and jet; angle; Events", 20, 0, 0.6);


  RooDataSet* data = NULL;
  RooDataSet* unwdata =NULL;
  RooArgSet* set = NULL;
  RooArgSet* unwset = NULL;
  RooDataSet* gendata = NULL;
  RooArgSet* genset = NULL;
  RooRealVar* mass = new RooRealVar("invMass","#mu#mu mass", 2.6, 3.4, "GeV/c^{2}");
  RooRealVar* genmass = new RooRealVar("geninvMass","#mu#mu mass", 2.6, 3.4, "GeV/c^{2}");
  RooRealVar* ctau = new RooRealVar("ctau","c_{#tau}", -100000.0, 100000.0, "mm");
  RooRealVar* zed = new RooRealVar("z", "z_{J/#psi}", 0, 1);
  zed->setBins(10);
  RooRealVar* r = new RooRealVar("r", "#DeltaR", 0, 5);
  RooRealVar* genzed = new RooRealVar("genz", "genz_{J/#psi}", 0, 1);
  RooRealVar* weight  = new RooRealVar("weight","MC weight", 0.0, 10000000.0, "");
  genzed->setBins(10);
  mass->setBins(40);
  RooRealVar* genr = new RooRealVar("genr", "#genDeltaR", 0, 5);


  set = new RooArgSet(*mass, *ctau, *zed, *r, *weight);
  unwset = new RooArgSet(*mass, *ctau, *zed, *r);
  genset = new RooArgSet (*genmass, *genzed, *genr);
  data = new RooDataSet ("data", "data for reconstructed J/#psi", *set, WeightVar(*weight));
  unwdata = new RooDataSet ("unwdata", "unweighted data", *unwset);
  gendata = new RooDataSet ("gendata", "data for generated J/#psi", *genset);

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
      
	  if (HLT_HIL1DoubleMu0ForPPRef_v1)// && pPAprimaryVertexFilter)
	    {
	      for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++) 
		{
		  TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
		  TLorentzVector *GenQQmupl4mom = (TLorentzVector*) Gen_QQ_mupl_4mom->At(iQQ);
		  TLorentzVector *GenQQmumi4mom = (TLorentzVector*) Gen_QQ_mumi_4mom->At(iQQ);
		  jpsi_pt = GenQQ4mom->Pt();
		  jpsi_eta =GenQQ4mom->Eta();
		  if (
		      jpsi_pt > 6 &&
		      (areGenMuonsInAcceptance2015(iQQ))
		      //(isMatchedGenDiMuon(iQQ))		   
		      )
		    {
		      genmass->setVal(GenQQ4mom->M());
		      if (jpsi_eta>(-2.4) && jpsi_eta < 2.4) 
			{
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
				      //if (drmin < 0.6 )
				      z= jpsi_pt/genpt[ijet];
				    }
				}
			    }
			  genzed->setVal(z);
			  genr->setVal(drmin);
			  if (drmin != 10000)
			    {
			      ptg->Fill(jpsi_pt);
			      zz->Fill(z);
			      etag->Fill(abs(jpsi_eta));
			      ptetag->Fill(jpsi_eta, jpsi_pt);
			      if (isMatchedGenDiMuon(iQQ))
				{
				  ptgm->Fill(jpsi_pt);
				  zmg->Fill(z);
				  ptetamg->Fill(jpsi_eta, jpsi_pt);
				  etagm->Fill(abs(jpsi_eta));
				}
			    }
			}
		    }
		}
	      gendata->add(*genset);

	      TEfficiency* gptef = new TEfficiency("gptef", "reconstruction efficiency fct of pt", 6, ptbins);
	      if(TEfficiency::CheckConsistency(*ptgm,*ptg))
		gptef = new TEfficiency (*ptgm,*ptg);
	      else 
		cout<<endl<<"non consistent";

	      TEfficiency* gzef = new TEfficiency("gzef", "reconstruction efficiency fct of z",10, 0, 10);
	      if(TEfficiency::CheckConsistency(*zmg,*zz))
		gzef = new TEfficiency (*zmg,*zz);
	      else 
		cout<<endl<<"non consistent";

	      TEfficiency* getaef = new TEfficiency("getaef", "reconstruction efficiency fct of eta", 6, 0, 2.4);
	      if(TEfficiency::CheckConsistency(*etagm,*etag))
		getaef = new TEfficiency (*etagm,*etag);
	      else 
		cout<<endl<<"non consistent";

	      TEfficiency* gptetaef = new TEfficiency("gptetaef", "reconstruction efficiency fct of pt and eta; eta; pt; eff", 12, etabins, 6, ptbins);
	      if(TEfficiency::CheckConsistency(*ptetamg,*ptetag))
		gptetaef = new TEfficiency (*ptetamg,*ptetag);
	      else 
		cout<<endl<<"non consistent";

	      TFile ef ("npefficiencies.root", "RECREATE");
	      gptef->Write();
	      gzef->Write();
	      gptetaef->Write();
	      getaef->Write();
	      //gptetaef= (TEfficiency*) ef.FindObjectAny("ptetag_clone");
	      ef.Close();
		

	      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) 
		{
		  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
		  TLorentzVector *RecoQQmupl4mom = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iQQ);
		  TLorentzVector *RecoQQmumi4mom = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iQQ);
		  mass->setVal(RecoQQ4mom->M());
		  jpsi_pt = RecoQQ4mom->Pt();
		  jpsi_eta = RecoQQ4mom->Eta();
		  if (
		      jpsi_pt > 6  &&
		      (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		      (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		      (isTriggerMatch(iQQ, triggerIndex_PP)) &&// if it matches the trigger 
		      ( isMatchedRecoDiMuon(iQQ))
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
				      //if (drmin < 0.6 )
				      z= jpsi_pt/jtpt[ijet];
				    }
				}
			    }
			  zmr->Fill(z);
			  zed->setVal(z);
			  r->setVal(drmin);
			  //weight->setVal(2);
			  weight->setVal(1/(gptetaef->GetEfficiency(gptetaef->GetGlobalBin(jpsi_eta, jpsi_pt))));
			  if (drmin<0.4)
			    opan->Fill(drmin);
			}
		    }
		}
	    }	    
	  data->add(*set);
	  unwdata->add(*unwset);
	}
	
      TFile *datafile = TFile::Open("MCsample","RECREATE");
      datafile->cd();
      data->Write("MCdata");
      unwdata->Write("unwdata");
      gendata->Write("MCgendata");
      datafile->Write();
      datafile->Close();
      delete datafile;

      //TFile fsa("mg.root","RECREATE");
      //zmg->Write();
      //fsa.Close();
      //TFile fsb ("mr.root", "RECREATE");
      //zmr->Write();
      //fsb.Close();
      //TFile npr("nonprompt.root","RECREATE");
      //opan->Write();
      //npr.Close();
	
    }

  else 
    {
      TFile datafile ("MCsample");
      data= (RooDataSet*) datafile.Get("MCdata");
      unwdata= (RooDataSet*) datafile.Get("unwdata");
      gendata= (RooDataSet*) datafile.Get("MCgendata");
	
      RooDataSet* d = (RooDataSet*) data->reduce("invMass>2.6 && invMass<3.4");
      RooDataSet* unwd = (RooDataSet*) unwdata->reduce("invMass>2.6 && invMass<3.4");
      RooDataSet* gend = (RooDataSet*) gendata->reduce("geninvMass>2.6 && geninvMass<3.4");
      RooPlot* mframe = mass->frame();
      RooPlot* genframe = genmass->frame();
      RooPlot* comp = mass->frame();
      RooPlot* zframe =zed->frame();
      RooPlot* rframe =r->frame();
      RooPlot* genzframe= genzed->frame();
      RooRealVar mean ("mean", "mean", 3.096, 3, 3.2);
      RooRealVar sigma ("sigma", "sigma", 0.03, 0, 0.1);
      RooRealVar alpha ("alpha", "alpha", 1.94, 0, 3);
      RooRealVar n ("n", "n", 1.64, 0, 3);
      RooRealVar sigma1 ("sigma1", "sigma1", 0.03, 0, 0.1);


      RooCBShape masspeak ("masspeak", "fit of invariant mass", *mass, mean, sigma, alpha, n);
      RooCBShape cb ("cb", "cb", *mass, mean, sigma, alpha, n);
      RooCBShape cb1 ("cb1", "cb1", *mass, mean, sigma1, alpha, n);

      RooCBShape gencb ("cb", "cb", *genmass, mean, sigma, alpha, n);
      RooCBShape gencb1 ("cb1", "cb1", *genmass, mean, sigma1, alpha, n);


      RooRealVar c0("c0","coefficient #0", 0.5,-1,1) ;
      RooRealVar c1("c1","coefficient #1", 0.1,-1,1);
      RooRealVar c2("c2","coefficient #2",-0.1,-10,10);
      RooChebychev bkg("bkg","background p.d.f.",*mass,RooArgList(c0,c1,c2)) ;
      RooRealVar fb ("fb", "background nb", 5000, 0, data->numEntries());
      RooRealVar fs ("fs", "signal nb", 5000, 0,  data->numEntries());
      RooRealVar fcb ("fcb", "cb", 5000, 0, data->numEntries());
      RooRealVar fcb1 ("fcb1", "cb1", 5000, 0, data->numEntries());

      RooRealVar genfcb ("genfcb", "gencb", 5000, 0, gendata->numEntries());
      RooRealVar genfcb1 ("genfcb1", "gencb1", 5000, 0, gendata->numEntries());

      mass->setRange("window",2.9, 3.2) ;
      genmass->setRange("window1", 2.9, 3.2);
      RooExtendPdf es ("es", "es", masspeak, fs, "window");
      RooExtendPdf eb ("eb", "eb", bkg, fb);
      RooExtendPdf ecb ("ecb", "ecb", cb, fcb);
      RooExtendPdf ecb1 ("ecb1", "ecb1", cb1, fcb1);

      RooExtendPdf genecb ("genecb", "genecb", gencb, genfcb);
      RooExtendPdf genecb1 ("genecb1", "genecb1", gencb1, genfcb1);
 
      RooAddPdf model ("model", "model", RooArgList(ecb, ecb1));
      RooAddPdf genmodel ("genmodel", "genmodel", RooArgList(genecb, genecb1));


      d->plotOn(mframe,DataError(RooAbsData::SumW2));
      d->plotOn(comp,DataError(RooAbsData::SumW2));
      gend->plotOn(genframe);
      unwd->plotOn(comp);
      data->plotOn(zframe);
      gendata->plotOn(genzframe);
      data->plotOn(rframe);

      model.fitTo(*d, Extended(kTRUE), SumW2Error(kTRUE));
      model.plotOn(mframe,Normalization(d->numEntries(), RooAbsReal::NumEvent), Name("total"), LineColor(kRed));
      model.paramOn(mframe);

      model.plotOn(comp ,Normalization(d->numEntries(), RooAbsReal::NumEvent), Name("total"), LineColor(kRed));
      model.fitTo(*unwd, Extended(kTRUE));
      model.plotOn(comp,Normalization(unwd->numEntries(), RooAbsReal::NumEvent), Name("tot"), LineColor(kBlue), LineStyle(kDashed));

      TFile fsave("file.root","RECREATE");
      mframe->Write("frame_invMass");
      rframe->Write("frame_r");


      for (int i =0; i<10; i++)
	{
	  RooPlot* mzframe = mass->frame();
	  RooWorkspace* w = new RooWorkspace(Form("w%d",i));
	  RooDataSet* dz = (RooDataSet*) d->reduce(Form("z >= %f && z<%f && invMass>2.6 && invMass<3.4", i*0.1, (i+1)*0.1));
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
	  w->pdf("model")->paramOn(mzframe);
	  mzframe->Write(Form("frame_invMass_%d%d",i,(i+1)));
	}
      zframe->Write("frame_z");
      genzframe->Write("genz");
      comp->Write("comparaison");
      genframe->Write("genmass");
      fsave.Close();
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

Bool_t myTree::areGenMuonsInAcceptance2015 (Int_t iGenQQ)
  {
    TLorentzVector *GenQQmupl = (TLorentzVector*) Gen_QQ_mupl_4mom->At(iGenQQ);
    TLorentzVector *GenQQmumi = (TLorentzVector*) Gen_QQ_mumi_4mom->At(iGenQQ);
    return ( isGlobalMuonInAccept2015(GenQQmupl) && isGlobalMuonInAccept2015(GenQQmumi) );
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
  
  bool isMatched(false);
  int iRecoMuon(0);
  while ( !isMatched && (iRecoMuon < Reco_QQ_size) )
  {
    TLorentzVector *RecoMuonpl = (TLorentzVector*)Reco_QQ_mupl_4mom->At(iRecoMuon);
    TLorentzVector *RecoMuonmi = (TLorentzVector*)Reco_QQ_mumi_4mom->At(iRecoMuon);
    double dRpl = deltaR(GenMuonpl,RecoMuonpl);
    double dRmi = deltaR(GenMuonmi,RecoMuonmi);
    if ( (dRpl < maxDeltaR) && (dRmi < maxDeltaR) && (areMuonsInAcceptance2015(iRecoMuon))&& (isTriggerMatch(iRecoMuon, 0))) isMatched = true;
    iRecoMuon++;
  }
  
  return isMatched;
};
