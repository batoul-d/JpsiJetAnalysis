//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar 23 15:36:01 2017 by ROOT version 6.02/13
// from TTree myTree/My TTree of dimuons
// found on file: /data_CMS/cms/mnguyen/jPsiJet/HiForestAOD.root
//////////////////////////////////////////////////////////

// /grid_mnt/vol__vol_U__u/llr/cms/mnguyen/prod/forest/compositeJet/CMSSW_7_5_8_patch5/src/HeavyIonsAnalysis/JetAnalysis/test/HiForestAODWithJPsi.root

// /grid_mnt/vol__vol_U__u/llr/cms/mnguyen/prod/forest/compositeJet/CMSSW_7_5_8_patch5/src/HeavyIonsAnalysis/JetAnalysis/test/HiForestAODWithoutJPsi.root

// /data_CMS/cms/mnguyen/jPsiJet/v2/merged_HiForestAOD.root


#ifndef myTree_h
#define myTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

using namespace std;

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"

class myTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          eventNb;
   UInt_t          runNb;
   UInt_t          LS;
   Float_t         zVtx;
   Float_t         nPV;
   Int_t           Centrality;
   Int_t           nTrig;
   Int_t           trigPrescale[15];   //[nTrig]
   ULong64_t       HLTriggers;
   Int_t           Npix;
   Int_t           NpixelTracks;
   Int_t           Ntracks;
   Float_t         SumET_HF;
   Float_t         SumET_HFplus;
   Float_t         SumET_HFminus;
   Float_t         SumET_HFplusEta4;
   Float_t         SumET_HFminusEta4;
   Float_t         SumET_ET;
   Float_t         SumET_EE;
   Float_t         SumET_EB;
   Float_t         SumET_EEplus;
   Float_t         SumET_EEminus;
   Float_t         SumET_ZDC;
   Float_t         SumET_ZDCplus;
   Float_t         SumET_ZDCminus;
   Int_t           Reco_QQ_size;
   Int_t           Reco_QQ_type[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_sign[15];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_4mom;
   TClonesArray    *Reco_QQ_mupl_4mom;
   TClonesArray    *Reco_QQ_mumi_4mom;
   ULong64_t       Reco_QQ_trig[15];   //[Reco_QQ_size]
   ULong64_t       Reco_QQ_mupl_trig[15];   //[Reco_QQ_size]
   ULong64_t       Reco_QQ_mumi_trig[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_isCowboy[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctau[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauErr[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctau3D[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauErr3D[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_VtxProb[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_dca[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_MassErr[15];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_vtx;
   Int_t           Reco_QQ_Ntrk[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_SelectionType[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_SelectionType[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mupl_isGoodMuon[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mumi_isGoodMuon[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mupl_highPurity[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mumi_highPurity[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mupl_TrkMuArb[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mumi_TrkMuArb[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mupl_TMOneStaTight[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mumi_TMOneStaTight[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nPixValHits[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nPixValHits[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nMuValHits[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nMuValHits[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nTrkHits[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nTrkHits[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_normChi2_inner[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_normChi2_inner[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_normChi2_global[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_normChi2_global[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nPixWMea[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nPixWMea[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nTrkWMea[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nTrkWMea[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_StationsMatched[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_StationsMatched[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dxy[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dxy[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dxyErr[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dxyErr[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dz[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dz[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dzErr[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dzErr[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_pt_inner[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_pt_inner[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_pt_global[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_pt_global[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_ptErr_inner[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_ptErr_inner[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_ptErr_global[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_ptErr_global[15];   //[Reco_QQ_size]
   Int_t           Reco_mu_size;
   Int_t           Reco_mu_type[6];   //[Reco_mu_size]
   Int_t           Reco_mu_SelectionType[6];   //[Reco_mu_size]
   Int_t           Reco_mu_charge[6];   //[Reco_mu_size]
   TClonesArray    *Reco_mu_4mom;
   ULong64_t       Reco_mu_trig[6];   //[Reco_mu_size]
   Bool_t          Reco_mu_isGoodMuon[6];   //[Reco_mu_size]
   Bool_t          Reco_mu_highPurity[6];   //[Reco_mu_size]
   Bool_t          Reco_mu_TrkMuArb[6];   //[Reco_mu_size]
   Bool_t          Reco_mu_TMOneStaTight[6];   //[Reco_mu_size]
   Int_t           Reco_mu_nPixValHits[6];   //[Reco_mu_size]
   Int_t           Reco_mu_nMuValHits[6];   //[Reco_mu_size]
   Int_t           Reco_mu_nTrkHits[6];   //[Reco_mu_size]
   Float_t         Reco_mu_normChi2_inner[6];   //[Reco_mu_size]
   Float_t         Reco_mu_normChi2_global[6];   //[Reco_mu_size]
   Int_t           Reco_mu_nPixWMea[6];   //[Reco_mu_size]
   Int_t           Reco_mu_nTrkWMea[6];   //[Reco_mu_size]
   Int_t           Reco_mu_StationsMatched[6];   //[Reco_mu_size]
   Float_t         Reco_mu_dxy[6];   //[Reco_mu_size]
   Float_t         Reco_mu_dxyErr[6];   //[Reco_mu_size]
   Float_t         Reco_mu_dz[6];   //[Reco_mu_size]
   Float_t         Reco_mu_dzErr[6];   //[Reco_mu_size]
   Float_t         Reco_mu_pt_inner[6];   //[Reco_mu_size]
   Float_t         Reco_mu_pt_global[6];   //[Reco_mu_size]
   Float_t         Reco_mu_ptErr_inner[6];   //[Reco_mu_size]
   Float_t         Reco_mu_ptErr_global[6];   //[Reco_mu_size]

   // Declaration of leaf types from HltTree.h
   Int_t           NL1IsolEm;
   Float_t         L1IsolEmEt[4];   //[NL1IsolEm]
   Float_t         L1IsolEmE[4];   //[NL1IsolEm]
   Float_t         L1IsolEmEta[4];   //[NL1IsolEm]
   Float_t         L1IsolEmPhi[4];   //[NL1IsolEm]
   Int_t           NL1NIsolEm;
   Float_t         L1NIsolEmEt[4];   //[NL1NIsolEm]
   Float_t         L1NIsolEmE[4];   //[NL1NIsolEm]
   Float_t         L1NIsolEmEta[4];   //[NL1NIsolEm]
   Float_t         L1NIsolEmPhi[4];   //[NL1NIsolEm]
   Int_t           NL1Mu;
   Float_t         L1MuPt[4];   //[NL1Mu]
   Float_t         L1MuE[4];   //[NL1Mu]
   Float_t         L1MuEta[4];   //[NL1Mu]
   Float_t         L1MuPhi[4];   //[NL1Mu]
   Int_t           L1MuIsol[4];   //[NL1Mu]
   Int_t           L1MuMip[4];   //[NL1Mu]
   Int_t           L1MuFor[4];   //[NL1Mu]
   Int_t           L1MuRPC[4];   //[NL1Mu]
   Int_t           L1MuQal[4];   //[NL1Mu]
   Int_t           L1MuChg[4];   //[NL1Mu]
   Int_t           NL1CenJet;
   Float_t         L1CenJetEt[4];   //[NL1CenJet]
   Float_t         L1CenJetE[4];   //[NL1CenJet]
   Float_t         L1CenJetEta[4];   //[NL1CenJet]
   Float_t         L1CenJetPhi[4];   //[NL1CenJet]
   Int_t           NL1ForJet;
   Float_t         L1ForJetEt[4];   //[NL1ForJet]
   Float_t         L1ForJetE[4];   //[NL1ForJet]
   Float_t         L1ForJetEta[4];   //[NL1ForJet]
   Float_t         L1ForJetPhi[4];   //[NL1ForJet]
   Int_t           NL1Tau;
   Float_t         L1TauEt[4];   //[NL1Tau]
   Float_t         L1TauE[4];   //[NL1Tau]
   Float_t         L1TauEta[4];   //[NL1Tau]
   Float_t         L1TauPhi[4];   //[NL1Tau]
   Float_t         L1Met;
   Float_t         L1MetPhi;
   Float_t         L1EtTot;
   Float_t         L1Mht;
   Float_t         L1MhtPhi;
   Float_t         L1EtHad;
   Int_t           L1HfRing1EtSumPositiveEta;
   Int_t           L1HfRing2EtSumPositiveEta;
   Int_t           L1HfRing1EtSumNegativeEta;
   Int_t           L1HfRing2EtSumNegativeEta;
   Int_t           L1HfTowerCountPositiveEtaRing1;
   Int_t           L1HfTowerCountNegativeEtaRing1;
   Int_t           L1HfTowerCountPositiveEtaRing2;
   Int_t           L1HfTowerCountNegativeEtaRing2;
   ULong64_t       Event;
   Int_t           LumiBlock;
   Int_t           Run;
   Int_t           Bx;
   Int_t           Orbit;
   Double_t        AvgInstDelLumi;
   Int_t           HLTriggerFirstPath;
   Int_t           HLTriggerFirstPath_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity60_v2;
   Int_t           HLT_PixelTracks_Multiplicity60_v2_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity85_v2;
   Int_t           HLT_PixelTracks_Multiplicity85_v2_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity110_v2;
   Int_t           HLT_PixelTracks_Multiplicity110_v2_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity135_v2;
   Int_t           HLT_PixelTracks_Multiplicity135_v2_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity160_v2;
   Int_t           HLT_PixelTracks_Multiplicity160_v2_Prescl;
   Int_t           HLT_Physics_v2;
   Int_t           HLT_Physics_v2_Prescl;
   Int_t           DST_Physics_v1;
   Int_t           DST_Physics_v1_Prescl;
   Int_t           HLT_Random_v1;
   Int_t           HLT_Random_v1_Prescl;
   Int_t           HLT_EcalCalibration_v1;
   Int_t           HLT_EcalCalibration_v1_Prescl;
   Int_t           HLT_HcalCalibration_v1;
   Int_t           HLT_HcalCalibration_v1_Prescl;
   Int_t           AlCa_EcalPhiSym_v3;
   Int_t           AlCa_EcalPhiSym_v3_Prescl;
   Int_t           HLT_L1Tech6_BPTX_MinusOnly_v1;
   Int_t           HLT_L1Tech6_BPTX_MinusOnly_v1_Prescl;
   Int_t           HLT_L1Tech5_BPTX_PlusOnly_v2;
   Int_t           HLT_L1Tech5_BPTX_PlusOnly_v2_Prescl;
   Int_t           HLT_L1Tech7_NoBPTX_v1;
   Int_t           HLT_L1Tech7_NoBPTX_v1_Prescl;
   Int_t           HLT_L1TOTEM1_MinBias_v1;
   Int_t           HLT_L1TOTEM1_MinBias_v1_Prescl;
   Int_t           HLT_L1TOTEM2_ZeroBias_v1;
   Int_t           HLT_L1TOTEM2_ZeroBias_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF2OR_v1;
   Int_t           HLT_L1MinimumBiasHF2OR_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1AND_v1;
   Int_t           HLT_L1MinimumBiasHF1AND_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF2AND_v1;
   Int_t           HLT_L1MinimumBiasHF2AND_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF2ORNoBptxGating_v1;
   Int_t           HLT_L1MinimumBiasHF2ORNoBptxGating_v1_Prescl;
   Int_t           AlCa_RPCMuonNoTriggers_v2;
   Int_t           AlCa_RPCMuonNoTriggers_v2_Prescl;
   Int_t           AlCa_RPCMuonNoHits_v2;
   Int_t           AlCa_RPCMuonNoHits_v2_Prescl;
   Int_t           AlCa_RPCMuonNormalisation_v2;
   Int_t           AlCa_RPCMuonNormalisation_v2_Prescl;
   Int_t           AlCa_LumiPixels_Random_v1;
   Int_t           AlCa_LumiPixels_Random_v1_Prescl;
   Int_t           AlCa_LumiPixels_ZeroBias_v2;
   Int_t           AlCa_LumiPixels_ZeroBias_v2_Prescl;
   Int_t           HLT_ZeroBias_v2;
   Int_t           HLT_ZeroBias_v2_Prescl;
   Int_t           HLT_ZeroBias_part0_v2;
   Int_t           HLT_ZeroBias_part0_v2_Prescl;
   Int_t           HLT_ZeroBias_part1_v2;
   Int_t           HLT_ZeroBias_part1_v2_Prescl;
   Int_t           HLT_ZeroBias_part2_v2;
   Int_t           HLT_ZeroBias_part2_v2_Prescl;
   Int_t           HLT_ZeroBias_part3_v2;
   Int_t           HLT_ZeroBias_part3_v2_Prescl;
   Int_t           HLT_ZeroBias_part4_v2;
   Int_t           HLT_ZeroBias_part4_v2_Prescl;
   Int_t           HLT_ZeroBias_part5_v2;
   Int_t           HLT_ZeroBias_part5_v2_Prescl;
   Int_t           HLT_ZeroBias_part6_v2;
   Int_t           HLT_ZeroBias_part6_v2_Prescl;
   Int_t           HLT_ZeroBias_part7_v2;
   Int_t           HLT_ZeroBias_part7_v2_Prescl;
   Int_t           HLT_ZeroBias_part8_v2;
   Int_t           HLT_ZeroBias_part8_v2_Prescl;
   Int_t           HLT_ZeroBias_part9_v2;
   Int_t           HLT_ZeroBias_part9_v2_Prescl;
   Int_t           HLT_ZeroBias_part10_v2;
   Int_t           HLT_ZeroBias_part10_v2_Prescl;
   Int_t           HLT_ZeroBias_part11_v2;
   Int_t           HLT_ZeroBias_part11_v2_Prescl;
   Int_t           HLT_ZeroBias_part12_v2;
   Int_t           HLT_ZeroBias_part12_v2_Prescl;
   Int_t           HLT_ZeroBias_part13_v2;
   Int_t           HLT_ZeroBias_part13_v2_Prescl;
   Int_t           HLT_ZeroBias_part14_v2;
   Int_t           HLT_ZeroBias_part14_v2_Prescl;
   Int_t           HLT_ZeroBias_part15_v2;
   Int_t           HLT_ZeroBias_part15_v2_Prescl;
   Int_t           HLT_ZeroBias_part16_v2;
   Int_t           HLT_ZeroBias_part16_v2_Prescl;
   Int_t           HLT_ZeroBias_part17_v2;
   Int_t           HLT_ZeroBias_part17_v2_Prescl;
   Int_t           HLT_ZeroBias_part18_v2;
   Int_t           HLT_ZeroBias_part18_v2_Prescl;
   Int_t           HLT_ZeroBias_part19_v2;
   Int_t           HLT_ZeroBias_part19_v2_Prescl;
   Int_t           HLT_AK4CaloJet40_Eta5p1_v1;
   Int_t           HLT_AK4CaloJet40_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet60_Eta5p1_v1;
   Int_t           HLT_AK4CaloJet60_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet80_Eta5p1_v1;
   Int_t           HLT_AK4CaloJet80_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet100_Eta5p1_v1;
   Int_t           HLT_AK4CaloJet100_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet110_Eta5p1_v1;
   Int_t           HLT_AK4CaloJet110_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet120_Eta5p1_v1;
   Int_t           HLT_AK4CaloJet120_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet150_v1;
   Int_t           HLT_AK4CaloJet150_v1_Prescl;
   Int_t           HLT_AK4PFJet40_Eta5p1_v1;
   Int_t           HLT_AK4PFJet40_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4PFJet60_Eta5p1_v1;
   Int_t           HLT_AK4PFJet60_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4PFJet80_Eta5p1_v1;
   Int_t           HLT_AK4PFJet80_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4PFJet100_Eta5p1_v1;
   Int_t           HLT_AK4PFJet100_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4PFJet110_Eta5p1_v1;
   Int_t           HLT_AK4PFJet110_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4PFJet120_Eta5p1_v1;
   Int_t           HLT_AK4PFJet120_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet80_Jet35_Eta1p1_v1;
   Int_t           HLT_AK4CaloJet80_Jet35_Eta1p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet80_Jet35_Eta0p7_v1;
   Int_t           HLT_AK4CaloJet80_Jet35_Eta0p7_v1_Prescl;
   Int_t           HLT_AK4CaloJet100_Jet35_Eta1p1_v1;
   Int_t           HLT_AK4CaloJet100_Jet35_Eta1p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet100_Jet35_Eta0p7_v1;
   Int_t           HLT_AK4CaloJet100_Jet35_Eta0p7_v1_Prescl;
   Int_t           HLT_AK4CaloJet80_45_45_Eta2p1_v1;
   Int_t           HLT_AK4CaloJet80_45_45_Eta2p1_v1_Prescl;
   Int_t           HLT_HISinglePhoton10_Eta1p5_v1;
   Int_t           HLT_HISinglePhoton10_Eta1p5_v1_Prescl;
   Int_t           HLT_HISinglePhoton15_Eta1p5_v1;
   Int_t           HLT_HISinglePhoton15_Eta1p5_v1_Prescl;
   Int_t           HLT_HISinglePhoton20_Eta1p5_v1;
   Int_t           HLT_HISinglePhoton20_Eta1p5_v1_Prescl;
   Int_t           HLT_HISinglePhoton30_Eta1p5_v1;
   Int_t           HLT_HISinglePhoton30_Eta1p5_v1_Prescl;
   Int_t           HLT_HISinglePhoton40_Eta1p5_v1;
   Int_t           HLT_HISinglePhoton40_Eta1p5_v1_Prescl;
   Int_t           HLT_HISinglePhoton50_Eta1p5_v1;
   Int_t           HLT_HISinglePhoton50_Eta1p5_v1_Prescl;
   Int_t           HLT_HISinglePhoton60_Eta1p5_v1;
   Int_t           HLT_HISinglePhoton60_Eta1p5_v1_Prescl;
   Int_t           HLT_HISinglePhoton10_Eta3p1_v1;
   Int_t           HLT_HISinglePhoton10_Eta3p1_v1_Prescl;
   Int_t           HLT_HISinglePhoton15_Eta3p1_v1;
   Int_t           HLT_HISinglePhoton15_Eta3p1_v1_Prescl;
   Int_t           HLT_HISinglePhoton20_Eta3p1_v1;
   Int_t           HLT_HISinglePhoton20_Eta3p1_v1_Prescl;
   Int_t           HLT_HISinglePhoton30_Eta3p1_v1;
   Int_t           HLT_HISinglePhoton30_Eta3p1_v1_Prescl;
   Int_t           HLT_HISinglePhoton40_Eta3p1_v1;
   Int_t           HLT_HISinglePhoton40_Eta3p1_v1_Prescl;
   Int_t           HLT_HISinglePhoton50_Eta3p1_v1;
   Int_t           HLT_HISinglePhoton50_Eta3p1_v1_Prescl;
   Int_t           HLT_HISinglePhoton60_Eta3p1_v1;
   Int_t           HLT_HISinglePhoton60_Eta3p1_v1_Prescl;
   Int_t           HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_v1;
   Int_t           HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_v1_Prescl;
   Int_t           HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_R9HECut_v1;
   Int_t           HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_R9HECut_v1_Prescl;
   Int_t           HLT_HIDoublePhoton15_Eta2p1_Mass50_1000_R9Cut_v1;
   Int_t           HLT_HIDoublePhoton15_Eta2p1_Mass50_1000_R9Cut_v1_Prescl;
   Int_t           HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1;
   Int_t           HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1_Prescl;
   Int_t           HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1_Prescl;
   Int_t           HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1_Prescl;
   Int_t           HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1_Prescl;
   Int_t           HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1_Prescl;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1_Prescl;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1_Prescl;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1_Prescl;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1_Prescl;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1_Prescl;
   Int_t           HLT_HIL1DoubleMu0_v1;
   Int_t           HLT_HIL1DoubleMu0_v1_Prescl;
   Int_t           HLT_HIL1DoubleMu10_v1;
   Int_t           HLT_HIL1DoubleMu10_v1_Prescl;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v1;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v1_Prescl;
   Int_t           HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1;
   Int_t           HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1_Prescl;
   Int_t           HLT_HIL3DoubleMu0_OS_m7to14_v1;
   Int_t           HLT_HIL3DoubleMu0_OS_m7to14_v1_Prescl;
   Int_t           HLT_HIL2Mu3_NHitQ10_v1;
   Int_t           HLT_HIL2Mu3_NHitQ10_v1_Prescl;
   Int_t           HLT_HIL3Mu3_NHitQ15_v1;
   Int_t           HLT_HIL3Mu3_NHitQ15_v1_Prescl;
   Int_t           HLT_HIL2Mu5_NHitQ10_v1;
   Int_t           HLT_HIL2Mu5_NHitQ10_v1_Prescl;
   Int_t           HLT_HIL3Mu5_NHitQ15_v1;
   Int_t           HLT_HIL3Mu5_NHitQ15_v1_Prescl;
   Int_t           HLT_HIL2Mu7_NHitQ10_v1;
   Int_t           HLT_HIL2Mu7_NHitQ10_v1_Prescl;
   Int_t           HLT_HIL3Mu7_NHitQ15_v1;
   Int_t           HLT_HIL3Mu7_NHitQ15_v1_Prescl;
   Int_t           HLT_HIL2Mu15_v1;
   Int_t           HLT_HIL2Mu15_v1_Prescl;
   Int_t           HLT_HIL3Mu15_v1;
   Int_t           HLT_HIL3Mu15_v1_Prescl;
   Int_t           HLT_HIL2Mu20_v1;
   Int_t           HLT_HIL2Mu20_v1_Prescl;
   Int_t           HLT_HIL3Mu20_v1;
   Int_t           HLT_HIL3Mu20_v1_Prescl;
   Int_t           HLT_FullTrack18ForPPRef_v1;
   Int_t           HLT_FullTrack18ForPPRef_v1_Prescl;
   Int_t           HLT_FullTrack24ForPPRef_v1;
   Int_t           HLT_FullTrack24ForPPRef_v1_Prescl;
   Int_t           HLT_FullTrack34ForPPRef_v1;
   Int_t           HLT_FullTrack34ForPPRef_v1_Prescl;
   Int_t           HLT_FullTrack45ForPPRef_v1;
   Int_t           HLT_FullTrack45ForPPRef_v1_Prescl;
   Int_t           HLT_HIUPCL1DoubleMuOpenNotHF2_v1;
   Int_t           HLT_HIUPCL1DoubleMuOpenNotHF2_v1_Prescl;
   Int_t           HLT_HIUPCDoubleMuNotHF2Pixel_SingleTrack_v1;
   Int_t           HLT_HIUPCDoubleMuNotHF2Pixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIUPCL1MuOpen_NotMinimumBiasHF2_AND_v1;
   Int_t           HLT_HIUPCL1MuOpen_NotMinimumBiasHF2_AND_v1_Prescl;
   Int_t           HLT_HIUPCMuOpen_NotMinimumBiasHF2_ANDPixel_SingleTrack_v1;
   Int_t           HLT_HIUPCMuOpen_NotMinimumBiasHF2_ANDPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIUPCL1NotMinimumBiasHF2_AND_v1;
   Int_t           HLT_HIUPCL1NotMinimumBiasHF2_AND_v1_Prescl;
   Int_t           HLT_HIUPCNotMinimumBiasHF2_ANDPixel_SingleTrack_v1;
   Int_t           HLT_HIUPCNotMinimumBiasHF2_ANDPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIUPCL1ZdcOR_BptxAND_v1;
   Int_t           HLT_HIUPCL1ZdcOR_BptxAND_v1_Prescl;
   Int_t           HLT_HIUPCZdcOR_BptxANDPixel_SingleTrack_v1;
   Int_t           HLT_HIUPCZdcOR_BptxANDPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIUPCL1ZdcXOR_BptxAND_v1;
   Int_t           HLT_HIUPCL1ZdcXOR_BptxAND_v1_Prescl;
   Int_t           HLT_HIUPCZdcXOR_BptxANDPixel_SingleTrack_v1;
   Int_t           HLT_HIUPCZdcXOR_BptxANDPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIUPCL1NotZdcOR_BptxAND_v1;
   Int_t           HLT_HIUPCL1NotZdcOR_BptxAND_v1_Prescl;
   Int_t           HLT_HIUPCNotZdcOR_BptxANDPixel_SingleTrack_v1;
   Int_t           HLT_HIUPCNotZdcOR_BptxANDPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIL1CastorMediumJet_v1;
   Int_t           HLT_HIL1CastorMediumJet_v1_Prescl;
   Int_t           HLT_HICastorMediumJetPixel_SingleTrack_v1;
   Int_t           HLT_HICastorMediumJetPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt15_v1;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt15_v1_Prescl;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt20_v1;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt20_v1_Prescl;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt30_v1;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt30_v1_Prescl;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt40_v1;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt40_v1_Prescl;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt50_v1;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt50_v1_Prescl;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt60_v1;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt60_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part0_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part0_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part1_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part1_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part2_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part2_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part3_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part3_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part4_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part4_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part5_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part5_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part6_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part6_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part7_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part7_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part8_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part8_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part9_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part9_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part10_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part10_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part11_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part11_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part12_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part12_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part13_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part13_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part14_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part14_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part15_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part15_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part16_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part16_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part17_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part17_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part18_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part18_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_part19_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_part19_v1_Prescl;
   Int_t           HLT_AK4PFBJetBCSV60_Eta2p1_v1;
   Int_t           HLT_AK4PFBJetBCSV60_Eta2p1_v1_Prescl;
   Int_t           HLT_AK4PFBJetBCSV80_Eta2p1_v1;
   Int_t           HLT_AK4PFBJetBCSV80_Eta2p1_v1_Prescl;
   Int_t           HLT_AK4PFDJet60_Eta2p1_v1;
   Int_t           HLT_AK4PFDJet60_Eta2p1_v1_Prescl;
   Int_t           HLT_AK4PFDJet80_Eta2p1_v1;
   Int_t           HLT_AK4PFDJet80_Eta2p1_v1_Prescl;
   Int_t           HLT_AK4PFBJetBSSV60_Eta2p1_v1;
   Int_t           HLT_AK4PFBJetBSSV60_Eta2p1_v1_Prescl;
   Int_t           HLT_AK4PFBJetBSSV80_Eta2p1_v1;
   Int_t           HLT_AK4PFBJetBSSV80_Eta2p1_v1_Prescl;
   Int_t           HLTriggerFinalPath;
   Int_t           HLTriggerFinalPath_Prescl;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt8_v1;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt8_v1_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity60_v1;
   Int_t           HLT_PixelTracks_Multiplicity60_v1_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity60_v3;
   Int_t           HLT_PixelTracks_Multiplicity60_v3_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity85_v1;
   Int_t           HLT_PixelTracks_Multiplicity85_v1_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity110_v1;
   Int_t           HLT_PixelTracks_Multiplicity110_v1_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity135_v1;
   Int_t           HLT_PixelTracks_Multiplicity135_v1_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity160_v1;
   Int_t           HLT_PixelTracks_Multiplicity160_v1_Prescl;
   Int_t           HLT_Physics_v1;
   Int_t           HLT_Physics_v1_Prescl;
   Int_t           HLT_L1Tech5_BPTX_PlusOnly_v1;
   Int_t           HLT_L1Tech5_BPTX_PlusOnly_v1_Prescl;
   Int_t           HLT_ZeroBias_v1;
   Int_t           HLT_ZeroBias_v1_Prescl;
   Int_t           HLT_ZeroBias_part0_v1;
   Int_t           HLT_ZeroBias_part0_v1_Prescl;
   Int_t           HLT_ZeroBias_part1_v1;
   Int_t           HLT_ZeroBias_part1_v1_Prescl;
   Int_t           HLT_ZeroBias_part2_v1;
   Int_t           HLT_ZeroBias_part2_v1_Prescl;
   Int_t           HLT_ZeroBias_part3_v1;
   Int_t           HLT_ZeroBias_part3_v1_Prescl;
   Int_t           HLT_ZeroBias_part4_v1;
   Int_t           HLT_ZeroBias_part4_v1_Prescl;
   Int_t           HLT_ZeroBias_part5_v1;
   Int_t           HLT_ZeroBias_part5_v1_Prescl;
   Int_t           HLT_ZeroBias_part6_v1;
   Int_t           HLT_ZeroBias_part6_v1_Prescl;
   Int_t           HLT_ZeroBias_part7_v1;
   Int_t           HLT_ZeroBias_part7_v1_Prescl;
   Int_t           HLT_ZeroBias_part8_v1;
   Int_t           HLT_ZeroBias_part8_v1_Prescl;
   Int_t           HLT_ZeroBias_part9_v1;
   Int_t           HLT_ZeroBias_part9_v1_Prescl;
   Int_t           HLT_ZeroBias_part10_v1;
   Int_t           HLT_ZeroBias_part10_v1_Prescl;
   Int_t           HLT_ZeroBias_part11_v1;
   Int_t           HLT_ZeroBias_part11_v1_Prescl;
   Int_t           HLT_ZeroBias_part12_v1;
   Int_t           HLT_ZeroBias_part12_v1_Prescl;
   Int_t           HLT_ZeroBias_part13_v1;
   Int_t           HLT_ZeroBias_part13_v1_Prescl;
   Int_t           HLT_ZeroBias_part14_v1;
   Int_t           HLT_ZeroBias_part14_v1_Prescl;
   Int_t           HLT_ZeroBias_part15_v1;
   Int_t           HLT_ZeroBias_part15_v1_Prescl;
   Int_t           HLT_ZeroBias_part16_v1;
   Int_t           HLT_ZeroBias_part16_v1_Prescl;
   Int_t           HLT_ZeroBias_part17_v1;
   Int_t           HLT_ZeroBias_part17_v1_Prescl;
   Int_t           HLT_ZeroBias_part18_v1;
   Int_t           HLT_ZeroBias_part18_v1_Prescl;
   Int_t           HLT_ZeroBias_part19_v1;
   Int_t           HLT_ZeroBias_part19_v1_Prescl;
   Int_t           HLT_FullTrack18ForPPRef_v2;
   Int_t           HLT_FullTrack18ForPPRef_v2_Prescl;
   Int_t           HLT_FullTrack18ForPPRef_v3;
   Int_t           HLT_FullTrack18ForPPRef_v3_Prescl;
   Int_t           HLT_FullTrack24ForPPRef_v2;
   Int_t           HLT_FullTrack24ForPPRef_v2_Prescl;
   Int_t           HLT_FullTrack24ForPPRef_v3;
   Int_t           HLT_FullTrack24ForPPRef_v3_Prescl;
   Int_t           HLT_FullTrack34ForPPRef_v2;
   Int_t           HLT_FullTrack34ForPPRef_v2_Prescl;
   Int_t           HLT_FullTrack34ForPPRef_v3;
   Int_t           HLT_FullTrack34ForPPRef_v3_Prescl;
   Int_t           HLT_FullTrack34ForPPRef_v4;
   Int_t           HLT_FullTrack34ForPPRef_v4_Prescl;
   Int_t           HLT_FullTrack45ForPPRef_v2;
   Int_t           HLT_FullTrack45ForPPRef_v2_Prescl;
   Int_t           HLT_FullTrack45ForPPRef_v3;
   Int_t           HLT_FullTrack45ForPPRef_v3_Prescl;
   Int_t           HLT_FullTrack53ForPPRef_v1;
   Int_t           HLT_FullTrack53ForPPRef_v1_Prescl;
   Int_t           HLT_FullTrack53ForPPRef_v2;
   Int_t           HLT_FullTrack53ForPPRef_v2_Prescl;
   Int_t           HLT_FullTrack53ForPPRef_v3;
   Int_t           HLT_FullTrack53ForPPRef_v3_Prescl;
   Int_t           L1_AlwaysTrue;
   Int_t           L1_AlwaysTrue_Prescl;
   Int_t           L1_CastorHighJet_BptxAND;
   Int_t           L1_CastorHighJet_BptxAND_Prescl;
   Int_t           L1_CastorMediumJet_BptxAND;
   Int_t           L1_CastorMediumJet_BptxAND_Prescl;
   Int_t           L1_DoubleMu0_BptxAND;
   Int_t           L1_DoubleMu0_BptxAND_Prescl;
   Int_t           L1_DoubleMu0_MinimumBiasHF1_AND;
   Int_t           L1_DoubleMu0_MinimumBiasHF1_AND_Prescl;
   Int_t           L1_DoubleMu3_BptxAND;
   Int_t           L1_DoubleMu3_BptxAND_Prescl;
   Int_t           L1_DoubleMuOpen_BptxAND;
   Int_t           L1_DoubleMuOpen_BptxAND_Prescl;
   Int_t           L1_DoubleMuOpen_NotMinimumBiasHF2_AND;
   Int_t           L1_DoubleMuOpen_NotMinimumBiasHF2_AND_Prescl;
   Int_t           L1_DoubleJet20;
   Int_t           L1_DoubleJet20_Prescl;
   Int_t           L1_DoubleJet28_BptxAND;
   Int_t           L1_DoubleJet28_BptxAND_Prescl;
   Int_t           L1_DoubleJet32;
   Int_t           L1_DoubleJet32_Prescl;
   Int_t           L1_ETT130;
   Int_t           L1_ETT130_Prescl;
   Int_t           L1_ETT15_BptxAND;
   Int_t           L1_ETT15_BptxAND_Prescl;
   Int_t           L1_ETT20_BptxAND;
   Int_t           L1_ETT20_BptxAND_Prescl;
   Int_t           L1_ETT30_BptxAND;
   Int_t           L1_ETT30_BptxAND_Prescl;
   Int_t           L1_ETT40;
   Int_t           L1_ETT40_Prescl;
   Int_t           L1_ETT40_BptxAND;
   Int_t           L1_ETT40_BptxAND_Prescl;
   Int_t           L1_ETT45_BptxAND;
   Int_t           L1_ETT45_BptxAND_Prescl;
   Int_t           L1_ETT50_BptxAND;
   Int_t           L1_ETT50_BptxAND_Prescl;
   Int_t           L1_ETT55_BptxAND;
   Int_t           L1_ETT55_BptxAND_Prescl;
   Int_t           L1_ETT60_BptxAND;
   Int_t           L1_ETT60_BptxAND_Prescl;
   Int_t           L1_ETT70_BptxAND;
   Int_t           L1_ETT70_BptxAND_Prescl;
   Int_t           L1_ETT90_BptxAND;
   Int_t           L1_ETT90_BptxAND_Prescl;
   Int_t           L1_MinimumBiasHF1_AND;
   Int_t           L1_MinimumBiasHF1_AND_Prescl;
   Int_t           L1_MinimumBiasHF1_OR;
   Int_t           L1_MinimumBiasHF1_OR_Prescl;
   Int_t           L1_MinimumBiasHF1_XOR;
   Int_t           L1_MinimumBiasHF1_XOR_Prescl;
   Int_t           L1_MinimumBiasHF2_AND;
   Int_t           L1_MinimumBiasHF2_AND_Prescl;
   Int_t           L1_MinimumBiasHF2_OR;
   Int_t           L1_MinimumBiasHF2_OR_Prescl;
   Int_t           L1_MinimumBiasHF2_OR_NoBptxGating;
   Int_t           L1_MinimumBiasHF2_OR_NoBptxGating_Prescl;
   Int_t           L1_MuOpen_NotMinimumBiasHF2_AND;
   Int_t           L1_MuOpen_NotMinimumBiasHF2_AND_Prescl;
   Int_t           L1_NotMinimumBiasHF1_OR;
   Int_t           L1_NotMinimumBiasHF1_OR_Prescl;
   Int_t           L1_NotMinimumBiasHF2_AND;
   Int_t           L1_NotMinimumBiasHF2_AND_Prescl;
   Int_t           L1_NotZdcOR_BptxAND;
   Int_t           L1_NotZdcOR_BptxAND_Prescl;
   Int_t           L1_SingleEG12_BptxAND;
   Int_t           L1_SingleEG12_BptxAND_Prescl;
   Int_t           L1_SingleEG15_BptxAND;
   Int_t           L1_SingleEG15_BptxAND_Prescl;
   Int_t           L1_SingleEG20;
   Int_t           L1_SingleEG20_Prescl;
   Int_t           L1_SingleEG20_BptxAND;
   Int_t           L1_SingleEG20_BptxAND_Prescl;
   Int_t           L1_SingleEG2_BptxAND;
   Int_t           L1_SingleEG2_BptxAND_Prescl;
   Int_t           L1_SingleEG30_BptxAND;
   Int_t           L1_SingleEG30_BptxAND_Prescl;
   Int_t           L1_SingleEG5;
   Int_t           L1_SingleEG5_Prescl;
   Int_t           L1_SingleEG5_BptxAND;
   Int_t           L1_SingleEG5_BptxAND_Prescl;
   Int_t           L1_SingleEG7_BptxAND;
   Int_t           L1_SingleEG7_BptxAND_Prescl;
   Int_t           L1_SingleMu12_BptxAND;
   Int_t           L1_SingleMu12_BptxAND_Prescl;
   Int_t           L1_SingleMu16_BptxAND;
   Int_t           L1_SingleMu16_BptxAND_Prescl;
   Int_t           L1_SingleMu16_MinimumBiasHF1_AND;
   Int_t           L1_SingleMu16_MinimumBiasHF1_AND_Prescl;
   Int_t           L1_SingleMu3_BptxAND;
   Int_t           L1_SingleMu3_BptxAND_Prescl;
   Int_t           L1_SingleMu3_MinimumBiasHF1_AND;
   Int_t           L1_SingleMu3_MinimumBiasHF1_AND_Prescl;
   Int_t           L1_SingleMu3p5;
   Int_t           L1_SingleMu3p5_Prescl;
   Int_t           L1_SingleMu5_BptxAND;
   Int_t           L1_SingleMu5_BptxAND_Prescl;
   Int_t           L1_SingleMu7_BptxAND;
   Int_t           L1_SingleMu7_BptxAND_Prescl;
   Int_t           L1_SingleMuBeamHalo;
   Int_t           L1_SingleMuBeamHalo_Prescl;
   Int_t           L1_SingleMuOpen;
   Int_t           L1_SingleMuOpen_Prescl;
   Int_t           L1_SingleMuOpen_BptxAND;
   Int_t           L1_SingleMuOpen_BptxAND_Prescl;
   Int_t           L1_SingleMuOpen_NotBptxOR;
   Int_t           L1_SingleMuOpen_NotBptxOR_Prescl;
   Int_t           L1_SingleJet12_BptxAND;
   Int_t           L1_SingleJet12_BptxAND_Prescl;
   Int_t           L1_SingleJet16;
   Int_t           L1_SingleJet16_Prescl;
   Int_t           L1_SingleJet16_BptxAND;
   Int_t           L1_SingleJet16_BptxAND_Prescl;
   Int_t           L1_SingleJet200;
   Int_t           L1_SingleJet200_Prescl;
   Int_t           L1_SingleJet20_BptxAND;
   Int_t           L1_SingleJet20_BptxAND_Prescl;
   Int_t           L1_SingleJet24_BptxAND;
   Int_t           L1_SingleJet24_BptxAND_Prescl;
   Int_t           L1_SingleJet28_BptxAND;
   Int_t           L1_SingleJet28_BptxAND_Prescl;
   Int_t           L1_SingleJet32_BptxAND;
   Int_t           L1_SingleJet32_BptxAND_Prescl;
   Int_t           L1_SingleJet36;
   Int_t           L1_SingleJet36_Prescl;
   Int_t           L1_SingleJet36_BptxAND;
   Int_t           L1_SingleJet36_BptxAND_Prescl;
   Int_t           L1_SingleJet40_BptxAND;
   Int_t           L1_SingleJet40_BptxAND_Prescl;
   Int_t           L1_SingleJet44_BptxAND;
   Int_t           L1_SingleJet44_BptxAND_Prescl;
   Int_t           L1_SingleJet48_BptxAND;
   Int_t           L1_SingleJet48_BptxAND_Prescl;
   Int_t           L1_SingleJet52_BptxAND;
   Int_t           L1_SingleJet52_BptxAND_Prescl;
   Int_t           L1_SingleJet60_BptxAND;
   Int_t           L1_SingleJet60_BptxAND_Prescl;
   Int_t           L1_SingleJet68;
   Int_t           L1_SingleJet68_Prescl;
   Int_t           L1_SingleJet68_BptxAND;
   Int_t           L1_SingleJet68_BptxAND_Prescl;
   Int_t           L1_SingleJet8_BptxAND;
   Int_t           L1_SingleJet8_BptxAND_Prescl;
   Int_t           L1_SingleJetC20_NotBptxOR;
   Int_t           L1_SingleJetC20_NotBptxOR_Prescl;
   Int_t           L1_SingleJetC32_NotBptxOR;
   Int_t           L1_SingleJetC32_NotBptxOR_Prescl;
   Int_t           L1_TOTEM_0;
   Int_t           L1_TOTEM_0_Prescl;
   Int_t           L1_TOTEM_1;
   Int_t           L1_TOTEM_1_Prescl;
   Int_t           L1_TOTEM_2;
   Int_t           L1_TOTEM_2_Prescl;
   Int_t           L1_TOTEM_3;
   Int_t           L1_TOTEM_3_Prescl;
   Int_t           L1_ZdcOR;
   Int_t           L1_ZdcOR_Prescl;
   Int_t           L1_ZdcOR_BptxAND;
   Int_t           L1_ZdcOR_BptxAND_Prescl;
   Int_t           L1_ZdcXOR;
   Int_t           L1_ZdcXOR_Prescl;
   Int_t           L1_ZdcXOR_BptxAND;
   Int_t           L1_ZdcXOR_BptxAND_Prescl;
   Int_t           L1_ZeroBias;
   Int_t           L1_ZeroBias_Prescl;
   Int_t           L1Tech_BPTX_PreBPTX_v0;
   Int_t           L1Tech_BPTX_PreBPTX_v0_Prescl;
   Int_t           L1Tech_BPTX_minus_v0;
   Int_t           L1Tech_BPTX_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_minus_AND_not_plus_v0;
   Int_t           L1Tech_BPTX_minus_AND_not_plus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_v0;
   Int_t           L1Tech_BPTX_plus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_NOT_minus_v0;
   Int_t           L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_minus_v0;
   Int_t           L1Tech_BPTX_plus_AND_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_minus_instance1_v0;
   Int_t           L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_OR_minus_v0;
   Int_t           L1Tech_BPTX_plus_OR_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_quiet_v0;
   Int_t           L1Tech_BPTX_quiet_v0_Prescl;
   Int_t           L1Tech_BRIL_bit28;
   Int_t           L1Tech_BRIL_bit28_Prescl;
   Int_t           L1Tech_BRIL_bit29;
   Int_t           L1Tech_BRIL_bit29_Prescl;
   Int_t           L1Tech_BRIL_bit30;
   Int_t           L1Tech_BRIL_bit30_Prescl;
   Int_t           L1Tech_BRIL_bit31;
   Int_t           L1Tech_BRIL_bit31_Prescl;
   Int_t           L1Tech_BRIL_bit32;
   Int_t           L1Tech_BRIL_bit32_Prescl;
   Int_t           L1Tech_BRIL_bit33;
   Int_t           L1Tech_BRIL_bit33_Prescl;
   Int_t           L1Tech_BRIL_bit34;
   Int_t           L1Tech_BRIL_bit34_Prescl;
   Int_t           L1Tech_BRIL_bit35;
   Int_t           L1Tech_BRIL_bit35_Prescl;
   Int_t           L1Tech_BRIL_bit36;
   Int_t           L1Tech_BRIL_bit36_Prescl;
   Int_t           L1Tech_BRIL_bit37;
   Int_t           L1Tech_BRIL_bit37_Prescl;
   Int_t           L1Tech_BRIL_bit38;
   Int_t           L1Tech_BRIL_bit38_Prescl;
   Int_t           L1Tech_BRIL_bit39;
   Int_t           L1Tech_BRIL_bit39_Prescl;
   Int_t           L1Tech_BRIL_bit40;
   Int_t           L1Tech_BRIL_bit40_Prescl;
   Int_t           L1Tech_BRIL_bit41;
   Int_t           L1Tech_BRIL_bit41_Prescl;
   Int_t           L1Tech_BRIL_bit42;
   Int_t           L1Tech_BRIL_bit42_Prescl;
   Int_t           L1Tech_BRIL_bit43;
   Int_t           L1Tech_BRIL_bit43_Prescl;
   Int_t           L1Tech_CASTOR_Gap_v0;
   Int_t           L1Tech_CASTOR_Gap_v0_Prescl;
   Int_t           L1Tech_CASTOR_HaloMuon_v0;
   Int_t           L1Tech_CASTOR_HaloMuon_v0_Prescl;
   Int_t           L1Tech_CASTOR_HighJet_v0;
   Int_t           L1Tech_CASTOR_HighJet_v0_Prescl;
   Int_t           L1Tech_CASTOR_MediumJet_v0;
   Int_t           L1Tech_CASTOR_MediumJet_v0_Prescl;
   Int_t           L1Tech_DT_GlobalOR_v0;
   Int_t           L1Tech_DT_GlobalOR_v0_Prescl;
   Int_t           L1Tech_HCAL_HBHE_totalOR_v0;
   Int_t           L1Tech_HCAL_HBHE_totalOR_v0_Prescl;
   Int_t           L1Tech_HCAL_HF_MMP_or_MPP_v1;
   Int_t           L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl;
   Int_t           L1Tech_HCAL_HF_coincidence_PM_v2;
   Int_t           L1Tech_HCAL_HF_coincidence_PM_v2_Prescl;
   Int_t           L1Tech_HCAL_HF_single_channel_v0;
   Int_t           L1Tech_HCAL_HF_single_channel_v0_Prescl;
   Int_t           L1Tech_HCAL_HO_totalOR_v0;
   Int_t           L1Tech_HCAL_HO_totalOR_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBplus1_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBplus2_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_barrel_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_pointing_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl;
   Int_t           L1Tech_TOTEM_0;
   Int_t           L1Tech_TOTEM_0_Prescl;
   Int_t           L1Tech_TOTEM_1;
   Int_t           L1Tech_TOTEM_1_Prescl;
   Int_t           L1Tech_TOTEM_2;
   Int_t           L1Tech_TOTEM_2_Prescl;
   Int_t           L1Tech_TOTEM_3;
   Int_t           L1Tech_TOTEM_3_Prescl;
   Int_t           L1Tech_ZDC_minus;
   Int_t           L1Tech_ZDC_minus_Prescl;
   Int_t           L1Tech_ZDC_plus;
   Int_t           L1Tech_ZDC_plus_Prescl;

   // Declaration of leaf types from HltTree1.h
   Int_t           Onia2MuMuPAT;
   Int_t           ana_step;
   Int_t           pHBHENoiseFilterResultProducer;
   Int_t           HBHENoiseFilterResult;
   Int_t           HBHENoiseFilterResultRun1;
   Int_t           HBHENoiseFilterResultRun2Loose;
   Int_t           HBHENoiseFilterResultRun2Tight;
   Int_t           HBHEIsoNoiseFilterResult;
   Int_t           pPAprimaryVertexFilter;
   Int_t           pBeamScrapingFilter;
   Int_t           pVertexFilterCutG;
   Int_t           pVertexFilterCutGloose;
   Int_t           pVertexFilterCutGtight;
   Int_t           pVertexFilterCutGplus;
   Int_t           pVertexFilterCutE;
   Int_t           pVertexFilterCutEandG;

   // Declaration of leaf types from HiTree.h
   UInt_t          run;
   //ULong64_t       evt;
   UInt_t          lumi;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Int_t           hiBin;
   Float_t         hiHF;
   Float_t         hiHFplus;
   Float_t         hiHFminus;
   Float_t         hiHFplusEta4;
   Float_t         hiHFminusEta4;
   Float_t         hiZDC;
   Float_t         hiZDCplus;
   Float_t         hiZDCminus;
   Float_t         hiHFhit;
   Float_t         hiHFhitPlus;
   Float_t         hiHFhitMinus;
   Float_t         hiET;
   Float_t         hiEE;
   Float_t         hiEB;
   Float_t         hiEEplus;
   Float_t         hiEEminus;
   Int_t           hiNpix;
   Int_t           hiNpixelTracks;
   Int_t           hiNtracks;
   Int_t           hiNtracksPtCut;
   Int_t           hiNtracksEtaCut;
   Int_t           hiNtracksEtaPtCut;

   // Declaration of leaf types from t.h
   Int_t           evt;
   Float_t         b;
   Int_t           nref;
   Float_t         rawpt[56];   //[nref]
   Float_t         jtpt[56];   //[nref]
   Float_t         jteta[56];   //[nref]
   Float_t         jty[56];   //[nref]
   Float_t         jtphi[56];   //[nref]
   Float_t         jtpu[56];   //[nref]
   Float_t         jtm[56];   //[nref]
   Float_t         jtarea[56];   //[nref]
   Float_t         jtPfCHF[56];   //[nref]
   Float_t         jtPfNHF[56];   //[nref]
   Float_t         jtPfCEF[56];   //[nref]
   Float_t         jtPfNEF[56];   //[nref]
   Float_t         jtPfMUF[56];   //[nref]
   Int_t           jtPfCHM[56];   //[nref]
   Int_t           jtPfNHM[56];   //[nref]
   Int_t           jtPfCEM[56];   //[nref]
   Int_t           jtPfNEM[56];   //[nref]
   Int_t           jtPfMUM[56];   //[nref]
   Float_t         jttau1[56];   //[nref]
   Float_t         jttau2[56];   //[nref]
   Float_t         jttau3[56];   //[nref]
   Float_t         discr_jetID_cuts[56];   //[nref]
   Float_t         discr_jetID_bdt[56];   //[nref]
   Float_t         discr_fr01[56];   //[nref]
   Float_t         trackMax[56];   //[nref]
   Float_t         trackSum[56];   //[nref]
   Int_t           trackN[56];   //[nref]
   Float_t         trackHardSum[56];   //[nref]
   Int_t           trackHardN[56];   //[nref]
   Float_t         chargedMax[56];   //[nref]
   Float_t         chargedSum[56];   //[nref]
   Int_t           chargedN[56];   //[nref]
   Float_t         chargedHardSum[56];   //[nref]
   Int_t           chargedHardN[56];   //[nref]
   Float_t         photonMax[56];   //[nref]
   Float_t         photonSum[56];   //[nref]
   Int_t           photonN[56];   //[nref]
   Float_t         photonHardSum[56];   //[nref]
   Int_t           photonHardN[56];   //[nref]
   Float_t         neutralMax[56];   //[nref]
   Float_t         neutralSum[56];   //[nref]
   Int_t           neutralN[56];   //[nref]
   Float_t         hcalSum[56];   //[nref]
   Float_t         ecalSum[56];   //[nref]
   Float_t         eMax[56];   //[nref]
   Float_t         eSum[56];   //[nref]
   Int_t           eN[56];   //[nref]
   Float_t         muMax[56];   //[nref]
   Float_t         muSum[56];   //[nref]
   Int_t           muN[56];   //[nref]
   Float_t         discr_ssvHighEff[56];   //[nref]
   Float_t         discr_ssvHighPur[56];   //[nref]
   Float_t         discr_csvV1[56];   //[nref]
   Float_t         discr_csvV2[56];   //[nref]
   Float_t         discr_muByIp3[56];   //[nref]
   Float_t         discr_muByPt[56];   //[nref]
   Float_t         discr_prob[56];   //[nref]
   Float_t         discr_probb[56];   //[nref]
   Float_t         discr_tcHighEff[56];   //[nref]
   Float_t         discr_tcHighPur[56];   //[nref]
   Float_t         ndiscr_ssvHighEff[56];   //[nref]
   Float_t         ndiscr_ssvHighPur[56];   //[nref]
   Float_t         ndiscr_csvV1[56];   //[nref]
   Float_t         ndiscr_csvV2[56];   //[nref]
   Float_t         ndiscr_muByPt[56];   //[nref]
   Float_t         pdiscr_csvV1[56];   //[nref]
   Float_t         pdiscr_csvV2[56];   //[nref]
   Int_t           nsvtx[56];   //[nref]
   Int_t           svtxntrk[56];   //[nref]
   Float_t         svtxdl[56];   //[nref]
   Float_t         svtxdls[56];   //[nref]
   Float_t         svtxdl2d[56];   //[nref]
   Float_t         svtxdls2d[56];   //[nref]
   Float_t         svtxm[56];   //[nref]
   Float_t         svtxpt[56];   //[nref]
   Float_t         svtxmcorr[56];   //[nref]
   Int_t           nIPtrk[56];   //[nref]
   Int_t           nselIPtrk[56];   //[nref]
   Float_t         mue[56];   //[nref]
   Float_t         mupt[56];   //[nref]
   Float_t         mueta[56];   //[nref]
   Float_t         muphi[56];   //[nref]
   Float_t         mudr[56];   //[nref]
   Float_t         muptrel[56];   //[nref]
   Int_t           muchg[56];   //[nref]

   // List of branches
   TBranch        *b_eventNb;   //!
   TBranch        *b_runNb;   //!
   TBranch        *b_LS;   //!
   TBranch        *b_zVtx;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_Centrality;   //!
   TBranch        *b_nTrig;   //!
   TBranch        *b_trigPrescale;   //!
   TBranch        *b_HLTriggers;   //!
   TBranch        *b_Npix;   //!
   TBranch        *b_NpixelTracks;   //!
   TBranch        *b_Ntracks;   //!
   TBranch        *b_SumET_HF;   //!
   TBranch        *b_SumET_HFplus;   //!
   TBranch        *b_SumET_HFminus;   //!
   TBranch        *b_SumET_HFplusEta4;   //!
   TBranch        *b_SumET_HFminusEta4;   //!
   TBranch        *b_SumET_ET;   //!
   TBranch        *b_SumET_EE;   //!
   TBranch        *b_SumET_EB;   //!
   TBranch        *b_SumET_EEplus;   //!
   TBranch        *b_SumET_EEminus;   //!
   TBranch        *b_SumET_ZDC;   //!
   TBranch        *b_SumET_ZDCplus;   //!
   TBranch        *b_SumET_ZDCminus;   //!
   TBranch        *b_Reco_QQ_size;   //!
   TBranch        *b_Reco_QQ_type;   //!
   TBranch        *b_Reco_QQ_sign;   //!
   TBranch        *b_Reco_QQ_4mom;   //!
   TBranch        *b_Reco_QQ_mupl_4mom;   //!
   TBranch        *b_Reco_QQ_mumi_4mom;   //!
   TBranch        *b_Reco_QQ_trig;   //!
   TBranch        *b_Reco_QQ_mupl_trig;   //!
   TBranch        *b_Reco_QQ_mumi_trig;   //!
   TBranch        *b_Reco_QQ_isCowboy;   //!
   TBranch        *b_Reco_QQ_ctau;   //!
   TBranch        *b_Reco_QQ_ctauErr;   //!
   TBranch        *b_Reco_QQ_ctau3D;   //!
   TBranch        *b_Reco_QQ_ctauErr3D;   //!
   TBranch        *b_Reco_QQ_VtxProb;   //!
   TBranch        *b_Reco_QQ_dca;   //!
   TBranch        *b_Reco_QQ_MassErr;   //!
   TBranch        *b_Reco_QQ_vtx;   //!
   TBranch        *b_Reco_QQ_Ntrk;   //!
   TBranch        *b_Reco_QQ_mupl_SelectionType;   //!
   TBranch        *b_Reco_QQ_mumi_SelectionType;   //!
   TBranch        *b_Reco_QQ_mupl_isGoodMuon;   //!
   TBranch        *b_Reco_QQ_mumi_isGoodMuon;   //!
   TBranch        *b_Reco_QQ_mupl_highPurity;   //!
   TBranch        *b_Reco_QQ_mumi_highPurity;   //!
   TBranch        *b_Reco_QQ_mupl_TrkMuArb;   //!
   TBranch        *b_Reco_QQ_mumi_TrkMuArb;   //!
   TBranch        *b_Reco_QQ_mupl_TMOneStaTight;   //!
   TBranch        *b_Reco_QQ_mumi_TMOneStaTight;   //!
   TBranch        *b_Reco_QQ_mupl_nPixValHits;   //!
   TBranch        *b_Reco_QQ_mumi_nPixValHits;   //!
   TBranch        *b_Reco_QQ_mupl_nMuValHits;   //!
   TBranch        *b_Reco_QQ_mumi_nMuValHits;   //!
   TBranch        *b_Reco_QQ_mupl_nTrkHits;   //!
   TBranch        *b_Reco_QQ_mumi_nTrkHits;   //!
   TBranch        *b_Reco_QQ_mupl_normChi2_inner;   //!
   TBranch        *b_Reco_QQ_mumi_normChi2_inner;   //!
   TBranch        *b_Reco_QQ_mupl_normChi2_global;   //!
   TBranch        *b_Reco_QQ_mumi_normChi2_global;   //!
   TBranch        *b_Reco_QQ_mupl_nPixWMea;   //!
   TBranch        *b_Reco_QQ_mumi_nPixWMea;   //!
   TBranch        *b_Reco_QQ_mupl_nTrkWMea;   //!
   TBranch        *b_Reco_QQ_mumi_nTrkWMea;   //!
   TBranch        *b_Reco_QQ_mupl_StationsMatched;   //!
   TBranch        *b_Reco_QQ_mumi_StationsMatched;   //!
   TBranch        *b_Reco_QQ_mupl_dxy;   //!
   TBranch        *b_Reco_QQ_mumi_dxy;   //!
   TBranch        *b_Reco_QQ_mupl_dxyErr;   //!
   TBranch        *b_Reco_QQ_mumi_dxyErr;   //!
   TBranch        *b_Reco_QQ_mupl_dz;   //!
   TBranch        *b_Reco_QQ_mumi_dz;   //!
   TBranch        *b_Reco_QQ_mupl_dzErr;   //!
   TBranch        *b_Reco_QQ_mumi_dzErr;   //!
   TBranch        *b_Reco_QQ_mupl_pt_inner;   //!
   TBranch        *b_Reco_QQ_mumi_pt_inner;   //!
   TBranch        *b_Reco_QQ_mupl_pt_global;   //!
   TBranch        *b_Reco_QQ_mumi_pt_global;   //!
   TBranch        *b_Reco_QQ_mupl_ptErr_inner;   //!
   TBranch        *b_Reco_QQ_mumi_ptErr_inner;   //!
   TBranch        *b_Reco_QQ_mupl_ptErr_global;   //!
   TBranch        *b_Reco_QQ_mumi_ptErr_global;   //!
   TBranch        *b_Reco_mu_size;   //!
   TBranch        *b_Reco_mu_type;   //!
   TBranch        *b_Reco_mu_SelectionType;   //!
   TBranch        *b_Reco_mu_charge;   //!
   TBranch        *b_Reco_mu_4mom;   //!
   TBranch        *b_Reco_mu_trig;   //!
   TBranch        *b_Reco_mu_isGoodMuon;   //!
   TBranch        *b_Reco_mu_highPurity;   //!
   TBranch        *b_Reco_mu_TrkMuArb;   //!
   TBranch        *b_Reco_mu_TMOneStaTight;   //!
   TBranch        *b_Reco_mu_nPixValHits;   //!
   TBranch        *b_Reco_mu_nMuValHits;   //!
   TBranch        *b_Reco_mu_nTrkHits;   //!
   TBranch        *b_Reco_mu_normChi2_inner;   //!
   TBranch        *b_Reco_mu_normChi2_global;   //!
   TBranch        *b_Reco_mu_nPixWMea;   //!
   TBranch        *b_Reco_mu_nTrkWMea;   //!
   TBranch        *b_Reco_mu_StationsMatched;   //!
   TBranch        *b_Reco_mu_dxy;   //!
   TBranch        *b_Reco_mu_dxyErr;   //!
   TBranch        *b_Reco_mu_dz;   //!
   TBranch        *b_Reco_mu_dzErr;   //!
   TBranch        *b_Reco_mu_pt_inner;   //!
   TBranch        *b_Reco_mu_pt_global;   //!
   TBranch        *b_Reco_mu_ptErr_inner;   //!
   TBranch        *b_Reco_mu_ptErr_global;   //!

   //List of branshes from HltTree.h
   TBranch        *b_NL1IsolEm;   //!
   TBranch        *b_L1IsolEmEt;   //!
   TBranch        *b_L1IsolEmE;   //!
   TBranch        *b_L1IsolEmEta;   //!
   TBranch        *b_L1IsolEmPhi;   //!
   TBranch        *b_NL1NIsolEm;   //!
   TBranch        *b_L1NIsolEmEt;   //!
   TBranch        *b_L1NIsolEmE;   //!
   TBranch        *b_L1NIsolEmEta;   //!
   TBranch        *b_L1NIsolEmPhi;   //!
   TBranch        *b_NL1Mu;   //!
   TBranch        *b_L1MuPt;   //!
   TBranch        *b_L1MuE;   //!
   TBranch        *b_L1MuEta;   //!
   TBranch        *b_L1MuPhi;   //!
   TBranch        *b_L1MuIsol;   //!
   TBranch        *b_L1MuMip;   //!
   TBranch        *b_L1MuFor;   //!
   TBranch        *b_L1MuRPC;   //!
   TBranch        *b_L1MuQal;   //!
   TBranch        *b_L1MuChg;   //!
   TBranch        *b_NL1CenJet;   //!
   TBranch        *b_L1CenJetEt;   //!
   TBranch        *b_L1CenJetE;   //!
   TBranch        *b_L1CenJetEta;   //!
   TBranch        *b_L1CenJetPhi;   //!
   TBranch        *b_NL1ForJet;   //!
   TBranch        *b_L1ForJetEt;   //!
   TBranch        *b_L1ForJetE;   //!
   TBranch        *b_L1ForJetEta;   //!
   TBranch        *b_L1ForJetPhi;   //!
   TBranch        *b_NL1Tau;   //!
   TBranch        *b_L1TauEt;   //!
   TBranch        *b_L1TauE;   //!
   TBranch        *b_L1TauEta;   //!
   TBranch        *b_L1TauPhi;   //!
   TBranch        *b_L1Met;   //!
   TBranch        *b_L1MetPhi;   //!
   TBranch        *b_L1EtTot;   //!
   TBranch        *b_L1Mht;   //!
   TBranch        *b_L1MhtPhi;   //!
   TBranch        *b_L1EtHad;   //!
   TBranch        *b_L1HfRing1EtSumPositiveEta;   //!
   TBranch        *b_L1HfRing2EtSumPositiveEta;   //!
   TBranch        *b_L1HfRing1EtSumNegativeEta;   //!
   TBranch        *b_L1HfRing2EtSumNegativeEta;   //!
   TBranch        *b_L1HfTowerCountPositiveEtaRing1;   //!
   TBranch        *b_L1HfTowerCountNegativeEtaRing1;   //!
   TBranch        *b_L1HfTowerCountPositiveEtaRing2;   //!
   TBranch        *b_L1HfTowerCountNegativeEtaRing2;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Bx;   //!
   TBranch        *b_Orbit;   //!
   TBranch        *b_AvgInstDelLumi;   //!
   TBranch        *b_HLTriggerFirstPath;   //!
   TBranch        *b_HLTriggerFirstPath_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity60_v2;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity60_v2_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity85_v2;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity85_v2_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity110_v2;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity110_v2_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity135_v2;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity135_v2_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity160_v2;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity160_v2_Prescl;   //!
   TBranch        *b_HLT_Physics_v2;   //!
   TBranch        *b_HLT_Physics_v2_Prescl;   //!
   TBranch        *b_DST_Physics_v1;   //!
   TBranch        *b_DST_Physics_v1_Prescl;   //!
   TBranch        *b_HLT_Random_v1;   //!
   TBranch        *b_HLT_Random_v1_Prescl;   //!
   TBranch        *b_HLT_EcalCalibration_v1;   //!
   TBranch        *b_HLT_EcalCalibration_v1_Prescl;   //!
   TBranch        *b_HLT_HcalCalibration_v1;   //!
   TBranch        *b_HLT_HcalCalibration_v1_Prescl;   //!
   TBranch        *b_AlCa_EcalPhiSym_v3;   //!
   TBranch        *b_AlCa_EcalPhiSym_v3_Prescl;   //!
   TBranch        *b_HLT_L1Tech6_BPTX_MinusOnly_v1;   //!
   TBranch        *b_HLT_L1Tech6_BPTX_MinusOnly_v1_Prescl;   //!
   TBranch        *b_HLT_L1Tech5_BPTX_PlusOnly_v2;   //!
   TBranch        *b_HLT_L1Tech5_BPTX_PlusOnly_v2_Prescl;   //!
   TBranch        *b_HLT_L1Tech7_NoBPTX_v1;   //!
   TBranch        *b_HLT_L1Tech7_NoBPTX_v1_Prescl;   //!
   TBranch        *b_HLT_L1TOTEM1_MinBias_v1;   //!
   TBranch        *b_HLT_L1TOTEM1_MinBias_v1_Prescl;   //!
   TBranch        *b_HLT_L1TOTEM2_ZeroBias_v1;   //!
   TBranch        *b_HLT_L1TOTEM2_ZeroBias_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF2OR_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF2OR_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1AND_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1AND_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF2AND_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF2AND_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF2ORNoBptxGating_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF2ORNoBptxGating_v1_Prescl;   //!
   TBranch        *b_AlCa_RPCMuonNoTriggers_v2;   //!
   TBranch        *b_AlCa_RPCMuonNoTriggers_v2_Prescl;   //!
   TBranch        *b_AlCa_RPCMuonNoHits_v2;   //!
   TBranch        *b_AlCa_RPCMuonNoHits_v2_Prescl;   //!
   TBranch        *b_AlCa_RPCMuonNormalisation_v2;   //!
   TBranch        *b_AlCa_RPCMuonNormalisation_v2_Prescl;   //!
   TBranch        *b_AlCa_LumiPixels_Random_v1;   //!
   TBranch        *b_AlCa_LumiPixels_Random_v1_Prescl;   //!
   TBranch        *b_AlCa_LumiPixels_ZeroBias_v2;   //!
   TBranch        *b_AlCa_LumiPixels_ZeroBias_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_v2;   //!
   TBranch        *b_HLT_ZeroBias_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part0_v2;   //!
   TBranch        *b_HLT_ZeroBias_part0_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part1_v2;   //!
   TBranch        *b_HLT_ZeroBias_part1_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part2_v2;   //!
   TBranch        *b_HLT_ZeroBias_part2_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part3_v2;   //!
   TBranch        *b_HLT_ZeroBias_part3_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part4_v2;   //!
   TBranch        *b_HLT_ZeroBias_part4_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part5_v2;   //!
   TBranch        *b_HLT_ZeroBias_part5_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part6_v2;   //!
   TBranch        *b_HLT_ZeroBias_part6_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part7_v2;   //!
   TBranch        *b_HLT_ZeroBias_part7_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part8_v2;   //!
   TBranch        *b_HLT_ZeroBias_part8_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part9_v2;   //!
   TBranch        *b_HLT_ZeroBias_part9_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part10_v2;   //!
   TBranch        *b_HLT_ZeroBias_part10_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part11_v2;   //!
   TBranch        *b_HLT_ZeroBias_part11_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part12_v2;   //!
   TBranch        *b_HLT_ZeroBias_part12_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part13_v2;   //!
   TBranch        *b_HLT_ZeroBias_part13_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part14_v2;   //!
   TBranch        *b_HLT_ZeroBias_part14_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part15_v2;   //!
   TBranch        *b_HLT_ZeroBias_part15_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part16_v2;   //!
   TBranch        *b_HLT_ZeroBias_part16_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part17_v2;   //!
   TBranch        *b_HLT_ZeroBias_part17_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part18_v2;   //!
   TBranch        *b_HLT_ZeroBias_part18_v2_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part19_v2;   //!
   TBranch        *b_HLT_ZeroBias_part19_v2_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet40_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet40_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet60_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet60_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet80_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet80_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet100_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet100_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet110_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet110_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet120_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet120_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet150_v1;   //!
   TBranch        *b_HLT_AK4CaloJet150_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFJet40_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4PFJet40_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFJet60_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4PFJet60_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFJet80_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4PFJet80_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFJet100_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4PFJet100_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFJet110_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4PFJet110_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFJet120_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4PFJet120_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet80_Jet35_Eta1p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet80_Jet35_Eta1p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet80_Jet35_Eta0p7_v1;   //!
   TBranch        *b_HLT_AK4CaloJet80_Jet35_Eta0p7_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet100_Jet35_Eta1p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet100_Jet35_Eta1p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet100_Jet35_Eta0p7_v1;   //!
   TBranch        *b_HLT_AK4CaloJet100_Jet35_Eta0p7_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet80_45_45_Eta2p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet80_45_45_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton10_Eta1p5_v1;   //!
   TBranch        *b_HLT_HISinglePhoton10_Eta1p5_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton15_Eta1p5_v1;   //!
   TBranch        *b_HLT_HISinglePhoton15_Eta1p5_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton20_Eta1p5_v1;   //!
   TBranch        *b_HLT_HISinglePhoton20_Eta1p5_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton30_Eta1p5_v1;   //!
   TBranch        *b_HLT_HISinglePhoton30_Eta1p5_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton40_Eta1p5_v1;   //!
   TBranch        *b_HLT_HISinglePhoton40_Eta1p5_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton50_Eta1p5_v1;   //!
   TBranch        *b_HLT_HISinglePhoton50_Eta1p5_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton60_Eta1p5_v1;   //!
   TBranch        *b_HLT_HISinglePhoton60_Eta1p5_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton10_Eta3p1_v1;   //!
   TBranch        *b_HLT_HISinglePhoton10_Eta3p1_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton15_Eta3p1_v1;   //!
   TBranch        *b_HLT_HISinglePhoton15_Eta3p1_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton20_Eta3p1_v1;   //!
   TBranch        *b_HLT_HISinglePhoton20_Eta3p1_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton30_Eta3p1_v1;   //!
   TBranch        *b_HLT_HISinglePhoton30_Eta3p1_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton40_Eta3p1_v1;   //!
   TBranch        *b_HLT_HISinglePhoton40_Eta3p1_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton50_Eta3p1_v1;   //!
   TBranch        *b_HLT_HISinglePhoton50_Eta3p1_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton60_Eta3p1_v1;   //!
   TBranch        *b_HLT_HISinglePhoton60_Eta3p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_v1;   //!
   TBranch        *b_HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_v1_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_R9HECut_v1;   //!
   TBranch        *b_HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_R9HECut_v1_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton15_Eta2p1_Mass50_1000_R9Cut_v1;   //!
   TBranch        *b_HLT_HIDoublePhoton15_Eta2p1_Mass50_1000_R9Cut_v1_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1;   //!
   TBranch        *b_HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_v1;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_v1_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMu10_v1;   //!
   TBranch        *b_HLT_HIL1DoubleMu10_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v1;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1;   //!
   TBranch        *b_HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMu0_OS_m7to14_v1;   //!
   TBranch        *b_HLT_HIL3DoubleMu0_OS_m7to14_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ10_v1;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ10_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu3_NHitQ15_v1;   //!
   TBranch        *b_HLT_HIL3Mu3_NHitQ15_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu5_NHitQ10_v1;   //!
   TBranch        *b_HLT_HIL2Mu5_NHitQ10_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu5_NHitQ15_v1;   //!
   TBranch        *b_HLT_HIL3Mu5_NHitQ15_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu7_NHitQ10_v1;   //!
   TBranch        *b_HLT_HIL2Mu7_NHitQ10_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu7_NHitQ15_v1;   //!
   TBranch        *b_HLT_HIL3Mu7_NHitQ15_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu15_v1;   //!
   TBranch        *b_HLT_HIL2Mu15_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu15_v1;   //!
   TBranch        *b_HLT_HIL3Mu15_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu20_v1;   //!
   TBranch        *b_HLT_HIL2Mu20_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu20_v1;   //!
   TBranch        *b_HLT_HIL3Mu20_v1_Prescl;   //!
   TBranch        *b_HLT_FullTrack18ForPPRef_v1;   //!
   TBranch        *b_HLT_FullTrack18ForPPRef_v1_Prescl;   //!
   TBranch        *b_HLT_FullTrack24ForPPRef_v1;   //!
   TBranch        *b_HLT_FullTrack24ForPPRef_v1_Prescl;   //!
   TBranch        *b_HLT_FullTrack34ForPPRef_v1;   //!
   TBranch        *b_HLT_FullTrack34ForPPRef_v1_Prescl;   //!
   TBranch        *b_HLT_FullTrack45ForPPRef_v1;   //!
   TBranch        *b_HLT_FullTrack45ForPPRef_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCL1DoubleMuOpenNotHF2_v1;   //!
   TBranch        *b_HLT_HIUPCL1DoubleMuOpenNotHF2_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCDoubleMuNotHF2Pixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCDoubleMuNotHF2Pixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCL1MuOpen_NotMinimumBiasHF2_AND_v1;   //!
   TBranch        *b_HLT_HIUPCL1MuOpen_NotMinimumBiasHF2_AND_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCMuOpen_NotMinimumBiasHF2_ANDPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCMuOpen_NotMinimumBiasHF2_ANDPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCL1NotMinimumBiasHF2_AND_v1;   //!
   TBranch        *b_HLT_HIUPCL1NotMinimumBiasHF2_AND_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCNotMinimumBiasHF2_ANDPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCNotMinimumBiasHF2_ANDPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCL1ZdcOR_BptxAND_v1;   //!
   TBranch        *b_HLT_HIUPCL1ZdcOR_BptxAND_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCZdcOR_BptxANDPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCZdcOR_BptxANDPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCL1ZdcXOR_BptxAND_v1;   //!
   TBranch        *b_HLT_HIUPCL1ZdcXOR_BptxAND_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCZdcXOR_BptxANDPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCZdcXOR_BptxANDPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCL1NotZdcOR_BptxAND_v1;   //!
   TBranch        *b_HLT_HIUPCL1NotZdcOR_BptxAND_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCNotZdcOR_BptxANDPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCNotZdcOR_BptxANDPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIL1CastorMediumJet_v1;   //!
   TBranch        *b_HLT_HIL1CastorMediumJet_v1_Prescl;   //!
   TBranch        *b_HLT_HICastorMediumJetPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HICastorMediumJetPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt15_v1;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt15_v1_Prescl;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt20_v1;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt20_v1_Prescl;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt30_v1;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt30_v1_Prescl;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt40_v1;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt40_v1_Prescl;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt50_v1;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt50_v1_Prescl;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt60_v1;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt60_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part0_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part0_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part1_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part1_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part2_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part2_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part3_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part3_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part4_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part4_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part5_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part5_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part6_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part6_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part7_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part7_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part8_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part8_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part9_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part9_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part10_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part10_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part11_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part11_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part12_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part12_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part13_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part13_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part14_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part14_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part15_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part15_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part16_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part16_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part17_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part17_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part18_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part18_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part19_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_part19_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFBJetBCSV60_Eta2p1_v1;   //!
   TBranch        *b_HLT_AK4PFBJetBCSV60_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFBJetBCSV80_Eta2p1_v1;   //!
   TBranch        *b_HLT_AK4PFBJetBCSV80_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFDJet60_Eta2p1_v1;   //!
   TBranch        *b_HLT_AK4PFDJet60_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFDJet80_Eta2p1_v1;   //!
   TBranch        *b_HLT_AK4PFDJet80_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFBJetBSSV60_Eta2p1_v1;   //!
   TBranch        *b_HLT_AK4PFBJetBSSV60_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFBJetBSSV80_Eta2p1_v1;   //!
   TBranch        *b_HLT_AK4PFBJetBSSV80_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLTriggerFinalPath;   //!
   TBranch        *b_HLTriggerFinalPath_Prescl;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt8_v1;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt8_v1_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity60_v1;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity60_v1_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity60_v3;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity60_v3_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity85_v1;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity85_v1_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity110_v1;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity110_v1_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity135_v1;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity135_v1_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity160_v1;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity160_v1_Prescl;   //!
   TBranch        *b_HLT_Physics_v1;   //!
   TBranch        *b_HLT_Physics_v1_Prescl;   //!
   TBranch        *b_HLT_L1Tech5_BPTX_PlusOnly_v1;   //!
   TBranch        *b_HLT_L1Tech5_BPTX_PlusOnly_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_v1;   //!
   TBranch        *b_HLT_ZeroBias_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part0_v1;   //!
   TBranch        *b_HLT_ZeroBias_part0_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part1_v1;   //!
   TBranch        *b_HLT_ZeroBias_part1_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part2_v1;   //!
   TBranch        *b_HLT_ZeroBias_part2_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part3_v1;   //!
   TBranch        *b_HLT_ZeroBias_part3_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part4_v1;   //!
   TBranch        *b_HLT_ZeroBias_part4_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part5_v1;   //!
   TBranch        *b_HLT_ZeroBias_part5_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part6_v1;   //!
   TBranch        *b_HLT_ZeroBias_part6_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part7_v1;   //!
   TBranch        *b_HLT_ZeroBias_part7_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part8_v1;   //!
   TBranch        *b_HLT_ZeroBias_part8_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part9_v1;   //!
   TBranch        *b_HLT_ZeroBias_part9_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part10_v1;   //!
   TBranch        *b_HLT_ZeroBias_part10_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part11_v1;   //!
   TBranch        *b_HLT_ZeroBias_part11_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part12_v1;   //!
   TBranch        *b_HLT_ZeroBias_part12_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part13_v1;   //!
   TBranch        *b_HLT_ZeroBias_part13_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part14_v1;   //!
   TBranch        *b_HLT_ZeroBias_part14_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part15_v1;   //!
   TBranch        *b_HLT_ZeroBias_part15_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part16_v1;   //!
   TBranch        *b_HLT_ZeroBias_part16_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part17_v1;   //!
   TBranch        *b_HLT_ZeroBias_part17_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part18_v1;   //!
   TBranch        *b_HLT_ZeroBias_part18_v1_Prescl;   //!
   TBranch        *b_HLT_ZeroBias_part19_v1;   //!
   TBranch        *b_HLT_ZeroBias_part19_v1_Prescl;   //!
   TBranch        *b_HLT_FullTrack18ForPPRef_v2;   //!
   TBranch        *b_HLT_FullTrack18ForPPRef_v2_Prescl;   //!
   TBranch        *b_HLT_FullTrack18ForPPRef_v3;   //!
   TBranch        *b_HLT_FullTrack18ForPPRef_v3_Prescl;   //!
   TBranch        *b_HLT_FullTrack24ForPPRef_v2;   //!
   TBranch        *b_HLT_FullTrack24ForPPRef_v2_Prescl;   //!
   TBranch        *b_HLT_FullTrack24ForPPRef_v3;   //!
   TBranch        *b_HLT_FullTrack24ForPPRef_v3_Prescl;   //!
   TBranch        *b_HLT_FullTrack34ForPPRef_v2;   //!
   TBranch        *b_HLT_FullTrack34ForPPRef_v2_Prescl;   //!
   TBranch        *b_HLT_FullTrack34ForPPRef_v3;   //!
   TBranch        *b_HLT_FullTrack34ForPPRef_v3_Prescl;   //!
   TBranch        *b_HLT_FullTrack34ForPPRef_v4;   //!
   TBranch        *b_HLT_FullTrack34ForPPRef_v4_Prescl;   //!
   TBranch        *b_HLT_FullTrack45ForPPRef_v2;   //!
   TBranch        *b_HLT_FullTrack45ForPPRef_v2_Prescl;   //!
   TBranch        *b_HLT_FullTrack45ForPPRef_v3;   //!
   TBranch        *b_HLT_FullTrack45ForPPRef_v3_Prescl;   //!
   TBranch        *b_HLT_FullTrack53ForPPRef_v1;   //!
   TBranch        *b_HLT_FullTrack53ForPPRef_v1_Prescl;   //!
   TBranch        *b_HLT_FullTrack53ForPPRef_v2;   //!
   TBranch        *b_HLT_FullTrack53ForPPRef_v2_Prescl;   //!
   TBranch        *b_HLT_FullTrack53ForPPRef_v3;   //!
   TBranch        *b_HLT_FullTrack53ForPPRef_v3_Prescl;   //!
   TBranch        *b_L1_AlwaysTrue;   //!
   TBranch        *b_L1_AlwaysTrue_Prescl;   //!
   TBranch        *b_L1_CastorHighJet_BptxAND;   //!
   TBranch        *b_L1_CastorHighJet_BptxAND_Prescl;   //!
   TBranch        *b_L1_CastorMediumJet_BptxAND;   //!
   TBranch        *b_L1_CastorMediumJet_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleMu0_BptxAND;   //!
   TBranch        *b_L1_DoubleMu0_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleMu0_MinimumBiasHF1_AND;   //!
   TBranch        *b_L1_DoubleMu0_MinimumBiasHF1_AND_Prescl;   //!
   TBranch        *b_L1_DoubleMu3_BptxAND;   //!
   TBranch        *b_L1_DoubleMu3_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleMuOpen_BptxAND;   //!
   TBranch        *b_L1_DoubleMuOpen_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleMuOpen_NotMinimumBiasHF2_AND;   //!
   TBranch        *b_L1_DoubleMuOpen_NotMinimumBiasHF2_AND_Prescl;   //!
   TBranch        *b_L1_DoubleJet20;   //!
   TBranch        *b_L1_DoubleJet20_Prescl;   //!
   TBranch        *b_L1_DoubleJet28_BptxAND;   //!
   TBranch        *b_L1_DoubleJet28_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleJet32;   //!
   TBranch        *b_L1_DoubleJet32_Prescl;   //!
   TBranch        *b_L1_ETT130;   //!
   TBranch        *b_L1_ETT130_Prescl;   //!
   TBranch        *b_L1_ETT15_BptxAND;   //!
   TBranch        *b_L1_ETT15_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT20_BptxAND;   //!
   TBranch        *b_L1_ETT20_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT30_BptxAND;   //!
   TBranch        *b_L1_ETT30_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT40;   //!
   TBranch        *b_L1_ETT40_Prescl;   //!
   TBranch        *b_L1_ETT40_BptxAND;   //!
   TBranch        *b_L1_ETT40_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT45_BptxAND;   //!
   TBranch        *b_L1_ETT45_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT50_BptxAND;   //!
   TBranch        *b_L1_ETT50_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT55_BptxAND;   //!
   TBranch        *b_L1_ETT55_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT60_BptxAND;   //!
   TBranch        *b_L1_ETT60_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT70_BptxAND;   //!
   TBranch        *b_L1_ETT70_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT90_BptxAND;   //!
   TBranch        *b_L1_ETT90_BptxAND_Prescl;   //!
   TBranch        *b_L1_MinimumBiasHF1_AND;   //!
   TBranch        *b_L1_MinimumBiasHF1_AND_Prescl;   //!
   TBranch        *b_L1_MinimumBiasHF1_OR;   //!
   TBranch        *b_L1_MinimumBiasHF1_OR_Prescl;   //!
   TBranch        *b_L1_MinimumBiasHF1_XOR;   //!
   TBranch        *b_L1_MinimumBiasHF1_XOR_Prescl;   //!
   TBranch        *b_L1_MinimumBiasHF2_AND;   //!
   TBranch        *b_L1_MinimumBiasHF2_AND_Prescl;   //!
   TBranch        *b_L1_MinimumBiasHF2_OR;   //!
   TBranch        *b_L1_MinimumBiasHF2_OR_Prescl;   //!
   TBranch        *b_L1_MinimumBiasHF2_OR_NoBptxGating;   //!
   TBranch        *b_L1_MinimumBiasHF2_OR_NoBptxGating_Prescl;   //!
   TBranch        *b_L1_MuOpen_NotMinimumBiasHF2_AND;   //!
   TBranch        *b_L1_MuOpen_NotMinimumBiasHF2_AND_Prescl;   //!
   TBranch        *b_L1_NotMinimumBiasHF1_OR;   //!
   TBranch        *b_L1_NotMinimumBiasHF1_OR_Prescl;   //!
   TBranch        *b_L1_NotMinimumBiasHF2_AND;   //!
   TBranch        *b_L1_NotMinimumBiasHF2_AND_Prescl;   //!
   TBranch        *b_L1_NotZdcOR_BptxAND;   //!
   TBranch        *b_L1_NotZdcOR_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG12_BptxAND;   //!
   TBranch        *b_L1_SingleEG12_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG15_BptxAND;   //!
   TBranch        *b_L1_SingleEG15_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG20;   //!
   TBranch        *b_L1_SingleEG20_Prescl;   //!
   TBranch        *b_L1_SingleEG20_BptxAND;   //!
   TBranch        *b_L1_SingleEG20_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG2_BptxAND;   //!
   TBranch        *b_L1_SingleEG2_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG30_BptxAND;   //!
   TBranch        *b_L1_SingleEG30_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG5;   //!
   TBranch        *b_L1_SingleEG5_Prescl;   //!
   TBranch        *b_L1_SingleEG5_BptxAND;   //!
   TBranch        *b_L1_SingleEG5_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG7_BptxAND;   //!
   TBranch        *b_L1_SingleEG7_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu12_BptxAND;   //!
   TBranch        *b_L1_SingleMu12_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu16_BptxAND;   //!
   TBranch        *b_L1_SingleMu16_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu16_MinimumBiasHF1_AND;   //!
   TBranch        *b_L1_SingleMu16_MinimumBiasHF1_AND_Prescl;   //!
   TBranch        *b_L1_SingleMu3_BptxAND;   //!
   TBranch        *b_L1_SingleMu3_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu3_MinimumBiasHF1_AND;   //!
   TBranch        *b_L1_SingleMu3_MinimumBiasHF1_AND_Prescl;   //!
   TBranch        *b_L1_SingleMu3p5;   //!
   TBranch        *b_L1_SingleMu3p5_Prescl;   //!
   TBranch        *b_L1_SingleMu5_BptxAND;   //!
   TBranch        *b_L1_SingleMu5_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu7_BptxAND;   //!
   TBranch        *b_L1_SingleMu7_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMuBeamHalo;   //!
   TBranch        *b_L1_SingleMuBeamHalo_Prescl;   //!
   TBranch        *b_L1_SingleMuOpen;   //!
   TBranch        *b_L1_SingleMuOpen_Prescl;   //!
   TBranch        *b_L1_SingleMuOpen_BptxAND;   //!
   TBranch        *b_L1_SingleMuOpen_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMuOpen_NotBptxOR;   //!
   TBranch        *b_L1_SingleMuOpen_NotBptxOR_Prescl;   //!
   TBranch        *b_L1_SingleJet12_BptxAND;   //!
   TBranch        *b_L1_SingleJet12_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet16;   //!
   TBranch        *b_L1_SingleJet16_Prescl;   //!
   TBranch        *b_L1_SingleJet16_BptxAND;   //!
   TBranch        *b_L1_SingleJet16_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet200;   //!
   TBranch        *b_L1_SingleJet200_Prescl;   //!
   TBranch        *b_L1_SingleJet20_BptxAND;   //!
   TBranch        *b_L1_SingleJet20_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet24_BptxAND;   //!
   TBranch        *b_L1_SingleJet24_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet28_BptxAND;   //!
   TBranch        *b_L1_SingleJet28_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet32_BptxAND;   //!
   TBranch        *b_L1_SingleJet32_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet36;   //!
   TBranch        *b_L1_SingleJet36_Prescl;   //!
   TBranch        *b_L1_SingleJet36_BptxAND;   //!
   TBranch        *b_L1_SingleJet36_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet40_BptxAND;   //!
   TBranch        *b_L1_SingleJet40_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet44_BptxAND;   //!
   TBranch        *b_L1_SingleJet44_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet48_BptxAND;   //!
   TBranch        *b_L1_SingleJet48_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet52_BptxAND;   //!
   TBranch        *b_L1_SingleJet52_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet60_BptxAND;   //!
   TBranch        *b_L1_SingleJet60_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet68;   //!
   TBranch        *b_L1_SingleJet68_Prescl;   //!
   TBranch        *b_L1_SingleJet68_BptxAND;   //!
   TBranch        *b_L1_SingleJet68_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet8_BptxAND;   //!
   TBranch        *b_L1_SingleJet8_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJetC20_NotBptxOR;   //!
   TBranch        *b_L1_SingleJetC20_NotBptxOR_Prescl;   //!
   TBranch        *b_L1_SingleJetC32_NotBptxOR;   //!
   TBranch        *b_L1_SingleJetC32_NotBptxOR_Prescl;   //!
   TBranch        *b_L1_TOTEM_0;   //!
   TBranch        *b_L1_TOTEM_0_Prescl;   //!
   TBranch        *b_L1_TOTEM_1;   //!
   TBranch        *b_L1_TOTEM_1_Prescl;   //!
   TBranch        *b_L1_TOTEM_2;   //!
   TBranch        *b_L1_TOTEM_2_Prescl;   //!
   TBranch        *b_L1_TOTEM_3;   //!
   TBranch        *b_L1_TOTEM_3_Prescl;   //!
   TBranch        *b_L1_ZdcOR;   //!
   TBranch        *b_L1_ZdcOR_Prescl;   //!
   TBranch        *b_L1_ZdcOR_BptxAND;   //!
   TBranch        *b_L1_ZdcOR_BptxAND_Prescl;   //!
   TBranch        *b_L1_ZdcXOR;   //!
   TBranch        *b_L1_ZdcXOR_Prescl;   //!
   TBranch        *b_L1_ZdcXOR_BptxAND;   //!
   TBranch        *b_L1_ZdcXOR_BptxAND_Prescl;   //!
   TBranch        *b_L1_ZeroBias;   //!
   TBranch        *b_L1_ZeroBias_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_PreBPTX_v0;   //!
   TBranch        *b_L1Tech_BPTX_PreBPTX_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_minus_AND_not_plus_v0;   //!
   TBranch        *b_L1Tech_BPTX_minus_AND_not_plus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_NOT_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_instance1_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_OR_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_OR_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_quiet_v0;   //!
   TBranch        *b_L1Tech_BPTX_quiet_v0_Prescl;   //!
   TBranch        *b_L1Tech_BRIL_bit28;   //!
   TBranch        *b_L1Tech_BRIL_bit28_Prescl;   //!
   TBranch        *b_L1Tech_BRIL_bit29;   //!
   TBranch        *b_L1Tech_BRIL_bit29_Prescl;   //!
   TBranch        *b_L1Tech_BRIL_bit30;   //!
   TBranch        *b_L1Tech_BRIL_bit30_Prescl;   //!
   TBranch        *b_L1Tech_BRIL_bit31;   //!
   TBranch        *b_L1Tech_BRIL_bit31_Prescl;   //!
   TBranch        *b_L1Tech_BRIL_bit32;   //!
   TBranch        *b_L1Tech_BRIL_bit32_Prescl;   //!
   TBranch        *b_L1Tech_BRIL_bit33;   //!
   TBranch        *b_L1Tech_BRIL_bit33_Prescl;   //!
   TBranch        *b_L1Tech_BRIL_bit34;   //!
   TBranch        *b_L1Tech_BRIL_bit34_Prescl;   //!
   TBranch        *b_L1Tech_BRIL_bit35;   //!
   TBranch        *b_L1Tech_BRIL_bit35_Prescl;   //!
   TBranch        *b_L1Tech_BRIL_bit36;   //!
   TBranch        *b_L1Tech_BRIL_bit36_Prescl;   //!
   TBranch        *b_L1Tech_BRIL_bit37;   //!
   TBranch        *b_L1Tech_BRIL_bit37_Prescl;   //!
   TBranch        *b_L1Tech_BRIL_bit38;   //!
   TBranch        *b_L1Tech_BRIL_bit38_Prescl;   //!
   TBranch        *b_L1Tech_BRIL_bit39;   //!
   TBranch        *b_L1Tech_BRIL_bit39_Prescl;   //!
   TBranch        *b_L1Tech_BRIL_bit40;   //!
   TBranch        *b_L1Tech_BRIL_bit40_Prescl;   //!
   TBranch        *b_L1Tech_BRIL_bit41;   //!
   TBranch        *b_L1Tech_BRIL_bit41_Prescl;   //!
   TBranch        *b_L1Tech_BRIL_bit42;   //!
   TBranch        *b_L1Tech_BRIL_bit42_Prescl;   //!
   TBranch        *b_L1Tech_BRIL_bit43;   //!
   TBranch        *b_L1Tech_BRIL_bit43_Prescl;   //!
   TBranch        *b_L1Tech_CASTOR_Gap_v0;   //!
   TBranch        *b_L1Tech_CASTOR_Gap_v0_Prescl;   //!
   TBranch        *b_L1Tech_CASTOR_HaloMuon_v0;   //!
   TBranch        *b_L1Tech_CASTOR_HaloMuon_v0_Prescl;   //!
   TBranch        *b_L1Tech_CASTOR_HighJet_v0;   //!
   TBranch        *b_L1Tech_CASTOR_HighJet_v0_Prescl;   //!
   TBranch        *b_L1Tech_CASTOR_MediumJet_v0;   //!
   TBranch        *b_L1Tech_CASTOR_MediumJet_v0_Prescl;   //!
   TBranch        *b_L1Tech_DT_GlobalOR_v0;   //!
   TBranch        *b_L1Tech_DT_GlobalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HBHE_totalOR_v0;   //!
   TBranch        *b_L1Tech_HCAL_HBHE_totalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_MMP_or_MPP_v1;   //!
   TBranch        *b_L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_coincidence_PM_v2;   //!
   TBranch        *b_L1Tech_HCAL_HF_coincidence_PM_v2_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_single_channel_v0;   //!
   TBranch        *b_L1Tech_HCAL_HF_single_channel_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HO_totalOR_v0;   //!
   TBranch        *b_L1Tech_HCAL_HO_totalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_barrel_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_pointing_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_TOTEM_0;   //!
   TBranch        *b_L1Tech_TOTEM_0_Prescl;   //!
   TBranch        *b_L1Tech_TOTEM_1;   //!
   TBranch        *b_L1Tech_TOTEM_1_Prescl;   //!
   TBranch        *b_L1Tech_TOTEM_2;   //!
   TBranch        *b_L1Tech_TOTEM_2_Prescl;   //!
   TBranch        *b_L1Tech_TOTEM_3;   //!
   TBranch        *b_L1Tech_TOTEM_3_Prescl;   //!
   TBranch        *b_L1Tech_ZDC_minus;   //!
   TBranch        *b_L1Tech_ZDC_minus_Prescl;   //!
   TBranch        *b_L1Tech_ZDC_plus;   //!
   TBranch        *b_L1Tech_ZDC_plus_Prescl;   //!

   //List of branshes from HltTree1.h
   TBranch        *b_Onia2MuMuPAT;   //!
   TBranch        *b_ana_step;   //!
   TBranch        *b_pHBHENoiseFilterResultProducer;   //!
   TBranch        *b_HBHENoiseFilterResult;   //!
   TBranch        *b_HBHENoiseFilterResultRun1;   //!
   TBranch        *b_HBHENoiseFilterResultRun2Loose;   //!
   TBranch        *b_HBHENoiseFilterResultRun2Tight;   //!
   TBranch        *b_HBHEIsoNoiseFilterResult;   //!
   TBranch        *b_pPAprimaryVertexFilter;   //!
   TBranch        *b_pBeamScrapingFilter;   //!
   TBranch        *b_pVertexFilterCutG;   //!
   TBranch        *b_pVertexFilterCutGloose;   //!
   TBranch        *b_pVertexFilterCutGtight;   //!
   TBranch        *b_pVertexFilterCutGplus;   //!
   TBranch        *b_pVertexFilterCutE;   //!
   TBranch        *b_pVertexFilterCutEandG;   //!

   // List of branches from HiTree.h
   TBranch        *b_run;   //!
   //   TBranch        *b_evt;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_hiBin;   //!
   TBranch        *b_hiHF;   //!
   TBranch        *b_hiHFplus;   //!
   TBranch        *b_hiHFminus;   //!
   TBranch        *b_hiHFplusEta4;   //!
   TBranch        *b_hiHFminusEta4;   //!
   TBranch        *b_hiZDC;   //!
   TBranch        *b_hiZDCplus;   //!
   TBranch        *b_hiZDCminus;   //!
   TBranch        *b_hiHFhit;   //!
   TBranch        *b_hiHFhitPlus;   //!
   TBranch        *b_hiHFhitMinus;   //!
   TBranch        *b_hiET;   //!
   TBranch        *b_hiEE;   //!
   TBranch        *b_hiEB;   //!
   TBranch        *b_hiEEplus;   //!
   TBranch        *b_hiEEminus;   //!
   TBranch        *b_hiNpix;   //!
   TBranch        *b_hiNpixelTracks;   //!
   TBranch        *b_hiNtracks;   //!
   TBranch        *b_hiNtracksPtCut;   //!
   TBranch        *b_hiNtracksEtaCut;   //!
   TBranch        *b_hiNtracksEtaPtCut;   //!

   // List of branches from t.h
   TBranch        *b_evt;   //!
   TBranch        *b_b;   //!
   TBranch        *b_nref;   //!
   TBranch        *b_rawpt;   //!
   TBranch        *b_jtpt;   //!
   TBranch        *b_jteta;   //!
   TBranch        *b_jty;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_jtpu;   //!
   TBranch        *b_jtm;   //!
   TBranch        *b_jtarea;   //!
   TBranch        *b_jtPfCHF;   //!
   TBranch        *b_jtPfNHF;   //!
   TBranch        *b_jtPfCEF;   //!
   TBranch        *b_jtPfNEF;   //!
   TBranch        *b_jtPfMUF;   //!
   TBranch        *b_jtPfCHM;   //!
   TBranch        *b_jtPfNHM;   //!
   TBranch        *b_jtPfCEM;   //!
   TBranch        *b_jtPfNEM;   //!
   TBranch        *b_jtPfMUM;   //!
   TBranch        *b_jttau1;   //!
   TBranch        *b_jttau2;   //!
   TBranch        *b_jttau3;   //!
   TBranch        *b_discr_jetID_cuts;   //!
   TBranch        *b_discr_jetID_bdt;   //!
   TBranch        *b_discr_fr01;   //!
   TBranch        *b_trackMax;   //!
   TBranch        *b_trackSum;   //!
   TBranch        *b_trackN;   //!
   TBranch        *b_trackHardSum;   //!
   TBranch        *b_trackHardN;   //!
   TBranch        *b_chargedMax;   //!
   TBranch        *b_chargedSum;   //!
   TBranch        *b_chargedN;   //!
   TBranch        *b_chargedHardSum;   //!
   TBranch        *b_chargedHardN;   //!
   TBranch        *b_photonMax;   //!
   TBranch        *b_photonSum;   //!
   TBranch        *b_photonN;   //!
   TBranch        *b_photonHardSum;   //!
   TBranch        *b_photonHardN;   //!
   TBranch        *b_neutralMax;   //!
   TBranch        *b_neutralSum;   //!
   TBranch        *b_neutralN;   //!
   TBranch        *b_hcalSum;   //!
   TBranch        *b_ecalSum;   //!
   TBranch        *b_eMax;   //!
   TBranch        *b_eSum;   //!
   TBranch        *b_eN;   //!
   TBranch        *b_muMax;   //!
   TBranch        *b_muSum;   //!
   TBranch        *b_muN;   //!
   TBranch        *b_discr_ssvHighEff;   //!
   TBranch        *b_discr_ssvHighPur;   //!
   TBranch        *b_discr_csvV1;   //!
   TBranch        *b_discr_csvV2;   //!
   TBranch        *b_discr_muByIp3;   //!
   TBranch        *b_discr_muByPt;   //!
   TBranch        *b_discr_prob;   //!
   TBranch        *b_discr_probb;   //!
   TBranch        *b_discr_tcHighEff;   //!
   TBranch        *b_discr_tcHighPur;   //!
   TBranch        *b_ndiscr_ssvHighEff;   //!
   TBranch        *b_ndiscr_ssvHighPur;   //!
   TBranch        *b_ndiscr_csvV1;   //!
   TBranch        *b_ndiscr_csvV2;   //!
   TBranch        *b_ndiscr_muByPt;   //!
   TBranch        *b_pdiscr_csvV1;   //!
   TBranch        *b_pdiscr_csvV2;   //!
   TBranch        *b_nsvtx;   //!
   TBranch        *b_svtxntrk;   //!
   TBranch        *b_svtxdl;   //!
   TBranch        *b_svtxdls;   //!
   TBranch        *b_svtxdl2d;   //!
   TBranch        *b_svtxdls2d;   //!
   TBranch        *b_svtxm;   //!
   TBranch        *b_svtxpt;   //!
   TBranch        *b_svtxmcorr;   //!
   TBranch        *b_nIPtrk;   //!
   TBranch        *b_nselIPtrk;   //!
   TBranch        *b_mue;   //!
   TBranch        *b_mupt;   //!
   TBranch        *b_mueta;   //!
   TBranch        *b_muphi;   //!
   TBranch        *b_mudr;   //!
   TBranch        *b_muptrel;   //!
   TBranch        *b_muchg;   //!

   //myTree(const char * input_file = "/data_CMS/cms/mnguyen/jPsiJet/HiForestAOD.root");

   //myTree(const char * input_file = "/grid_mnt/vol__vol_U__u/llr/cms/mnguyen/prod/forest/compositeJet/CMSSW_7_5_8_patch5/src/HeavyIonsAnalysis/JetAnalysis/test/HiForestAODWithJPsi.root");

   myTree(const char * input_file =" /data_CMS/cms/mnguyen/jPsiJet/v2/merged_HiForestAOD.root");

   virtual ~myTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int q=0);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);


   virtual Bool_t isTriggerMatch (Int_t iRecoQQ, Int_t TriggerBit);
   virtual Bool_t isGlobalMuonInAccept2015 (TLorentzVector* Muon);
   virtual Bool_t areMuonsInAcceptance2015 (Int_t iRecoQQ);
   virtual Bool_t passQualityCuts2015 (Int_t iRecoQQ);
};

#endif

#ifdef myTree_cxx
myTree::myTree(const char * input_file) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  TFile* f(0x0);
if (!strcmp(input_file,"")) {

  //f = (TFile*)gROOT->GetListOfFiles()->FindObject("/grid_mnt/vol__vol_U__u/llr/cms/mnguyen/prod/forest/compositeJet/CMSSW_7_5_8_patch5/src/HeavyIonsAnalysis/JetAnalysis/test/HiForestAODWithJPsi.root");

  //f = (TFile*)gROOT->GetListOfFiles()->FindObject(" /data_CMS/cms/mnguyen/jPsiJet/HiForestAOD.root");

  f = (TFile*)gROOT->GetListOfFiles()->FindObject(" /data_CMS/cms/mnguyen/jPsiJet/v2/merged_HiForestAOD.root");

 }
      if (!f || !f->IsOpen()) {
         f = new TFile(input_file);
      }
      TTree * tree (0x0);
      TTree *hitree (0x0);
      TTree *hlttree(0x0);
      TTree *hlttree1(0x0);
      TTree *tr(0x0);
      TDirectory * dir = (TDirectory*)f->Get("hionia");
      dir->GetObject("myTree",tree);
      if(!tree)
	cout <<"error in the onia tree"<<endl;

      TDirectory * dirhi = (TDirectory*)f->Get("hiEvtAnalyzer");
      dirhi ->GetObject("HiTree",hitree);
 if(!hitree)
   {cout <<"error in the hitree"<<endl;}

      TDirectory * dirhlt = (TDirectory*)f->Get("hltanalysis");
      dirhlt ->GetObject("HltTree",hlttree);
 if(!hlttree)
   {cout <<"error in the hlttree"<<endl;}

      TDirectory * dirhlt1 = (TDirectory*)f->Get("skimanalysis");
      dirhlt1->GetObject("HltTree",hlttree1);
 if(!hlttree1)
   {cout <<"error in the hlttree1"<<endl;}
      
      TDirectory * dirt = (TDirectory*)f->Get("ak4PFJetAnalyzer");
      dirt->GetObject("t",tr);
 if(!tr)
   {cout <<"error in the t"<<endl;}

      tree -> AddFriend(hitree);
      tree -> AddFriend(hlttree);
      tree -> AddFriend(hlttree1);
      tree ->AddFriend(tr);


   Init(tree);
}

myTree::~myTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t myTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t myTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void myTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Reco_QQ_4mom = 0;
   Reco_QQ_mupl_4mom = 0;
   Reco_QQ_mumi_4mom = 0;
   Reco_QQ_vtx = 0;
   Reco_mu_4mom = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
   fChain->SetBranchAddress("runNb", &runNb, &b_runNb);
   fChain->SetBranchAddress("LS", &LS, &b_LS);
   fChain->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
   fChain->SetBranchAddress("nTrig", &nTrig, &b_nTrig);
   fChain->SetBranchAddress("trigPrescale", trigPrescale, &b_trigPrescale);
   fChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
   fChain->SetBranchAddress("Npix", &Npix, &b_Npix);
   fChain->SetBranchAddress("NpixelTracks", &NpixelTracks, &b_NpixelTracks);
   fChain->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
   fChain->SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
   fChain->SetBranchAddress("SumET_HFplus", &SumET_HFplus, &b_SumET_HFplus);
   fChain->SetBranchAddress("SumET_HFminus", &SumET_HFminus, &b_SumET_HFminus);
   fChain->SetBranchAddress("SumET_HFplusEta4", &SumET_HFplusEta4, &b_SumET_HFplusEta4);
   fChain->SetBranchAddress("SumET_HFminusEta4", &SumET_HFminusEta4, &b_SumET_HFminusEta4);
   fChain->SetBranchAddress("SumET_ET", &SumET_ET, &b_SumET_ET);
   fChain->SetBranchAddress("SumET_EE", &SumET_EE, &b_SumET_EE);
   fChain->SetBranchAddress("SumET_EB", &SumET_EB, &b_SumET_EB);
   fChain->SetBranchAddress("SumET_EEplus", &SumET_EEplus, &b_SumET_EEplus);
   fChain->SetBranchAddress("SumET_EEminus", &SumET_EEminus, &b_SumET_EEminus);
   fChain->SetBranchAddress("SumET_ZDC", &SumET_ZDC, &b_SumET_ZDC);
   fChain->SetBranchAddress("SumET_ZDCplus", &SumET_ZDCplus, &b_SumET_ZDCplus);
   fChain->SetBranchAddress("SumET_ZDCminus", &SumET_ZDCminus, &b_SumET_ZDCminus);
   fChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
   fChain->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
   fChain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
   fChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
   fChain->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
   fChain->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
   fChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
   fChain->SetBranchAddress("Reco_QQ_mupl_trig", Reco_QQ_mupl_trig, &b_Reco_QQ_mupl_trig);
   fChain->SetBranchAddress("Reco_QQ_mumi_trig", Reco_QQ_mumi_trig, &b_Reco_QQ_mumi_trig);
   fChain->SetBranchAddress("Reco_QQ_isCowboy", Reco_QQ_isCowboy, &b_Reco_QQ_isCowboy);
   fChain->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau, &b_Reco_QQ_ctau);
   fChain->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr, &b_Reco_QQ_ctauErr);
   fChain->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
   fChain->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);
   fChain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
   fChain->SetBranchAddress("Reco_QQ_dca", Reco_QQ_dca, &b_Reco_QQ_dca);
   fChain->SetBranchAddress("Reco_QQ_MassErr", Reco_QQ_MassErr, &b_Reco_QQ_MassErr);
   fChain->SetBranchAddress("Reco_QQ_vtx", &Reco_QQ_vtx, &b_Reco_QQ_vtx);
   fChain->SetBranchAddress("Reco_QQ_Ntrk", Reco_QQ_Ntrk, &b_Reco_QQ_Ntrk);
   fChain->SetBranchAddress("Reco_QQ_mupl_SelectionType", Reco_QQ_mupl_SelectionType, &b_Reco_QQ_mupl_SelectionType);
   fChain->SetBranchAddress("Reco_QQ_mumi_SelectionType", Reco_QQ_mumi_SelectionType, &b_Reco_QQ_mumi_SelectionType);
   fChain->SetBranchAddress("Reco_QQ_mupl_isGoodMuon", Reco_QQ_mupl_isGoodMuon, &b_Reco_QQ_mupl_isGoodMuon);
   fChain->SetBranchAddress("Reco_QQ_mumi_isGoodMuon", Reco_QQ_mumi_isGoodMuon, &b_Reco_QQ_mumi_isGoodMuon);
   fChain->SetBranchAddress("Reco_QQ_mupl_highPurity", Reco_QQ_mupl_highPurity, &b_Reco_QQ_mupl_highPurity);
   fChain->SetBranchAddress("Reco_QQ_mumi_highPurity", Reco_QQ_mumi_highPurity, &b_Reco_QQ_mumi_highPurity);
   fChain->SetBranchAddress("Reco_QQ_mupl_TrkMuArb", Reco_QQ_mupl_TrkMuArb, &b_Reco_QQ_mupl_TrkMuArb);
   fChain->SetBranchAddress("Reco_QQ_mumi_TrkMuArb", Reco_QQ_mumi_TrkMuArb, &b_Reco_QQ_mumi_TrkMuArb);
   fChain->SetBranchAddress("Reco_QQ_mupl_TMOneStaTight", Reco_QQ_mupl_TMOneStaTight, &b_Reco_QQ_mupl_TMOneStaTight);
   fChain->SetBranchAddress("Reco_QQ_mumi_TMOneStaTight", Reco_QQ_mumi_TMOneStaTight, &b_Reco_QQ_mumi_TMOneStaTight);
   fChain->SetBranchAddress("Reco_QQ_mupl_nPixValHits", Reco_QQ_mupl_nPixValHits, &b_Reco_QQ_mupl_nPixValHits);
   fChain->SetBranchAddress("Reco_QQ_mumi_nPixValHits", Reco_QQ_mumi_nPixValHits, &b_Reco_QQ_mumi_nPixValHits);
   fChain->SetBranchAddress("Reco_QQ_mupl_nMuValHits", Reco_QQ_mupl_nMuValHits, &b_Reco_QQ_mupl_nMuValHits);
   fChain->SetBranchAddress("Reco_QQ_mumi_nMuValHits", Reco_QQ_mumi_nMuValHits, &b_Reco_QQ_mumi_nMuValHits);
   fChain->SetBranchAddress("Reco_QQ_mupl_nTrkHits", Reco_QQ_mupl_nTrkHits, &b_Reco_QQ_mupl_nTrkHits);
   fChain->SetBranchAddress("Reco_QQ_mumi_nTrkHits", Reco_QQ_mumi_nTrkHits, &b_Reco_QQ_mumi_nTrkHits);
   fChain->SetBranchAddress("Reco_QQ_mupl_normChi2_inner", Reco_QQ_mupl_normChi2_inner, &b_Reco_QQ_mupl_normChi2_inner);
   fChain->SetBranchAddress("Reco_QQ_mumi_normChi2_inner", Reco_QQ_mumi_normChi2_inner, &b_Reco_QQ_mumi_normChi2_inner);
   fChain->SetBranchAddress("Reco_QQ_mupl_normChi2_global", Reco_QQ_mupl_normChi2_global, &b_Reco_QQ_mupl_normChi2_global);
   fChain->SetBranchAddress("Reco_QQ_mumi_normChi2_global", Reco_QQ_mumi_normChi2_global, &b_Reco_QQ_mumi_normChi2_global);
   fChain->SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea, &b_Reco_QQ_mupl_nPixWMea);
   fChain->SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea, &b_Reco_QQ_mumi_nPixWMea);
   fChain->SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea, &b_Reco_QQ_mupl_nTrkWMea);
   fChain->SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea, &b_Reco_QQ_mumi_nTrkWMea);
   fChain->SetBranchAddress("Reco_QQ_mupl_StationsMatched", Reco_QQ_mupl_StationsMatched, &b_Reco_QQ_mupl_StationsMatched);
   fChain->SetBranchAddress("Reco_QQ_mumi_StationsMatched", Reco_QQ_mumi_StationsMatched, &b_Reco_QQ_mumi_StationsMatched);
   fChain->SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
   fChain->SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
   fChain->SetBranchAddress("Reco_QQ_mupl_dxyErr", Reco_QQ_mupl_dxyErr, &b_Reco_QQ_mupl_dxyErr);
   fChain->SetBranchAddress("Reco_QQ_mumi_dxyErr", Reco_QQ_mumi_dxyErr, &b_Reco_QQ_mumi_dxyErr);
   fChain->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
   fChain->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);
   fChain->SetBranchAddress("Reco_QQ_mupl_dzErr", Reco_QQ_mupl_dzErr, &b_Reco_QQ_mupl_dzErr);
   fChain->SetBranchAddress("Reco_QQ_mumi_dzErr", Reco_QQ_mumi_dzErr, &b_Reco_QQ_mumi_dzErr);
   fChain->SetBranchAddress("Reco_QQ_mupl_pt_inner", Reco_QQ_mupl_pt_inner, &b_Reco_QQ_mupl_pt_inner);
   fChain->SetBranchAddress("Reco_QQ_mumi_pt_inner", Reco_QQ_mumi_pt_inner, &b_Reco_QQ_mumi_pt_inner);
   fChain->SetBranchAddress("Reco_QQ_mupl_pt_global", Reco_QQ_mupl_pt_global, &b_Reco_QQ_mupl_pt_global);
   fChain->SetBranchAddress("Reco_QQ_mumi_pt_global", Reco_QQ_mumi_pt_global, &b_Reco_QQ_mumi_pt_global);
   fChain->SetBranchAddress("Reco_QQ_mupl_ptErr_inner", Reco_QQ_mupl_ptErr_inner, &b_Reco_QQ_mupl_ptErr_inner);
   fChain->SetBranchAddress("Reco_QQ_mumi_ptErr_inner", Reco_QQ_mumi_ptErr_inner, &b_Reco_QQ_mumi_ptErr_inner);
   fChain->SetBranchAddress("Reco_QQ_mupl_ptErr_global", Reco_QQ_mupl_ptErr_global, &b_Reco_QQ_mupl_ptErr_global);
   fChain->SetBranchAddress("Reco_QQ_mumi_ptErr_global", Reco_QQ_mumi_ptErr_global, &b_Reco_QQ_mumi_ptErr_global);
   fChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
   fChain->SetBranchAddress("Reco_mu_type", Reco_mu_type, &b_Reco_mu_type);
   fChain->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);
   fChain->SetBranchAddress("Reco_mu_charge", Reco_mu_charge, &b_Reco_mu_charge);
   fChain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
   fChain->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
   fChain->SetBranchAddress("Reco_mu_isGoodMuon", Reco_mu_isGoodMuon, &b_Reco_mu_isGoodMuon);
   fChain->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
   fChain->SetBranchAddress("Reco_mu_TrkMuArb", Reco_mu_TrkMuArb, &b_Reco_mu_TrkMuArb);
   fChain->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
   fChain->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
   fChain->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
   fChain->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
   fChain->SetBranchAddress("Reco_mu_normChi2_inner", Reco_mu_normChi2_inner, &b_Reco_mu_normChi2_inner);
   fChain->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
   fChain->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
   fChain->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
   fChain->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
   fChain->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
   fChain->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
   fChain->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
   fChain->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
   fChain->SetBranchAddress("Reco_mu_pt_inner", Reco_mu_pt_inner, &b_Reco_mu_pt_inner);
   fChain->SetBranchAddress("Reco_mu_pt_global", Reco_mu_pt_global, &b_Reco_mu_pt_global);
   fChain->SetBranchAddress("Reco_mu_ptErr_inner", Reco_mu_ptErr_inner, &b_Reco_mu_ptErr_inner);
   fChain->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
   fChain->SetBranchAddress("hiHF", &hiHF, &b_hiHF);
   fChain->SetBranchAddress("hiHFplus", &hiHFplus, &b_hiHFplus);
   fChain->SetBranchAddress("hiHFminus", &hiHFminus, &b_hiHFminus);
   fChain->SetBranchAddress("hiHFplusEta4", &hiHFplusEta4, &b_hiHFplusEta4);
   fChain->SetBranchAddress("hiHFminusEta4", &hiHFminusEta4, &b_hiHFminusEta4);
   fChain->SetBranchAddress("hiZDC", &hiZDC, &b_hiZDC);
   fChain->SetBranchAddress("hiZDCplus", &hiZDCplus, &b_hiZDCplus);
   fChain->SetBranchAddress("hiZDCminus", &hiZDCminus, &b_hiZDCminus);
   fChain->SetBranchAddress("hiHFhit", &hiHFhit, &b_hiHFhit);
   fChain->SetBranchAddress("hiHFhitPlus", &hiHFhitPlus, &b_hiHFhitPlus);
   fChain->SetBranchAddress("hiHFhitMinus", &hiHFhitMinus, &b_hiHFhitMinus);
   fChain->SetBranchAddress("hiET", &hiET, &b_hiET);
   fChain->SetBranchAddress("hiEE", &hiEE, &b_hiEE);
   fChain->SetBranchAddress("hiEB", &hiEB, &b_hiEB);
   fChain->SetBranchAddress("hiEEplus", &hiEEplus, &b_hiEEplus);
   fChain->SetBranchAddress("hiEEminus", &hiEEminus, &b_hiEEminus);
   fChain->SetBranchAddress("hiNpix", &hiNpix, &b_hiNpix);
   fChain->SetBranchAddress("hiNpixelTracks", &hiNpixelTracks, &b_hiNpixelTracks);
   fChain->SetBranchAddress("hiNtracks", &hiNtracks, &b_hiNtracks);
   fChain->SetBranchAddress("hiNtracksPtCut", &hiNtracksPtCut, &b_hiNtracksPtCut);
   fChain->SetBranchAddress("hiNtracksEtaCut", &hiNtracksEtaCut, &b_hiNtracksEtaCut);
   fChain->SetBranchAddress("hiNtracksEtaPtCut", &hiNtracksEtaPtCut, &b_hiNtracksEtaPtCut);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
   fChain->SetBranchAddress("hiHF", &hiHF, &b_hiHF);
   fChain->SetBranchAddress("hiHFplus", &hiHFplus, &b_hiHFplus);
   fChain->SetBranchAddress("hiHFminus", &hiHFminus, &b_hiHFminus);
   fChain->SetBranchAddress("hiHFplusEta4", &hiHFplusEta4, &b_hiHFplusEta4);
   fChain->SetBranchAddress("hiHFminusEta4", &hiHFminusEta4, &b_hiHFminusEta4);
   fChain->SetBranchAddress("hiZDC", &hiZDC, &b_hiZDC);
   fChain->SetBranchAddress("hiZDCplus", &hiZDCplus, &b_hiZDCplus);
   fChain->SetBranchAddress("hiZDCminus", &hiZDCminus, &b_hiZDCminus);
   fChain->SetBranchAddress("hiHFhit", &hiHFhit, &b_hiHFhit);
   fChain->SetBranchAddress("hiHFhitPlus", &hiHFhitPlus, &b_hiHFhitPlus);
   fChain->SetBranchAddress("hiHFhitMinus", &hiHFhitMinus, &b_hiHFhitMinus);
   fChain->SetBranchAddress("hiET", &hiET, &b_hiET);
   fChain->SetBranchAddress("hiEE", &hiEE, &b_hiEE);
   fChain->SetBranchAddress("hiEB", &hiEB, &b_hiEB);
   fChain->SetBranchAddress("hiEEplus", &hiEEplus, &b_hiEEplus);
   fChain->SetBranchAddress("hiEEminus", &hiEEminus, &b_hiEEminus);
   fChain->SetBranchAddress("hiNpix", &hiNpix, &b_hiNpix);
   fChain->SetBranchAddress("hiNpixelTracks", &hiNpixelTracks, &b_hiNpixelTracks);
   fChain->SetBranchAddress("hiNtracks", &hiNtracks, &b_hiNtracks);
   fChain->SetBranchAddress("hiNtracksPtCut", &hiNtracksPtCut, &b_hiNtracksPtCut);
   fChain->SetBranchAddress("hiNtracksEtaCut", &hiNtracksEtaCut, &b_hiNtracksEtaCut);
   fChain->SetBranchAddress("hiNtracksEtaPtCut", &hiNtracksEtaPtCut, &b_hiNtracksEtaPtCut);

   fChain->SetBranchAddress("Onia2MuMuPAT", &Onia2MuMuPAT, &b_Onia2MuMuPAT);
   fChain->SetBranchAddress("ana_step", &ana_step, &b_ana_step);
   fChain->SetBranchAddress("pHBHENoiseFilterResultProducer", &pHBHENoiseFilterResultProducer, &b_pHBHENoiseFilterResultProducer);
   fChain->SetBranchAddress("HBHENoiseFilterResult", &HBHENoiseFilterResult, &b_HBHENoiseFilterResult);
   fChain->SetBranchAddress("HBHENoiseFilterResultRun1", &HBHENoiseFilterResultRun1, &b_HBHENoiseFilterResultRun1);
   fChain->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose, &b_HBHENoiseFilterResultRun2Loose);
   fChain->SetBranchAddress("HBHENoiseFilterResultRun2Tight", &HBHENoiseFilterResultRun2Tight, &b_HBHENoiseFilterResultRun2Tight);
   fChain->SetBranchAddress("HBHEIsoNoiseFilterResult", &HBHEIsoNoiseFilterResult, &b_HBHEIsoNoiseFilterResult);
   fChain->SetBranchAddress("pPAprimaryVertexFilter", &pPAprimaryVertexFilter, &b_pPAprimaryVertexFilter);
   fChain->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter, &b_pBeamScrapingFilter);
   fChain->SetBranchAddress("pVertexFilterCutG", &pVertexFilterCutG, &b_pVertexFilterCutG);
   fChain->SetBranchAddress("pVertexFilterCutGloose", &pVertexFilterCutGloose, &b_pVertexFilterCutGloose);
   fChain->SetBranchAddress("pVertexFilterCutGtight", &pVertexFilterCutGtight, &b_pVertexFilterCutGtight);
   fChain->SetBranchAddress("pVertexFilterCutGplus", &pVertexFilterCutGplus, &b_pVertexFilterCutGplus);
   fChain->SetBranchAddress("pVertexFilterCutE", &pVertexFilterCutE, &b_pVertexFilterCutE);
   fChain->SetBranchAddress("pVertexFilterCutEandG", &pVertexFilterCutEandG, &b_pVertexFilterCutEandG);

   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("b", &b, &b_b);
   fChain->SetBranchAddress("nref", &nref, &b_nref);
   fChain->SetBranchAddress("rawpt", rawpt, &b_rawpt);
   fChain->SetBranchAddress("jtpt", jtpt, &b_jtpt);
   fChain->SetBranchAddress("jteta", jteta, &b_jteta);
   fChain->SetBranchAddress("jty", jty, &b_jty);
   fChain->SetBranchAddress("jtphi", jtphi, &b_jtphi);
   fChain->SetBranchAddress("jtpu", jtpu, &b_jtpu);
   fChain->SetBranchAddress("jtm", jtm, &b_jtm);
   fChain->SetBranchAddress("jtarea", jtarea, &b_jtarea);
   fChain->SetBranchAddress("jtPfCHF", jtPfCHF, &b_jtPfCHF);
   fChain->SetBranchAddress("jtPfNHF", jtPfNHF, &b_jtPfNHF);
   fChain->SetBranchAddress("jtPfCEF", jtPfCEF, &b_jtPfCEF);
   fChain->SetBranchAddress("jtPfNEF", jtPfNEF, &b_jtPfNEF);
   fChain->SetBranchAddress("jtPfMUF", jtPfMUF, &b_jtPfMUF);
   fChain->SetBranchAddress("jtPfCHM", jtPfCHM, &b_jtPfCHM);
   fChain->SetBranchAddress("jtPfNHM", jtPfNHM, &b_jtPfNHM);
   fChain->SetBranchAddress("jtPfCEM", jtPfCEM, &b_jtPfCEM);
   fChain->SetBranchAddress("jtPfNEM", jtPfNEM, &b_jtPfNEM);
   fChain->SetBranchAddress("jtPfMUM", jtPfMUM, &b_jtPfMUM);
   fChain->SetBranchAddress("jttau1", jttau1, &b_jttau1);
   fChain->SetBranchAddress("jttau2", jttau2, &b_jttau2);
   fChain->SetBranchAddress("jttau3", jttau3, &b_jttau3);
   fChain->SetBranchAddress("discr_jetID_cuts", discr_jetID_cuts, &b_discr_jetID_cuts);
   fChain->SetBranchAddress("discr_jetID_bdt", discr_jetID_bdt, &b_discr_jetID_bdt);
   fChain->SetBranchAddress("discr_fr01", discr_fr01, &b_discr_fr01);
   fChain->SetBranchAddress("trackMax", trackMax, &b_trackMax);
   fChain->SetBranchAddress("trackSum", trackSum, &b_trackSum);
   fChain->SetBranchAddress("trackN", trackN, &b_trackN);
   fChain->SetBranchAddress("trackHardSum", trackHardSum, &b_trackHardSum);
   fChain->SetBranchAddress("trackHardN", trackHardN, &b_trackHardN);
   fChain->SetBranchAddress("chargedMax", chargedMax, &b_chargedMax);
   fChain->SetBranchAddress("chargedSum", chargedSum, &b_chargedSum);
   fChain->SetBranchAddress("chargedN", chargedN, &b_chargedN);
   fChain->SetBranchAddress("chargedHardSum", chargedHardSum, &b_chargedHardSum);
   fChain->SetBranchAddress("chargedHardN", chargedHardN, &b_chargedHardN);
   fChain->SetBranchAddress("photonMax", photonMax, &b_photonMax);
   fChain->SetBranchAddress("photonSum", photonSum, &b_photonSum);
   fChain->SetBranchAddress("photonN", photonN, &b_photonN);
   fChain->SetBranchAddress("photonHardSum", photonHardSum, &b_photonHardSum);
   fChain->SetBranchAddress("photonHardN", photonHardN, &b_photonHardN);
   fChain->SetBranchAddress("neutralMax", neutralMax, &b_neutralMax);
   fChain->SetBranchAddress("neutralSum", neutralSum, &b_neutralSum);
   fChain->SetBranchAddress("neutralN", neutralN, &b_neutralN);
   fChain->SetBranchAddress("hcalSum", hcalSum, &b_hcalSum);
   fChain->SetBranchAddress("ecalSum", ecalSum, &b_ecalSum);
   fChain->SetBranchAddress("eMax", eMax, &b_eMax);
   fChain->SetBranchAddress("eSum", eSum, &b_eSum);
   fChain->SetBranchAddress("eN", eN, &b_eN);
   fChain->SetBranchAddress("muMax", muMax, &b_muMax);
   fChain->SetBranchAddress("muSum", muSum, &b_muSum);
   fChain->SetBranchAddress("muN", muN, &b_muN);
   fChain->SetBranchAddress("discr_ssvHighEff", discr_ssvHighEff, &b_discr_ssvHighEff);
   fChain->SetBranchAddress("discr_ssvHighPur", discr_ssvHighPur, &b_discr_ssvHighPur);
   fChain->SetBranchAddress("discr_csvV1", discr_csvV1, &b_discr_csvV1);
   fChain->SetBranchAddress("discr_csvV2", discr_csvV2, &b_discr_csvV2);
   fChain->SetBranchAddress("discr_muByIp3", discr_muByIp3, &b_discr_muByIp3);
   fChain->SetBranchAddress("discr_muByPt", discr_muByPt, &b_discr_muByPt);
   fChain->SetBranchAddress("discr_prob", discr_prob, &b_discr_prob);
   fChain->SetBranchAddress("discr_probb", discr_probb, &b_discr_probb);
   fChain->SetBranchAddress("discr_tcHighEff", discr_tcHighEff, &b_discr_tcHighEff);
   fChain->SetBranchAddress("discr_tcHighPur", discr_tcHighPur, &b_discr_tcHighPur);
   fChain->SetBranchAddress("ndiscr_ssvHighEff", ndiscr_ssvHighEff, &b_ndiscr_ssvHighEff);
   fChain->SetBranchAddress("ndiscr_ssvHighPur", ndiscr_ssvHighPur, &b_ndiscr_ssvHighPur);
   fChain->SetBranchAddress("ndiscr_csvV1", ndiscr_csvV1, &b_ndiscr_csvV1);
   fChain->SetBranchAddress("ndiscr_csvV2", ndiscr_csvV2, &b_ndiscr_csvV2);
   fChain->SetBranchAddress("ndiscr_muByPt", ndiscr_muByPt, &b_ndiscr_muByPt);
   fChain->SetBranchAddress("pdiscr_csvV1", pdiscr_csvV1, &b_pdiscr_csvV1);
   fChain->SetBranchAddress("pdiscr_csvV2", pdiscr_csvV2, &b_pdiscr_csvV2);
   fChain->SetBranchAddress("nsvtx", nsvtx, &b_nsvtx);
   fChain->SetBranchAddress("svtxntrk", svtxntrk, &b_svtxntrk);
   fChain->SetBranchAddress("svtxdl", svtxdl, &b_svtxdl);
   fChain->SetBranchAddress("svtxdls", svtxdls, &b_svtxdls);
   fChain->SetBranchAddress("svtxdl2d", svtxdl2d, &b_svtxdl2d);
   fChain->SetBranchAddress("svtxdls2d", svtxdls2d, &b_svtxdls2d);
   fChain->SetBranchAddress("svtxm", svtxm, &b_svtxm);
   fChain->SetBranchAddress("svtxpt", svtxpt, &b_svtxpt);
   fChain->SetBranchAddress("svtxmcorr", svtxmcorr, &b_svtxmcorr);
   fChain->SetBranchAddress("nIPtrk", nIPtrk, &b_nIPtrk);
   fChain->SetBranchAddress("nselIPtrk", nselIPtrk, &b_nselIPtrk);
   fChain->SetBranchAddress("mue", mue, &b_mue);
   fChain->SetBranchAddress("mupt", mupt, &b_mupt);
   fChain->SetBranchAddress("mueta", mueta, &b_mueta);
   fChain->SetBranchAddress("muphi", muphi, &b_muphi);
   fChain->SetBranchAddress("mudr", mudr, &b_mudr);
   fChain->SetBranchAddress("muptrel", muptrel, &b_muptrel);
   fChain->SetBranchAddress("muchg", muchg, &b_muchg);
   fChain->SetBranchAddress("NL1IsolEm", &NL1IsolEm, &b_NL1IsolEm);
   fChain->SetBranchAddress("L1IsolEmEt", L1IsolEmEt, &b_L1IsolEmEt);
   fChain->SetBranchAddress("L1IsolEmE", L1IsolEmE, &b_L1IsolEmE);
   fChain->SetBranchAddress("L1IsolEmEta", L1IsolEmEta, &b_L1IsolEmEta);
   fChain->SetBranchAddress("L1IsolEmPhi", L1IsolEmPhi, &b_L1IsolEmPhi);
   fChain->SetBranchAddress("NL1NIsolEm", &NL1NIsolEm, &b_NL1NIsolEm);
   fChain->SetBranchAddress("L1NIsolEmEt", L1NIsolEmEt, &b_L1NIsolEmEt);
   fChain->SetBranchAddress("L1NIsolEmE", L1NIsolEmE, &b_L1NIsolEmE);
   fChain->SetBranchAddress("L1NIsolEmEta", L1NIsolEmEta, &b_L1NIsolEmEta);
   fChain->SetBranchAddress("L1NIsolEmPhi", L1NIsolEmPhi, &b_L1NIsolEmPhi);
   fChain->SetBranchAddress("NL1Mu", &NL1Mu, &b_NL1Mu);
   fChain->SetBranchAddress("L1MuPt", L1MuPt, &b_L1MuPt);
   fChain->SetBranchAddress("L1MuE", L1MuE, &b_L1MuE);
   fChain->SetBranchAddress("L1MuEta", L1MuEta, &b_L1MuEta);
   fChain->SetBranchAddress("L1MuPhi", L1MuPhi, &b_L1MuPhi);
   fChain->SetBranchAddress("L1MuIsol", L1MuIsol, &b_L1MuIsol);
   fChain->SetBranchAddress("L1MuMip", L1MuMip, &b_L1MuMip);
   fChain->SetBranchAddress("L1MuFor", L1MuFor, &b_L1MuFor);
   fChain->SetBranchAddress("L1MuRPC", L1MuRPC, &b_L1MuRPC);
   fChain->SetBranchAddress("L1MuQal", L1MuQal, &b_L1MuQal);
   fChain->SetBranchAddress("L1MuChg", L1MuChg, &b_L1MuChg);
   fChain->SetBranchAddress("NL1CenJet", &NL1CenJet, &b_NL1CenJet);
   fChain->SetBranchAddress("L1CenJetEt", L1CenJetEt, &b_L1CenJetEt);
   fChain->SetBranchAddress("L1CenJetE", L1CenJetE, &b_L1CenJetE);
   fChain->SetBranchAddress("L1CenJetEta", L1CenJetEta, &b_L1CenJetEta);
   fChain->SetBranchAddress("L1CenJetPhi", L1CenJetPhi, &b_L1CenJetPhi);
   fChain->SetBranchAddress("NL1ForJet", &NL1ForJet, &b_NL1ForJet);
   fChain->SetBranchAddress("L1ForJetEt", L1ForJetEt, &b_L1ForJetEt);
   fChain->SetBranchAddress("L1ForJetE", L1ForJetE, &b_L1ForJetE);
   fChain->SetBranchAddress("L1ForJetEta", L1ForJetEta, &b_L1ForJetEta);
   fChain->SetBranchAddress("L1ForJetPhi", L1ForJetPhi, &b_L1ForJetPhi);
   fChain->SetBranchAddress("NL1Tau", &NL1Tau, &b_NL1Tau);
   fChain->SetBranchAddress("L1TauEt", L1TauEt, &b_L1TauEt);
   fChain->SetBranchAddress("L1TauE", L1TauE, &b_L1TauE);
   fChain->SetBranchAddress("L1TauEta", L1TauEta, &b_L1TauEta);
   fChain->SetBranchAddress("L1TauPhi", L1TauPhi, &b_L1TauPhi);
   fChain->SetBranchAddress("L1Met", &L1Met, &b_L1Met);
   fChain->SetBranchAddress("L1MetPhi", &L1MetPhi, &b_L1MetPhi);
   fChain->SetBranchAddress("L1EtTot", &L1EtTot, &b_L1EtTot);
   fChain->SetBranchAddress("L1Mht", &L1Mht, &b_L1Mht);
   fChain->SetBranchAddress("L1MhtPhi", &L1MhtPhi, &b_L1MhtPhi);
   fChain->SetBranchAddress("L1EtHad", &L1EtHad, &b_L1EtHad);
   fChain->SetBranchAddress("L1HfRing1EtSumPositiveEta", &L1HfRing1EtSumPositiveEta, &b_L1HfRing1EtSumPositiveEta);
   fChain->SetBranchAddress("L1HfRing2EtSumPositiveEta", &L1HfRing2EtSumPositiveEta, &b_L1HfRing2EtSumPositiveEta);
   fChain->SetBranchAddress("L1HfRing1EtSumNegativeEta", &L1HfRing1EtSumNegativeEta, &b_L1HfRing1EtSumNegativeEta);
   fChain->SetBranchAddress("L1HfRing2EtSumNegativeEta", &L1HfRing2EtSumNegativeEta, &b_L1HfRing2EtSumNegativeEta);
   fChain->SetBranchAddress("L1HfTowerCountPositiveEtaRing1", &L1HfTowerCountPositiveEtaRing1, &b_L1HfTowerCountPositiveEtaRing1);
   fChain->SetBranchAddress("L1HfTowerCountNegativeEtaRing1", &L1HfTowerCountNegativeEtaRing1, &b_L1HfTowerCountNegativeEtaRing1);
   fChain->SetBranchAddress("L1HfTowerCountPositiveEtaRing2", &L1HfTowerCountPositiveEtaRing2, &b_L1HfTowerCountPositiveEtaRing2);
   fChain->SetBranchAddress("L1HfTowerCountNegativeEtaRing2", &L1HfTowerCountNegativeEtaRing2, &b_L1HfTowerCountNegativeEtaRing2);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Bx", &Bx, &b_Bx);
   fChain->SetBranchAddress("Orbit", &Orbit, &b_Orbit);
   fChain->SetBranchAddress("AvgInstDelLumi", &AvgInstDelLumi, &b_AvgInstDelLumi);
   fChain->SetBranchAddress("HLTriggerFirstPath", &HLTriggerFirstPath, &b_HLTriggerFirstPath);
   fChain->SetBranchAddress("HLTriggerFirstPath_Prescl", &HLTriggerFirstPath_Prescl, &b_HLTriggerFirstPath_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity60_v2", &HLT_PixelTracks_Multiplicity60_v2, &b_HLT_PixelTracks_Multiplicity60_v2);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity60_v2_Prescl", &HLT_PixelTracks_Multiplicity60_v2_Prescl, &b_HLT_PixelTracks_Multiplicity60_v2_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity85_v2", &HLT_PixelTracks_Multiplicity85_v2, &b_HLT_PixelTracks_Multiplicity85_v2);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity85_v2_Prescl", &HLT_PixelTracks_Multiplicity85_v2_Prescl, &b_HLT_PixelTracks_Multiplicity85_v2_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity110_v2", &HLT_PixelTracks_Multiplicity110_v2, &b_HLT_PixelTracks_Multiplicity110_v2);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity110_v2_Prescl", &HLT_PixelTracks_Multiplicity110_v2_Prescl, &b_HLT_PixelTracks_Multiplicity110_v2_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity135_v2", &HLT_PixelTracks_Multiplicity135_v2, &b_HLT_PixelTracks_Multiplicity135_v2);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity135_v2_Prescl", &HLT_PixelTracks_Multiplicity135_v2_Prescl, &b_HLT_PixelTracks_Multiplicity135_v2_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity160_v2", &HLT_PixelTracks_Multiplicity160_v2, &b_HLT_PixelTracks_Multiplicity160_v2);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity160_v2_Prescl", &HLT_PixelTracks_Multiplicity160_v2_Prescl, &b_HLT_PixelTracks_Multiplicity160_v2_Prescl);
   fChain->SetBranchAddress("HLT_Physics_v2", &HLT_Physics_v2, &b_HLT_Physics_v2);
   fChain->SetBranchAddress("HLT_Physics_v2_Prescl", &HLT_Physics_v2_Prescl, &b_HLT_Physics_v2_Prescl);
   fChain->SetBranchAddress("DST_Physics_v1", &DST_Physics_v1, &b_DST_Physics_v1);
   fChain->SetBranchAddress("DST_Physics_v1_Prescl", &DST_Physics_v1_Prescl, &b_DST_Physics_v1_Prescl);
   fChain->SetBranchAddress("HLT_Random_v1", &HLT_Random_v1, &b_HLT_Random_v1);
   fChain->SetBranchAddress("HLT_Random_v1_Prescl", &HLT_Random_v1_Prescl, &b_HLT_Random_v1_Prescl);
   fChain->SetBranchAddress("HLT_EcalCalibration_v1", &HLT_EcalCalibration_v1, &b_HLT_EcalCalibration_v1);
   fChain->SetBranchAddress("HLT_EcalCalibration_v1_Prescl", &HLT_EcalCalibration_v1_Prescl, &b_HLT_EcalCalibration_v1_Prescl);
   fChain->SetBranchAddress("HLT_HcalCalibration_v1", &HLT_HcalCalibration_v1, &b_HLT_HcalCalibration_v1);
   fChain->SetBranchAddress("HLT_HcalCalibration_v1_Prescl", &HLT_HcalCalibration_v1_Prescl, &b_HLT_HcalCalibration_v1_Prescl);
   fChain->SetBranchAddress("AlCa_EcalPhiSym_v3", &AlCa_EcalPhiSym_v3, &b_AlCa_EcalPhiSym_v3);
   fChain->SetBranchAddress("AlCa_EcalPhiSym_v3_Prescl", &AlCa_EcalPhiSym_v3_Prescl, &b_AlCa_EcalPhiSym_v3_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech6_BPTX_MinusOnly_v1", &HLT_L1Tech6_BPTX_MinusOnly_v1, &b_HLT_L1Tech6_BPTX_MinusOnly_v1);
   fChain->SetBranchAddress("HLT_L1Tech6_BPTX_MinusOnly_v1_Prescl", &HLT_L1Tech6_BPTX_MinusOnly_v1_Prescl, &b_HLT_L1Tech6_BPTX_MinusOnly_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech5_BPTX_PlusOnly_v2", &HLT_L1Tech5_BPTX_PlusOnly_v2, &b_HLT_L1Tech5_BPTX_PlusOnly_v2);
   fChain->SetBranchAddress("HLT_L1Tech5_BPTX_PlusOnly_v2_Prescl", &HLT_L1Tech5_BPTX_PlusOnly_v2_Prescl, &b_HLT_L1Tech5_BPTX_PlusOnly_v2_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech7_NoBPTX_v1", &HLT_L1Tech7_NoBPTX_v1, &b_HLT_L1Tech7_NoBPTX_v1);
   fChain->SetBranchAddress("HLT_L1Tech7_NoBPTX_v1_Prescl", &HLT_L1Tech7_NoBPTX_v1_Prescl, &b_HLT_L1Tech7_NoBPTX_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1TOTEM1_MinBias_v1", &HLT_L1TOTEM1_MinBias_v1, &b_HLT_L1TOTEM1_MinBias_v1);
   fChain->SetBranchAddress("HLT_L1TOTEM1_MinBias_v1_Prescl", &HLT_L1TOTEM1_MinBias_v1_Prescl, &b_HLT_L1TOTEM1_MinBias_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1TOTEM2_ZeroBias_v1", &HLT_L1TOTEM2_ZeroBias_v1, &b_HLT_L1TOTEM2_ZeroBias_v1);
   fChain->SetBranchAddress("HLT_L1TOTEM2_ZeroBias_v1_Prescl", &HLT_L1TOTEM2_ZeroBias_v1_Prescl, &b_HLT_L1TOTEM2_ZeroBias_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_v1", &HLT_L1MinimumBiasHF1OR_v1, &b_HLT_L1MinimumBiasHF1OR_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_v1_Prescl", &HLT_L1MinimumBiasHF1OR_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF2OR_v1", &HLT_L1MinimumBiasHF2OR_v1, &b_HLT_L1MinimumBiasHF2OR_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF2OR_v1_Prescl", &HLT_L1MinimumBiasHF2OR_v1_Prescl, &b_HLT_L1MinimumBiasHF2OR_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1AND_v1", &HLT_L1MinimumBiasHF1AND_v1, &b_HLT_L1MinimumBiasHF1AND_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1AND_v1_Prescl", &HLT_L1MinimumBiasHF1AND_v1_Prescl, &b_HLT_L1MinimumBiasHF1AND_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF2AND_v1", &HLT_L1MinimumBiasHF2AND_v1, &b_HLT_L1MinimumBiasHF2AND_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF2AND_v1_Prescl", &HLT_L1MinimumBiasHF2AND_v1_Prescl, &b_HLT_L1MinimumBiasHF2AND_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF2ORNoBptxGating_v1", &HLT_L1MinimumBiasHF2ORNoBptxGating_v1, &b_HLT_L1MinimumBiasHF2ORNoBptxGating_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF2ORNoBptxGating_v1_Prescl", &HLT_L1MinimumBiasHF2ORNoBptxGating_v1_Prescl, &b_HLT_L1MinimumBiasHF2ORNoBptxGating_v1_Prescl);
   fChain->SetBranchAddress("AlCa_RPCMuonNoTriggers_v2", &AlCa_RPCMuonNoTriggers_v2, &b_AlCa_RPCMuonNoTriggers_v2);
   fChain->SetBranchAddress("AlCa_RPCMuonNoTriggers_v2_Prescl", &AlCa_RPCMuonNoTriggers_v2_Prescl, &b_AlCa_RPCMuonNoTriggers_v2_Prescl);
   fChain->SetBranchAddress("AlCa_RPCMuonNoHits_v2", &AlCa_RPCMuonNoHits_v2, &b_AlCa_RPCMuonNoHits_v2);
   fChain->SetBranchAddress("AlCa_RPCMuonNoHits_v2_Prescl", &AlCa_RPCMuonNoHits_v2_Prescl, &b_AlCa_RPCMuonNoHits_v2_Prescl);
   fChain->SetBranchAddress("AlCa_RPCMuonNormalisation_v2", &AlCa_RPCMuonNormalisation_v2, &b_AlCa_RPCMuonNormalisation_v2);
   fChain->SetBranchAddress("AlCa_RPCMuonNormalisation_v2_Prescl", &AlCa_RPCMuonNormalisation_v2_Prescl, &b_AlCa_RPCMuonNormalisation_v2_Prescl);
   fChain->SetBranchAddress("AlCa_LumiPixels_Random_v1", &AlCa_LumiPixels_Random_v1, &b_AlCa_LumiPixels_Random_v1);
   fChain->SetBranchAddress("AlCa_LumiPixels_Random_v1_Prescl", &AlCa_LumiPixels_Random_v1_Prescl, &b_AlCa_LumiPixels_Random_v1_Prescl);
   fChain->SetBranchAddress("AlCa_LumiPixels_ZeroBias_v2", &AlCa_LumiPixels_ZeroBias_v2, &b_AlCa_LumiPixels_ZeroBias_v2);
   fChain->SetBranchAddress("AlCa_LumiPixels_ZeroBias_v2_Prescl", &AlCa_LumiPixels_ZeroBias_v2_Prescl, &b_AlCa_LumiPixels_ZeroBias_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_v2", &HLT_ZeroBias_v2, &b_HLT_ZeroBias_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_v2_Prescl", &HLT_ZeroBias_v2_Prescl, &b_HLT_ZeroBias_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part0_v2", &HLT_ZeroBias_part0_v2, &b_HLT_ZeroBias_part0_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part0_v2_Prescl", &HLT_ZeroBias_part0_v2_Prescl, &b_HLT_ZeroBias_part0_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part1_v2", &HLT_ZeroBias_part1_v2, &b_HLT_ZeroBias_part1_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part1_v2_Prescl", &HLT_ZeroBias_part1_v2_Prescl, &b_HLT_ZeroBias_part1_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part2_v2", &HLT_ZeroBias_part2_v2, &b_HLT_ZeroBias_part2_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part2_v2_Prescl", &HLT_ZeroBias_part2_v2_Prescl, &b_HLT_ZeroBias_part2_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part3_v2", &HLT_ZeroBias_part3_v2, &b_HLT_ZeroBias_part3_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part3_v2_Prescl", &HLT_ZeroBias_part3_v2_Prescl, &b_HLT_ZeroBias_part3_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part4_v2", &HLT_ZeroBias_part4_v2, &b_HLT_ZeroBias_part4_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part4_v2_Prescl", &HLT_ZeroBias_part4_v2_Prescl, &b_HLT_ZeroBias_part4_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part5_v2", &HLT_ZeroBias_part5_v2, &b_HLT_ZeroBias_part5_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part5_v2_Prescl", &HLT_ZeroBias_part5_v2_Prescl, &b_HLT_ZeroBias_part5_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part6_v2", &HLT_ZeroBias_part6_v2, &b_HLT_ZeroBias_part6_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part6_v2_Prescl", &HLT_ZeroBias_part6_v2_Prescl, &b_HLT_ZeroBias_part6_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part7_v2", &HLT_ZeroBias_part7_v2, &b_HLT_ZeroBias_part7_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part7_v2_Prescl", &HLT_ZeroBias_part7_v2_Prescl, &b_HLT_ZeroBias_part7_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part8_v2", &HLT_ZeroBias_part8_v2, &b_HLT_ZeroBias_part8_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part8_v2_Prescl", &HLT_ZeroBias_part8_v2_Prescl, &b_HLT_ZeroBias_part8_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part9_v2", &HLT_ZeroBias_part9_v2, &b_HLT_ZeroBias_part9_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part9_v2_Prescl", &HLT_ZeroBias_part9_v2_Prescl, &b_HLT_ZeroBias_part9_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part10_v2", &HLT_ZeroBias_part10_v2, &b_HLT_ZeroBias_part10_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part10_v2_Prescl", &HLT_ZeroBias_part10_v2_Prescl, &b_HLT_ZeroBias_part10_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part11_v2", &HLT_ZeroBias_part11_v2, &b_HLT_ZeroBias_part11_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part11_v2_Prescl", &HLT_ZeroBias_part11_v2_Prescl, &b_HLT_ZeroBias_part11_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part12_v2", &HLT_ZeroBias_part12_v2, &b_HLT_ZeroBias_part12_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part12_v2_Prescl", &HLT_ZeroBias_part12_v2_Prescl, &b_HLT_ZeroBias_part12_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part13_v2", &HLT_ZeroBias_part13_v2, &b_HLT_ZeroBias_part13_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part13_v2_Prescl", &HLT_ZeroBias_part13_v2_Prescl, &b_HLT_ZeroBias_part13_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part14_v2", &HLT_ZeroBias_part14_v2, &b_HLT_ZeroBias_part14_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part14_v2_Prescl", &HLT_ZeroBias_part14_v2_Prescl, &b_HLT_ZeroBias_part14_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part15_v2", &HLT_ZeroBias_part15_v2, &b_HLT_ZeroBias_part15_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part15_v2_Prescl", &HLT_ZeroBias_part15_v2_Prescl, &b_HLT_ZeroBias_part15_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part16_v2", &HLT_ZeroBias_part16_v2, &b_HLT_ZeroBias_part16_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part16_v2_Prescl", &HLT_ZeroBias_part16_v2_Prescl, &b_HLT_ZeroBias_part16_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part17_v2", &HLT_ZeroBias_part17_v2, &b_HLT_ZeroBias_part17_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part17_v2_Prescl", &HLT_ZeroBias_part17_v2_Prescl, &b_HLT_ZeroBias_part17_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part18_v2", &HLT_ZeroBias_part18_v2, &b_HLT_ZeroBias_part18_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part18_v2_Prescl", &HLT_ZeroBias_part18_v2_Prescl, &b_HLT_ZeroBias_part18_v2_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part19_v2", &HLT_ZeroBias_part19_v2, &b_HLT_ZeroBias_part19_v2);
   fChain->SetBranchAddress("HLT_ZeroBias_part19_v2_Prescl", &HLT_ZeroBias_part19_v2_Prescl, &b_HLT_ZeroBias_part19_v2_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet40_Eta5p1_v1", &HLT_AK4CaloJet40_Eta5p1_v1, &b_HLT_AK4CaloJet40_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet40_Eta5p1_v1_Prescl", &HLT_AK4CaloJet40_Eta5p1_v1_Prescl, &b_HLT_AK4CaloJet40_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet60_Eta5p1_v1", &HLT_AK4CaloJet60_Eta5p1_v1, &b_HLT_AK4CaloJet60_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet60_Eta5p1_v1_Prescl", &HLT_AK4CaloJet60_Eta5p1_v1_Prescl, &b_HLT_AK4CaloJet60_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1", &HLT_AK4CaloJet80_Eta5p1_v1, &b_HLT_AK4CaloJet80_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1_Prescl", &HLT_AK4CaloJet80_Eta5p1_v1_Prescl, &b_HLT_AK4CaloJet80_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet100_Eta5p1_v1", &HLT_AK4CaloJet100_Eta5p1_v1, &b_HLT_AK4CaloJet100_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet100_Eta5p1_v1_Prescl", &HLT_AK4CaloJet100_Eta5p1_v1_Prescl, &b_HLT_AK4CaloJet100_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet110_Eta5p1_v1", &HLT_AK4CaloJet110_Eta5p1_v1, &b_HLT_AK4CaloJet110_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet110_Eta5p1_v1_Prescl", &HLT_AK4CaloJet110_Eta5p1_v1_Prescl, &b_HLT_AK4CaloJet110_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet120_Eta5p1_v1", &HLT_AK4CaloJet120_Eta5p1_v1, &b_HLT_AK4CaloJet120_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet120_Eta5p1_v1_Prescl", &HLT_AK4CaloJet120_Eta5p1_v1_Prescl, &b_HLT_AK4CaloJet120_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet150_v1", &HLT_AK4CaloJet150_v1, &b_HLT_AK4CaloJet150_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet150_v1_Prescl", &HLT_AK4CaloJet150_v1_Prescl, &b_HLT_AK4CaloJet150_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFJet40_Eta5p1_v1", &HLT_AK4PFJet40_Eta5p1_v1, &b_HLT_AK4PFJet40_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFJet40_Eta5p1_v1_Prescl", &HLT_AK4PFJet40_Eta5p1_v1_Prescl, &b_HLT_AK4PFJet40_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFJet60_Eta5p1_v1", &HLT_AK4PFJet60_Eta5p1_v1, &b_HLT_AK4PFJet60_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFJet60_Eta5p1_v1_Prescl", &HLT_AK4PFJet60_Eta5p1_v1_Prescl, &b_HLT_AK4PFJet60_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFJet80_Eta5p1_v1", &HLT_AK4PFJet80_Eta5p1_v1, &b_HLT_AK4PFJet80_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFJet80_Eta5p1_v1_Prescl", &HLT_AK4PFJet80_Eta5p1_v1_Prescl, &b_HLT_AK4PFJet80_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFJet100_Eta5p1_v1", &HLT_AK4PFJet100_Eta5p1_v1, &b_HLT_AK4PFJet100_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFJet100_Eta5p1_v1_Prescl", &HLT_AK4PFJet100_Eta5p1_v1_Prescl, &b_HLT_AK4PFJet100_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFJet110_Eta5p1_v1", &HLT_AK4PFJet110_Eta5p1_v1, &b_HLT_AK4PFJet110_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFJet110_Eta5p1_v1_Prescl", &HLT_AK4PFJet110_Eta5p1_v1_Prescl, &b_HLT_AK4PFJet110_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFJet120_Eta5p1_v1", &HLT_AK4PFJet120_Eta5p1_v1, &b_HLT_AK4PFJet120_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFJet120_Eta5p1_v1_Prescl", &HLT_AK4PFJet120_Eta5p1_v1_Prescl, &b_HLT_AK4PFJet120_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet80_Jet35_Eta1p1_v1", &HLT_AK4CaloJet80_Jet35_Eta1p1_v1, &b_HLT_AK4CaloJet80_Jet35_Eta1p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet80_Jet35_Eta1p1_v1_Prescl", &HLT_AK4CaloJet80_Jet35_Eta1p1_v1_Prescl, &b_HLT_AK4CaloJet80_Jet35_Eta1p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet80_Jet35_Eta0p7_v1", &HLT_AK4CaloJet80_Jet35_Eta0p7_v1, &b_HLT_AK4CaloJet80_Jet35_Eta0p7_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet80_Jet35_Eta0p7_v1_Prescl", &HLT_AK4CaloJet80_Jet35_Eta0p7_v1_Prescl, &b_HLT_AK4CaloJet80_Jet35_Eta0p7_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet100_Jet35_Eta1p1_v1", &HLT_AK4CaloJet100_Jet35_Eta1p1_v1, &b_HLT_AK4CaloJet100_Jet35_Eta1p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet100_Jet35_Eta1p1_v1_Prescl", &HLT_AK4CaloJet100_Jet35_Eta1p1_v1_Prescl, &b_HLT_AK4CaloJet100_Jet35_Eta1p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet100_Jet35_Eta0p7_v1", &HLT_AK4CaloJet100_Jet35_Eta0p7_v1, &b_HLT_AK4CaloJet100_Jet35_Eta0p7_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet100_Jet35_Eta0p7_v1_Prescl", &HLT_AK4CaloJet100_Jet35_Eta0p7_v1_Prescl, &b_HLT_AK4CaloJet100_Jet35_Eta0p7_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet80_45_45_Eta2p1_v1", &HLT_AK4CaloJet80_45_45_Eta2p1_v1, &b_HLT_AK4CaloJet80_45_45_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet80_45_45_Eta2p1_v1_Prescl", &HLT_AK4CaloJet80_45_45_Eta2p1_v1_Prescl, &b_HLT_AK4CaloJet80_45_45_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton10_Eta1p5_v1", &HLT_HISinglePhoton10_Eta1p5_v1, &b_HLT_HISinglePhoton10_Eta1p5_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton10_Eta1p5_v1_Prescl", &HLT_HISinglePhoton10_Eta1p5_v1_Prescl, &b_HLT_HISinglePhoton10_Eta1p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton15_Eta1p5_v1", &HLT_HISinglePhoton15_Eta1p5_v1, &b_HLT_HISinglePhoton15_Eta1p5_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton15_Eta1p5_v1_Prescl", &HLT_HISinglePhoton15_Eta1p5_v1_Prescl, &b_HLT_HISinglePhoton15_Eta1p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton20_Eta1p5_v1", &HLT_HISinglePhoton20_Eta1p5_v1, &b_HLT_HISinglePhoton20_Eta1p5_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton20_Eta1p5_v1_Prescl", &HLT_HISinglePhoton20_Eta1p5_v1_Prescl, &b_HLT_HISinglePhoton20_Eta1p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton30_Eta1p5_v1", &HLT_HISinglePhoton30_Eta1p5_v1, &b_HLT_HISinglePhoton30_Eta1p5_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton30_Eta1p5_v1_Prescl", &HLT_HISinglePhoton30_Eta1p5_v1_Prescl, &b_HLT_HISinglePhoton30_Eta1p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton40_Eta1p5_v1", &HLT_HISinglePhoton40_Eta1p5_v1, &b_HLT_HISinglePhoton40_Eta1p5_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton40_Eta1p5_v1_Prescl", &HLT_HISinglePhoton40_Eta1p5_v1_Prescl, &b_HLT_HISinglePhoton40_Eta1p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton50_Eta1p5_v1", &HLT_HISinglePhoton50_Eta1p5_v1, &b_HLT_HISinglePhoton50_Eta1p5_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton50_Eta1p5_v1_Prescl", &HLT_HISinglePhoton50_Eta1p5_v1_Prescl, &b_HLT_HISinglePhoton50_Eta1p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton60_Eta1p5_v1", &HLT_HISinglePhoton60_Eta1p5_v1, &b_HLT_HISinglePhoton60_Eta1p5_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton60_Eta1p5_v1_Prescl", &HLT_HISinglePhoton60_Eta1p5_v1_Prescl, &b_HLT_HISinglePhoton60_Eta1p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton10_Eta3p1_v1", &HLT_HISinglePhoton10_Eta3p1_v1, &b_HLT_HISinglePhoton10_Eta3p1_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton10_Eta3p1_v1_Prescl", &HLT_HISinglePhoton10_Eta3p1_v1_Prescl, &b_HLT_HISinglePhoton10_Eta3p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton15_Eta3p1_v1", &HLT_HISinglePhoton15_Eta3p1_v1, &b_HLT_HISinglePhoton15_Eta3p1_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton15_Eta3p1_v1_Prescl", &HLT_HISinglePhoton15_Eta3p1_v1_Prescl, &b_HLT_HISinglePhoton15_Eta3p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton20_Eta3p1_v1", &HLT_HISinglePhoton20_Eta3p1_v1, &b_HLT_HISinglePhoton20_Eta3p1_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton20_Eta3p1_v1_Prescl", &HLT_HISinglePhoton20_Eta3p1_v1_Prescl, &b_HLT_HISinglePhoton20_Eta3p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton30_Eta3p1_v1", &HLT_HISinglePhoton30_Eta3p1_v1, &b_HLT_HISinglePhoton30_Eta3p1_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton30_Eta3p1_v1_Prescl", &HLT_HISinglePhoton30_Eta3p1_v1_Prescl, &b_HLT_HISinglePhoton30_Eta3p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton40_Eta3p1_v1", &HLT_HISinglePhoton40_Eta3p1_v1, &b_HLT_HISinglePhoton40_Eta3p1_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton40_Eta3p1_v1_Prescl", &HLT_HISinglePhoton40_Eta3p1_v1_Prescl, &b_HLT_HISinglePhoton40_Eta3p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton50_Eta3p1_v1", &HLT_HISinglePhoton50_Eta3p1_v1, &b_HLT_HISinglePhoton50_Eta3p1_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton50_Eta3p1_v1_Prescl", &HLT_HISinglePhoton50_Eta3p1_v1_Prescl, &b_HLT_HISinglePhoton50_Eta3p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton60_Eta3p1_v1", &HLT_HISinglePhoton60_Eta3p1_v1, &b_HLT_HISinglePhoton60_Eta3p1_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton60_Eta3p1_v1_Prescl", &HLT_HISinglePhoton60_Eta3p1_v1_Prescl, &b_HLT_HISinglePhoton60_Eta3p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_v1", &HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_v1, &b_HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_v1);
   fChain->SetBranchAddress("HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_v1_Prescl", &HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_v1_Prescl, &b_HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_R9HECut_v1", &HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_R9HECut_v1, &b_HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_R9HECut_v1);
   fChain->SetBranchAddress("HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_R9HECut_v1_Prescl", &HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_R9HECut_v1_Prescl, &b_HLT_HIDoublePhoton15_Eta1p5_Mass50_1000_R9HECut_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIDoublePhoton15_Eta2p1_Mass50_1000_R9Cut_v1", &HLT_HIDoublePhoton15_Eta2p1_Mass50_1000_R9Cut_v1, &b_HLT_HIDoublePhoton15_Eta2p1_Mass50_1000_R9Cut_v1);
   fChain->SetBranchAddress("HLT_HIDoublePhoton15_Eta2p1_Mass50_1000_R9Cut_v1_Prescl", &HLT_HIDoublePhoton15_Eta2p1_Mass50_1000_R9Cut_v1_Prescl, &b_HLT_HIDoublePhoton15_Eta2p1_Mass50_1000_R9Cut_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1", &HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1, &b_HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1);
   fChain->SetBranchAddress("HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1_Prescl", &HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1_Prescl, &b_HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1", &HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1, &b_HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1_Prescl", &HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1", &HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1, &b_HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1_Prescl", &HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1", &HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1, &b_HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1_Prescl", &HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1", &HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1, &b_HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1_Prescl", &HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1", &HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1, &b_HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1_Prescl", &HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1", &HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1, &b_HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1_Prescl", &HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1", &HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1, &b_HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1_Prescl", &HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1", &HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1, &b_HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1_Prescl", &HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1", &HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1, &b_HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1_Prescl", &HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL1DoubleMu0_v1", &HLT_HIL1DoubleMu0_v1, &b_HLT_HIL1DoubleMu0_v1);
   fChain->SetBranchAddress("HLT_HIL1DoubleMu0_v1_Prescl", &HLT_HIL1DoubleMu0_v1_Prescl, &b_HLT_HIL1DoubleMu0_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL1DoubleMu10_v1", &HLT_HIL1DoubleMu10_v1, &b_HLT_HIL1DoubleMu10_v1);
   fChain->SetBranchAddress("HLT_HIL1DoubleMu10_v1_Prescl", &HLT_HIL1DoubleMu10_v1_Prescl, &b_HLT_HIL1DoubleMu10_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v1", &HLT_HIL2DoubleMu0_NHitQ_v1, &b_HLT_HIL2DoubleMu0_NHitQ_v1);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v1_Prescl", &HLT_HIL2DoubleMu0_NHitQ_v1_Prescl, &b_HLT_HIL2DoubleMu0_NHitQ_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1", &HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1, &b_HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1);
   fChain->SetBranchAddress("HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1_Prescl", &HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1_Prescl, &b_HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMu0_OS_m7to14_v1", &HLT_HIL3DoubleMu0_OS_m7to14_v1, &b_HLT_HIL3DoubleMu0_OS_m7to14_v1);
   fChain->SetBranchAddress("HLT_HIL3DoubleMu0_OS_m7to14_v1_Prescl", &HLT_HIL3DoubleMu0_OS_m7to14_v1_Prescl, &b_HLT_HIL3DoubleMu0_OS_m7to14_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3_NHitQ10_v1", &HLT_HIL2Mu3_NHitQ10_v1, &b_HLT_HIL2Mu3_NHitQ10_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3_NHitQ10_v1_Prescl", &HLT_HIL2Mu3_NHitQ10_v1_Prescl, &b_HLT_HIL2Mu3_NHitQ10_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3Mu3_NHitQ15_v1", &HLT_HIL3Mu3_NHitQ15_v1, &b_HLT_HIL3Mu3_NHitQ15_v1);
   fChain->SetBranchAddress("HLT_HIL3Mu3_NHitQ15_v1_Prescl", &HLT_HIL3Mu3_NHitQ15_v1_Prescl, &b_HLT_HIL3Mu3_NHitQ15_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu5_NHitQ10_v1", &HLT_HIL2Mu5_NHitQ10_v1, &b_HLT_HIL2Mu5_NHitQ10_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu5_NHitQ10_v1_Prescl", &HLT_HIL2Mu5_NHitQ10_v1_Prescl, &b_HLT_HIL2Mu5_NHitQ10_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3Mu5_NHitQ15_v1", &HLT_HIL3Mu5_NHitQ15_v1, &b_HLT_HIL3Mu5_NHitQ15_v1);
   fChain->SetBranchAddress("HLT_HIL3Mu5_NHitQ15_v1_Prescl", &HLT_HIL3Mu5_NHitQ15_v1_Prescl, &b_HLT_HIL3Mu5_NHitQ15_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu7_NHitQ10_v1", &HLT_HIL2Mu7_NHitQ10_v1, &b_HLT_HIL2Mu7_NHitQ10_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu7_NHitQ10_v1_Prescl", &HLT_HIL2Mu7_NHitQ10_v1_Prescl, &b_HLT_HIL2Mu7_NHitQ10_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3Mu7_NHitQ15_v1", &HLT_HIL3Mu7_NHitQ15_v1, &b_HLT_HIL3Mu7_NHitQ15_v1);
   fChain->SetBranchAddress("HLT_HIL3Mu7_NHitQ15_v1_Prescl", &HLT_HIL3Mu7_NHitQ15_v1_Prescl, &b_HLT_HIL3Mu7_NHitQ15_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu15_v1", &HLT_HIL2Mu15_v1, &b_HLT_HIL2Mu15_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu15_v1_Prescl", &HLT_HIL2Mu15_v1_Prescl, &b_HLT_HIL2Mu15_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3Mu15_v1", &HLT_HIL3Mu15_v1, &b_HLT_HIL3Mu15_v1);
   fChain->SetBranchAddress("HLT_HIL3Mu15_v1_Prescl", &HLT_HIL3Mu15_v1_Prescl, &b_HLT_HIL3Mu15_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu20_v1", &HLT_HIL2Mu20_v1, &b_HLT_HIL2Mu20_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu20_v1_Prescl", &HLT_HIL2Mu20_v1_Prescl, &b_HLT_HIL2Mu20_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3Mu20_v1", &HLT_HIL3Mu20_v1, &b_HLT_HIL3Mu20_v1);
   fChain->SetBranchAddress("HLT_HIL3Mu20_v1_Prescl", &HLT_HIL3Mu20_v1_Prescl, &b_HLT_HIL3Mu20_v1_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack18ForPPRef_v1", &HLT_FullTrack18ForPPRef_v1, &b_HLT_FullTrack18ForPPRef_v1);
   fChain->SetBranchAddress("HLT_FullTrack18ForPPRef_v1_Prescl", &HLT_FullTrack18ForPPRef_v1_Prescl, &b_HLT_FullTrack18ForPPRef_v1_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack24ForPPRef_v1", &HLT_FullTrack24ForPPRef_v1, &b_HLT_FullTrack24ForPPRef_v1);
   fChain->SetBranchAddress("HLT_FullTrack24ForPPRef_v1_Prescl", &HLT_FullTrack24ForPPRef_v1_Prescl, &b_HLT_FullTrack24ForPPRef_v1_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack34ForPPRef_v1", &HLT_FullTrack34ForPPRef_v1, &b_HLT_FullTrack34ForPPRef_v1);
   fChain->SetBranchAddress("HLT_FullTrack34ForPPRef_v1_Prescl", &HLT_FullTrack34ForPPRef_v1_Prescl, &b_HLT_FullTrack34ForPPRef_v1_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack45ForPPRef_v1", &HLT_FullTrack45ForPPRef_v1, &b_HLT_FullTrack45ForPPRef_v1);
   fChain->SetBranchAddress("HLT_FullTrack45ForPPRef_v1_Prescl", &HLT_FullTrack45ForPPRef_v1_Prescl, &b_HLT_FullTrack45ForPPRef_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCL1DoubleMuOpenNotHF2_v1", &HLT_HIUPCL1DoubleMuOpenNotHF2_v1, &b_HLT_HIUPCL1DoubleMuOpenNotHF2_v1);
   fChain->SetBranchAddress("HLT_HIUPCL1DoubleMuOpenNotHF2_v1_Prescl", &HLT_HIUPCL1DoubleMuOpenNotHF2_v1_Prescl, &b_HLT_HIUPCL1DoubleMuOpenNotHF2_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCDoubleMuNotHF2Pixel_SingleTrack_v1", &HLT_HIUPCDoubleMuNotHF2Pixel_SingleTrack_v1, &b_HLT_HIUPCDoubleMuNotHF2Pixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HIUPCDoubleMuNotHF2Pixel_SingleTrack_v1_Prescl", &HLT_HIUPCDoubleMuNotHF2Pixel_SingleTrack_v1_Prescl, &b_HLT_HIUPCDoubleMuNotHF2Pixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCL1MuOpen_NotMinimumBiasHF2_AND_v1", &HLT_HIUPCL1MuOpen_NotMinimumBiasHF2_AND_v1, &b_HLT_HIUPCL1MuOpen_NotMinimumBiasHF2_AND_v1);
   fChain->SetBranchAddress("HLT_HIUPCL1MuOpen_NotMinimumBiasHF2_AND_v1_Prescl", &HLT_HIUPCL1MuOpen_NotMinimumBiasHF2_AND_v1_Prescl, &b_HLT_HIUPCL1MuOpen_NotMinimumBiasHF2_AND_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCMuOpen_NotMinimumBiasHF2_ANDPixel_SingleTrack_v1", &HLT_HIUPCMuOpen_NotMinimumBiasHF2_ANDPixel_SingleTrack_v1, &b_HLT_HIUPCMuOpen_NotMinimumBiasHF2_ANDPixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HIUPCMuOpen_NotMinimumBiasHF2_ANDPixel_SingleTrack_v1_Prescl", &HLT_HIUPCMuOpen_NotMinimumBiasHF2_ANDPixel_SingleTrack_v1_Prescl, &b_HLT_HIUPCMuOpen_NotMinimumBiasHF2_ANDPixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCL1NotMinimumBiasHF2_AND_v1", &HLT_HIUPCL1NotMinimumBiasHF2_AND_v1, &b_HLT_HIUPCL1NotMinimumBiasHF2_AND_v1);
   fChain->SetBranchAddress("HLT_HIUPCL1NotMinimumBiasHF2_AND_v1_Prescl", &HLT_HIUPCL1NotMinimumBiasHF2_AND_v1_Prescl, &b_HLT_HIUPCL1NotMinimumBiasHF2_AND_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNotMinimumBiasHF2_ANDPixel_SingleTrack_v1", &HLT_HIUPCNotMinimumBiasHF2_ANDPixel_SingleTrack_v1, &b_HLT_HIUPCNotMinimumBiasHF2_ANDPixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HIUPCNotMinimumBiasHF2_ANDPixel_SingleTrack_v1_Prescl", &HLT_HIUPCNotMinimumBiasHF2_ANDPixel_SingleTrack_v1_Prescl, &b_HLT_HIUPCNotMinimumBiasHF2_ANDPixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCL1ZdcOR_BptxAND_v1", &HLT_HIUPCL1ZdcOR_BptxAND_v1, &b_HLT_HIUPCL1ZdcOR_BptxAND_v1);
   fChain->SetBranchAddress("HLT_HIUPCL1ZdcOR_BptxAND_v1_Prescl", &HLT_HIUPCL1ZdcOR_BptxAND_v1_Prescl, &b_HLT_HIUPCL1ZdcOR_BptxAND_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCZdcOR_BptxANDPixel_SingleTrack_v1", &HLT_HIUPCZdcOR_BptxANDPixel_SingleTrack_v1, &b_HLT_HIUPCZdcOR_BptxANDPixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HIUPCZdcOR_BptxANDPixel_SingleTrack_v1_Prescl", &HLT_HIUPCZdcOR_BptxANDPixel_SingleTrack_v1_Prescl, &b_HLT_HIUPCZdcOR_BptxANDPixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCL1ZdcXOR_BptxAND_v1", &HLT_HIUPCL1ZdcXOR_BptxAND_v1, &b_HLT_HIUPCL1ZdcXOR_BptxAND_v1);
   fChain->SetBranchAddress("HLT_HIUPCL1ZdcXOR_BptxAND_v1_Prescl", &HLT_HIUPCL1ZdcXOR_BptxAND_v1_Prescl, &b_HLT_HIUPCL1ZdcXOR_BptxAND_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCZdcXOR_BptxANDPixel_SingleTrack_v1", &HLT_HIUPCZdcXOR_BptxANDPixel_SingleTrack_v1, &b_HLT_HIUPCZdcXOR_BptxANDPixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HIUPCZdcXOR_BptxANDPixel_SingleTrack_v1_Prescl", &HLT_HIUPCZdcXOR_BptxANDPixel_SingleTrack_v1_Prescl, &b_HLT_HIUPCZdcXOR_BptxANDPixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCL1NotZdcOR_BptxAND_v1", &HLT_HIUPCL1NotZdcOR_BptxAND_v1, &b_HLT_HIUPCL1NotZdcOR_BptxAND_v1);
   fChain->SetBranchAddress("HLT_HIUPCL1NotZdcOR_BptxAND_v1_Prescl", &HLT_HIUPCL1NotZdcOR_BptxAND_v1_Prescl, &b_HLT_HIUPCL1NotZdcOR_BptxAND_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNotZdcOR_BptxANDPixel_SingleTrack_v1", &HLT_HIUPCNotZdcOR_BptxANDPixel_SingleTrack_v1, &b_HLT_HIUPCNotZdcOR_BptxANDPixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HIUPCNotZdcOR_BptxANDPixel_SingleTrack_v1_Prescl", &HLT_HIUPCNotZdcOR_BptxANDPixel_SingleTrack_v1_Prescl, &b_HLT_HIUPCNotZdcOR_BptxANDPixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL1CastorMediumJet_v1", &HLT_HIL1CastorMediumJet_v1, &b_HLT_HIL1CastorMediumJet_v1);
   fChain->SetBranchAddress("HLT_HIL1CastorMediumJet_v1_Prescl", &HLT_HIL1CastorMediumJet_v1_Prescl, &b_HLT_HIL1CastorMediumJet_v1_Prescl);
   fChain->SetBranchAddress("HLT_HICastorMediumJetPixel_SingleTrack_v1", &HLT_HICastorMediumJetPixel_SingleTrack_v1, &b_HLT_HICastorMediumJetPixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HICastorMediumJetPixel_SingleTrack_v1_Prescl", &HLT_HICastorMediumJetPixel_SingleTrack_v1_Prescl, &b_HLT_HICastorMediumJetPixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt15_v1", &HLT_DmesonPPTrackingGlobal_Dpt15_v1, &b_HLT_DmesonPPTrackingGlobal_Dpt15_v1);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt15_v1_Prescl", &HLT_DmesonPPTrackingGlobal_Dpt15_v1_Prescl, &b_HLT_DmesonPPTrackingGlobal_Dpt15_v1_Prescl);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt20_v1", &HLT_DmesonPPTrackingGlobal_Dpt20_v1, &b_HLT_DmesonPPTrackingGlobal_Dpt20_v1);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt20_v1_Prescl", &HLT_DmesonPPTrackingGlobal_Dpt20_v1_Prescl, &b_HLT_DmesonPPTrackingGlobal_Dpt20_v1_Prescl);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt30_v1", &HLT_DmesonPPTrackingGlobal_Dpt30_v1, &b_HLT_DmesonPPTrackingGlobal_Dpt30_v1);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt30_v1_Prescl", &HLT_DmesonPPTrackingGlobal_Dpt30_v1_Prescl, &b_HLT_DmesonPPTrackingGlobal_Dpt30_v1_Prescl);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt40_v1", &HLT_DmesonPPTrackingGlobal_Dpt40_v1, &b_HLT_DmesonPPTrackingGlobal_Dpt40_v1);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt40_v1_Prescl", &HLT_DmesonPPTrackingGlobal_Dpt40_v1_Prescl, &b_HLT_DmesonPPTrackingGlobal_Dpt40_v1_Prescl);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt50_v1", &HLT_DmesonPPTrackingGlobal_Dpt50_v1, &b_HLT_DmesonPPTrackingGlobal_Dpt50_v1);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt50_v1_Prescl", &HLT_DmesonPPTrackingGlobal_Dpt50_v1_Prescl, &b_HLT_DmesonPPTrackingGlobal_Dpt50_v1_Prescl);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt60_v1", &HLT_DmesonPPTrackingGlobal_Dpt60_v1, &b_HLT_DmesonPPTrackingGlobal_Dpt60_v1);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt60_v1_Prescl", &HLT_DmesonPPTrackingGlobal_Dpt60_v1_Prescl, &b_HLT_DmesonPPTrackingGlobal_Dpt60_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part0_v1", &HLT_L1MinimumBiasHF1OR_part0_v1, &b_HLT_L1MinimumBiasHF1OR_part0_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part0_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part0_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part0_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part1_v1", &HLT_L1MinimumBiasHF1OR_part1_v1, &b_HLT_L1MinimumBiasHF1OR_part1_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part1_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part1_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part1_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part2_v1", &HLT_L1MinimumBiasHF1OR_part2_v1, &b_HLT_L1MinimumBiasHF1OR_part2_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part2_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part2_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part2_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part3_v1", &HLT_L1MinimumBiasHF1OR_part3_v1, &b_HLT_L1MinimumBiasHF1OR_part3_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part3_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part3_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part3_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part4_v1", &HLT_L1MinimumBiasHF1OR_part4_v1, &b_HLT_L1MinimumBiasHF1OR_part4_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part4_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part4_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part4_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part5_v1", &HLT_L1MinimumBiasHF1OR_part5_v1, &b_HLT_L1MinimumBiasHF1OR_part5_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part5_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part5_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part5_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part6_v1", &HLT_L1MinimumBiasHF1OR_part6_v1, &b_HLT_L1MinimumBiasHF1OR_part6_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part6_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part6_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part6_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part7_v1", &HLT_L1MinimumBiasHF1OR_part7_v1, &b_HLT_L1MinimumBiasHF1OR_part7_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part7_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part7_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part7_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part8_v1", &HLT_L1MinimumBiasHF1OR_part8_v1, &b_HLT_L1MinimumBiasHF1OR_part8_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part8_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part8_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part8_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part9_v1", &HLT_L1MinimumBiasHF1OR_part9_v1, &b_HLT_L1MinimumBiasHF1OR_part9_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part9_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part9_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part9_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part10_v1", &HLT_L1MinimumBiasHF1OR_part10_v1, &b_HLT_L1MinimumBiasHF1OR_part10_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part10_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part10_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part10_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part11_v1", &HLT_L1MinimumBiasHF1OR_part11_v1, &b_HLT_L1MinimumBiasHF1OR_part11_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part11_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part11_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part11_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part12_v1", &HLT_L1MinimumBiasHF1OR_part12_v1, &b_HLT_L1MinimumBiasHF1OR_part12_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part12_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part12_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part12_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part13_v1", &HLT_L1MinimumBiasHF1OR_part13_v1, &b_HLT_L1MinimumBiasHF1OR_part13_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part13_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part13_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part13_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part14_v1", &HLT_L1MinimumBiasHF1OR_part14_v1, &b_HLT_L1MinimumBiasHF1OR_part14_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part14_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part14_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part14_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part15_v1", &HLT_L1MinimumBiasHF1OR_part15_v1, &b_HLT_L1MinimumBiasHF1OR_part15_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part15_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part15_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part15_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part16_v1", &HLT_L1MinimumBiasHF1OR_part16_v1, &b_HLT_L1MinimumBiasHF1OR_part16_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part16_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part16_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part16_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part17_v1", &HLT_L1MinimumBiasHF1OR_part17_v1, &b_HLT_L1MinimumBiasHF1OR_part17_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part17_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part17_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part17_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part18_v1", &HLT_L1MinimumBiasHF1OR_part18_v1, &b_HLT_L1MinimumBiasHF1OR_part18_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part18_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part18_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part18_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part19_v1", &HLT_L1MinimumBiasHF1OR_part19_v1, &b_HLT_L1MinimumBiasHF1OR_part19_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part19_v1_Prescl", &HLT_L1MinimumBiasHF1OR_part19_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_part19_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFBJetBCSV60_Eta2p1_v1", &HLT_AK4PFBJetBCSV60_Eta2p1_v1, &b_HLT_AK4PFBJetBCSV60_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFBJetBCSV60_Eta2p1_v1_Prescl", &HLT_AK4PFBJetBCSV60_Eta2p1_v1_Prescl, &b_HLT_AK4PFBJetBCSV60_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFBJetBCSV80_Eta2p1_v1", &HLT_AK4PFBJetBCSV80_Eta2p1_v1, &b_HLT_AK4PFBJetBCSV80_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFBJetBCSV80_Eta2p1_v1_Prescl", &HLT_AK4PFBJetBCSV80_Eta2p1_v1_Prescl, &b_HLT_AK4PFBJetBCSV80_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFDJet60_Eta2p1_v1", &HLT_AK4PFDJet60_Eta2p1_v1, &b_HLT_AK4PFDJet60_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFDJet60_Eta2p1_v1_Prescl", &HLT_AK4PFDJet60_Eta2p1_v1_Prescl, &b_HLT_AK4PFDJet60_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFDJet80_Eta2p1_v1", &HLT_AK4PFDJet80_Eta2p1_v1, &b_HLT_AK4PFDJet80_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFDJet80_Eta2p1_v1_Prescl", &HLT_AK4PFDJet80_Eta2p1_v1_Prescl, &b_HLT_AK4PFDJet80_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFBJetBSSV60_Eta2p1_v1", &HLT_AK4PFBJetBSSV60_Eta2p1_v1, &b_HLT_AK4PFBJetBSSV60_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFBJetBSSV60_Eta2p1_v1_Prescl", &HLT_AK4PFBJetBSSV60_Eta2p1_v1_Prescl, &b_HLT_AK4PFBJetBSSV60_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFBJetBSSV80_Eta2p1_v1", &HLT_AK4PFBJetBSSV80_Eta2p1_v1, &b_HLT_AK4PFBJetBSSV80_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFBJetBSSV80_Eta2p1_v1_Prescl", &HLT_AK4PFBJetBSSV80_Eta2p1_v1_Prescl, &b_HLT_AK4PFBJetBSSV80_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath, &b_HLTriggerFinalPath);
   fChain->SetBranchAddress("HLTriggerFinalPath_Prescl", &HLTriggerFinalPath_Prescl, &b_HLTriggerFinalPath_Prescl);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt8_v1", &HLT_DmesonPPTrackingGlobal_Dpt8_v1, &b_HLT_DmesonPPTrackingGlobal_Dpt8_v1);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt8_v1_Prescl", &HLT_DmesonPPTrackingGlobal_Dpt8_v1_Prescl, &b_HLT_DmesonPPTrackingGlobal_Dpt8_v1_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity60_v1", &HLT_PixelTracks_Multiplicity60_v1, &b_HLT_PixelTracks_Multiplicity60_v1);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity60_v1_Prescl", &HLT_PixelTracks_Multiplicity60_v1_Prescl, &b_HLT_PixelTracks_Multiplicity60_v1_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity60_v3", &HLT_PixelTracks_Multiplicity60_v3, &b_HLT_PixelTracks_Multiplicity60_v3);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity60_v3_Prescl", &HLT_PixelTracks_Multiplicity60_v3_Prescl, &b_HLT_PixelTracks_Multiplicity60_v3_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity85_v1", &HLT_PixelTracks_Multiplicity85_v1, &b_HLT_PixelTracks_Multiplicity85_v1);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity85_v1_Prescl", &HLT_PixelTracks_Multiplicity85_v1_Prescl, &b_HLT_PixelTracks_Multiplicity85_v1_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity110_v1", &HLT_PixelTracks_Multiplicity110_v1, &b_HLT_PixelTracks_Multiplicity110_v1);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity110_v1_Prescl", &HLT_PixelTracks_Multiplicity110_v1_Prescl, &b_HLT_PixelTracks_Multiplicity110_v1_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity135_v1", &HLT_PixelTracks_Multiplicity135_v1, &b_HLT_PixelTracks_Multiplicity135_v1);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity135_v1_Prescl", &HLT_PixelTracks_Multiplicity135_v1_Prescl, &b_HLT_PixelTracks_Multiplicity135_v1_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity160_v1", &HLT_PixelTracks_Multiplicity160_v1, &b_HLT_PixelTracks_Multiplicity160_v1);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity160_v1_Prescl", &HLT_PixelTracks_Multiplicity160_v1_Prescl, &b_HLT_PixelTracks_Multiplicity160_v1_Prescl);
   fChain->SetBranchAddress("HLT_Physics_v1", &HLT_Physics_v1, &b_HLT_Physics_v1);
   fChain->SetBranchAddress("HLT_Physics_v1_Prescl", &HLT_Physics_v1_Prescl, &b_HLT_Physics_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech5_BPTX_PlusOnly_v1", &HLT_L1Tech5_BPTX_PlusOnly_v1, &b_HLT_L1Tech5_BPTX_PlusOnly_v1);
   fChain->SetBranchAddress("HLT_L1Tech5_BPTX_PlusOnly_v1_Prescl", &HLT_L1Tech5_BPTX_PlusOnly_v1_Prescl, &b_HLT_L1Tech5_BPTX_PlusOnly_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_v1", &HLT_ZeroBias_v1, &b_HLT_ZeroBias_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_v1_Prescl", &HLT_ZeroBias_v1_Prescl, &b_HLT_ZeroBias_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part0_v1", &HLT_ZeroBias_part0_v1, &b_HLT_ZeroBias_part0_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part0_v1_Prescl", &HLT_ZeroBias_part0_v1_Prescl, &b_HLT_ZeroBias_part0_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part1_v1", &HLT_ZeroBias_part1_v1, &b_HLT_ZeroBias_part1_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part1_v1_Prescl", &HLT_ZeroBias_part1_v1_Prescl, &b_HLT_ZeroBias_part1_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part2_v1", &HLT_ZeroBias_part2_v1, &b_HLT_ZeroBias_part2_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part2_v1_Prescl", &HLT_ZeroBias_part2_v1_Prescl, &b_HLT_ZeroBias_part2_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part3_v1", &HLT_ZeroBias_part3_v1, &b_HLT_ZeroBias_part3_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part3_v1_Prescl", &HLT_ZeroBias_part3_v1_Prescl, &b_HLT_ZeroBias_part3_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part4_v1", &HLT_ZeroBias_part4_v1, &b_HLT_ZeroBias_part4_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part4_v1_Prescl", &HLT_ZeroBias_part4_v1_Prescl, &b_HLT_ZeroBias_part4_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part5_v1", &HLT_ZeroBias_part5_v1, &b_HLT_ZeroBias_part5_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part5_v1_Prescl", &HLT_ZeroBias_part5_v1_Prescl, &b_HLT_ZeroBias_part5_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part6_v1", &HLT_ZeroBias_part6_v1, &b_HLT_ZeroBias_part6_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part6_v1_Prescl", &HLT_ZeroBias_part6_v1_Prescl, &b_HLT_ZeroBias_part6_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part7_v1", &HLT_ZeroBias_part7_v1, &b_HLT_ZeroBias_part7_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part7_v1_Prescl", &HLT_ZeroBias_part7_v1_Prescl, &b_HLT_ZeroBias_part7_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part8_v1", &HLT_ZeroBias_part8_v1, &b_HLT_ZeroBias_part8_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part8_v1_Prescl", &HLT_ZeroBias_part8_v1_Prescl, &b_HLT_ZeroBias_part8_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part9_v1", &HLT_ZeroBias_part9_v1, &b_HLT_ZeroBias_part9_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part9_v1_Prescl", &HLT_ZeroBias_part9_v1_Prescl, &b_HLT_ZeroBias_part9_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part10_v1", &HLT_ZeroBias_part10_v1, &b_HLT_ZeroBias_part10_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part10_v1_Prescl", &HLT_ZeroBias_part10_v1_Prescl, &b_HLT_ZeroBias_part10_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part11_v1", &HLT_ZeroBias_part11_v1, &b_HLT_ZeroBias_part11_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part11_v1_Prescl", &HLT_ZeroBias_part11_v1_Prescl, &b_HLT_ZeroBias_part11_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part12_v1", &HLT_ZeroBias_part12_v1, &b_HLT_ZeroBias_part12_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part12_v1_Prescl", &HLT_ZeroBias_part12_v1_Prescl, &b_HLT_ZeroBias_part12_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part13_v1", &HLT_ZeroBias_part13_v1, &b_HLT_ZeroBias_part13_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part13_v1_Prescl", &HLT_ZeroBias_part13_v1_Prescl, &b_HLT_ZeroBias_part13_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part14_v1", &HLT_ZeroBias_part14_v1, &b_HLT_ZeroBias_part14_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part14_v1_Prescl", &HLT_ZeroBias_part14_v1_Prescl, &b_HLT_ZeroBias_part14_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part15_v1", &HLT_ZeroBias_part15_v1, &b_HLT_ZeroBias_part15_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part15_v1_Prescl", &HLT_ZeroBias_part15_v1_Prescl, &b_HLT_ZeroBias_part15_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part16_v1", &HLT_ZeroBias_part16_v1, &b_HLT_ZeroBias_part16_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part16_v1_Prescl", &HLT_ZeroBias_part16_v1_Prescl, &b_HLT_ZeroBias_part16_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part17_v1", &HLT_ZeroBias_part17_v1, &b_HLT_ZeroBias_part17_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part17_v1_Prescl", &HLT_ZeroBias_part17_v1_Prescl, &b_HLT_ZeroBias_part17_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part18_v1", &HLT_ZeroBias_part18_v1, &b_HLT_ZeroBias_part18_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part18_v1_Prescl", &HLT_ZeroBias_part18_v1_Prescl, &b_HLT_ZeroBias_part18_v1_Prescl);
   fChain->SetBranchAddress("HLT_ZeroBias_part19_v1", &HLT_ZeroBias_part19_v1, &b_HLT_ZeroBias_part19_v1);
   fChain->SetBranchAddress("HLT_ZeroBias_part19_v1_Prescl", &HLT_ZeroBias_part19_v1_Prescl, &b_HLT_ZeroBias_part19_v1_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack18ForPPRef_v2", &HLT_FullTrack18ForPPRef_v2, &b_HLT_FullTrack18ForPPRef_v2);
   fChain->SetBranchAddress("HLT_FullTrack18ForPPRef_v2_Prescl", &HLT_FullTrack18ForPPRef_v2_Prescl, &b_HLT_FullTrack18ForPPRef_v2_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack18ForPPRef_v3", &HLT_FullTrack18ForPPRef_v3, &b_HLT_FullTrack18ForPPRef_v3);
   fChain->SetBranchAddress("HLT_FullTrack18ForPPRef_v3_Prescl", &HLT_FullTrack18ForPPRef_v3_Prescl, &b_HLT_FullTrack18ForPPRef_v3_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack24ForPPRef_v2", &HLT_FullTrack24ForPPRef_v2, &b_HLT_FullTrack24ForPPRef_v2);
   fChain->SetBranchAddress("HLT_FullTrack24ForPPRef_v2_Prescl", &HLT_FullTrack24ForPPRef_v2_Prescl, &b_HLT_FullTrack24ForPPRef_v2_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack24ForPPRef_v3", &HLT_FullTrack24ForPPRef_v3, &b_HLT_FullTrack24ForPPRef_v3);
   fChain->SetBranchAddress("HLT_FullTrack24ForPPRef_v3_Prescl", &HLT_FullTrack24ForPPRef_v3_Prescl, &b_HLT_FullTrack24ForPPRef_v3_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack34ForPPRef_v2", &HLT_FullTrack34ForPPRef_v2, &b_HLT_FullTrack34ForPPRef_v2);
   fChain->SetBranchAddress("HLT_FullTrack34ForPPRef_v2_Prescl", &HLT_FullTrack34ForPPRef_v2_Prescl, &b_HLT_FullTrack34ForPPRef_v2_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack34ForPPRef_v3", &HLT_FullTrack34ForPPRef_v3, &b_HLT_FullTrack34ForPPRef_v3);
   fChain->SetBranchAddress("HLT_FullTrack34ForPPRef_v3_Prescl", &HLT_FullTrack34ForPPRef_v3_Prescl, &b_HLT_FullTrack34ForPPRef_v3_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack34ForPPRef_v4", &HLT_FullTrack34ForPPRef_v4, &b_HLT_FullTrack34ForPPRef_v4);
   fChain->SetBranchAddress("HLT_FullTrack34ForPPRef_v4_Prescl", &HLT_FullTrack34ForPPRef_v4_Prescl, &b_HLT_FullTrack34ForPPRef_v4_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack45ForPPRef_v2", &HLT_FullTrack45ForPPRef_v2, &b_HLT_FullTrack45ForPPRef_v2);
   fChain->SetBranchAddress("HLT_FullTrack45ForPPRef_v2_Prescl", &HLT_FullTrack45ForPPRef_v2_Prescl, &b_HLT_FullTrack45ForPPRef_v2_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack45ForPPRef_v3", &HLT_FullTrack45ForPPRef_v3, &b_HLT_FullTrack45ForPPRef_v3);
   fChain->SetBranchAddress("HLT_FullTrack45ForPPRef_v3_Prescl", &HLT_FullTrack45ForPPRef_v3_Prescl, &b_HLT_FullTrack45ForPPRef_v3_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack53ForPPRef_v1", &HLT_FullTrack53ForPPRef_v1, &b_HLT_FullTrack53ForPPRef_v1);
   fChain->SetBranchAddress("HLT_FullTrack53ForPPRef_v1_Prescl", &HLT_FullTrack53ForPPRef_v1_Prescl, &b_HLT_FullTrack53ForPPRef_v1_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack53ForPPRef_v2", &HLT_FullTrack53ForPPRef_v2, &b_HLT_FullTrack53ForPPRef_v2);
   fChain->SetBranchAddress("HLT_FullTrack53ForPPRef_v2_Prescl", &HLT_FullTrack53ForPPRef_v2_Prescl, &b_HLT_FullTrack53ForPPRef_v2_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack53ForPPRef_v3", &HLT_FullTrack53ForPPRef_v3, &b_HLT_FullTrack53ForPPRef_v3);
   fChain->SetBranchAddress("HLT_FullTrack53ForPPRef_v3_Prescl", &HLT_FullTrack53ForPPRef_v3_Prescl, &b_HLT_FullTrack53ForPPRef_v3_Prescl);
   fChain->SetBranchAddress("L1_AlwaysTrue", &L1_AlwaysTrue, &b_L1_AlwaysTrue);
   fChain->SetBranchAddress("L1_AlwaysTrue_Prescl", &L1_AlwaysTrue_Prescl, &b_L1_AlwaysTrue_Prescl);
   fChain->SetBranchAddress("L1_CastorHighJet_BptxAND", &L1_CastorHighJet_BptxAND, &b_L1_CastorHighJet_BptxAND);
   fChain->SetBranchAddress("L1_CastorHighJet_BptxAND_Prescl", &L1_CastorHighJet_BptxAND_Prescl, &b_L1_CastorHighJet_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_CastorMediumJet_BptxAND", &L1_CastorMediumJet_BptxAND, &b_L1_CastorMediumJet_BptxAND);
   fChain->SetBranchAddress("L1_CastorMediumJet_BptxAND_Prescl", &L1_CastorMediumJet_BptxAND_Prescl, &b_L1_CastorMediumJet_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu0_BptxAND", &L1_DoubleMu0_BptxAND, &b_L1_DoubleMu0_BptxAND);
   fChain->SetBranchAddress("L1_DoubleMu0_BptxAND_Prescl", &L1_DoubleMu0_BptxAND_Prescl, &b_L1_DoubleMu0_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu0_MinimumBiasHF1_AND", &L1_DoubleMu0_MinimumBiasHF1_AND, &b_L1_DoubleMu0_MinimumBiasHF1_AND);
   fChain->SetBranchAddress("L1_DoubleMu0_MinimumBiasHF1_AND_Prescl", &L1_DoubleMu0_MinimumBiasHF1_AND_Prescl, &b_L1_DoubleMu0_MinimumBiasHF1_AND_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu3_BptxAND", &L1_DoubleMu3_BptxAND, &b_L1_DoubleMu3_BptxAND);
   fChain->SetBranchAddress("L1_DoubleMu3_BptxAND_Prescl", &L1_DoubleMu3_BptxAND_Prescl, &b_L1_DoubleMu3_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_DoubleMuOpen_BptxAND", &L1_DoubleMuOpen_BptxAND, &b_L1_DoubleMuOpen_BptxAND);
   fChain->SetBranchAddress("L1_DoubleMuOpen_BptxAND_Prescl", &L1_DoubleMuOpen_BptxAND_Prescl, &b_L1_DoubleMuOpen_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_DoubleMuOpen_NotMinimumBiasHF2_AND", &L1_DoubleMuOpen_NotMinimumBiasHF2_AND, &b_L1_DoubleMuOpen_NotMinimumBiasHF2_AND);
   fChain->SetBranchAddress("L1_DoubleMuOpen_NotMinimumBiasHF2_AND_Prescl", &L1_DoubleMuOpen_NotMinimumBiasHF2_AND_Prescl, &b_L1_DoubleMuOpen_NotMinimumBiasHF2_AND_Prescl);
   fChain->SetBranchAddress("L1_DoubleJet20", &L1_DoubleJet20, &b_L1_DoubleJet20);
   fChain->SetBranchAddress("L1_DoubleJet20_Prescl", &L1_DoubleJet20_Prescl, &b_L1_DoubleJet20_Prescl);
   fChain->SetBranchAddress("L1_DoubleJet28_BptxAND", &L1_DoubleJet28_BptxAND, &b_L1_DoubleJet28_BptxAND);
   fChain->SetBranchAddress("L1_DoubleJet28_BptxAND_Prescl", &L1_DoubleJet28_BptxAND_Prescl, &b_L1_DoubleJet28_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_DoubleJet32", &L1_DoubleJet32, &b_L1_DoubleJet32);
   fChain->SetBranchAddress("L1_DoubleJet32_Prescl", &L1_DoubleJet32_Prescl, &b_L1_DoubleJet32_Prescl);
   fChain->SetBranchAddress("L1_ETT130", &L1_ETT130, &b_L1_ETT130);
   fChain->SetBranchAddress("L1_ETT130_Prescl", &L1_ETT130_Prescl, &b_L1_ETT130_Prescl);
   fChain->SetBranchAddress("L1_ETT15_BptxAND", &L1_ETT15_BptxAND, &b_L1_ETT15_BptxAND);
   fChain->SetBranchAddress("L1_ETT15_BptxAND_Prescl", &L1_ETT15_BptxAND_Prescl, &b_L1_ETT15_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT20_BptxAND", &L1_ETT20_BptxAND, &b_L1_ETT20_BptxAND);
   fChain->SetBranchAddress("L1_ETT20_BptxAND_Prescl", &L1_ETT20_BptxAND_Prescl, &b_L1_ETT20_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT30_BptxAND", &L1_ETT30_BptxAND, &b_L1_ETT30_BptxAND);
   fChain->SetBranchAddress("L1_ETT30_BptxAND_Prescl", &L1_ETT30_BptxAND_Prescl, &b_L1_ETT30_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT40", &L1_ETT40, &b_L1_ETT40);
   fChain->SetBranchAddress("L1_ETT40_Prescl", &L1_ETT40_Prescl, &b_L1_ETT40_Prescl);
   fChain->SetBranchAddress("L1_ETT40_BptxAND", &L1_ETT40_BptxAND, &b_L1_ETT40_BptxAND);
   fChain->SetBranchAddress("L1_ETT40_BptxAND_Prescl", &L1_ETT40_BptxAND_Prescl, &b_L1_ETT40_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT45_BptxAND", &L1_ETT45_BptxAND, &b_L1_ETT45_BptxAND);
   fChain->SetBranchAddress("L1_ETT45_BptxAND_Prescl", &L1_ETT45_BptxAND_Prescl, &b_L1_ETT45_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT50_BptxAND", &L1_ETT50_BptxAND, &b_L1_ETT50_BptxAND);
   fChain->SetBranchAddress("L1_ETT50_BptxAND_Prescl", &L1_ETT50_BptxAND_Prescl, &b_L1_ETT50_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT55_BptxAND", &L1_ETT55_BptxAND, &b_L1_ETT55_BptxAND);
   fChain->SetBranchAddress("L1_ETT55_BptxAND_Prescl", &L1_ETT55_BptxAND_Prescl, &b_L1_ETT55_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT60_BptxAND", &L1_ETT60_BptxAND, &b_L1_ETT60_BptxAND);
   fChain->SetBranchAddress("L1_ETT60_BptxAND_Prescl", &L1_ETT60_BptxAND_Prescl, &b_L1_ETT60_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT70_BptxAND", &L1_ETT70_BptxAND, &b_L1_ETT70_BptxAND);
   fChain->SetBranchAddress("L1_ETT70_BptxAND_Prescl", &L1_ETT70_BptxAND_Prescl, &b_L1_ETT70_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT90_BptxAND", &L1_ETT90_BptxAND, &b_L1_ETT90_BptxAND);
   fChain->SetBranchAddress("L1_ETT90_BptxAND_Prescl", &L1_ETT90_BptxAND_Prescl, &b_L1_ETT90_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_MinimumBiasHF1_AND", &L1_MinimumBiasHF1_AND, &b_L1_MinimumBiasHF1_AND);
   fChain->SetBranchAddress("L1_MinimumBiasHF1_AND_Prescl", &L1_MinimumBiasHF1_AND_Prescl, &b_L1_MinimumBiasHF1_AND_Prescl);
   fChain->SetBranchAddress("L1_MinimumBiasHF1_OR", &L1_MinimumBiasHF1_OR, &b_L1_MinimumBiasHF1_OR);
   fChain->SetBranchAddress("L1_MinimumBiasHF1_OR_Prescl", &L1_MinimumBiasHF1_OR_Prescl, &b_L1_MinimumBiasHF1_OR_Prescl);
   fChain->SetBranchAddress("L1_MinimumBiasHF1_XOR", &L1_MinimumBiasHF1_XOR, &b_L1_MinimumBiasHF1_XOR);
   fChain->SetBranchAddress("L1_MinimumBiasHF1_XOR_Prescl", &L1_MinimumBiasHF1_XOR_Prescl, &b_L1_MinimumBiasHF1_XOR_Prescl);
   fChain->SetBranchAddress("L1_MinimumBiasHF2_AND", &L1_MinimumBiasHF2_AND, &b_L1_MinimumBiasHF2_AND);
   fChain->SetBranchAddress("L1_MinimumBiasHF2_AND_Prescl", &L1_MinimumBiasHF2_AND_Prescl, &b_L1_MinimumBiasHF2_AND_Prescl);
   fChain->SetBranchAddress("L1_MinimumBiasHF2_OR", &L1_MinimumBiasHF2_OR, &b_L1_MinimumBiasHF2_OR);
   fChain->SetBranchAddress("L1_MinimumBiasHF2_OR_Prescl", &L1_MinimumBiasHF2_OR_Prescl, &b_L1_MinimumBiasHF2_OR_Prescl);
   fChain->SetBranchAddress("L1_MinimumBiasHF2_OR_NoBptxGating", &L1_MinimumBiasHF2_OR_NoBptxGating, &b_L1_MinimumBiasHF2_OR_NoBptxGating);
   fChain->SetBranchAddress("L1_MinimumBiasHF2_OR_NoBptxGating_Prescl", &L1_MinimumBiasHF2_OR_NoBptxGating_Prescl, &b_L1_MinimumBiasHF2_OR_NoBptxGating_Prescl);
   fChain->SetBranchAddress("L1_MuOpen_NotMinimumBiasHF2_AND", &L1_MuOpen_NotMinimumBiasHF2_AND, &b_L1_MuOpen_NotMinimumBiasHF2_AND);
   fChain->SetBranchAddress("L1_MuOpen_NotMinimumBiasHF2_AND_Prescl", &L1_MuOpen_NotMinimumBiasHF2_AND_Prescl, &b_L1_MuOpen_NotMinimumBiasHF2_AND_Prescl);
   fChain->SetBranchAddress("L1_NotMinimumBiasHF1_OR", &L1_NotMinimumBiasHF1_OR, &b_L1_NotMinimumBiasHF1_OR);
   fChain->SetBranchAddress("L1_NotMinimumBiasHF1_OR_Prescl", &L1_NotMinimumBiasHF1_OR_Prescl, &b_L1_NotMinimumBiasHF1_OR_Prescl);
   fChain->SetBranchAddress("L1_NotMinimumBiasHF2_AND", &L1_NotMinimumBiasHF2_AND, &b_L1_NotMinimumBiasHF2_AND);
   fChain->SetBranchAddress("L1_NotMinimumBiasHF2_AND_Prescl", &L1_NotMinimumBiasHF2_AND_Prescl, &b_L1_NotMinimumBiasHF2_AND_Prescl);
   fChain->SetBranchAddress("L1_NotZdcOR_BptxAND", &L1_NotZdcOR_BptxAND, &b_L1_NotZdcOR_BptxAND);
   fChain->SetBranchAddress("L1_NotZdcOR_BptxAND_Prescl", &L1_NotZdcOR_BptxAND_Prescl, &b_L1_NotZdcOR_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleEG12_BptxAND", &L1_SingleEG12_BptxAND, &b_L1_SingleEG12_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG12_BptxAND_Prescl", &L1_SingleEG12_BptxAND_Prescl, &b_L1_SingleEG12_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleEG15_BptxAND", &L1_SingleEG15_BptxAND, &b_L1_SingleEG15_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG15_BptxAND_Prescl", &L1_SingleEG15_BptxAND_Prescl, &b_L1_SingleEG15_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleEG20", &L1_SingleEG20, &b_L1_SingleEG20);
   fChain->SetBranchAddress("L1_SingleEG20_Prescl", &L1_SingleEG20_Prescl, &b_L1_SingleEG20_Prescl);
   fChain->SetBranchAddress("L1_SingleEG20_BptxAND", &L1_SingleEG20_BptxAND, &b_L1_SingleEG20_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG20_BptxAND_Prescl", &L1_SingleEG20_BptxAND_Prescl, &b_L1_SingleEG20_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleEG2_BptxAND", &L1_SingleEG2_BptxAND, &b_L1_SingleEG2_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG2_BptxAND_Prescl", &L1_SingleEG2_BptxAND_Prescl, &b_L1_SingleEG2_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleEG30_BptxAND", &L1_SingleEG30_BptxAND, &b_L1_SingleEG30_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG30_BptxAND_Prescl", &L1_SingleEG30_BptxAND_Prescl, &b_L1_SingleEG30_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleEG5", &L1_SingleEG5, &b_L1_SingleEG5);
   fChain->SetBranchAddress("L1_SingleEG5_Prescl", &L1_SingleEG5_Prescl, &b_L1_SingleEG5_Prescl);
   fChain->SetBranchAddress("L1_SingleEG5_BptxAND", &L1_SingleEG5_BptxAND, &b_L1_SingleEG5_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG5_BptxAND_Prescl", &L1_SingleEG5_BptxAND_Prescl, &b_L1_SingleEG5_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleEG7_BptxAND", &L1_SingleEG7_BptxAND, &b_L1_SingleEG7_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG7_BptxAND_Prescl", &L1_SingleEG7_BptxAND_Prescl, &b_L1_SingleEG7_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu12_BptxAND", &L1_SingleMu12_BptxAND, &b_L1_SingleMu12_BptxAND);
   fChain->SetBranchAddress("L1_SingleMu12_BptxAND_Prescl", &L1_SingleMu12_BptxAND_Prescl, &b_L1_SingleMu12_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu16_BptxAND", &L1_SingleMu16_BptxAND, &b_L1_SingleMu16_BptxAND);
   fChain->SetBranchAddress("L1_SingleMu16_BptxAND_Prescl", &L1_SingleMu16_BptxAND_Prescl, &b_L1_SingleMu16_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu16_MinimumBiasHF1_AND", &L1_SingleMu16_MinimumBiasHF1_AND, &b_L1_SingleMu16_MinimumBiasHF1_AND);
   fChain->SetBranchAddress("L1_SingleMu16_MinimumBiasHF1_AND_Prescl", &L1_SingleMu16_MinimumBiasHF1_AND_Prescl, &b_L1_SingleMu16_MinimumBiasHF1_AND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu3_BptxAND", &L1_SingleMu3_BptxAND, &b_L1_SingleMu3_BptxAND);
   fChain->SetBranchAddress("L1_SingleMu3_BptxAND_Prescl", &L1_SingleMu3_BptxAND_Prescl, &b_L1_SingleMu3_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu3_MinimumBiasHF1_AND", &L1_SingleMu3_MinimumBiasHF1_AND, &b_L1_SingleMu3_MinimumBiasHF1_AND);
   fChain->SetBranchAddress("L1_SingleMu3_MinimumBiasHF1_AND_Prescl", &L1_SingleMu3_MinimumBiasHF1_AND_Prescl, &b_L1_SingleMu3_MinimumBiasHF1_AND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu3p5", &L1_SingleMu3p5, &b_L1_SingleMu3p5);
   fChain->SetBranchAddress("L1_SingleMu3p5_Prescl", &L1_SingleMu3p5_Prescl, &b_L1_SingleMu3p5_Prescl);
   fChain->SetBranchAddress("L1_SingleMu5_BptxAND", &L1_SingleMu5_BptxAND, &b_L1_SingleMu5_BptxAND);
   fChain->SetBranchAddress("L1_SingleMu5_BptxAND_Prescl", &L1_SingleMu5_BptxAND_Prescl, &b_L1_SingleMu5_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu7_BptxAND", &L1_SingleMu7_BptxAND, &b_L1_SingleMu7_BptxAND);
   fChain->SetBranchAddress("L1_SingleMu7_BptxAND_Prescl", &L1_SingleMu7_BptxAND_Prescl, &b_L1_SingleMu7_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMuBeamHalo", &L1_SingleMuBeamHalo, &b_L1_SingleMuBeamHalo);
   fChain->SetBranchAddress("L1_SingleMuBeamHalo_Prescl", &L1_SingleMuBeamHalo_Prescl, &b_L1_SingleMuBeamHalo_Prescl);
   fChain->SetBranchAddress("L1_SingleMuOpen", &L1_SingleMuOpen, &b_L1_SingleMuOpen);
   fChain->SetBranchAddress("L1_SingleMuOpen_Prescl", &L1_SingleMuOpen_Prescl, &b_L1_SingleMuOpen_Prescl);
   fChain->SetBranchAddress("L1_SingleMuOpen_BptxAND", &L1_SingleMuOpen_BptxAND, &b_L1_SingleMuOpen_BptxAND);
   fChain->SetBranchAddress("L1_SingleMuOpen_BptxAND_Prescl", &L1_SingleMuOpen_BptxAND_Prescl, &b_L1_SingleMuOpen_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMuOpen_NotBptxOR", &L1_SingleMuOpen_NotBptxOR, &b_L1_SingleMuOpen_NotBptxOR);
   fChain->SetBranchAddress("L1_SingleMuOpen_NotBptxOR_Prescl", &L1_SingleMuOpen_NotBptxOR_Prescl, &b_L1_SingleMuOpen_NotBptxOR_Prescl);
   fChain->SetBranchAddress("L1_SingleJet12_BptxAND", &L1_SingleJet12_BptxAND, &b_L1_SingleJet12_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet12_BptxAND_Prescl", &L1_SingleJet12_BptxAND_Prescl, &b_L1_SingleJet12_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet16", &L1_SingleJet16, &b_L1_SingleJet16);
   fChain->SetBranchAddress("L1_SingleJet16_Prescl", &L1_SingleJet16_Prescl, &b_L1_SingleJet16_Prescl);
   fChain->SetBranchAddress("L1_SingleJet16_BptxAND", &L1_SingleJet16_BptxAND, &b_L1_SingleJet16_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet16_BptxAND_Prescl", &L1_SingleJet16_BptxAND_Prescl, &b_L1_SingleJet16_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet200", &L1_SingleJet200, &b_L1_SingleJet200);
   fChain->SetBranchAddress("L1_SingleJet200_Prescl", &L1_SingleJet200_Prescl, &b_L1_SingleJet200_Prescl);
   fChain->SetBranchAddress("L1_SingleJet20_BptxAND", &L1_SingleJet20_BptxAND, &b_L1_SingleJet20_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet20_BptxAND_Prescl", &L1_SingleJet20_BptxAND_Prescl, &b_L1_SingleJet20_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet24_BptxAND", &L1_SingleJet24_BptxAND, &b_L1_SingleJet24_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet24_BptxAND_Prescl", &L1_SingleJet24_BptxAND_Prescl, &b_L1_SingleJet24_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet28_BptxAND", &L1_SingleJet28_BptxAND, &b_L1_SingleJet28_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet28_BptxAND_Prescl", &L1_SingleJet28_BptxAND_Prescl, &b_L1_SingleJet28_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet32_BptxAND", &L1_SingleJet32_BptxAND, &b_L1_SingleJet32_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet32_BptxAND_Prescl", &L1_SingleJet32_BptxAND_Prescl, &b_L1_SingleJet32_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet36", &L1_SingleJet36, &b_L1_SingleJet36);
   fChain->SetBranchAddress("L1_SingleJet36_Prescl", &L1_SingleJet36_Prescl, &b_L1_SingleJet36_Prescl);
   fChain->SetBranchAddress("L1_SingleJet36_BptxAND", &L1_SingleJet36_BptxAND, &b_L1_SingleJet36_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl", &L1_SingleJet36_BptxAND_Prescl, &b_L1_SingleJet36_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet40_BptxAND", &L1_SingleJet40_BptxAND, &b_L1_SingleJet40_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet40_BptxAND_Prescl", &L1_SingleJet40_BptxAND_Prescl, &b_L1_SingleJet40_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet44_BptxAND", &L1_SingleJet44_BptxAND, &b_L1_SingleJet44_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet44_BptxAND_Prescl", &L1_SingleJet44_BptxAND_Prescl, &b_L1_SingleJet44_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet48_BptxAND", &L1_SingleJet48_BptxAND, &b_L1_SingleJet48_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet48_BptxAND_Prescl", &L1_SingleJet48_BptxAND_Prescl, &b_L1_SingleJet48_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet52_BptxAND", &L1_SingleJet52_BptxAND, &b_L1_SingleJet52_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl", &L1_SingleJet52_BptxAND_Prescl, &b_L1_SingleJet52_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet60_BptxAND", &L1_SingleJet60_BptxAND, &b_L1_SingleJet60_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet60_BptxAND_Prescl", &L1_SingleJet60_BptxAND_Prescl, &b_L1_SingleJet60_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet68", &L1_SingleJet68, &b_L1_SingleJet68);
   fChain->SetBranchAddress("L1_SingleJet68_Prescl", &L1_SingleJet68_Prescl, &b_L1_SingleJet68_Prescl);
   fChain->SetBranchAddress("L1_SingleJet68_BptxAND", &L1_SingleJet68_BptxAND, &b_L1_SingleJet68_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet68_BptxAND_Prescl", &L1_SingleJet68_BptxAND_Prescl, &b_L1_SingleJet68_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet8_BptxAND", &L1_SingleJet8_BptxAND, &b_L1_SingleJet8_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet8_BptxAND_Prescl", &L1_SingleJet8_BptxAND_Prescl, &b_L1_SingleJet8_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJetC20_NotBptxOR", &L1_SingleJetC20_NotBptxOR, &b_L1_SingleJetC20_NotBptxOR);
   fChain->SetBranchAddress("L1_SingleJetC20_NotBptxOR_Prescl", &L1_SingleJetC20_NotBptxOR_Prescl, &b_L1_SingleJetC20_NotBptxOR_Prescl);
   fChain->SetBranchAddress("L1_SingleJetC32_NotBptxOR", &L1_SingleJetC32_NotBptxOR, &b_L1_SingleJetC32_NotBptxOR);
   fChain->SetBranchAddress("L1_SingleJetC32_NotBptxOR_Prescl", &L1_SingleJetC32_NotBptxOR_Prescl, &b_L1_SingleJetC32_NotBptxOR_Prescl);
   fChain->SetBranchAddress("L1_TOTEM_0", &L1_TOTEM_0, &b_L1_TOTEM_0);
   fChain->SetBranchAddress("L1_TOTEM_0_Prescl", &L1_TOTEM_0_Prescl, &b_L1_TOTEM_0_Prescl);
   fChain->SetBranchAddress("L1_TOTEM_1", &L1_TOTEM_1, &b_L1_TOTEM_1);
   fChain->SetBranchAddress("L1_TOTEM_1_Prescl", &L1_TOTEM_1_Prescl, &b_L1_TOTEM_1_Prescl);
   fChain->SetBranchAddress("L1_TOTEM_2", &L1_TOTEM_2, &b_L1_TOTEM_2);
   fChain->SetBranchAddress("L1_TOTEM_2_Prescl", &L1_TOTEM_2_Prescl, &b_L1_TOTEM_2_Prescl);
   fChain->SetBranchAddress("L1_TOTEM_3", &L1_TOTEM_3, &b_L1_TOTEM_3);
   fChain->SetBranchAddress("L1_TOTEM_3_Prescl", &L1_TOTEM_3_Prescl, &b_L1_TOTEM_3_Prescl);
   fChain->SetBranchAddress("L1_ZdcOR", &L1_ZdcOR, &b_L1_ZdcOR);
   fChain->SetBranchAddress("L1_ZdcOR_Prescl", &L1_ZdcOR_Prescl, &b_L1_ZdcOR_Prescl);
   fChain->SetBranchAddress("L1_ZdcOR_BptxAND", &L1_ZdcOR_BptxAND, &b_L1_ZdcOR_BptxAND);
   fChain->SetBranchAddress("L1_ZdcOR_BptxAND_Prescl", &L1_ZdcOR_BptxAND_Prescl, &b_L1_ZdcOR_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ZdcXOR", &L1_ZdcXOR, &b_L1_ZdcXOR);
   fChain->SetBranchAddress("L1_ZdcXOR_Prescl", &L1_ZdcXOR_Prescl, &b_L1_ZdcXOR_Prescl);
   fChain->SetBranchAddress("L1_ZdcXOR_BptxAND", &L1_ZdcXOR_BptxAND, &b_L1_ZdcXOR_BptxAND);
   fChain->SetBranchAddress("L1_ZdcXOR_BptxAND_Prescl", &L1_ZdcXOR_BptxAND_Prescl, &b_L1_ZdcXOR_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ZeroBias", &L1_ZeroBias, &b_L1_ZeroBias);
   fChain->SetBranchAddress("L1_ZeroBias_Prescl", &L1_ZeroBias_Prescl, &b_L1_ZeroBias_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_PreBPTX.v0", &L1Tech_BPTX_PreBPTX_v0, &b_L1Tech_BPTX_PreBPTX_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_PreBPTX.v0_Prescl", &L1Tech_BPTX_PreBPTX_v0_Prescl, &b_L1Tech_BPTX_PreBPTX_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_minus.v0", &L1Tech_BPTX_minus_v0, &b_L1Tech_BPTX_minus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_minus.v0_Prescl", &L1Tech_BPTX_minus_v0_Prescl, &b_L1Tech_BPTX_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_minus_AND_not_plus.v0", &L1Tech_BPTX_minus_AND_not_plus_v0, &b_L1Tech_BPTX_minus_AND_not_plus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_minus_AND_not_plus.v0_Prescl", &L1Tech_BPTX_minus_AND_not_plus_v0_Prescl, &b_L1Tech_BPTX_minus_AND_not_plus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus.v0", &L1Tech_BPTX_plus_v0, &b_L1Tech_BPTX_plus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus.v0_Prescl", &L1Tech_BPTX_plus_v0_Prescl, &b_L1Tech_BPTX_plus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_NOT_minus.v0", &L1Tech_BPTX_plus_AND_NOT_minus_v0, &b_L1Tech_BPTX_plus_AND_NOT_minus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_NOT_minus.v0_Prescl", &L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl, &b_L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0", &L1Tech_BPTX_plus_AND_minus_v0, &b_L1Tech_BPTX_plus_AND_minus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0_Prescl", &L1Tech_BPTX_plus_AND_minus_v0_Prescl, &b_L1Tech_BPTX_plus_AND_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus_instance1.v0", &L1Tech_BPTX_plus_AND_minus_instance1_v0, &b_L1Tech_BPTX_plus_AND_minus_instance1_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus_instance1.v0_Prescl", &L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl, &b_L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_OR_minus.v0", &L1Tech_BPTX_plus_OR_minus_v0, &b_L1Tech_BPTX_plus_OR_minus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_OR_minus.v0_Prescl", &L1Tech_BPTX_plus_OR_minus_v0_Prescl, &b_L1Tech_BPTX_plus_OR_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_quiet.v0", &L1Tech_BPTX_quiet_v0, &b_L1Tech_BPTX_quiet_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_quiet.v0_Prescl", &L1Tech_BPTX_quiet_v0_Prescl, &b_L1Tech_BPTX_quiet_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BRIL_bit28", &L1Tech_BRIL_bit28, &b_L1Tech_BRIL_bit28);
   fChain->SetBranchAddress("L1Tech_BRIL_bit28_Prescl", &L1Tech_BRIL_bit28_Prescl, &b_L1Tech_BRIL_bit28_Prescl);
   fChain->SetBranchAddress("L1Tech_BRIL_bit29", &L1Tech_BRIL_bit29, &b_L1Tech_BRIL_bit29);
   fChain->SetBranchAddress("L1Tech_BRIL_bit29_Prescl", &L1Tech_BRIL_bit29_Prescl, &b_L1Tech_BRIL_bit29_Prescl);
   fChain->SetBranchAddress("L1Tech_BRIL_bit30", &L1Tech_BRIL_bit30, &b_L1Tech_BRIL_bit30);
   fChain->SetBranchAddress("L1Tech_BRIL_bit30_Prescl", &L1Tech_BRIL_bit30_Prescl, &b_L1Tech_BRIL_bit30_Prescl);
   fChain->SetBranchAddress("L1Tech_BRIL_bit31", &L1Tech_BRIL_bit31, &b_L1Tech_BRIL_bit31);
   fChain->SetBranchAddress("L1Tech_BRIL_bit31_Prescl", &L1Tech_BRIL_bit31_Prescl, &b_L1Tech_BRIL_bit31_Prescl);
   fChain->SetBranchAddress("L1Tech_BRIL_bit32", &L1Tech_BRIL_bit32, &b_L1Tech_BRIL_bit32);
   fChain->SetBranchAddress("L1Tech_BRIL_bit32_Prescl", &L1Tech_BRIL_bit32_Prescl, &b_L1Tech_BRIL_bit32_Prescl);
   fChain->SetBranchAddress("L1Tech_BRIL_bit33", &L1Tech_BRIL_bit33, &b_L1Tech_BRIL_bit33);
   fChain->SetBranchAddress("L1Tech_BRIL_bit33_Prescl", &L1Tech_BRIL_bit33_Prescl, &b_L1Tech_BRIL_bit33_Prescl);
   fChain->SetBranchAddress("L1Tech_BRIL_bit34", &L1Tech_BRIL_bit34, &b_L1Tech_BRIL_bit34);
   fChain->SetBranchAddress("L1Tech_BRIL_bit34_Prescl", &L1Tech_BRIL_bit34_Prescl, &b_L1Tech_BRIL_bit34_Prescl);
   fChain->SetBranchAddress("L1Tech_BRIL_bit35", &L1Tech_BRIL_bit35, &b_L1Tech_BRIL_bit35);
   fChain->SetBranchAddress("L1Tech_BRIL_bit35_Prescl", &L1Tech_BRIL_bit35_Prescl, &b_L1Tech_BRIL_bit35_Prescl);
   fChain->SetBranchAddress("L1Tech_BRIL_bit36", &L1Tech_BRIL_bit36, &b_L1Tech_BRIL_bit36);
   fChain->SetBranchAddress("L1Tech_BRIL_bit36_Prescl", &L1Tech_BRIL_bit36_Prescl, &b_L1Tech_BRIL_bit36_Prescl);
   fChain->SetBranchAddress("L1Tech_BRIL_bit37", &L1Tech_BRIL_bit37, &b_L1Tech_BRIL_bit37);
   fChain->SetBranchAddress("L1Tech_BRIL_bit37_Prescl", &L1Tech_BRIL_bit37_Prescl, &b_L1Tech_BRIL_bit37_Prescl);
   fChain->SetBranchAddress("L1Tech_BRIL_bit38", &L1Tech_BRIL_bit38, &b_L1Tech_BRIL_bit38);
   fChain->SetBranchAddress("L1Tech_BRIL_bit38_Prescl", &L1Tech_BRIL_bit38_Prescl, &b_L1Tech_BRIL_bit38_Prescl);
   fChain->SetBranchAddress("L1Tech_BRIL_bit39", &L1Tech_BRIL_bit39, &b_L1Tech_BRIL_bit39);
   fChain->SetBranchAddress("L1Tech_BRIL_bit39_Prescl", &L1Tech_BRIL_bit39_Prescl, &b_L1Tech_BRIL_bit39_Prescl);
   fChain->SetBranchAddress("L1Tech_BRIL_bit40", &L1Tech_BRIL_bit40, &b_L1Tech_BRIL_bit40);
   fChain->SetBranchAddress("L1Tech_BRIL_bit40_Prescl", &L1Tech_BRIL_bit40_Prescl, &b_L1Tech_BRIL_bit40_Prescl);
   fChain->SetBranchAddress("L1Tech_BRIL_bit41", &L1Tech_BRIL_bit41, &b_L1Tech_BRIL_bit41);
   fChain->SetBranchAddress("L1Tech_BRIL_bit41_Prescl", &L1Tech_BRIL_bit41_Prescl, &b_L1Tech_BRIL_bit41_Prescl);
   fChain->SetBranchAddress("L1Tech_BRIL_bit42", &L1Tech_BRIL_bit42, &b_L1Tech_BRIL_bit42);
   fChain->SetBranchAddress("L1Tech_BRIL_bit42_Prescl", &L1Tech_BRIL_bit42_Prescl, &b_L1Tech_BRIL_bit42_Prescl);
   fChain->SetBranchAddress("L1Tech_BRIL_bit43", &L1Tech_BRIL_bit43, &b_L1Tech_BRIL_bit43);
   fChain->SetBranchAddress("L1Tech_BRIL_bit43_Prescl", &L1Tech_BRIL_bit43_Prescl, &b_L1Tech_BRIL_bit43_Prescl);
   fChain->SetBranchAddress("L1Tech_CASTOR_Gap.v0", &L1Tech_CASTOR_Gap_v0, &b_L1Tech_CASTOR_Gap_v0);
   fChain->SetBranchAddress("L1Tech_CASTOR_Gap.v0_Prescl", &L1Tech_CASTOR_Gap_v0_Prescl, &b_L1Tech_CASTOR_Gap_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_CASTOR_HaloMuon.v0", &L1Tech_CASTOR_HaloMuon_v0, &b_L1Tech_CASTOR_HaloMuon_v0);
   fChain->SetBranchAddress("L1Tech_CASTOR_HaloMuon.v0_Prescl", &L1Tech_CASTOR_HaloMuon_v0_Prescl, &b_L1Tech_CASTOR_HaloMuon_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_CASTOR_HighJet.v0", &L1Tech_CASTOR_HighJet_v0, &b_L1Tech_CASTOR_HighJet_v0);
   fChain->SetBranchAddress("L1Tech_CASTOR_HighJet.v0_Prescl", &L1Tech_CASTOR_HighJet_v0_Prescl, &b_L1Tech_CASTOR_HighJet_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_CASTOR_MediumJet.v0", &L1Tech_CASTOR_MediumJet_v0, &b_L1Tech_CASTOR_MediumJet_v0);
   fChain->SetBranchAddress("L1Tech_CASTOR_MediumJet.v0_Prescl", &L1Tech_CASTOR_MediumJet_v0_Prescl, &b_L1Tech_CASTOR_MediumJet_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_DT_GlobalOR.v0", &L1Tech_DT_GlobalOR_v0, &b_L1Tech_DT_GlobalOR_v0);
   fChain->SetBranchAddress("L1Tech_DT_GlobalOR.v0_Prescl", &L1Tech_DT_GlobalOR_v0_Prescl, &b_L1Tech_DT_GlobalOR_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HBHE_totalOR.v0", &L1Tech_HCAL_HBHE_totalOR_v0, &b_L1Tech_HCAL_HBHE_totalOR_v0);
   fChain->SetBranchAddress("L1Tech_HCAL_HBHE_totalOR.v0_Prescl", &L1Tech_HCAL_HBHE_totalOR_v0_Prescl, &b_L1Tech_HCAL_HBHE_totalOR_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_MMP_or_MPP.v1", &L1Tech_HCAL_HF_MMP_or_MPP_v1, &b_L1Tech_HCAL_HF_MMP_or_MPP_v1);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_MMP_or_MPP.v1_Prescl", &L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl, &b_L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_coincidence_PM.v2", &L1Tech_HCAL_HF_coincidence_PM_v2, &b_L1Tech_HCAL_HF_coincidence_PM_v2);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_coincidence_PM.v2_Prescl", &L1Tech_HCAL_HF_coincidence_PM_v2_Prescl, &b_L1Tech_HCAL_HF_coincidence_PM_v2_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_single_channel.v0", &L1Tech_HCAL_HF_single_channel_v0, &b_L1Tech_HCAL_HF_single_channel_v0);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_single_channel.v0_Prescl", &L1Tech_HCAL_HF_single_channel_v0_Prescl, &b_L1Tech_HCAL_HF_single_channel_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HO_totalOR.v0", &L1Tech_HCAL_HO_totalOR_v0, &b_L1Tech_HCAL_HO_totalOR_v0);
   fChain->SetBranchAddress("L1Tech_HCAL_HO_totalOR.v0_Prescl", &L1Tech_HCAL_HO_totalOR_v0_Prescl, &b_L1Tech_HCAL_HO_totalOR_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus1_Cosmics.v0", &L1Tech_RPC_TTU_RBplus1_Cosmics_v0, &b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus1_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus2_Cosmics.v0", &L1Tech_RPC_TTU_RBplus2_Cosmics_v0, &b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus2_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_barrel_Cosmics.v0", &L1Tech_RPC_TTU_barrel_Cosmics_v0, &b_L1Tech_RPC_TTU_barrel_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_barrel_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_pointing_Cosmics.v0", &L1Tech_RPC_TTU_pointing_Cosmics_v0, &b_L1Tech_RPC_TTU_pointing_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_pointing_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_TOTEM_0", &L1Tech_TOTEM_0, &b_L1Tech_TOTEM_0);
   fChain->SetBranchAddress("L1Tech_TOTEM_0_Prescl", &L1Tech_TOTEM_0_Prescl, &b_L1Tech_TOTEM_0_Prescl);
   fChain->SetBranchAddress("L1Tech_TOTEM_1", &L1Tech_TOTEM_1, &b_L1Tech_TOTEM_1);
   fChain->SetBranchAddress("L1Tech_TOTEM_1_Prescl", &L1Tech_TOTEM_1_Prescl, &b_L1Tech_TOTEM_1_Prescl);
   fChain->SetBranchAddress("L1Tech_TOTEM_2", &L1Tech_TOTEM_2, &b_L1Tech_TOTEM_2);
   fChain->SetBranchAddress("L1Tech_TOTEM_2_Prescl", &L1Tech_TOTEM_2_Prescl, &b_L1Tech_TOTEM_2_Prescl);
   fChain->SetBranchAddress("L1Tech_TOTEM_3", &L1Tech_TOTEM_3, &b_L1Tech_TOTEM_3);
   fChain->SetBranchAddress("L1Tech_TOTEM_3_Prescl", &L1Tech_TOTEM_3_Prescl, &b_L1Tech_TOTEM_3_Prescl);
   fChain->SetBranchAddress("L1Tech_ZDC_minus", &L1Tech_ZDC_minus, &b_L1Tech_ZDC_minus);
   fChain->SetBranchAddress("L1Tech_ZDC_minus_Prescl", &L1Tech_ZDC_minus_Prescl, &b_L1Tech_ZDC_minus_Prescl);
   fChain->SetBranchAddress("L1Tech_ZDC_plus", &L1Tech_ZDC_plus, &b_L1Tech_ZDC_plus);
   fChain->SetBranchAddress("L1Tech_ZDC_plus_Prescl", &L1Tech_ZDC_plus_Prescl, &b_L1Tech_ZDC_plus_Prescl);

   Notify();
}

Bool_t myTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void myTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t myTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#endif // #ifdef myTree_cxx
