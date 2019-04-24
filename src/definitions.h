#include "TRandom.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"
#include "TEfficiency.h"
#include "TVector2.h"
// ==============================================
// class for loading and retrieving GJets NLO
// weights
// ----------------------------------------------
enum NLO_weight_type {kLO,kNLO,kUpVar,kDnVar};

class NLO_weights{
 public:
    TFile* inputFile;
    TH1F *LO,*NLO,*UpVar,*DnVar;
    vector<TH1F*> histo;

    NLO_weights(){
        inputFile = new TFile("../data/kfactors.root");
        LO = (TH1F*) inputFile->Get("GJets_LO/inv_pt_G");
        NLO = (TH1F*) inputFile->Get("GJets_1j_NLO/nominal_G");
        UpVar = (TH1F*) inputFile->Get("GJets_1j_NLO/ren_up_G");
        DnVar = (TH1F*) inputFile->Get("GJets_1j_NLO/ren_down_G");
        histo.push_back(LO);
        histo.push_back(NLO);
        histo.push_back(UpVar);
        histo.push_back(DnVar);
    }

    double get(double pt, NLO_weight_type weightType=kLO){
        return histo[weightType]->GetBinContent(histo[weightType]->FindBin(pt));
    }
};

// ==============================================

// constants
// ==============================================
//double bbtagCut = 0.4;
//TFile* puWeightFile = new TFile("../data/PileupHistograms_0121_69p2mb_pm4p6.root");
//TH1F* puWeightHist = (TH1F*) puWeightFile->Get("pu_weights_down");
//NLO_weights NLOw;
// ==============================================

double CalcdPhi( double phi1 , double phi2 ){

  double dPhi = phi1-phi2;
  if( dPhi < -TMath::Pi() )
    dPhi += 2*TMath::Pi() ;
  if( dPhi > TMath::Pi() )
    dPhi -= 2*TMath::Pi() ;
  return fabs(dPhi);

}

template<typename ntupleType>void ntupleBranchStatus(ntupleType* ntuple){
  ntuple->fChain->SetBranchStatus("*",0);
  ntuple->fChain->SetBranchStatus("NMuons",1);
  ntuple->fChain->SetBranchStatus("NElectrons",1);
  ntuple->fChain->SetBranchStatus("iso*Tracks",1);
  ntuple->fChain->SetBranchStatus("Jets",1);
  ntuple->fChain->SetBranchStatus("DeltaPhi*",1);
  ntuple->fChain->SetBranchStatus("HT*",1);
  ntuple->fChain->SetBranchStatus("NJets",1);
  ntuple->fChain->SetBranchStatus("BTags*",1);
  ntuple->fChain->SetBranchStatus("MHT",1);
  ntuple->fChain->SetBranchStatus("MHTPhi",1);
  ntuple->fChain->SetBranchStatus("Weight",1);
  ntuple->fChain->SetBranchStatus("pu*",1);
  ntuple->fChain->SetBranchStatus("TriggerPass",1);
  ntuple->fChain->SetBranchStatus("Photons*",1);
  ntuple->fChain->SetBranchStatus("madMinPhotonDeltaR",1);
  ntuple->fChain->SetBranchStatus("GenParticles*",1);
  ntuple->fChain->SetBranchStatus("madHT",1);
  ntuple->fChain->SetBranchStatus("*Filter",1);
  ntuple->fChain->SetBranchStatus("PFCaloMETRatio",1);
  ntuple->fChain->SetBranchStatus("TrueNumInteractions",1);
  ntuple->fChain->SetBranchStatus("NVtx",1);
  ntuple->fChain->SetBranchStatus("JetID",1);
  ntuple->fChain->SetBranchStatus("ZCandidates",1);
  ntuple->fChain->SetBranchStatus("NonPrefiringProb",1);
  ntuple->fChain->SetBranchStatus("NonPrefiringProbUp",1);
  ntuple->fChain->SetBranchStatus("NonPrefiringProbDn",1);
  ntuple->fChain->SetBranchStatus("HTRatioDPhiFilter",1);
  ntuple->fChain->SetBranchStatus("RunNum",1);
  ntuple->fChain->SetBranchStatus("LumiBlockNum",1);
  ntuple->fChain->SetBranchStatus("EvtNum",1);
  ntuple->fChain->SetBranchStatus("Jets_jecFactor",1);
}
/******************************************************************/
/* - - - - - - - - - - - - cut flow function - - - - - - - - - -  */
/******************************************************************/
template<typename ntupleType> bool cutFlow_none(ntupleType* ntuple){
    return true;
}

template<typename ntupleType> bool cutFlow_HT300(ntupleType* ntuple){
    return ntuple->HT>300.;
}

template<typename ntupleType> bool cutFlow_MHT300(ntupleType* ntuple){
    return ntuple->MHT>300.;
}

template<typename ntupleType> bool cutFlow_NJets2plus(ntupleType* ntuple){
    return ntuple->NJets>=2;
}

template<typename ntupleType> bool cutFlow_leptonVeto(ntupleType* ntuple){
  return (ntuple->NMuons==0 &&
	  ntuple->NElectrons==0 &&
	  ntuple->isoElectronTracks==0 &&
	  ntuple->isoMuonTracks==0 && 
	  ntuple->isoPionTracks==0);
}

template<typename ntupleType> bool cutFlow_onePhoton(ntupleType* ntuple){
  return ( ntuple->Photons->size()==1
	   && isPromptPhoton(ntuple)
	   && ntuple->Photons_fullID->at(0)==1
	   && ntuple->Photons_hasPixelSeed->at(0)==0 );
}

template<typename ntupleType> bool onePhoton(ntupleType* ntuple){
  return ( ntuple->Photons->size()==1
	   && ntuple->Photons_hasPixelSeed->at(0)==0 );
}

template<typename ntupleType> bool cutFlow_onePhoton_data(ntupleType* ntuple){
  return ( ntuple->Photons->size()==1
	   && ntuple->Photons_fullID->at(0)==1
	   && ntuple->Photons_hasPixelSeed->at(0)==0 );
}

template<typename ntupleType> bool cutFlow_trigger(ntupleType* ntuple){
  return ntuple->TriggerPass->at(140)==1 || ntuple->TriggerPass->at(141)==1 ;
}

template<typename ntupleType> bool cutFlow_photonPt200(ntupleType* ntuple){
    return ntuple->Photons->at(0).Pt()>=200.;
}
template<typename ntupleType> bool cutFlow_minDR(ntupleType* ntuple){
    return ntuple->madMinPhotonDeltaR>0.4;
}
template<typename ntupleType> bool cutFlow_HTbin1(ntupleType* ntuple){
    return ntuple->HT>=300. && ntuple->HT<500.;
}
template<typename ntupleType> bool cutFlow_MHTbin1(ntupleType* ntuple){
    return ntuple->MHT>=300. && ntuple->MHT<350;
}
template<typename ntupleType> bool cutFlow_nJetsTwo(ntupleType* ntuple){
    return ntuple->NJets==2;
}
template<typename ntupleType> bool cutFlow_btagsZero(ntupleType* ntuple){
    return ntuple->BTagsDeepCSV==0;
}
template<typename ntupleType> bool cutFlow_filters(ntupleType* ntuple){
    return ( ntuple->HTRatioDPhiFilter == 1
	     && ntuple->PFCaloMETRatio < 5. 
	     && ntuple->globalSuperTightHalo2016Filter==1
	     && ntuple->HBHENoiseFilter==1
	     && ntuple->HBHEIsoNoiseFilter==1
	     && ntuple->eeBadScFilter==1
	     && ntuple->EcalDeadCellTriggerPrimitiveFilter == 1
	     && ntuple->BadChargedCandidateFilter == 1
	     && ntuple->BadPFMuonFilter == 1
	     && ntuple->NVtx>0
	     && ntuple->JetID == 1
	     );
}

/******************************************/
/* custom weights			   */
/******************************************/
//template<typename ntupleType> double customPUweights(ntupleType* ntuple){
//    int nVtx = ntuple->TrueNumInteractions;
//    return puWeightHist->GetBinContent(puWeightHist->GetXaxis()->FindBin(min(ntuple->TrueNumInteractions,puWeightHist->GetBinLowEdge(puWeightHist->GetNbinsX()+1))));
//}


template<typename ntupleType> double dRweights(ntupleType* ntuple){
    return 1. / ( (min(ntuple->HT, 900.0) *(0.0001906)) + 0.8477);
}

template<typename ntupleType> double GJets0p4Weights(ntupleType* ntuple){
    if( ntuple->madHT > 100. && ntuple->madHT < 200. )
        return 5391./5000.;
    else if( ntuple->madHT > 200. && ntuple->madHT < 400. )
        return 1168./1079.;
    else if( ntuple->madHT > 400. && ntuple->madHT < 600. )
        return 132.5/125.9;
    else if( ntuple->madHT > 600. )
        return 44.05/43.36;
    else
        return 1.;
}

TRandom randGen;

double truncGauss(double mean, double errUp, double errDown){
    double result=randGen.Gaus()*(errUp/2.+errDown/2.)+mean;
    if( result > 1. )
        return 1. ;
    else if( result < 0. )
        return 0. ;
    else
        return result;
}

template<typename ntupleType> double photonTriggerWeightRand( ntupleType* ntuple ){

  if( ntuple->Photons->size() == 0 ) return -9999.;

  double EBtrigger[5]={0.991,0.988,0.988,0.986,0.978};
  double ECtrigger[5]={0.989,0.985,0.986,0.991,1.000};
  double EBtriggerErrDown[5]={0.001,0.001,0.001,0.001,0.005};
  double ECtriggerErrDown[5]={0.001,0.002,0.001,0.002,0.014};
  double EBtriggerErrUp[5]={0.001,0.001,0.001,0.001,0.004};
  double ECtriggerErrUp[5]={0.001,0.001,0.001,0.002,0.000};

  double MHT = ntuple->MHT;
  if( ntuple->Photons_isEB->at(0) ){
    if( MHT > 250. && MHT < 300. ){
      return truncGauss(EBtrigger[0],EBtriggerErrUp[0],EBtriggerErrDown[0]);
    }else if( MHT > 300. && MHT < 350. ){
        return truncGauss(EBtrigger[1],EBtriggerErrUp[1],EBtriggerErrDown[1]);
    }else if( MHT > 350.  && MHT < 500. ){
        return truncGauss(EBtrigger[2],EBtriggerErrUp[2],EBtriggerErrDown[2]);
    }else if( MHT > 500.  && MHT < 750. ){
        return truncGauss(EBtrigger[3],EBtriggerErrUp[3],EBtriggerErrDown[3]);
    }else if( MHT > 750. ){
        return truncGauss(EBtrigger[4],EBtriggerErrUp[4],EBtriggerErrDown[4]);
    }else
        return 1.;
  }else{
    if( MHT > 250. && MHT < 300. ){
        return truncGauss(ECtrigger[0],ECtriggerErrUp[0],ECtriggerErrDown[0]);
    }else if( MHT > 300. && MHT < 350. ){
        return truncGauss(ECtrigger[1],ECtriggerErrUp[1],ECtriggerErrDown[1]);
    }else if( MHT > 350. && MHT < 500. ){
        return truncGauss(ECtrigger[2],ECtriggerErrUp[2],ECtriggerErrDown[2]);
    }else if( MHT > 500.  && MHT < 750. ){
        return truncGauss(ECtrigger[3],ECtriggerErrUp[3],ECtriggerErrDown[3]);
    }else if( MHT > 750. ){
        return truncGauss(ECtrigger[4],ECtriggerErrUp[4],ECtriggerErrDown[4]);
    }else
      return 1.;
  }

}

template<typename ntupleType> double photonTriggerWeight( ntupleType* ntuple ){

  if( ntuple->Photons->size() == 0 ) return -9999.;

  double EBtrigger[5]={0.991,0.988,0.988,0.986,0.978};
  double ECtrigger[5]={0.989,0.985,0.986,0.991,1.000};

  double MHT = ntuple->MHT;
  if( ntuple->Photons_isEB->at(0) ){
    if( MHT > 250. && MHT < 300. ){
      return EBtrigger[0];
    }else if( MHT > 300. && MHT < 350. ){
      return EBtrigger[1];
    }else if( MHT > 350.  && MHT < 500. ){
      return EBtrigger[2];
    }else if( MHT > 500.  && MHT < 750. ){
      return EBtrigger[3];
    }else if( MHT > 750. ){
      return EBtrigger[4];
    }else
      return 1.;
  }else{
    if( MHT > 250. && MHT < 300. ){
      return ECtrigger[0];
    }else if( MHT > 300. && MHT < 350. ){
      return ECtrigger[1];
    }else if( MHT > 350. && MHT < 500. ){
      return ECtrigger[2];
    }else if( MHT > 500.  && MHT < 750. ){
      return ECtrigger[3];
    }else if( MHT > 750. ){
      return ECtrigger[4];
    }else
      return 1.;
  }

}

//////////////////////
//////////////////////
//////////////////////
//   DY functions   //
//////////////////////
//////////////////////
//////////////////////
template<typename ntupleType> bool isDY(ntupleType* ntuple){
  for( int i = 0 ; i < ntuple->ZCandidates->size() ; i++ ){
    if( fabs( ntuple->ZCandidates->at(i).M() - 91.19 ) < 10. ) return true;
  }
  return false;
}

//////////////////////
//////////////////////
//////////////////////
// Photon functions //
//////////////////////
//////////////////////
//////////////////////
template<typename ntupleType> int isPromptPhoton(ntupleType* ntuple){
  return ntuple->Photons_nonPrompt->at(0)==0;
}

template<typename ntupleType> double photonPt(ntupleType* ntuple){
  return ntuple->Photons->at(0).Pt();
}

template<typename ntupleType> double photonEta(ntupleType* ntuple){
  return ntuple->Photons->at(0).Eta();
}

template<typename ntupleType> double photonPhi(ntupleType* ntuple){
  return ntuple->Photons->at(0).Phi();
}



template<typename ntupleType> double jetPt(ntupleType* ntuple){
  return ntuple->Jets->at(0).Pt();
}

template<typename ntupleType> double jetPhi(ntupleType* ntuple){
  return ntuple->Jets->at(0).Phi();
}

template<typename ntupleType> double jetEta(ntupleType* ntuple){
  return ntuple->Jets->at(0).Eta();
}




template<typename ntupleType> double jetSubleadPt(ntupleType* ntuple){
  return ntuple->Jets->at(1).Pt();
}

template<typename ntupleType> double jetSubleadPhi(ntupleType* ntuple){
  return ntuple->Jets->at(1).Phi();
}

template<typename ntupleType> double jetSubleadEta(ntupleType* ntuple){
  return ntuple->Jets->at(1).Eta();
}

template<typename ntupleType> double photonSieie(ntupleType* ntuple){
  return ntuple->Photons_sigmaIetaIeta->at(0);
}

template<typename ntupleType> double photonIsoChrg(ntupleType* ntuple){
  return ntuple->Photons_pfChargedIsoRhoCorr->at(0);
}

template<typename ntupleType> double photonIsoGam(ntupleType* ntuple){
  return ntuple->Photons_pfGammaIsoRhoCorr->at(0);
}

template<typename ntupleType> double photonIsoNeu(ntupleType* ntuple){
  return ntuple->Photons_pfNeutralIsoRhoCorr->at(0);
}

template<typename ntupleType> double photonHoverE(ntupleType* ntuple){
  return ntuple->Photons_hadTowOverEM->at(0);
}

//////////////////////
//////////////////////
//////////////////////
// Lepton functions //
//////////////////////
//////////////////////
//////////////////////


///////////////////////////
// RECO Muon definitions //
///////////////////////////
template<typename ntupleType> int numMuons(ntupleType* ntuple){
  return ntuple->Muons->size();
}

template<typename ntupleType> double fillNumMuons(ntupleType* ntuple){
  return double(ntuple->Muons->size());
}

template<typename ntupleType> int numElectrons(ntupleType* ntuple){
  return ntuple->Electrons->size();
}

template<typename ntupleType> double fillLepPt(ntupleType* ntuple){
  if( ntuple->Electrons->size() > 0 && ntuple->Muons->size() == 0 )
    return ntuple->Electrons->at(0)->Pt();
  else if( ntuple->Electrons->size() ==0 && ntuple->Muons->size() > 0 )
    return ntuple->Muons->at(0)->Pt();
  else
    return -9999.;
}

template<typename ntupleType> double fillLepActivity(ntupleType* ntuple){
  if( ntuple->Electrons->size() > 0 && ntuple->Muons->size() == 0 )
    return ntuple->Electrons_MT2Activity->at(0);
  else if( ntuple->Electrons->size() ==0 && ntuple->Muons->size() > 0 )
    return ntuple->Muons_MT2Activity->at(0);
  else
    return -9999.;
}

//////////////////////
// END LEPTON STUFF //
//////////////////////

/**********************************************/
/* - - - - - - - PHOTON STUFF - - - - - - - - */
/**********************************************/

template<typename ntupleType> double fillPhotonChargedIso(ntupleType* ntuple){
  return 0.;
}
/**********************************************/
/* - - - - - - END PHOTON STUFF - - - - - - - */
/**********************************************/

////////////////////////////////////////////////////////////
// - - - - - - - - EVENT LEVEL VARIABLES - - - - - - - -  //
////////////////////////////////////////////////////////////
template<typename ntupleType> double fillNumVertices(ntupleType* ntuple){
    return ntuple->NVtx;
}

template<typename ntupleType> double fillMadMinPhotonDeltaR(ntupleType* ntuple){
  return ntuple->madMinPhotonDeltaR;
}

template<typename ntupleType> double fillPythiaMinPhotonDeltaR(ntupleType* ntuple){

  if( ntuple->Photons->size() == 0 ) return 999.;
  //cout << "num reco photons: " << ntuple->Photons->size() << endl;

  double minDR_genRecoPhotons = 99999.;
  int genPhotonIndex = -1;
  double dRtemp;
  for( int i = 0 ; i < ntuple->GenParticles->size() ; i++){
    if( ntuple->GenParticles_PdgId->at(i) != 22 ) continue;
    dRtemp = ntuple->GenParticles->at(i).DeltaR(ntuple->Photons->at(0));
    if( dRtemp < minDR_genRecoPhotons ){
      minDR_genRecoPhotons = dRtemp;
      genPhotonIndex = i;
    }
  }
  if( genPhotonIndex == -1 ) return 9999.;

  double minDR_partonPhoton = 999999.;

  for( int i = 0 ; i < ntuple->GenParticles->size() ; i++){
    if( abs(ntuple->GenParticles_PdgId->at(i)) > 6 && abs(ntuple->GenParticles_PdgId->at(i)) != 21 ) continue;
    if( ntuple->GenParticles->at(i).Pt() == 0 ) continue;
    //cout << "gen pdg id: " << ntuple->GenParticles_PdgId->at(i) << endl;
    //cout << "gen status: " << ntuple->GenParticles_Status->at(i) << endl;
    dRtemp = ntuple->GenParticles->at(i).DeltaR(ntuple->GenParticles->at(genPhotonIndex));
    //cout << "dRtemp: " << dRtemp << endl;
    if( dRtemp < minDR_partonPhoton )
      minDR_partonPhoton = dRtemp;
  }

  return minDR_partonPhoton;
}

template<typename ntupleType> double fillDeltaPhi1(ntupleType* ntuple){
  return ntuple->DeltaPhi1;
}

template<typename ntupleType> double fillRecoPhotonDeltaR(ntupleType* ntuple){
    //int minIndex=-1;
    double minDR=999999.;
    if( ntuple->Photons->size()<1 ) return -9999.;
    for( int iJet = 0 ; iJet < ntuple->Jets->size() ; iJet++ ){
        //cout << "temp DR: " << ntuple->Jets->at(iJet).DeltaR(ntuple->Photons->at(0)) << endl;
        if( ntuple->Jets->at(iJet).DeltaR(ntuple->Photons->at(0)) < minDR ){
            minDR = ntuple->Jets->at(iJet).DeltaR(ntuple->Photons->at(0));
            //minIndex = iJet;
        }
    }
    //cout << "minDR: " << minDR << endl;
    return minDR;
}

template<typename ntupleType> double fillDeltaPhi2(ntupleType* ntuple){
  return ntuple->DeltaPhi2;
}

template<typename ntupleType> double fillDeltaPhi3(ntupleType* ntuple){
  return ntuple->DeltaPhi3;
}

template<typename ntupleType> double fillDeltaPhi4(ntupleType* ntuple){
  return ntuple->DeltaPhi4;
}

template<typename ntupleType> double fillHT(ntupleType* ntuple){
  return ntuple->HT;
}

template<typename ntupleType> double fillMHT(ntupleType* ntuple){
  return ntuple->MHT;
}

template<typename ntupleType> double fillOne(ntupleType* ntuple){
  return 1.;
}

template<typename ntupleType> double fillNJets(ntupleType* ntuple){
  return ntuple->NJets;
}

template<typename ntupleType> double fillBTags(ntupleType* ntuple){
  return ntuple->BTagsDeepCSV;
}

  /////////////////
 // OTHER STUFF //
/////////////////
template<typename ntupleType> double fillAnalysisBins(ntupleType* ntuple){
  double MHT = ntuple->MHT;
  double HT = ntuple->HT;

  if( MHT > 300. && MHT < 600. ){
    if( HT > 300. && HT < 1000. ){
      return 1.;
    }else if( HT > 1000. && HT < 2000. ){
      return 2.;
    }else if( HT > 2000. ){
      return 3.;
    }else
      return -1.;
  }else if( MHT > 600. && MHT < 1000. ){
    if( HT > 600. && HT < 1000. ){
      return 4.;
    }else if( HT > 1000. && HT < 2000. ){
      return 5.;
    }else if( HT > 2000. ){
      return 6.;
    }else
      return -1.;
  }else if( MHT > 1000. ){
    if( HT > 1000. && HT < 2000. ){
      return 7.;
    }else if( HT > 2000. ){
      return 8.;
    }else
      return -1.;
  }else
    return -1.;
}

template<typename ntupleType> double fillRA2b10Bins(ntupleType* ntuple){
  double MHT = ntuple->MHT;
  double HT = ntuple->HT;
  double HT5 = ntuple->HT5;
  double DeltaPhi1 = ntuple->DeltaPhi1;

if( HT >= MHT && DeltaPhi1 >= ( (1.025*(HT5/HT)) - 0.5875) ){

  if( MHT > 300. && MHT < 350. ){
    if( HT > 300. && HT < 600. ){
      if(ntuple->NJets>=8)
        return -999999.;
      else
        return 1.;
    }else if( HT > 600. && HT < 1200. ){
      if(ntuple->NJets>=8)
        return 1.;
      else
        return 2.;
    }else if( HT > 1200. ){
      if(ntuple->NJets>=8)
        return 2.;
      else
        return 3.;
    }else
      return -999999.;
  }else if( MHT > 350. && MHT < 600. ){
    if( HT > 350. && HT < 600. ){
      if(ntuple->NJets>=8)
        return -999999.;
      else
        return 4.;
    }else if( HT > 600. && HT < 1200. ){
      if(ntuple->NJets>=8)
        return 3.;
      else
        return 5.;
    }else if( HT > 1200. ){
      if(ntuple->NJets>=8)
        return 4.;
      else
        return 6.;
    }else
      return -999999.;
  }else if( MHT > 600. && MHT < 850. ){
    if( HT >600. && HT < 1200. ){
      if(ntuple->NJets>=8)
        return 5.;
      else
        return 7.;
    }else if( HT > 1200. ){
      if(ntuple->NJets>=8)
        return 6.;
      else
        return 8.;
    }else
      return -999999.;
  }else if( MHT > 850. ){
    if( HT > 850. && HT < 1700. ){
      if(ntuple->NJets>=8)
        return 7.;
      else
        return 9.;
    }else if( HT > 1700. ){
      if(ntuple->NJets>=8)
        return 8.;
      else
        return 10.;
    }else
      return -999999.;
  }else
    return -999999.;
}
else
   return -999999.;
}

template<typename ntupleType> double fillRA2b46Bins( ntupleType* ntuple ){

  int BTags = int(ntuple->BTagsDeepCSV);
  if( BTags != 0 ) return -999999.;
  int NJets = int(ntuple->NJets);

  if( NJets >= 2 && NJets <= 3 ){
    return fillRA2b10Bins(ntuple);
  }else if( NJets >= 4 && NJets <= 5 ){
    return 10.+fillRA2b10Bins(ntuple);
  }else if( NJets >= 6 && NJets <= 7 ){
    return 20.+fillRA2b10Bins(ntuple);
  }else if( NJets >= 8 && NJets <= 9 ){
    return 30.+fillRA2b10Bins(ntuple);
  }else if( NJets >= 10 ){
    return 38.+fillRA2b10Bins(ntuple);
  }else
    return -999999.;
}

template<typename ntupleType> double fillRA2bNJet2Bins( ntupleType* ntuple ){

  int BTags = int(ntuple->BTagsDeepCSV);
  int NJets = int(ntuple->NJets);

  if( NJets == 2 ){
    if( BTags == 0 )
      return fillRA2b10Bins(ntuple);
    else if( BTags == 1 )
      return 10.+fillRA2b10Bins(ntuple);
    else if( BTags >= 2 )
      return 20.+fillRA2b10Bins(ntuple);
  }else
    return -999999.;
}


template<typename ntupleType> double fillRA2b13Bins(ntupleType* ntuple){
  double MHT = ntuple->MHT;
  double HT = ntuple->HT;
  double HT5 = ntuple->HT5;
  double DeltaPhi1 = ntuple->DeltaPhi1;

if( HT >= MHT && DeltaPhi1 >= ( (1.025*(HT5/HT)) - 0.5875) ){

  if( MHT > 250. && MHT < 300. ){
    if( HT > 300. && HT < 600. )
      if( ntuple->NJets>=8 )
        return -999999.;
      else
        return 1.;
    else if( HT > 600. && HT < 1200. )
      if( ntuple->NJets>=8 )
        return 1.;
      else
        return 2.;
    else if( HT > 1200. )
      if( ntuple->NJets>=8 )
        return 2.;
      else
        return 3.;
  }else if( MHT > 300. && MHT < 350. ){
    if( HT > 300. && HT < 600. ){
      if(ntuple->NJets>=8)
        return -999999.;
      else
        return 4.;
    }else if( HT > 600. && HT < 1200. ){
      if(ntuple->NJets>=8)
        return 3.;
      else
        return 5.;
    }else if( HT > 1200. ){
      if(ntuple->NJets>=8)
        return 4.;
      else
        return 6.;
    }else
      return -999999.;
  }else if( MHT > 350. && MHT < 600. ){
    if( HT > 350. && HT < 600. ){
      if(ntuple->NJets>=8)
        return -999999.;
      else
        return 7.;
    }else if( HT > 600. && HT < 1200. ){
      if(ntuple->NJets>=8)
        return 5.;
      else
        return 8.;
    }else if( HT > 1200. ){
      if(ntuple->NJets>=8)
        return 6.;
      else
        return 9.;
    }else
      return -999999.;
  }else if( MHT > 600. && MHT < 850. ){
    if( HT > 600. && HT < 1200. ){
      if(ntuple->NJets>=8)
        return 7.;
      else
        return 10.;
    }else if( HT > 1200. ){
      if(ntuple->NJets>=8)
        return 8.;
      else
        return 11.;
    }else
      return -999999.;
  }else if( MHT > 850. ){
    if( HT > 850. && HT < 1700. ){
      if(ntuple->NJets>=8)
        return 9.;
      else
        return 12.;
    }else if( HT > 1700. ){
      if(ntuple->NJets>=8)
        return 10.;
      else
        return 13.;
    }else
      return -999999.;
  }else
    return -999999.;
 }else
   return -999999.;

}
template<typename ntupleType> double fillRA2b59Bins( ntupleType* ntuple ){

  int BTags = int(ntuple->BTagsDeepCSV);
  if( BTags != 0 ) return -999999.;
  int NJets = int(ntuple->NJets);

  if(NJets >= 2 && NJets <= 3){
      return fillRA2b13Bins(ntuple);
  }else if( NJets >= 4 && NJets <= 5 ){
    return 13.+fillRA2b13Bins(ntuple);
  }else if( NJets >= 6 && NJets <= 7 ){
    return 26.+fillRA2b13Bins(ntuple);
  }else if( NJets >= 8 && NJets <= 9 ){
    return 39.+fillRA2b13Bins(ntuple);
  }else if( NJets >= 10 ){
    return 49.+fillRA2b13Bins(ntuple);
  }else
    return -999999.;
}

template<typename ntupleType> double fillRA2b174Bins( ntupleType* ntuple ){

  int BTags = int(ntuple->BTagsDeepCSV);
  int NJets = int(ntuple->NJets);

  if( NJets >= 2 && NJets <= 3 ){
    if( BTags == 0 ){
      return fillRA2b10Bins(ntuple);
    }else if( BTags == 1 ){
      return 10.+fillRA2b10Bins(ntuple);
    }else if( BTags>=2 ){
      return 20.+fillRA2b10Bins(ntuple);
    }else
      return -999999.;
  }else if(NJets >= 4 && NJets <= 5 ){
    if( BTags == 0 )
      return 30.+fillRA2b10Bins(ntuple);
    else if( BTags == 1 )
      return 40.+fillRA2b10Bins(ntuple);
    else if( BTags == 2 )
      return 50.+fillRA2b10Bins(ntuple);
    else if( BTags >= 3 )
      return 60.+fillRA2b10Bins(ntuple);
  }else if( NJets >= 6 && NJets <= 7 ){
    if( BTags == 0 )
      return 70.+fillRA2b10Bins(ntuple);
    else if( BTags == 1 )
      return 80.+fillRA2b10Bins(ntuple);
    else if( BTags == 2 )
      return 90.+fillRA2b10Bins(ntuple);
    else if( BTags >= 3 )
      return 100.+fillRA2b10Bins(ntuple);
  }else if( NJets >= 8 && NJets <= 9 ){
    if( BTags == 0 )
      return 110.+fillRA2b10Bins(ntuple);
    else if( BTags == 1 )
      return 118.+fillRA2b10Bins(ntuple);
    else if( BTags == 2 )
      return 126.+fillRA2b10Bins(ntuple);
    else if( BTags >= 3 )
      return 134.+fillRA2b10Bins(ntuple);
  }else if( NJets >= 10 ){
    if( BTags == 0 )
      return 142.+fillRA2b10Bins(ntuple);
    else if( BTags == 1 )
      return 150.+fillRA2b10Bins(ntuple);
    else if( BTags == 2 )
      return 158.+fillRA2b10Bins(ntuple);
    else if( BTags >= 3 )
      return 166.+fillRA2b10Bins(ntuple);
  }else
    return -999999.;
}


template<typename ntupleType> bool ptBinCut(double pt , int ithBin){
  if( ithBin > 5 ) return false;
  double ptCut[6] = {300.,400.,500.,700.,1000.,999999.};
  return pt>ptCut[ithBin] && pt<ptCut[ithBin+1];
}

template<typename ntupleType> bool RA2bBaselineCut(ntupleType* ntuple){

  double DeltaPhi1 = ntuple->DeltaPhi1;
  double DeltaPhi2 = ntuple->DeltaPhi2;
  double DeltaPhi3 = ntuple->DeltaPhi3;
  double DeltaPhi4 = ntuple->DeltaPhi4;

  double HT = ntuple->HT;
  double MHT = ntuple->MHT;
  int NJets = ntuple->NJets;

  return ( (NJets==2 && DeltaPhi1>0.5 && DeltaPhi2>0.5)
           || (NJets == 3 && DeltaPhi1 > 0.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3)
           || (NJets > 3 && DeltaPhi1 > 0.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3 && DeltaPhi4 > 0.3 ) )
    && MHT>300. && HT>300.
    && cutFlow_leptonVeto(ntuple)
    && cutFlow_filters(ntuple)
      ;

}

template<typename ntupleType> bool RA2bBaselinePhotonCut(ntupleType* ntuple){

  double DeltaPhi1 = ntuple->DeltaPhi1;
  double DeltaPhi2 = ntuple->DeltaPhi2;
  double DeltaPhi3 = ntuple->DeltaPhi3;
  double DeltaPhi4 = ntuple->DeltaPhi4;

  double HT = ntuple->HT;
  double MHT = ntuple->MHT;
  int NJets = ntuple->NJets;

  return ( (NJets==2 && DeltaPhi1>0.5 && DeltaPhi2>0.5)
           || (NJets == 3 && DeltaPhi1 > 0.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3)
           || (NJets > 3 && DeltaPhi1 > 0.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3 && DeltaPhi4 > 0.3 ) )
    && MHT>300. && HT>300.
    && cutFlow_leptonVeto(ntuple)
    && onePhoton(ntuple)
    && cutFlow_photonPt200(ntuple)
    && cutFlow_filters(ntuple)
      ;

}

template<typename ntupleType> bool RA2bBaselineWideCut(ntupleType* ntuple){

  double DeltaPhi1 = ntuple->DeltaPhi1;
  double DeltaPhi2 = ntuple->DeltaPhi2;
  double DeltaPhi3 = ntuple->DeltaPhi3;
  double DeltaPhi4 = ntuple->DeltaPhi4;

  double HT = ntuple->HT;
  double MHT = ntuple->MHT;
  int NJets = ntuple->NJets;

  return ( (NJets==2 && DeltaPhi1>0.5 && DeltaPhi2>0.5)
           || (NJets == 3 && DeltaPhi1 > 0.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3)
           || (NJets > 3 && DeltaPhi1 > 0.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3 && DeltaPhi4 > 0.3 ) )
      && MHT>250. && HT>300.
      && cutFlow_filters<ntupleType>(ntuple)
      ;

}

template<typename ntupleType> bool RA2bBaselinePhotonWideCut(ntupleType* ntuple){

  double DeltaPhi1 = ntuple->DeltaPhi1;
  double DeltaPhi2 = ntuple->DeltaPhi2;
  double DeltaPhi3 = ntuple->DeltaPhi3;
  double DeltaPhi4 = ntuple->DeltaPhi4;

  double HT = ntuple->HT;
  double MHT = ntuple->MHT;
  int NJets = ntuple->NJets;

  return ( (NJets==2 && DeltaPhi1>0.5 && DeltaPhi2>0.5)
           || (NJets == 3 && DeltaPhi1 > 0.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3)
           || (NJets > 3 && DeltaPhi1 > 0.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3 && DeltaPhi4 > 0.3 ) )
      && MHT>250. && HT>300.
    && cutFlow_leptonVeto(ntuple)
    && onePhoton(ntuple)
    && cutFlow_photonPt200(ntuple)
    && cutFlow_filters(ntuple)
      ;

}

template<typename ntupleType> bool RA2bLDPBaselinePhotonCut(ntupleType* ntuple){

  double DeltaPhi1 = ntuple->DeltaPhi1;
  double DeltaPhi2 = ntuple->DeltaPhi2;
  double DeltaPhi3 = ntuple->DeltaPhi3;
  double DeltaPhi4 = ntuple->DeltaPhi4;

  double HT = ntuple->HT;
  double MHT = ntuple->MHT;
  int NJets = ntuple->NJets;

  return NJets >= 2 && MHT > 250. && HT > 300.
      && (DeltaPhi1 < 0.5 || DeltaPhi2 < 0.5 || DeltaPhi3 < 0.3 || DeltaPhi4 < 0.3)
    && cutFlow_leptonVeto(ntuple)
    && onePhoton(ntuple)
    && cutFlow_photonPt200(ntuple)
    && cutFlow_filters(ntuple);
      ;

}

template<typename ntupleType> bool RA2bLDPBaselineCut(ntupleType* ntuple){

  double DeltaPhi1 = ntuple->DeltaPhi1;
  double DeltaPhi2 = ntuple->DeltaPhi2;
  double DeltaPhi3 = ntuple->DeltaPhi3;
  double DeltaPhi4 = ntuple->DeltaPhi4;

  double HT = ntuple->HT;
  double MHT = ntuple->MHT;
  int NJets = ntuple->NJets;

  return NJets >= 2 && MHT > 250. && HT > 300.
      && (DeltaPhi1 < 0.5 || DeltaPhi2 < 0.5 || DeltaPhi3 < 0.3 || DeltaPhi4 < 0.3)
    && cutFlow_leptonVeto(ntuple)
    && cutFlow_filters(ntuple);
      ;

}

TFile *f1 = new TFile("../data/L1PrefiringMaps_new.root");
TH2F* h_photon = (TH2F*)f1->Get("L1prefiring_photonptvseta_2017BtoF");
TH2F* h_jet = (TH2F*)f1->Get("L1prefiring_jetptvseta_2017BtoF");

/* .......................Prefiring Weight ..............................................*/
/*****************************************************************************************/
 double prefiring_weight_photon(RA2bTree* ntuple){
          return ( 1 - h_photon->GetBinContent(h_photon->GetXaxis()->FindBin(ntuple->Photons->at(0).Eta()),h_photon->GetYaxis()->FindBin(ntuple->Photons->at(0).Pt())));
   }


 double prefiring_weight_jet(RA2bTree* ntuple,int unsigned s ){
          return ( 1 - h_jet->GetBinContent(h_jet->GetXaxis()->FindBin(ntuple->Jets->at(s).Eta()),h_jet->GetYaxis()->FindBin(ntuple->Jets->at(s).Pt()))) ;
 }


 /*****************.>>>>>>>>>>>>>>>>>>>>   Photon Trigger Efficiency <<<<<<<<<<<<<<<***********/
 /*****................................................................................********/
 TFile *etrigger = new TFile("../data/trigger_efficiency_PhotonPt_2018.root","READ");
  std::vector<TEfficiency*> eTrigEff_;
  void trig_eff_func()
	{
          eTrigEff_.clear();
	  eTrigEff_.push_back((TEfficiency*)etrigger->Get("f_trig_eb"));
	  eTrigEff_.push_back((TEfficiency*)etrigger->Get("f_trig_ec"));

	}
  
 double trig_eff(RA2bTree* ntuple){
 if( ntuple->Photons_isEB->at(0) && eTrigEff_.at(0) != nullptr ) {
       TH1F* htot = (TH1F*) eTrigEff_.at(0)->GetTotalHistogram();    
       return  eTrigEff_.at(0)->GetEfficiency(min(htot->GetNbinsX(), htot->FindBin(ntuple->Photons->at(0).Pt())));  
  } 
     
  else if (!ntuple->Photons_isEB->at(0) && eTrigEff_.at(1) != nullptr ) {
       TH1F* htot = (TH1F*) eTrigEff_.at(1)->GetTotalHistogram();    
       return  eTrigEff_.at(1)->GetEfficiency(min(htot->GetNbinsX(), htot->FindBin(ntuple->Photons->at(0).Pt())));  
  }
  else
       return 1;
   
} 
 /***********************************HEM veto  **********************************************/
//////***************************************************************************************/
 
 int StartHEM = 319077;
 bool isMC_, Run_number;
 double MHTPhi,JetPhi,DeltaPhi;
 double deltaPhi(double JetPhi,double MHTPhi) {return TVector2::Phi_mpi_pi(MHTPhi - JetPhi);}

 /*bool testHEM() {
    if (!isMC_ && RunNum < StartHEM) return false;
    if ( isMC_ && EvtNum % 1000 < 1000*21.0/59.6) return false;
    return true;
 }
 */
 
 bool passHEMjetVeto(RA2bTree* ntuple, double ptThresh = 30) {
   
    if (!isMC_ && ntuple->RunNum < StartHEM) return true; 
    MHTPhi = ntuple->MHTPhi;     

   for (int p = 0; p < ntuple->Jets->size(); p++){
   JetPhi = ntuple->Jets->at(p).Phi();
   DeltaPhi = deltaPhi(JetPhi,MHTPhi);
   if (-3.2 <= ntuple->Jets->at(p).Eta() && ntuple->Jets->at(p).Eta() <= -1.2 && 
      -1.77 <= ntuple->Jets->at(p).Phi() && ntuple->Jets->at(p).Phi() <= -0.67 &&
       ntuple->Jets->at(p).Pt() > ptThresh && abs(DeltaPhi) < 0.5)
   
    return false;
    } 
    return true;
};


/***************************************>>>>>> ZPt Reweighting <<<<<<<<<<<<<***********************************************/
/**************************************************************************************************************************/
   double ZPtWt = 1; 
   Double_t ptBins[297] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,35.,36.,37.,38.,39.,40.,41.,42.,43.,44.,45.,46.,47.,48.,49.,50.,51.,52.,53.,54.,55.,56.,57.,58.,59.,60.,61.,62.,63.,64.,65.,66.,67.,68.,69.,70.,71.,72.,73.,74.,75.,76.,77.,78.,79.,80.,81.,82.,83.,84.,85.,86.,87.,88.,89.,90.,91.,92.,93.,94.,95.,96.,97.,98.,99.,100.,101.,102.,103.,104.,105.,106.,107.,108.,109.,110.,111.,112.,113.,114.,115.,116.,117.,118.,119.,120.,121.,122.,123.,124.,125.,126.,127.,128.,129.,130.,131.,132.,133.,134.,135.,136.,137.,138.,139.,140.,141.,142.,143.,144.,145.,146.,147.,148.,149.,150.,151.,152.,153.,154.,155.,156.,157.,158.,159.,160.,161.,162.,163.,164.,165.,166.,167.,168.,169.,170.,171.,172.,173.,174.,175.,176.,177.,178.,179.,180.,181.,182.,183.,184.,185.,186.,187.,188.,189.,190.,191.,192.,193.,194.,195.,196.,197.,198.,199.,200.,202.,204.,206.,208.,210.,212.,214.,216.,218.,220.,222.,224.,226.,228.,230.,232.,234.,236.,238.,240.,242.,244.,246.,248.,250.,252.,254.,256.,258.,260.,262.,264.,266.,268.,270.,272.,274.,276.,278.,280.,282.,284.,286.,288.,290.,292.,294.,296.,298.,300.,304.,308.,312.,316.,320.,324.,328.,332.,336.,340.,344.,348.,352.,356.,360.,364.,368.,372.,376.,380.,384.,388.,392.,396.,400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,500.,520.,540.,560.,580.,600.,650.,700.,750.,800.,900.,1000.};
  
   Double_t ptWgts[297] = {0.615, 0.626, 0.649, 0.688, 0.739, 0.796, 0.852, 0.899, 0.936, 0.969, 0.995, 1.015, 1.034, 1.05, 1.065, 1.079, 1.092, 1.103, 1.118, 1.126, 1.139, 1.146, 1.156, 1.16, 1.164, 1.164, 1.165, 1.165, 1.164, 1.164, 1.164, 1.161, 1.158, 1.156, 1.154, 1.146, 1.147, 1.14, 1.135, 1.134, 1.129, 1.128, 1.121, 1.115, 1.11, 1.109, 1.111, 1.103, 1.094, 1.092, 1.093, 1.089, 1.085, 1.084, 1.074, 1.074, 1.067, 1.068, 1.062, 1.066, 1.063, 1.055, 1.05, 1.054, 1.048, 1.044, 1.042, 1.043, 1.034, 1.032, 1.035, 1.033, 1.028, 1.03, 1.029, 1.03, 1.012, 1.012, 1.018, 1.013, 1.011, 1, 1.007, 1.015, 0.989, 1.001, 0.994, 0.99, 0.99, 0.986, 0.981, 0.986, 0.972, 0.971, 0.977, 0.974, 0.978, 0.966, 0.985, 0.978, 0.966, 0.972, 0.969, 0.971, 0.964, 0.962, 0.948, 0.952, 0.943, 0.973, 0.936, 0.945, 0.935, 0.944, 0.953, 0.943, 0.935, 0.926, 0.946, 0.94, 0.943, 0.933, 0.92, 0.912, 0.923, 0.904, 0.919, 0.937, 0.941, 0.929, 0.923, 0.902, 0.925, 0.927, 0.896, 0.926, 0.905, 0.92, 0.908, 0.889, 0.91, 0.913, 0.911, 0.889, 0.881, 0.888, 0.904, 0.886, 0.891, 0.915, 0.879, 0.884, 0.9, 0.874, 0.85, 0.874, 0.873, 0.872, 0.894, 0.88, 0.87, 0.871, 0.863, 0.883, 0.863, 0.892, 0.845, 0.863, 0.892, 0.851, 0.878, 0.847, 0.881, 0.809, 0.86, 0.829, 0.851, 0.871, 0.846, 0.825, 0.86, 0.853, 0.894, 0.807, 0.793, 0.863, 0.832, 0.829, 0.87, 0.85, 0.831, 0.82, 0.829, 0.839, 0.88, 0.831, 0.804, 0.836, 0.822, 0.796, 0.812, 0.831, 0.83, 0.827, 0.802, 0.781, 0.855, 0.798, 0.774, 0.79, 0.81, 0.825, 0.799, 0.79, 0.811, 0.778, 0.803, 0.779, 0.795, 0.768, 0.809, 0.762, 0.767, 0.752, 0.822, 0.777, 0.791, 0.799, 0.721, 0.776, 0.718, 0.798, 0.754, 0.749, 0.758, 0.847, 0.766, 0.775, 0.718, 0.756, 0.826, 0.792, 0.792, 0.767, 0.679, 0.694, 0.75, 0.738, 0.707, 0.709, 0.711, 0.754, 0.762, 0.717, 0.722, 0.692, 0.714, 0.726, 0.703, 0.693, 0.704, 0.704, 0.653, 0.718, 0.748, 0.713, 0.733, 0.741, 0.742, 0.652, 0.586, 0.644, 0.656, 0.702, 0.721, 0.682, 0.758, 0.662, 0.578, 0.654, 0.723, 0.634, 0.68, 0.656, 0.602, 0.59, 0.566, 0.623, 0.562, 0.589, 0.604, 0.554, 0.525, 0.552, 0.503, 0.563, 0.476};
  double  ZPtWeight(RA2bTree* ntuple){
      double ptZ = -1.0;
      int nGen =  ntuple->GenParticles_PdgId->size() ;
        for (int iGen = 0 ; iGen < nGen; ++iGen) {
          if (ntuple->GenParticles_PdgId->at(iGen)==23 && ntuple->GenParticles_Status->at(iGen) == 62) {
            ptZ = ntuple->GenParticles->at(iGen).Pt();
            break;
          }
        }
        ZPtWt = 1.0;
        if (ptZ > 0.0) {
          int iptbin;
          for (iptbin = 1; iptbin<297; iptbin++) {
            if (ptZ < ptBins[iptbin]) break;
          }
          ZPtWt *= ptWgts[iptbin-1];
        }

    return ZPtWt; 
   }

/******************************** MC weight correction *************************/
/////////////////////////////////////////////////////////////////////////////////

double MCwtCorr(RA2bTree* ntuple){

          if      (100 <= ntuple->madHT && ntuple->madHT <= 200) return 1.05713;
          else if (200 <= ntuple->madHT && ntuple->madHT <= 400) return 1.20695;
          else if (400 <= ntuple->madHT && ntuple->madHT <= 600) return 1.30533;
          else if (600 <= ntuple->madHT && ntuple->madHT <= 800) return 1.38453;
          else if (800 <= ntuple->madHT && ntuple->madHT <= 1200) return 1.40301;
          else if (1200<= ntuple->madHT && ntuple->madHT <= 2500) return 1.42145;
          else if (2500<= ntuple->madHT && ntuple->madHT <= 10000000) return 1.11697;

}

//********************EcalNoiseJetFilter****************************************//

 	 Double_t MHTPhiv2Recipe;
         TVector3 MHT3Vecv2Recipe,temp3Vec;
	 std::vector<int>  MHTminusHTJetsIdxv2Recipe,MHTJetsIdxv2Recipe;
	 Double_t        MHTminusHTDeltaPhi1v2Recipe;
	 Double_t        MHTminusHTDeltaPhi2v2Recipe;
	 Double_t        MHTminusHTDeltaPhi3v2Recipe;
	 Double_t        MHTminusHTDeltaPhi4v2Recipe;
	 const bool EENoiseCutbyAditee = true;

  bool NoiseJetfilter(RA2bTree* ntuple,int iEvt){
	   ntuple->GetEntry(iEvt);
	   MHTPhiv2Recipe=-99.;

	   MHTminusHTJetsIdxv2Recipe.clear();
           MHTJetsIdxv2Recipe.clear();

	   for(int j = 0; j < ntuple->Jets->size(); ++j){
                     if(ntuple->Jets->at(j).Pt()>30 && fabs(ntuple->Jets->at(j).Eta()) < 5.0)
				MHTJetsIdxv2Recipe.push_back(j);
		     if(ntuple->Jets->at(j).Pt()>30 && fabs(ntuple->Jets->at(j).Eta()) > 2.4 && fabs(ntuple->Jets->at(j).Eta()) < 5.0)
			        MHTminusHTJetsIdxv2Recipe.push_back(j);
           }
  

           for(int i=0; i<MHTJetsIdxv2Recipe.size(); i++){
              
                     int jetIdx=MHTJetsIdxv2Recipe[i];
                     temp3Vec.SetPtEtaPhi(ntuple->Jets->at(jetIdx).Pt(),ntuple->Jets->at(jetIdx).Eta(),ntuple->Jets->at(jetIdx).Phi());
                     MHT3Vecv2Recipe-=temp3Vec;
            }
  		
	    MHTPhiv2Recipe=MHT3Vecv2Recipe.Phi();

	  if(MHTminusHTJetsIdxv2Recipe.size()>0)
	    MHTminusHTDeltaPhi1v2Recipe=fabs(TVector2::Phi_mpi_pi(ntuple->Jets->at(MHTminusHTJetsIdxv2Recipe[0]).Phi() - MHTPhiv2Recipe ));
	  if(MHTminusHTJetsIdxv2Recipe.size()>1)
  	    MHTminusHTDeltaPhi2v2Recipe=fabs(TVector2::Phi_mpi_pi(ntuple->Jets->at(MHTminusHTJetsIdxv2Recipe[1]).Phi() - MHTPhiv2Recipe ));
	  if(MHTminusHTJetsIdxv2Recipe.size()>2)
	    MHTminusHTDeltaPhi3v2Recipe=fabs(TVector2::Phi_mpi_pi(ntuple->Jets->at(MHTminusHTJetsIdxv2Recipe[2]).Phi() - MHTPhiv2Recipe ));
 	  if(MHTminusHTJetsIdxv2Recipe.size()>3)
	    MHTminusHTDeltaPhi4v2Recipe=fabs(TVector2::Phi_mpi_pi(ntuple->Jets->at(MHTminusHTJetsIdxv2Recipe[3]).Phi() - MHTPhiv2Recipe ));

	
	  if(EENoiseCutbyAditee){
		    if((MHTminusHTJetsIdxv2Recipe.size()>0 && ntuple->Jets->at(MHTminusHTJetsIdxv2Recipe[0]).Pt()>250 && (MHTminusHTDeltaPhi1v2Recipe>2.6 || MHTminusHTDeltaPhi1v2Recipe<0.1)) || (MHTminusHTJetsIdxv2Recipe.size()>1 && ntuple->Jets->at(MHTminusHTJetsIdxv2Recipe[1]).Pt()>250 && (MHTminusHTDeltaPhi2v2Recipe>2.6 || MHTminusHTDeltaPhi2v2Recipe<0.1)))
       			     return true;

	            else
   		      	    return false;
	  }
 }

/***********************************  SF weights ************************************/
TFile *f_phoSF = new TFile("~/SF_rootFiles/Photons2018_SF_all.root");
TH2F* h_SFweight = (TH2F*)f_phoSF->Get("EGamma_SF2D");
float photon_pt = 0; float photon_eta = 0;

  double SFweights(RA2bTree* ntuple){
        photon_pt = ntuple->Photons->at(0).Pt();
        photon_eta = ntuple->Photons->at(0).Eta();
        if (photon_pt>500) photon_pt=499.9;
        return ( h_SFweight->GetBinContent(h_SFweight->GetXaxis()->FindBin(photon_eta),h_SFweight->GetYaxis()->FindBin(photon_pt)));
 }
