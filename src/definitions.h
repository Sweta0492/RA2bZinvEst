#include "TRandom.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"

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
    return 1. /( (min(ntuple->MHT, 900.0) - 399.6)*( -0.00040321 *2/3) + 0.8389); 
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


template<typename ntupleType> double fillRA2b13Bins(ntupleType* ntuple){
  double MHT = ntuple->MHT;
  double HT = ntuple->HT;

  if( MHT > 250. && MHT < 300. ){
    if( HT > 300. && HT < 500. )
      if( ntuple->NJets>=7 )
	     return -999999.;
      else
	     return 1.;
    else if( HT > 500. && HT < 1000. )
      if( ntuple->NJets>=7 )
	     return 1.;
      else
	     return 2.;
    else if( HT > 1000. )
      if( ntuple->NJets>=7 )
	     return 2.;
      else
	     return 3.;
  }else if( MHT > 300. && MHT < 350. ){
    if( HT > 300. && HT < 500. ){
      if(ntuple->NJets>=7)
	     return -999999.;
      else
	     return 4.;
    }else if( HT > 500. && HT < 1000. ){
      if(ntuple->NJets>=7)
	     return 3.;
      else
	     return 5.;
    }else if( HT > 1000. ){
      if(ntuple->NJets>=7)
	     return 4.;
      else
	     return 6.;
    }else
      return -999999.;
  }else if( MHT > 350. && MHT < 500. ){
    if( HT > 350. && HT < 500. ){
      if(ntuple->NJets>=7)
	     return -999999.;
      else
	     return 7.;
    }else if( HT > 500. && HT < 1000. ){
      if(ntuple->NJets>=7)
	     return 5.;
      else
	     return 8.;
    }else if( HT > 1000. ){
      if(ntuple->NJets>=7)
	     return 6.;
      else
	     return 9.;
    }else
      return -999999.;
  }else if( MHT > 500. && MHT < 750. ){
    if( HT > 500. && HT < 1000. ){
      if(ntuple->NJets>=7)
	     return 7.;
      else
	     return 10.;
    }else if( HT > 1000. ){
      if(ntuple->NJets>=7)
	     return 8.;
      else
	     return 11.;
    }else
      return -999999.;
  }else if( MHT > 750. ){
    if( HT > 750. && HT < 1500. ){
      if(ntuple->NJets>=7)
	     return 9.;
      else
	     return 12.;
    }else if( HT > 1500. ){
      if(ntuple->NJets>=7)
	     return 10.;
      else
	     return 13.;
    }else
      return -999999.;
  }else
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

template<typename ntupleType> double fillRA2b59Bins( ntupleType* ntuple ){

  int BTags = int(ntuple->BTagsDeepCSV);
  if( BTags != 0 ) return -999999.;
  int NJets = int(ntuple->NJets);

  if( NJets == 2 ){
      return fillRA2b13Bins(ntuple);
  }else if( NJets >= 3 && NJets <=4 ){
    return 13.+fillRA2b13Bins(ntuple);
  }else if( NJets >= 5 && NJets <= 6 ){
    return 26.+fillRA2b13Bins(ntuple);
  }else if( NJets >= 7 && NJets <= 8 ){
    return 39.+fillRA2b13Bins(ntuple);
  }else if( NJets >= 9 ){
    return 49.+fillRA2b13Bins(ntuple);
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

template<typename ntupleType> double fillRA2b174Bins( ntupleType* ntuple ){

  int BTags = int(ntuple->BTagsDeepCSV);
  int NJets = int(ntuple->NJets);

  if( NJets == 2 ){
    if( BTags == 0 ){
      return fillRA2b10Bins(ntuple);
    }else if( BTags == 1 ){
      return 10.+fillRA2b10Bins(ntuple);
    }else if( BTags>=2 ){
      return 20.+fillRA2b10Bins(ntuple);
    }else
      return -999999.;
  }else if( NJets >= 3 && NJets <=4 ){
    if( BTags == 0 )
      return 30.+fillRA2b10Bins(ntuple);
    else if( BTags == 1 )
      return 40.+fillRA2b10Bins(ntuple);
    else if( BTags == 2 )
      return 50.+fillRA2b10Bins(ntuple);
    else if( BTags >= 3 )
      return 60.+fillRA2b10Bins(ntuple);
  }else if( NJets >= 5 && NJets <= 6 ){
    if( BTags == 0 )
      return 70.+fillRA2b10Bins(ntuple);
    else if( BTags == 1 )
      return 80.+fillRA2b10Bins(ntuple);
    else if( BTags == 2 )
      return 90.+fillRA2b10Bins(ntuple);
    else if( BTags >= 3 )
      return 100.+fillRA2b10Bins(ntuple);
  }else if( NJets >= 7 && NJets <= 8 ){
    if( BTags == 0 )
      return 110.+fillRA2b10Bins(ntuple);
    else if( BTags == 1 )
      return 118.+fillRA2b10Bins(ntuple);
    else if( BTags == 2 )
      return 126.+fillRA2b10Bins(ntuple);
    else if( BTags >= 3 )
      return 134.+fillRA2b10Bins(ntuple);
  }else if( NJets >= 9 ){
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
 double prefiring_weight_photon(RA2bTree* ntuple,int iEvt ){
          ntuple->GetEntry(iEvt);
          return ( 1 - h_photon->GetBinContent(h_photon->GetXaxis()->FindBin(ntuple->Photons->at(0).Eta()),h_photon->GetYaxis()->FindBin(ntuple->Photons->at(0).Pt())));
   }


   double prefiring_weight_jet(RA2bTree* ntuple,int iEvt,int unsigned s ){
          ntuple->GetEntry(iEvt);
          return ( 1 - h_jet->GetBinContent(h_jet->GetXaxis()->FindBin(ntuple->Jets->at(s).Eta()),h_jet->GetYaxis()->FindBin(ntuple->Jets->at(s).Pt()))) ;
 }


 /*****************.>>>>>>>>>>>>>>>>>>>>   Photon Trigger Efficiency <<<<<<<<<<<<<<<***********/
 /*****................................................................................********/
  TFile *ftrigger = new TFile("../data/trigger_efficiency_PhotonPt_2018.root","READ");
  std::vector<TF1*> fTrigEff_;
  void trig_eff_func()
	{
          fTrigEff_.clear();
	  fTrigEff_.push_back((TF1*)ftrigger->Get("f_trig_eb"));
	  fTrigEff_.push_back((TF1*)ftrigger->Get("f_trig_ec"));

	}
  
 double trig_eff(RA2bTree* ntuple, int iEvt){
 ntuple->GetEntry(iEvt);
 if( ntuple->Photons_isEB->at(0) && fTrigEff_.at(0) != nullptr ) 
       return  fTrigEff_.at(0)->Eval(max(double(205),ntuple->Photons->at(0).Pt()));  
     
  else if (!ntuple->Photons_isEB->at(0) && fTrigEff_.at(1) != nullptr ) 
       return  fTrigEff_.at(1)->Eval(max(double(205),ntuple->Photons->at(0).Pt()));  
  else
       return 1;
   
 }

/************************************HEM veto ***************************************\\\\\\/
 * ...................................................................................****/



 int StartHEM = 319077;
 bool isMC_, Run_number;

 bool passHEMobjVeto(TLorentzVector& obj, double ptThresh = 0) {
    if (!isMC_ && Run_number < StartHEM) return true;
   // if (isMC_ && runBlock_.find("HEM") == std::string::npos) return true;
    if (-3.0 <= obj.Eta() && obj.Eta() <= -1.4 && 
	-1.57 <= obj.Phi() && obj.Phi() <= -0.87 &&
	obj.Pt() > ptThresh)
      return false;
    else return true;
  };
 
 bool passHEMjetVeto(RA2bTree* ntuple, int iEvt,double ptThresh = 30) {
    ntuple->GetEntry(iEvt);
    if (!isMC_ && ntuple->RunNum < StartHEM) return true;
   // if (isMC_ && runBlock_.find("HEM") == std::string::npos) return true;
    Run_number=ntuple->RunNum;    
    for (auto & jet : *ntuple->Jets)
      if (!passHEMobjVeto(jet, ptThresh)) return false;
    return true;
  };























