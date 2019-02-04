#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TROOT.h"
#include "THStack.h"
#include "TPad.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "plotterUtils.cc"
#include "skimSamples.cc"
#include "definitions.h"
#include "RA2bTree.cc"

using namespace std;

static int MAX_EVENTS=99999999;

void process(int region, string backgroundSample, string dataSample){

    skimSamples::region reg = static_cast<skimSamples::region>(region);

    if( reg != skimSamples::kPhotonLoose && reg != skimSamples::kPhotonLDPLoose && reg != skimSamples::kDYe ){
        cerr << "only photon regions are valid" << endl;
        return;
    }

    vector<TString> backgroundSamples;
    if( backgroundSample != "" )
      backgroundSamples.push_back(TString(backgroundSample));
    vector<TString> signalSamples;
    vector<TString> dataSamples;
    if( dataSample != "" )
      dataSamples.push_back(TString(dataSample));

    skimSamples skims(reg,backgroundSamples,signalSamples,dataSamples);
    typedef plot<RA2bTree> plot;

    cout << skims.skimType << endl;
    for( int i = 0 ; i < skims.ntuples.size() ; i++){
        cout << "sample name: " << skims.sampleName[i] << endl;
    }

    // Triggers that will be used for data
    vector<int> trigger_indices = {140,141};
    //vector<int> trigger_indices = {52};

    double MHTbinning[7]={200.,225.,250.,300.,350.,600.,2000.};
    plot MHTplot(*fillMHT<RA2bTree>,"MHT_"+skims.regionNames[region],"MHT [GeV]",6,MHTbinning);

    vector<plot> projections;
    projections.push_back(MHTplot);
    vector<TString> projLabel;
    projLabel.push_back("MHT");

    vector<vector<plot> > chargeIsoEEHighSieieVersus;
    vector<vector<plot> > chargeIsoEELowSieieVersus;
    vector<vector<plot> > chargeIsoEBHighSieieVersus;
    vector<vector<plot> > chargeIsoEBLowSieieVersus;
    for( int iProj = 0 ; iProj<projections.size(); iProj++ ){
        vector<plot> chargeIsoEEHighSieieTemp;
        vector<plot> chargeIsoEELowSieieTemp;
        vector<plot> chargeIsoEBHighSieieTemp;
        vector<plot> chargeIsoEBLowSieieTemp;
        for( int iBin = 0 ; iBin <= projections[iProj].nbins ; iBin++ ){
            TString binEdge="_";
            binEdge+=projections[iProj].binEdges[iBin];
            chargeIsoEEHighSieieTemp.push_back(plot(*photonIsoChrg<RA2bTree>,"photonIsoChrgHighSieie_EE_"+projLabel[iProj]+binEdge+"_"+skims.regionNames[region],"Iso_{chrg}",25,0.,10.));
            chargeIsoEELowSieieTemp.push_back(plot(*photonIsoChrg<RA2bTree>,"photonIsoChrgLowSieie_EE_"+projLabel[iProj]+binEdge+"_"+skims.regionNames[region],"Iso_{chrg}",25,0.,10.));

            chargeIsoEBHighSieieTemp.push_back(plot(*photonIsoChrg<RA2bTree>,"photonIsoChrgHighSieie_EB_"+projLabel[iProj]+binEdge+"_"+skims.regionNames[region],"Iso_{chrg}",25,0.,10.));
            chargeIsoEBLowSieieTemp.push_back(plot(*photonIsoChrg<RA2bTree>,"photonIsoChrgLowSieie_EB_"+projLabel[iProj]+binEdge+"_"+skims.regionNames[region],"Iso_{chrg}",25,0.,10.));
        }
        chargeIsoEEHighSieieVersus.push_back(chargeIsoEEHighSieieTemp);
        chargeIsoEELowSieieVersus.push_back(chargeIsoEELowSieieTemp);
        chargeIsoEBHighSieieVersus.push_back(chargeIsoEBHighSieieTemp);
        chargeIsoEBLowSieieVersus.push_back(chargeIsoEBLowSieieTemp);
    }

    plot SieieEBplot(*photonSieie<RA2bTree>,"photonSieie_EB_"+skims.regionNames[region],"#sigma_{i#etai#eta}",100,0.,.03);
    plot IsoChrgEBplot(*photonIsoChrg<RA2bTree>,"photonIsoChrg_EB_"+skims.regionNames[region],"Iso_{chrg}",100,0.,10.);
    plot IsoGamEBplot(*photonIsoGam<RA2bTree>,"photonIsoGam_EB_"+skims.regionNames[region],"Iso_{gam}",100,0.,10.);
    plot IsoNeuEBplot(*photonIsoNeu<RA2bTree>,"photonIsoNeu_EB_"+skims.regionNames[region],"Iso_{neu}",100,0.,10.);
    plot HoverEEBplot(*photonHoverE<RA2bTree>,"photonHoverE_EB_"+skims.regionNames[region],"H/E",100,0.,.1);
    plot IsoChrgHighSieiEBplot(*photonIsoChrg<RA2bTree>,"photonIsoChrgHighSieie_EB_"+skims.regionNames[region],"Iso_{chrg}",100,0.,10.);
    plot IsoChrgLowSieieEBplot(*photonIsoChrg<RA2bTree>,"photonIsoChrgLowSieie_EB_"+skims.regionNames[region],"Iso_{chrg}",100,0.,10.);

    plot SieieEEplot(*photonSieie<RA2bTree>,"photonSieie_EE_"+skims.regionNames[region],"#sigma_{i#etai#eta}",100,0.,.07);
    plot IsoChrgEEplot(*photonIsoChrg<RA2bTree>,"photonIsoChrg_EE_"+skims.regionNames[region],"Iso_{chrg}",100,0.,10.);
    plot IsoChrgHighSieiEEplot(*photonIsoChrg<RA2bTree>,"photonIsoChrgHighSieie_EE_"+skims.regionNames[region],"Iso_{chrg}",100,0.,10.);
    plot IsoChrgLowSieieEEplot(*photonIsoChrg<RA2bTree>,"photonIsoChrgLowSieie_EE_"+skims.regionNames[region],"Iso_{chrg}",100,0.,10.);
    plot HoverEEEplot(*photonHoverE<RA2bTree>,"photonHoverE_EE_"+skims.regionNames[region],"H/E",100,0.,.1);
    plot IsoGamEEplot(*photonIsoGam<RA2bTree>,"photonIsoGam_EE_"+skims.regionNames[region],"Iso_{gam}",100,0.,10.);
    plot IsoNeuEEplot(*photonIsoNeu<RA2bTree>,"photonIsoNeu_EE_"+skims.regionNames[region],"Iso_{neu}",100,0.,10.);

    // need this plotted for QCD, GJets, Data above sieie, and Data below sieie
    // then, we ultimately need these histograms in different regions (i.e. MHT regions)

    vector<plot> plotsEB;
    plotsEB.push_back(SieieEBplot);
    plotsEB.push_back(IsoChrgEBplot);
    plotsEB.push_back(IsoGamEBplot);
    plotsEB.push_back(IsoNeuEBplot);
    plotsEB.push_back(HoverEEBplot);
    plotsEB.push_back(IsoChrgHighSieiEBplot);
    plotsEB.push_back(IsoChrgLowSieieEBplot);

    vector<plot> plotsEE;
    plotsEE.push_back(SieieEEplot);
    plotsEE.push_back(IsoChrgEEplot);
    plotsEE.push_back(IsoGamEEplot);
    plotsEE.push_back(IsoNeuEEplot);
    plotsEE.push_back(HoverEEEplot);
    plotsEE.push_back(IsoChrgHighSieiEEplot);
    plotsEE.push_back(IsoChrgLowSieieEEplot);

    // output file
    TFile* outputFile = new TFile("plotPurityProperties_"+skims.regionNames[region]+"_"+backgroundSample+".root","RECREATE");

    // background MC samples
    for( int iSample = 0 ; iSample < skims.ntuples.size() ; iSample++){

        RA2bTree* ntuple = skims.ntuples[iSample];

        for( int iPlot = 0 ; iPlot < plotsEB.size() ; iPlot++){
            plotsEB[iPlot].addNtuple(ntuple,skims.sampleName[iSample]);
            plotsEE[iPlot].addNtuple(ntuple,skims.sampleName[iSample]);
        }

        for( int iProj = 0 ; iProj < projections.size() ; iProj++){
            projections[iProj].addNtuple(ntuple,skims.sampleName[iSample]);
            for( int iBin = 0 ; iBin < chargeIsoEBHighSieieVersus[iProj].size() ; iBin++){
                chargeIsoEBHighSieieVersus[iProj][iBin].addNtuple(ntuple,skims.sampleName[iSample]);
                chargeIsoEEHighSieieVersus[iProj][iBin].addNtuple(ntuple,skims.sampleName[iSample]);
                chargeIsoEBLowSieieVersus[iProj][iBin].addNtuple(ntuple,skims.sampleName[iSample]);
                chargeIsoEELowSieieVersus[iProj][iBin].addNtuple(ntuple,skims.sampleName[iSample]);
            }
        }
	
        int numEvents = ntuple->fChain->GetEntries();
	cout << "number of " << skims.sampleName[iSample] << " events: " << numEvents << endl;
        ntupleBranchStatus<RA2bTree>(ntuple);
        double weight = 1.0;
        for( int iEvt = 0 ; iEvt < min(MAX_EVENTS,numEvents) ; iEvt++ ){
            ntuple->GetEntry(iEvt);
            if( iEvt % 1000000 == 0 ) cout << skims.sampleName[iSample] << ": " << iEvt << "/" << numEvents << endl;

            //if( ntuple->Photons->size() == 0 || (ntuple->Photons->size() > 0 && ntuple->Photons->at(0).Pt()<100.)) continue;
	    
	    // for studying proton egamma objects in DY control region
	    if( skims.regionNames[region] == "DYe" && !isDY(ntuple) ) continue;
            if( skims.regionNames[region] == "photonLDPLoose" || skims.regionNames[region] == "photonLoose" || skims.regionNames[region] == "photonLDP" || skims.regionNames[region] == "photon" ){
	      if( skims.sampleName[iSample] == "QCD" && isPromptPhoton(ntuple) ) continue;
	      if( skims.sampleName[iSample] == "GJets" && !isPromptPhoton(ntuple) ) continue;
            }

            if( ((skims.regionNames[region] == "photonLDPLoose"||skims.regionNames[region] == "photonLDP")&&!RA2bLDPBaselinePhotonCut(ntuple)) ) continue;
	    if( ((skims.regionNames[region] == "photonLoose"||skims.regionNames[region] == "photon")&&!RA2bBaselinePhotonWideCut(ntuple)) ) continue;
	    

            for( int iProj = 0 ; iProj<projections.size() ; iProj++ ){

                int iBin = projections[iProj].fill(ntuple);
                if( iBin <= 0 ) continue;

                if( ntuple->Photons_isEB->at(0) ){
                    if(ntuple->Photons_sigmaIetaIeta->at(0)>.0106){
                        chargeIsoEBHighSieieVersus[iProj][iBin-1].fill(ntuple);
                    }else{
                        chargeIsoEBLowSieieVersus[iProj][iBin-1].fill(ntuple);
                    }
                }
                else{
                    if(ntuple->Photons_sigmaIetaIeta->at(0)>.0272){
                        chargeIsoEEHighSieieVersus[iProj][iBin-1].fill(ntuple);
                    }else{
                        chargeIsoEELowSieieVersus[iProj][iBin-1].fill(ntuple);
                    }
                }
            }

            for( int iPlot = 0 ; iPlot < plotsEB.size()-2 ; iPlot++ ){
                if( ntuple->Photons_isEB->at(0) )
                    plotsEB[iPlot].fill(ntuple);
                else
                    plotsEE[iPlot].fill(ntuple);
            }
            if( ntuple->Photons_isEB->at(0) ){
                if(ntuple->Photons_sigmaIetaIeta->at(0)>.0106)
                    plotsEB[plotsEB.size()-2].fill(ntuple);
                else
                    plotsEB[plotsEB.size()-1].fill(ntuple);
            }else{
                if(ntuple->Photons_sigmaIetaIeta->at(0)>.0272)
                    plotsEE[plotsEB.size()-2].fill(ntuple);
                else
                    plotsEE[plotsEB.size()-1].fill(ntuple);
            }	    
        }	

	// ---------- write background files -----------
        for( int iPlot = 0 ; iPlot < plotsEB.size() ; iPlot++){
	  plotsEB[iPlot].Write(ntuple);
	  plotsEE[iPlot].Write(ntuple);
        }	
        for( int iProj = 0 ; iProj < projections.size() ; iProj++){
	  projections[iProj].Write(ntuple);
	  for( int iBin = 0 ; iBin < chargeIsoEBHighSieieVersus[iProj].size() ; iBin++){
	    chargeIsoEBHighSieieVersus[iProj][iBin].Write(ntuple);
	    chargeIsoEEHighSieieVersus[iProj][iBin].Write(ntuple);
	    chargeIsoEBLowSieieVersus[iProj][iBin].Write(ntuple);
	    chargeIsoEELowSieieVersus[iProj][iBin].Write(ntuple);
	  }
        }
	// ------------------------------------------------
    }// end loop over background samples
    outputFile->Close();
    delete outputFile;


    // output file
    TFile* dataOutputFile = new TFile("plotPurityProperties_"+skims.regionNames[region]+"_"+dataSample+".root","RECREATE");

    // Data samples
    for( int iSample = 0 ; iSample < skims.dataNtuple.size() ; iSample++){

      RA2bTree* ntuple = skims.dataNtuple[iSample];
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      for( int iPlot = 0 ; iPlot < plotsEB.size() ; iPlot++ ){
	plotsEB[iPlot].addDataNtuple(ntuple,skims.dataSampleName[iSample]);
	plotsEE[iPlot].addDataNtuple(ntuple,skims.dataSampleName[iSample]);
      }

      for( int iProj = 0 ; iProj < projections.size() ; iProj++){
	projections[iProj].addDataNtuple(ntuple,skims.dataSampleName[iSample]);
	for( int iBin = 0 ; iBin < chargeIsoEBHighSieieVersus[iProj].size() ; iBin++){
	  chargeIsoEBHighSieieVersus[iProj][iBin].addDataNtuple(ntuple,skims.dataSampleName[iSample]);
	  chargeIsoEEHighSieieVersus[iProj][iBin].addDataNtuple(ntuple,skims.dataSampleName[iSample]);
	  chargeIsoEBLowSieieVersus[iProj][iBin].addDataNtuple(ntuple,skims.dataSampleName[iSample]);
	  chargeIsoEELowSieieVersus[iProj][iBin].addDataNtuple(ntuple,skims.dataSampleName[iSample]);
	}
      }    
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      int numEvents = ntuple->fChain->GetEntries();
      cout << "number of data events: " << numEvents << endl;
      ntupleBranchStatus<RA2bTree>(ntuple);
      for( int iEvt = 0 ; iEvt < min(MAX_EVENTS,numEvents) ; iEvt++ ){
        ntuple->GetEntry(iEvt);
        if( iEvt % 1000000 == 0 ) cout << skims.dataSampleName[iSample] << ": " << iEvt << "/" << numEvents << endl;

        bool pass_trigger=true;
        for( unsigned int itrig = 0 ; itrig < trigger_indices.size() ; itrig++ )
          pass_trigger|=(ntuple->TriggerPass->at(trigger_indices[itrig])==1);
        if( !pass_trigger ) continue;

        //if( ntuple->Photons->size() == 0 || (ntuple->Photons->size() > 0 && ntuple->Photons->at(0).Pt()<100.) ) continue;
	
	// for studying proton egamma objects in DY control region	
	if( skims.regionNames[region] == "DYe" && !isDY(ntuple) ) continue;
        if( ((skims.regionNames[region] == "photonLDPLoose"||skims.regionNames[region] == "photonLDP")&&!RA2bLDPBaselineCut(ntuple)) ) continue;
	if( ((skims.regionNames[region] == "photonLoose"||skims.regionNames[region] == "photon")&&!RA2bBaselineWideCut(ntuple)) ) continue;

	// ------------- fill data histograms --------------
        for( int iProj = 0 ; iProj<projections.size() ; iProj++ ){
	  int iBin = projections[iProj].fillData(ntuple);
	  if( iBin <= 0 ) continue;
	  if( ntuple->Photons_isEB->at(0) ){
	    if(ntuple->Photons_sigmaIetaIeta->at(0)>.0106){
	      chargeIsoEBHighSieieVersus[iProj][iBin-1].fillData(ntuple);
	    }else{
	      chargeIsoEBLowSieieVersus[iProj][iBin-1].fillData(ntuple);
	    }
	  }else{
	    if(ntuple->Photons_sigmaIetaIeta->at(0)>.0272){
	      chargeIsoEEHighSieieVersus[iProj][iBin-1].fillData(ntuple);
	    }else{
	      chargeIsoEELowSieieVersus[iProj][iBin-1].fillData(ntuple);
	    }
	  }
	}
	   
        for( int iPlot = 0 ; iPlot < plotsEB.size()-2 ; iPlot++ ){
	  if( ntuple->Photons_isEB->at(0) ){
	    plotsEB[iPlot].fillData(ntuple);
	  }else{
	    plotsEE[iPlot].fillData(ntuple);
	  }
	}
        if( ntuple->Photons_isEB->at(0) ){
	  if(ntuple->Photons_sigmaIetaIeta->at(0)>.0106){
	    plotsEB[plotsEB.size()-2].fillData(ntuple);
	  }else{
	    plotsEB[plotsEB.size()-1].fillData(ntuple);
	  }
	}else{
	  if(ntuple->Photons_sigmaIetaIeta->at(0)>.0272){
	    plotsEE[plotsEB.size()-2].fillData(ntuple);
	  }else{
	    plotsEE[plotsEB.size()-1].fillData(ntuple);
	  }
        }
	// --------------------------------------------------
      }
      // ------------- write data histograms --------------
      for( int iPlot = 0 ; iPlot < plotsEB.size() ; iPlot++ ){
	plotsEB[iPlot].WriteData();
	plotsEE[iPlot].WriteData();
      }
    
      for( int iProj = 0 ; iProj < projections.size() ; iProj++){
	projections[iProj].WriteData();
	for( int iBin = 0 ; iBin < chargeIsoEBHighSieieVersus[iProj].size() ; iBin++){
	  chargeIsoEBHighSieieVersus[iProj][iBin].WriteData();
	  chargeIsoEEHighSieieVersus[iProj][iBin].WriteData();
	  chargeIsoEBLowSieieVersus[iProj][iBin].WriteData();
	  chargeIsoEELowSieieVersus[iProj][iBin].WriteData();
	}
      }    
      // --------------------------------------------------
    }
    dataOutputFile->Close();
    delete dataOutputFile;
}



int main(int argc, char** argv){

    gROOT->ProcessLine(".L tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    if(argc != 4 ){
        cerr << "please provide a region index, a background sample, and a data sample" << endl;
    }

    process(atoi(argv[1]),argv[2],argv[3]);

}
