#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TROOT.h"
#include "THStack.h"
#include "TPad.h"

#include <vector>
#include <map>
#include <iostream>
#include <string> 
#include <fstream>

#include "plotterUtils.cc"
#include "skimSamples.cc"
#include "definitions.h"
#include "RA2bTree.cc"

using namespace std;

const static int MAX_EVENTS=999999999;

int main(int argc, char** argv){

    gROOT->ProcessLine(".L tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");
  
    // background MC samples
    vector<RA2bTree*> samples;
    vector<TString> sampleNames;

    TString skimType;
    skimType=BASE_DIR+"tree_signal/";
    
    vector<TString> ZJetsFileNames;
    ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-100to200_MC2016.root");
    ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-200to400_MC2016.root");
    ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-400to600_MC2016.root");
    ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-600to800_MC2016.root");
    ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-800to1200_MC2016.root");
    ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-1200to2500_MC2016.root");
    ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-2500toInf_MC2016.root");
    TChain* ZJets = new TChain("tree");
    for( int i = 0 ; i < ZJetsFileNames.size() ; i++ ){
        ZJets->Add(skimType+"/"+ZJetsFileNames[i]);
    }
    samples.push_back(new RA2bTree(ZJets));
    sampleNames.push_back("ZJets");

    skimType=BASE_DIR+"tree_GJet_CleanVars";

    vector<TString> GJetsFileNames;
    GJetsFileNames.push_back("tree_GJets_DR-0p4_HT-100to200_MC2016.root");
    GJetsFileNames.push_back("tree_GJets_DR-0p4_HT-200to400_MC2016.root");
    GJetsFileNames.push_back("tree_GJets_DR-0p4_HT-400to600_MC2016.root");
    GJetsFileNames.push_back("tree_GJets_DR-0p4_HT-600toInf_MC2016.root");
    
    TChain* GJets = new TChain("tree");
    for( int i = 0 ; i < GJetsFileNames.size() ; i++ ){
        GJets->Add(skimType+"/"+GJetsFileNames[i]);
    }
    //samples.push_back(new RA2bTree(GJets));
    //sampleNames.push_back("GJets");
    
    typedef plot<RA2bTree> plot;
    typedef bool(*cutFunc)(RA2bTree*);
    vector<cutFunc> cutFlow;
    vector<TString> cutName;
    cutFlow.push_back(*cutFlow_none<RA2bTree>);
    cutName.push_back("none");
    cutFlow.push_back(*cutFlow_HT300<RA2bTree>);
    cutName.push_back("HT300");
    cutFlow.push_back(*cutFlow_MHT300<RA2bTree>);
    cutName.push_back("MHT300");
    cutFlow.push_back(*cutFlow_NJets2plus<RA2bTree>);
    cutName.push_back("NJets2");
    cutFlow.push_back(*cutFlow_leptonVeto<RA2bTree>);
    cutName.push_back("LeptonVeto");
    // cutFlow.push_back(*cutFlow_onePhoton<RA2bTree>);
    // cutName.push_back("onePhoton");
    // cutFlow.push_back(*cutFlow_photonPt200<RA2bTree>);
    // cutName.push_back("photonPt200");
    // cutFlow.push_back(*cutFlow_minDR<RA2bTree>);
    // cutName.push_back("minDR");
    cutFlow.push_back(*cutFlow_filters<RA2bTree>);
    cutName.push_back("filters");
    cutFlow.push_back(*cutFlow_btagsZero<RA2bTree>);
    cutName.push_back("btagsZero");
    cutFlow.push_back(*cutFlow_HTbin1<RA2bTree>);
    cutName.push_back("HTbin1");
    cutFlow.push_back(*cutFlow_MHTbin1<RA2bTree>);
    cutName.push_back("MHTbin1");
    cutFlow.push_back(*cutFlow_nJetsTwo<RA2bTree>);
    cutName.push_back("nJetsTwo");

    //vector<vector<plot> > plots;
    vector<vector<plot> > plots;

    for( int iCut = 0 ; iCut < cutFlow.size() ; iCut++ ){
        vector<plot> tempPlots;
        tempPlots.push_back(plot(*fillMHT<RA2bTree>,"MHT_"+cutName[iCut],"MHT [GeV]",20,300.,6300.));
        plots.push_back(tempPlots);
    }

    // background MC samples
    for( int iSample = 0 ; iSample < samples.size() ; iSample++){

        RA2bTree* ntuple = samples[iSample];

        for( int iCut = 0 ; iCut < plots.size() ; iCut++){
            for( int iPlot = 0 ; iPlot < plots[iCut].size() ; iPlot++){
                plots[iCut][iPlot].addNtuple(ntuple,sampleNames[iSample]);
                plots[iCut][iPlot].dataHist = NULL;
            }
        }

        int numEvents = ntuple->fChain->GetEntries();
        ntupleBranchStatus<RA2bTree>(ntuple);
        double weight=1.0;
        for( int iEvt = 0 ; iEvt < min(MAX_EVENTS,numEvents) ; iEvt++ ){
            ntuple->GetEntry(iEvt);
            if( iEvt % 100000 == 0 ) cout << sampleNames[iSample] << ": " << iEvt << "/" << numEvents << endl;
            for( int iCut = 0 ; iCut < cutFlow.size() ; iCut++ ){
                if( ! cutFlow[iCut](ntuple) ) break;
                for( int iPlot = 0 ; iPlot < plots[iCut].size() ; iPlot++ ){
		  //weight = ntuple->Weight*lumi;
                    if( sampleNames[iSample] == "GJets" ){
		      weight*=1;// photonTriggerWeight(ntuple)*GJets0p4Weights(ntuple)*dRweights(ntuple);
                        plots[iCut][iPlot].fill(ntuple,weight);
                    }else 
                        plots[iCut][iPlot].fill(ntuple,weight);
                }
            }
        }
    }

    TFile* outputFile = new TFile("cutFlow_synchTroy_MC2016.root","RECREATE");

    for( int iCut = 0 ; iCut < cutFlow.size() ; iCut++ ){ 
        for( int iPlot = 0 ; iPlot < plots[iCut].size() ; iPlot++ ){ 
            plots[iCut][iPlot].Write();
        }
    }
    outputFile->Close();
}
