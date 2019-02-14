#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TROOT.h"
#include "THStack.h"
#include "TPad.h"

#include <vector>
#include <map>
#include <iostream>

#include "plotterUtils.cc"
#include "skimSamples.cc"
#include "definitions.h"
#include "RA2bTree.cc"

using namespace std;

static const int MAX_EVENTS=99999999;

void process(int region, string backgroundSample, string dataSample){

    vector<int> trigger_indices = {140,141};

    skimSamples::region reg = static_cast<skimSamples::region>(region);

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
    
    plot MHTplot(*fillMHT<RA2bTree>,"MHT_"+skims.regionNames[region]+"_baseline","H_{T}^{miss} [GeV]",64,200.,1800.);
    plot HTplot(*fillHT<RA2bTree>,"HT_"+skims.regionNames[region]+"_baseline","H_{T} [GeV]",50,300,2800.);
    plot NJetsplot(*fillNJets<RA2bTree>,"NJets_"+skims.regionNames[region]+"_baseline","N_{jet}",14,1.5,15.5);
    plot BTagsplot(*fillBTags<RA2bTree>,"BTags_"+skims.regionNames[region]+"_baseline","N_{b-jet}",6,-0.5,5.5);
    plot Binsplot(*fillRA2b174Bins<RA2bTree>,"AnalysisBins_"+skims.regionNames[region]+"_baseline","Search Bin",174,0.5,174.5);
    plot Bins46plot(*fillRA2b46Bins<RA2bTree>,"AnalysisBins_BTag0_"+skims.regionNames[region]+"_baseline","Search Bin",46,0.5,46.5);
    plot Bins59plot(*fillRA2b59Bins<RA2bTree>,"AnalysisBins_BTag0plusQCDCR_"+skims.regionNames[region]+"_baseline","Search Bin",59,0.5,59.5);
    plot BinsNJet2plot(*fillRA2bNJet2Bins<RA2bTree>,"AnalysisBins_NJet2_"+skims.regionNames[region]+"_baseline","Search Bin",30,0.5,30.5);

    plot DeltaPhi1plot(*fillDeltaPhi1<RA2bTree>,"DeltaPhi1_"+skims.regionNames[region]+"_baseline","#Delta#Phi_{1}",20,0.,3.1415);
    plot DeltaPhi2plot(*fillDeltaPhi2<RA2bTree>,"DeltaPhi2_"+skims.regionNames[region]+"_baseline","#Delta#Phi_{2}",20,0.,3.1415);
    plot DeltaPhi3plot(*fillDeltaPhi3<RA2bTree>,"DeltaPhi3_"+skims.regionNames[region]+"_baseline","#Delta#Phi_{3}",20,0.,3.1415);
    plot DeltaPhi4plot(*fillDeltaPhi4<RA2bTree>,"DeltaPhi4_"+skims.regionNames[region]+"_baseline","#Delta#Phi_{4}",20,0.,3.1415);

    plot PhotonPtplot(*photonPt<RA2bTree>,"PhotonPt_"+skims.regionNames[region]+"_baseline","p_{T,#gamma} [GeV]",50,200,1400);
    plot PhotonEtaplot(*photonEta<RA2bTree>,"PhotonEta_"+skims.regionNames[region]+"_baseline","#eta_{#gamma}",60,-3.,3.);
    plot PhotonPhiplot(*photonPhi<RA2bTree>,"PhotonPhi_"+skims.regionNames[region]+"_baseline","#phi_{#gamma}",70,-3.5,3.5);

    plot JetPtplot(*jetPt<RA2bTree>,"JetPt_"+skims.regionNames[region]+"_baseline","p_{T,leading_Jet} [GeV]",120,20,1400);
    plot JetEtaplot(*jetEta<RA2bTree>,"JetEta_"+skims.regionNames[region]+"_baseline","#eta_{leading_Jet}",60,-3.,3.);
    plot JetPhiplot(*jetPhi<RA2bTree>,"JetPhi_"+skims.regionNames[region]+"_baseline","#phi_{leading_Jet}",70,-3.5,3.5);
    plot JetSubleadPtplot(*jetSubleadPt<RA2bTree>,"JetSubleadPt_"+skims.regionNames[region]+"_baseline","p_{T,Sub-leading_Jet} [GeV]",120,20,1400);
    plot JetSubleadEtaplot(*jetSubleadEta<RA2bTree>,"JetSubleadEta_"+skims.regionNames[region]+"_baseline","#eta_{Sub-leading_Jet}",60,-3.,3.);
    plot JetSubleadPhiplot(*jetSubleadPhi<RA2bTree>,"JetSubleadPhi_"+skims.regionNames[region]+"_baseline","#phi_{Sub-leading_Jet}",70,-3.5,3.5);

    plot PhotonMinDeltaR(*fillRecoPhotonDeltaR<RA2bTree>,"PhotonMinDeltaR_"+skims.regionNames[region]+"_baseline","min#Delta R(jet,#gamma)",40,0,4);

    plot verticesplot(*fillNumVertices<RA2bTree>,"NumVertices_"+skims.regionNames[region]+"_baseline","n_{vtx}",40,0,80);

    vector<plot> plotsAllEvents;
    plotsAllEvents.push_back(MHTplot);
    plotsAllEvents.push_back(HTplot);
    plotsAllEvents.push_back(NJetsplot);
    plotsAllEvents.push_back(BTagsplot);
    plotsAllEvents.push_back(Binsplot);
    plotsAllEvents.push_back(Bins46plot);
    plotsAllEvents.push_back(Bins59plot);
    plotsAllEvents.push_back(BinsNJet2plot);
    plotsAllEvents.push_back(DeltaPhi1plot);
    plotsAllEvents.push_back(DeltaPhi2plot);
    plotsAllEvents.push_back(DeltaPhi3plot);
    plotsAllEvents.push_back(DeltaPhi4plot);
    plotsAllEvents.push_back(PhotonPtplot);
    plotsAllEvents.push_back(PhotonEtaplot);
    plotsAllEvents.push_back(PhotonPhiplot);

    plotsAllEvents.push_back(JetPtplot);
    plotsAllEvents.push_back(JetEtaplot);
    plotsAllEvents.push_back(JetPhiplot);

    plotsAllEvents.push_back(JetSubleadPtplot);
    plotsAllEvents.push_back(JetSubleadEtaplot);
    plotsAllEvents.push_back(JetSubleadPhiplot);

    plotsAllEvents.push_back(PhotonMinDeltaR);
    plotsAllEvents.push_back(verticesplot);

    plot MHTplotEB(*fillMHT<RA2bTree>,"MHT_"+skims.regionNames[region]+"_baseline_EB","H_{T}^{miss} [GeV]",64,200.,1800.);
    plot HTplotEB(*fillHT<RA2bTree>,"HT_"+skims.regionNames[region]+"_baseline_EB","H_{T} [GeV]",50,300,2800.);
    plot NJetsplotEB(*fillNJets<RA2bTree>,"NJets_"+skims.regionNames[region]+"_baseline_EB","N_{jet}",14,1.5,15.5);
    plot BTagsplotEB(*fillBTags<RA2bTree>,"BTags_"+skims.regionNames[region]+"_baseline_EB","N_{b-jet}",6,-0.5,5.5);
    plot BinsplotEB(*fillRA2b174Bins<RA2bTree>,"AnalysisBins_"+skims.regionNames[region]+"_baseline_EB","Search Bin",174,0.5,174.5);
    plot Bins46plotEB(*fillRA2b46Bins<RA2bTree>,"AnalysisBins_BTag0_"+skims.regionNames[region]+"_baseline_EB","Search Bin",46,0.5,46.5);
    plot Bins59plotEB(*fillRA2b59Bins<RA2bTree>,"AnalysisBins_BTag0plusQCDCR_"+skims.regionNames[region]+"_baseline_EB","Search Bin",59,0.5,59.5);
    plot BinsNJet2plotEB(*fillRA2bNJet2Bins<RA2bTree>,"AnalysisBins_NJet2_"+skims.regionNames[region]+"_baseline_EB","Search Bin",30,0.5,30.5);

    plot DeltaPhi1plotEB(*fillDeltaPhi1<RA2bTree>,"DeltaPhi1_"+skims.regionNames[region]+"_baseline_EB","#Delta#Phi_{1}",20,0.,3.1415);
    plot DeltaPhi2plotEB(*fillDeltaPhi2<RA2bTree>,"DeltaPhi2_"+skims.regionNames[region]+"_baseline_EB","#Delta#Phi_{2}",20,0.,3.1415);
    plot DeltaPhi3plotEB(*fillDeltaPhi3<RA2bTree>,"DeltaPhi3_"+skims.regionNames[region]+"_baseline_EB","#Delta#Phi_{3}",20,0.,3.1415);
    plot DeltaPhi4plotEB(*fillDeltaPhi4<RA2bTree>,"DeltaPhi4_"+skims.regionNames[region]+"_baseline_EB","#Delta#Phi_{4}",20,0.,3.1415);

    plot PhotonPtplotEB(*photonPt<RA2bTree>,"PhotonPt_"+skims.regionNames[region]+"_baseline_EB","p_{T,#gamma} [GeV]",50,200,1400);
    plot PhotonEtaplotEB(*photonEta<RA2bTree>,"PhotonEta_"+skims.regionNames[region]+"_baseline_EB","#eta_{#gamma}",60,-3.,3.);
    plot PhotonPhiplotEB(*photonPhi<RA2bTree>,"PhotonPhi_"+skims.regionNames[region]+"_baseline_EB","#phi_{#gamma}",70,-3.5,3.5);

    plot JetPtplotEB(*jetPt<RA2bTree>,"JetPt_"+skims.regionNames[region]+"_baseline_EB","p_{T,leading_Jet} [GeV]",120,20,1400);
    plot JetEtaplotEB(*jetEta<RA2bTree>,"JetEta_"+skims.regionNames[region]+"_baseline_EB","#eta_{leading_Jet}",60,-3.,3.);
    plot JetPhiplotEB(*jetPhi<RA2bTree>,"JetPhi_"+skims.regionNames[region]+"_baseline_EB","#phi_{leading_Jet}",70,-3.5,3.5);
    plot JetSubleadPtplotEB(*jetSubleadPt<RA2bTree>,"JetSubleadPt_"+skims.regionNames[region]+"_baseline_EB","p_{T,Sub-leading_Jet} [GeV]",120,20,1400);
    plot JetSubleadEtaplotEB(*jetSubleadEta<RA2bTree>,"JetSubleadEta_"+skims.regionNames[region]+"_baseline_EB","#eta_{Sub-leading_Jet}",60,-3.,3.);
    plot JetSubleadPhiplotEB(*jetSubleadPhi<RA2bTree>,"JetSubleadPhi_"+skims.regionNames[region]+"_baseline_EB","#phi_{Sub-leading_Jet}",70,-3.5,3.5);

    plot PhotonMinDeltaREB(*fillRecoPhotonDeltaR<RA2bTree>,"PhotonMinDeltaR_"+skims.regionNames[region]+"_baseline_EB","min#Delta R(jet,#gamma)",40,0,4);
  
    plot verticesplotEB(*fillNumVertices<RA2bTree>,"NumVertices_"+skims.regionNames[region]+"_baseline_EB","n_{vtx}",40,0,80);

    vector<plot> plotsEBevents;
    plotsEBevents.push_back(MHTplotEB);
    plotsEBevents.push_back(HTplotEB);
    plotsEBevents.push_back(NJetsplotEB);
    plotsEBevents.push_back(BTagsplotEB);
    plotsEBevents.push_back(BinsplotEB);
    plotsEBevents.push_back(Bins46plotEB);
    plotsEBevents.push_back(Bins59plotEB);
    plotsEBevents.push_back(BinsNJet2plotEB);
    plotsEBevents.push_back(DeltaPhi1plotEB);
    plotsEBevents.push_back(DeltaPhi2plotEB);
    plotsEBevents.push_back(DeltaPhi3plotEB);
    plotsEBevents.push_back(DeltaPhi4plotEB);
    plotsEBevents.push_back(PhotonPtplotEB);
    plotsEBevents.push_back(PhotonEtaplotEB);
    plotsEBevents.push_back(PhotonPhiplotEB);

    plotsEBevents.push_back(JetPtplotEB);
    plotsEBevents.push_back(JetEtaplotEB);
    plotsEBevents.push_back(JetPhiplotEB);

    plotsEBevents.push_back(JetSubleadPtplotEB);
    plotsEBevents.push_back(JetSubleadEtaplotEB);
    plotsEBevents.push_back(JetSubleadPhiplotEB);
    plotsEBevents.push_back(PhotonMinDeltaREB);
    plotsEBevents.push_back(verticesplotEB);

    plot MHTplotEE(*fillMHT<RA2bTree>,"MHT_"+skims.regionNames[region]+"_baseline_EE","H_{T}^{miss} [GeV]",64,200.,1800.);
    plot HTplotEE(*fillHT<RA2bTree>,"HT_"+skims.regionNames[region]+"_baseline_EE","H_{T} [GeV]",50,300,2800.);
    plot NJetsplotEE(*fillNJets<RA2bTree>,"NJets_"+skims.regionNames[region]+"_baseline_EE","N_{jet}",14,1.5,15.5);
    plot BTagsplotEE(*fillBTags<RA2bTree>,"BTags_"+skims.regionNames[region]+"_baseline_EE","N_{b-jet}",6,-0.5,5.5);
    plot BinsplotEE(*fillRA2b174Bins<RA2bTree>,"AnalysisBins_"+skims.regionNames[region]+"_baseline_EE","Search Bin",174,0.5,174.5);
    plot Bins46plotEE(*fillRA2b46Bins<RA2bTree>,"AnalysisBins_BTag0_"+skims.regionNames[region]+"_baseline_EE","Search Bin",46,0.5,46.5);
    plot Bins59plotEE(*fillRA2b59Bins<RA2bTree>,"AnalysisBins_BTag0plusQCDCR_"+skims.regionNames[region]+"_baseline_EE","Search Bin",59,0.5,59.5);
    plot BinsNJet2plotEE(*fillRA2bNJet2Bins<RA2bTree>,"AnalysisBins_NJet2_"+skims.regionNames[region]+"_baseline_EE","Search Bin",30,0.5,30.5);

    plot DeltaPhi1plotEE(*fillDeltaPhi1<RA2bTree>,"DeltaPhi1_"+skims.regionNames[region]+"_baseline_EE","#Delta#Phi_{1}",20,0.,3.1415);
    plot DeltaPhi2plotEE(*fillDeltaPhi2<RA2bTree>,"DeltaPhi2_"+skims.regionNames[region]+"_baseline_EE","#Delta#Phi_{2}",20,0.,3.1415);
    plot DeltaPhi3plotEE(*fillDeltaPhi3<RA2bTree>,"DeltaPhi3_"+skims.regionNames[region]+"_baseline_EE","#Delta#Phi_{3}",20,0.,3.1415);
    plot DeltaPhi4plotEE(*fillDeltaPhi4<RA2bTree>,"DeltaPhi4_"+skims.regionNames[region]+"_baseline_EE","#Delta#Phi_{4}",20,0.,3.1415);

    plot PhotonPtplotEE(*photonPt<RA2bTree>,"PhotonPt_"+skims.regionNames[region]+"_baseline_EE","p_{T,#gamma} [GeV]",50,200,1400);
    plot PhotonEtaplotEE(*photonEta<RA2bTree>,"PhotonEta_"+skims.regionNames[region]+"_baseline_EE","#eta_{#gamma}",60,-3.,3.);
    plot PhotonPhiplotEE(*photonPhi<RA2bTree>,"PhotonPhi_"+skims.regionNames[region]+"_baseline_EE","#phi_{#gamma}",70,-3.5,3.5);

 
    plot JetPtplotEE(*jetPt<RA2bTree>,"JetPt_"+skims.regionNames[region]+"_baseline_EE","p_{T,leading_Jet} [GeV]",120,20,1400);
    plot JetEtaplotEE(*jetEta<RA2bTree>,"JetEta_"+skims.regionNames[region]+"_baseline_EE","#eta_{leading_Jet}",60,-3.,3.);
    plot JetPhiplotEE(*jetPhi<RA2bTree>,"JetPhi_"+skims.regionNames[region]+"_baseline_EE","#phi_{leading_Jet}",70,-3.5,3.5);
    plot JetSubleadPtplotEE(*jetSubleadPt<RA2bTree>,"JetSubleadPt_"+skims.regionNames[region]+"_baseline_EE","p_{T,Sub-leading_Jet} [GeV]",120,20,1400);
    plot JetSubleadEtaplotEE(*jetSubleadEta<RA2bTree>,"JetSubleadEta_"+skims.regionNames[region]+"_baseline_EE","#eta_{Sub-leading_Jet}",60,-3.,3.);
    plot JetSubleadPhiplotEE(*jetSubleadPhi<RA2bTree>,"JetSubleadPhi_"+skims.regionNames[region]+"_baseline_EE","#phi_{Sub-leading_Jet}",70,-3.5,3.5);

    plot PhotonMinDeltaREE(*fillRecoPhotonDeltaR<RA2bTree>,"PhotonMinDeltaR_"+skims.regionNames[region]+"_baseline_EE","min#Delta R(jet,#gamma)",40,0,4);

    plot verticesplotEE(*fillNumVertices<RA2bTree>,"NumVertices_"+skims.regionNames[region]+"_baseline_EE","n_{vtx}",40,0,80);

    vector<plot> plotsEEevents;
    plotsEEevents.push_back(MHTplotEE);
    plotsEEevents.push_back(HTplotEE);
    plotsEEevents.push_back(NJetsplotEE);
    plotsEEevents.push_back(BTagsplotEE);
    plotsEEevents.push_back(BinsplotEE);
    plotsEEevents.push_back(Bins46plotEE);
    plotsEEevents.push_back(Bins59plotEE);
    plotsEEevents.push_back(BinsNJet2plotEE);
    plotsEEevents.push_back(DeltaPhi1plotEE);
    plotsEEevents.push_back(DeltaPhi2plotEE);
    plotsEEevents.push_back(DeltaPhi3plotEE);
    plotsEEevents.push_back(DeltaPhi4plotEE);
    plotsEEevents.push_back(PhotonPtplotEE);
    plotsEEevents.push_back(PhotonEtaplotEE);
    plotsEEevents.push_back(PhotonPhiplotEE);

    plotsEEevents.push_back(JetPtplotEE);
    plotsEEevents.push_back(JetEtaplotEE);
    plotsEEevents.push_back(JetPhiplotEE);
    plotsEEevents.push_back(JetSubleadPtplotEE);
    plotsEEevents.push_back(JetSubleadEtaplotEE);
    plotsEEevents.push_back(JetSubleadPhiplotEE);
    
    plotsEEevents.push_back(PhotonMinDeltaREE);
    plotsEEevents.push_back(verticesplotEE);

    // background MC samples
    for( int iSample = 0 ; iSample < skims.ntuples.size() ; iSample++){

      TFile* outputFile = new TFile("plotObs_"+skims.regionNames[region]+"_baseline_"+skims.sampleName[iSample]+".root","UPDATE");

      RA2bTree* ntuple = skims.ntuples[iSample];

      for( int iPlot = 0 ; iPlot < plotsAllEvents.size() ; iPlot++){
	plotsAllEvents[iPlot].addNtuple(ntuple,skims.sampleName[iSample]);
	plotsEBevents[iPlot].addNtuple(ntuple,skims.sampleName[iSample]);
	plotsEEevents[iPlot].addNtuple(ntuple,skims.sampleName[iSample]);
      }

      int numEvents = ntuple->fChain->GetEntries();
      ntupleBranchStatus<RA2bTree>(ntuple);
      double weight;
      for( int iEvt = 0 ; iEvt < min(MAX_EVENTS,numEvents) ; iEvt++ ){

	//cout << "sample: " << skims.sampleName[iSample] << endl;
	ntuple->GetEntry(iEvt);
        weight = 1.0;
	if( iEvt % 1000000 == 0 ) cout << skims.sampleName[iSample] << ": " << iEvt << "/" << numEvents << endl;

	if( ( reg == skimSamples::kPhoton || reg == skimSamples::kPhotonLoose ) && !RA2bBaselinePhotonCut(ntuple) ) continue;
	if( ( reg == skimSamples::kPhotonLDP ) && !RA2bLDPBaselinePhotonCut(ntuple) ) continue;

	if( skims.regionNames[region] == "photonLDP" || skims.regionNames[region] == "photon" || skims.regionNames[region] == "photonLoose" ){
	  if( skims.sampleName[iSample] == "QCD" && isPromptPhoton(ntuple) ) continue;
	  if( skims.sampleName[iSample] == "GJets" && ( !isPromptPhoton(ntuple) || ntuple->madMinPhotonDeltaR < 0.4 ) ) continue;
	}

	// ----------- weights -----------------
	weight = lumi*ntuple->Weight*ntuple->NonPrefiringProb;
        if( skims.sampleName[iSample] == "GJets" ) weight *= dRweights(ntuple)*trig_eff(ntuple,iEvt);
        
	// -------------------------------------

	for( int iPlot = 0 ; iPlot < plotsAllEvents.size() ; iPlot++ ){
	  plotsAllEvents[iPlot].fill(ntuple,weight);
	  if( ntuple->Photons_isEB->at(0) ){
	    plotsEBevents[iPlot].fill(ntuple,weight);
	  }else{
	    plotsEEevents[iPlot].fill(ntuple,weight);
	  }
	}// end loop over plots
      }// end loop over number of events

      for( int iPlot = 0 ; iPlot < plotsAllEvents.size() ; iPlot++){
	plotsAllEvents[iPlot].Write(ntuple);
	plotsEBevents[iPlot].Write(ntuple);
	plotsEEevents[iPlot].Write(ntuple);
      }
      
      outputFile->Close();
    }// end loop over samples

    
    // Data samples
    for( int iSample = 0 ; iSample < skims.dataNtuple.size() ; iSample++){
      
      TFile* outputFile = new TFile("plotObs_"+skims.regionNames[region]+"_baseline_"+skims.dataSampleName[iSample]+".root","UPDATE");
      RA2bTree* ntuple = skims.dataNtuple[iSample];
      for( int iPlot = 0 ; iPlot < plotsAllEvents.size() ; iPlot++){
        plotsAllEvents[iPlot].addDataNtuple(ntuple,skims.dataSampleName[iSample]);
        plotsEBevents[iPlot].addDataNtuple(ntuple,skims.dataSampleName[iSample]);
        plotsEEevents[iPlot].addDataNtuple(ntuple,skims.dataSampleName[iSample]);
      }// end loop over plots
  
      int numEvents = ntuple->fChain->GetEntries();
      ntupleBranchStatus<RA2bTree>(ntuple);
      for( int iEvt = 0 ; iEvt < min(MAX_EVENTS,numEvents) ; iEvt++ ){
        ntuple->GetEntry(iEvt);
        if( iEvt % 1000000 == 0 ) cout << "data: " << iEvt << "/" << numEvents << endl;
        if( ( reg == skimSamples::kPhotonLDP ) && !RA2bLDPBaselinePhotonCut(ntuple) ) continue;   
        if( ( reg == skimSamples::kPhoton || reg == skimSamples::kPhotonLoose ) && !RA2bBaselinePhotonCut(ntuple) ) continue;
        
         /*......................Trigger Weight ...............................*/
	bool pass_trigger=false;
	for( unsigned int itrig = 0 ; itrig < trigger_indices.size() ; itrig++ )
	  pass_trigger|=(ntuple->TriggerPass->at(trigger_indices[itrig])==1);
	if( !pass_trigger ) continue;

        for( int iPlot = 0 ; iPlot < plotsAllEvents.size() ; iPlot++){
	  plotsAllEvents[iPlot].fillData(ntuple);
	  if( ntuple->Photons_isEB->at(0) )
	    plotsEBevents[iPlot].fillData(ntuple);
	  else
	    plotsEEevents[iPlot].fillData(ntuple);
        }// end loop over plots
      }// end loop over number of events
      for( int iPlot = 0 ; iPlot < plotsAllEvents.size() ; iPlot++){
	plotsAllEvents[iPlot].WriteData();
	plotsEBevents[iPlot].WriteData();
	plotsEEevents[iPlot].WriteData();
      }// end loop over plots

      outputFile->Close();
    }// end loop over samples

}

int main(int argc, char** argv){

    gROOT->ProcessLine(".L tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    trig_eff_func(); // trigger Efficiency Function
    if(argc != 4 ){
      cerr << "please provide a region index, a background sample name, and a data sample name" << endl;
    }

    process(atoi(argv[1]),argv[2],argv[3]);

}
