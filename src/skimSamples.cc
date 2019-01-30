#include "TChain.h"
#include "RA2bTree.cc"
#include "TString.h"

#include <iostream>
#include <vector>

static const TString BASE_DIR = "root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV16/";
static const bool RE_MINIAOD = false;
 
class skimSamples{

public : 

    TChain *WJets, *ZJets, *QCD, *TT,*GJets,*GJets0p4,*Other,*DY; 
    TChain *mGluino1300, *mGluino1400, *mGluino1500, *mGluino1600, *mGluino1700;
    TChain *mHiggsino1000,*mHiggsino900,*mHiggsino800,*mHiggsino700,*mHiggsino600,*mHiggsino500,*mHiggsino400,*mHiggsino300;
    TChain *data;
    std::vector<RA2bTree*> ntuples,signalNtuples;
    RA2bTree* dataNtuple;
    std::vector<TString> sampleName, signalSampleName;
    std::vector<TString> dataSampleName; 
    std::vector<int> fillColor, lineColor;
    
    enum region {kSignal,kLDP,kPhoton,kPhotonLoose,kPhotonLDP,kPhotonLDPLoose,kDYe,kDYm,kDYeLDP,kDYmLDP,kNumRegions};
    TString regionNames[kNumRegions]={"signal","LDP","photon","photonLoose","photonLDP","photonLDPLoose","DYe","DYm","DYeLDP","DYmLDP"};

    TString skimType;

    skimSamples(region r=kPhoton){

        skimType="";

        if( r == kSignal )
            skimType=BASE_DIR+"tree_signal/";
        if( r == kLDP )
            skimType=BASE_DIR+"tree_LDP/";
        if( r == kPhoton )
            skimType=BASE_DIR+"tree_GJet_CleanVars/";
        if( r == kPhotonLDP )
            skimType=BASE_DIR+"tree_GJetLDP_CleanVars/";
        if( r == kPhotonLoose )
            skimType=BASE_DIR+"tree_GJetLoose_CleanVars/";
        if( r == kPhotonLDPLoose )
            skimType=BASE_DIR+"tree_GJetLooseLDP_CleanVars/";
        if( r == kDYm )
            skimType=BASE_DIR+"tree_DYm_CleanVars";
        if( r == kDYmLDP )
            skimType=BASE_DIR+"tree_DYmLDP_CleanVars";
        if( r == kDYe )
            skimType=BASE_DIR+"tree_DYe_CleanVars";
        if( r == kDYeLDP )
           skimType=BASE_DIR+"tree_DYeLDP_CleanVars";
        vector<TString> ZJetsFileNames;
        ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-100to200_MC2017.root");
        ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-200to400_MC2017.root");
        ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-400to600_MC2017.root");
        ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-600to800_MC2017.root");
        ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-800to1200_MC2017.root");
        ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-1200to2500_MC2017.root");
        ZJetsFileNames.push_back("tree_ZJetsToNuNu_HT-2500toInf_MC2017.root");
        ZJets = new TChain("tree");
        for( int i = 0 ; i < ZJetsFileNames.size() ; i++ ){
            ZJets->Add(skimType+"/"+ZJetsFileNames[i]);
        }
        if( r == kSignal || r == kLDP ){ 
            ntuples.push_back(new RA2bTree(ZJets));
            sampleName.push_back("ZJets");
            fillColor.push_back(kGreen+1);
        }

        vector<TString> QCDFileNames;
        QCDFileNames.push_back("tree_QCD_HT-200to300_MC2017.root");
        QCDFileNames.push_back("tree_QCD_HT-300to500_MC2017.root");
        QCDFileNames.push_back("tree_QCD_HT-500to700_MC2017.root");
        QCDFileNames.push_back("tree_QCD_HT-700to1000_MC2017.root");
        QCDFileNames.push_back("tree_QCD_HT-1000to1500_MC2017.root");
        QCDFileNames.push_back("tree_QCD_HT-1500to2000_MC2017.root");
        QCDFileNames.push_back("tree_QCD_HT-2000toInf_MC2017.root");
        QCD = new TChain("tree");
        for( int i = 0 ; i < QCDFileNames.size() ; i++ ){
            QCD->Add(skimType+"/"+QCDFileNames[i]);
        }
        if( r == kSignal || r == kLDP || r == kPhoton || r == kPhotonLDP || r == kPhotonLoose || r == kPhotonLDPLoose ){ 
            ntuples.push_back(new RA2bTree(QCD));
            sampleName.push_back("QCD");
            fillColor.push_back(kGray);
        }

        vector<TString> WJetsFileNames;
        WJetsFileNames.push_back("tree_WJetsToLNu_HT-100to200_MC2017.root");
        WJetsFileNames.push_back("tree_WJetsToLNu_HT-1200to2500_MC2017.root");
        WJetsFileNames.push_back("tree_WJetsToLNu_HT-200to400_MC2017.root");
        WJetsFileNames.push_back("tree_WJetsToLNu_HT-2500toInf_MC2017.root");
        WJetsFileNames.push_back("tree_WJetsToLNu_HT-400to600_MC2017.root");
        WJetsFileNames.push_back("tree_WJetsToLNu_HT-600to800_MC2017.root");
        WJetsFileNames.push_back("tree_WJetsToLNu_HT-800to1200_MC2017.root");
        WJets = new TChain("tree");
        for( int i = 0 ; i < WJetsFileNames.size() ; i++ ){
            WJets->Add(skimType+"/"+WJetsFileNames[i]);
        }
        if( r == kSignal || r == kLDP ){
            ntuples.push_back(new RA2bTree(WJets));
            sampleName.push_back("WJets");
            fillColor.push_back(kBlue);
        }

       vector<TString> TTFileNames;
        TTFileNames.push_back("tree_TTJets_HT-1200to2500_MC2017.root");
        TTFileNames.push_back("tree_TTJets_HT-2500toInf_MC2017.root");
        TTFileNames.push_back("tree_TTJets_HT-600to800_MC2017.root");
        TTFileNames.push_back("tree_TTJets_HT-800to1200_MC2017.root");
        TT = new TChain("tree");
        for( int i = 0 ; i < TTFileNames.size() ; i++ ){
            TT->Add(skimType+"/"+TTFileNames[i]);
        }
        ntuples.push_back(new RA2bTree(TT));
        sampleName.push_back("TT");
        fillColor.push_back(kCyan);


        vector<TString> OtherFileNames;
        OtherFileNames.push_back("tree_ST_s-channel_MC2017.root");
        OtherFileNames.push_back("tree_ST_t-channel_antitop_MC2017.root");
        OtherFileNames.push_back("tree_ST_t-channel_top_MC2017.root");
        OtherFileNames.push_back("tree_ST_tW_antitop_MC2017.root");
        OtherFileNames.push_back("tree_ST_tW_top_MC2017.root");
        OtherFileNames.push_back("tree_WWTo1L1Nu2Q_MC2017.root");
       // OtherFileNames.push_back("tree_WWTo2L2Nu_MC2017.root");
       // OtherFileNames.push_back("tree_WWZ_MC2017.root");
        OtherFileNames.push_back("tree_WZTo1L1Nu2Q_MC2017.root");
        OtherFileNames.push_back("tree_WZTo1L3Nu_MC2017.root");
        OtherFileNames.push_back("tree_WZZ_MC2017.root");
        OtherFileNames.push_back("tree_ZZTo2L2Q_MC2017.root");
        //OtherFileNames.push_back("tree_ZZTo2Q2Nu_MC2017.root");
        OtherFileNames.push_back("tree_ZZZ_MC2017.root");
        OtherFileNames.push_back("tree_TTTT_MC2017.root");
        OtherFileNames.push_back("tree_TTWJetsToLNu_MC2017.root");
        OtherFileNames.push_back("tree_TTWJetsToQQ_MC2017.root");
        OtherFileNames.push_back("tree_TTGJets_MC2017.root");
        OtherFileNames.push_back("tree_TTZToLLNuNu_MC2017.root");
        OtherFileNames.push_back("tree_TTZToQQ_MC2017.root");
        Other = new TChain("tree");
        for( int i = 0 ; i < OtherFileNames.size() ; i++ ){
            Other->Add(skimType+"/"+OtherFileNames[i]);
        }
        ntuples.push_back(new RA2bTree(Other));
        sampleName.push_back("Other");
        fillColor.push_back(kRed+1);

        vector<TString> DYFileNames;
        DYFileNames.push_back("tree_DYJetsToLL_M-50_HT-100to200_MC2017.root");
        DYFileNames.push_back("tree_DYJetsToLL_M-50_HT-200to400_MC2017.root");
        DYFileNames.push_back("tree_DYJetsToLL_M-50_HT-400to600_MC2017.root");
        DYFileNames.push_back("tree_DYJetsToLL_M-50_HT-600toInf_MC2017.root");
        DY = new TChain("tree");
        for( int i = 0 ; i < DYFileNames.size() ; i++ ){
            DY->Add(skimType+"/"+DYFileNames[i]);
        }
        if( r == kDYe || r == kDYm || r == kDYeLDP || r == kDYmLDP ){
            ntuples.push_back(new RA2bTree(DY));
            sampleName.push_back("DY");
            fillColor.push_back(kGreen);
        }

        vector<TString> GJets0p4FileNames;
        /*GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-100to200_MC2017.root");
        GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-200to400_MC2017.root");
        GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-400to600_MC2017.root");
        GJets0p4FileNames.push_back("tree_GJets_DR-0p4_HT-600toInf_MC2017.root");
      */  GJets0p4FileNames.push_back("tree_GJets_HT-100to200_MC2017.root");
        GJets0p4FileNames.push_back("tree_GJets_HT-200to400_MC2017.root");
        GJets0p4FileNames.push_back("tree_GJets_HT-400to600_MC2017.root");
        GJets0p4FileNames.push_back("tree_GJets_HT-600toInf_MC2017.root");
        GJets0p4 = new TChain("tree");
        for( int i = 0 ; i < GJets0p4FileNames.size() ; i++ ){
            GJets0p4->Add(skimType+"/"+GJets0p4FileNames[i]);
        }
        if( r == kPhoton || r== kPhotonLDP || r == kPhotonLoose || r == kPhotonLDPLoose ){
        ntuples.push_back(new RA2bTree(GJets0p4));
        sampleName.push_back("GJets");
        fillColor.push_back(kGreen);
        }

        vector<TString> GJetsFileNames;
        GJetsFileNames.push_back("tree_GJets_HT-100to200_MC2017.root");
        GJetsFileNames.push_back("tree_GJets_HT-200to400_MC2017.root");
        GJetsFileNames.push_back("tree_GJets_HT-400to600_MC2017.root");
        GJetsFileNames.push_back("tree_GJets_HT-600toInf_MC2017.root");
        GJets = new TChain("tree");
        for( int i = 0 ; i < GJetsFileNames.size() ; i++ ){
       GJets->Add(skimType+"/"+GJetsFileNames[i]);
}
        vector<TString> HTMHTFileNames;
        TString basename;
        if( RE_MINIAOD )
            basename="tree_HTMHT_re";
        else 
            basename="tree_HTMHT_";
        HTMHTFileNames.push_back(basename+"2017B.root");
        HTMHTFileNames.push_back(basename+"2017C.root");
        HTMHTFileNames.push_back(basename+"2017D.root");
        HTMHTFileNames.push_back(basename+"2017E.root");
        HTMHTFileNames.push_back(basename+"2017F.root");

        if( r == kSignal || r == kLDP ){
            data = new TChain("tree");
            for( int i = 0 ; i < HTMHTFileNames.size() ; i++ ){
                data->Add(skimType+"/"+HTMHTFileNames[i]);
            }    
            dataNtuple = new RA2bTree(data); 
        }

        vector<TString> SingleElectronNames;
        if( RE_MINIAOD )
            basename="tree_SingleElectron_re";
        else 
            basename="tree_SingleElectron_";
        SingleElectronNames.push_back(basename+"2017B.root");
        SingleElectronNames.push_back(basename+"2017C.root");
        SingleElectronNames.push_back(basename+"2017D.root");
        SingleElectronNames.push_back(basename+"2017E.root");
        SingleElectronNames.push_back(basename+"2017F.root");

        if( r == kDYe || r == kDYeLDP ){
            data = new TChain("tree");
            for( int i = 0 ; i < SingleElectronNames.size() ; i++ ){
                data->Add(skimType+"/"+SingleElectronNames[i]);
            }
            dataNtuple = new RA2bTree(data);
        }

        vector<TString> SingleMuonNames;
        if( RE_MINIAOD )
            basename="tree_SingleMuon_re";
        else 
            basename="tree_SingleMuon_";
        SingleMuonNames.push_back(basename+"2017B.root");
        SingleMuonNames.push_back(basename+"2017C.root");
        SingleMuonNames.push_back(basename+"2017D.root");
        SingleMuonNames.push_back(basename+"2017E.root");
        SingleMuonNames.push_back(basename+"2017F.root");

        if( r == kDYm || r == kDYmLDP ){
            data = new TChain("tree");
            for( int i = 0 ; i < SingleMuonNames.size() ; i++ ){
                data->Add(skimType+"/"+SingleMuonNames[i]);
            }
            dataNtuple = new RA2bTree(data);
        }

        vector<TString> SinglePhotonFileNames;
        if( RE_MINIAOD )
            basename="tree_SinglePhoton_re";
        else
            basename="tree_SinglePhoton_";

        SinglePhotonFileNames.push_back(basename+"2017B.root");
        SinglePhotonFileNames.push_back(basename+"2017C.root");
        SinglePhotonFileNames.push_back(basename+"2017D.root");
        SinglePhotonFileNames.push_back(basename+"2017E.root");
        SinglePhotonFileNames.push_back(basename+"2017F.root");
    
       if( r == kPhoton || r == kPhotonLDP || r == kPhotonLoose || r == kPhotonLDPLoose ){
            data = new TChain("tree");
            for( int i = 0 ; i < SinglePhotonFileNames.size() ; i++ ){
                data->Add(skimType+"/"+SinglePhotonFileNames[i]);
            }
            dataNtuple = new RA2bTree(data);
        }

    };

    void dumpRegions(){
        cout << "index  name" << endl;
        for( int iReg = 0 ; iReg < kNumRegions ; iReg++ ){
            cout << iReg << "      " << regionNames[iReg] << endl;
        }
    };

    RA2bTree* findNtuple(TString name){
        for( int iSam = 0 ; iSam < sampleName.size() ; iSam++ ){
            if( sampleName[iSam] == name )
                return ntuples[iSam] ;
        }
        for( int iSam = 0 ; iSam < signalSampleName.size() ; iSam++ ){
            if( signalSampleName[iSam] == name )
                return signalNtuples[iSam] ;
        }
        return NULL;
    };

};
