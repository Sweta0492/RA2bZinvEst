#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include "TH1D.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "CMS_lumi.cc"

#include <iostream>
#include <vector>
#include <map>

using namespace std;

float lumi=41525.897;

template <typename ntupleType> class plot{

 public: 

  plot( double (*fillerFunc_)(ntupleType*) , TString label_="var"){
    fillerFunc = fillerFunc_;
    label = label_;
    xlabel = "" ; 
    nbins = 40; lower = 200.; upper = 1200.;
    //binEdges = NULL;
    for( int ibin=1 ; ibin<nbins+1 ; ibin++ ){
      binEdges[ibin]=binEdges[ibin-1]+(upper-lower)/nbins;
    }
    stack=new THStack(label+"_stack",label+"_stack");
    dataHist=NULL;
  };

  plot( double (*fillerFunc_)(ntupleType*) , 
	TString label_="var",TString xlabel_="var",
	int nbins_=20 , double lower_=200. , double upper_=2200.){
    fillerFunc = fillerFunc_;
    label = label_;
    xlabel =xlabel_ ; 
    nbins = nbins_; lower = lower_; upper = upper_;
    //binEdges=NULL;
    binEdges.push_back(lower);
    for( int ibin=1 ; ibin<nbins+1 ; ibin++ ){
      binEdges.push_back(binEdges[ibin-1]+(upper-lower)/nbins);
    }
    stack=new THStack(label+"_stack",label+"_stack");
  };

  plot( double (*fillerFunc_)(ntupleType*) , 
	TString label_,TString xlabel_,
	int nbins_ , double* bins_){
    fillerFunc = fillerFunc_;
    label = label_;
    xlabel =xlabel_ ; 
    binEdges.assign(bins_,bins_+nbins_+1);
    nbins = nbins_; lower = binEdges[0]; upper = binEdges[nbins];
    stack=new THStack(label+"_stack",label+"_stack");
  };
  
  void setBinning(int nbins_ , double lower_ , double upper_ ){
    nbins = nbins_ ; lower = lower_ ; upper = upper_ ;
  };
  
  void addNtuple(ntupleType* ntuple_,TString tag="test"){ 
    //cout << "nbins: " << nbins << " lower: " << lower << " upper: " << upper << endl;
    tagMap[ntuple_] = tag;
    if( binEdges.size() )
      histoMap[ntuple_] = new TH1D(label+"_"+tag,xlabel,nbins,&binEdges[0]);
    else
      histoMap[ntuple_] = new TH1D(label+"_"+tag,xlabel,nbins,lower,upper);
    histoMap[ntuple_]->Sumw2();
  };

  void addSignalNtuple(ntupleType* ntuple_,TString tag="test"){ 
    //cout << "nbins: " << nbins << " lower: " << lower << " upper: " << upper << endl;
    tagMap[ntuple_] = tag;
    if( binEdges.size() )
      signalHistoMap[ntuple_] = new TH1D(label+"_"+tag,xlabel,nbins,&binEdges[0]);   
    else
      signalHistoMap[ntuple_] = new TH1D(label+"_"+tag,xlabel,nbins,lower,upper);

    signalHistoMap[ntuple_]->Sumw2();
  };
  
  void addDataNtuple(ntupleType* ntuple_,TString tag="test"){
    //cout << "nbins: " << nbins << " lower: " << lower << " upper: " << upper << endl;
    tagMap[ntuple_] = tag ;
    if( binEdges.size() ){
        dataHist = new TH1D(label+"_"+tag,xlabel,nbins,&binEdges[0]);
    }else
      dataHist = new TH1D(label+"_"+tag,xlabel,nbins,lower,upper);
    dataHist->SetMarkerStyle(8);
  };
  
  int fill(ntupleType* ntuple ){
    if( histoMap.find(ntuple) != histoMap.end() ){
      return histoMap[ntuple]->Fill(fillerFunc(ntuple),ntuple->Weight*lumi);
    }else{
      cout << "plot::fill - ERROR: key not found: " << ntuple << endl;
      return 0;
    }
  };

  int fillSignal(ntupleType* ntuple ){
    if( signalHistoMap.find(ntuple) != signalHistoMap.end() ){
      return signalHistoMap[ntuple]->Fill(fillerFunc(ntuple),ntuple->Weight*lumi);
    }else{
      cout << "plot::fillSignal - ERROR: key not found, " << ntuple << endl;
      return 0;
    }
  };

  int fill(ntupleType* ntuple , float customWeight ){
    if( histoMap.find(ntuple) != histoMap.end() ){
      return histoMap[ntuple]->Fill(fillerFunc(ntuple),customWeight);
    }else{
      cout << "plot::fill - ERROR: key not found: " << ntuple << endl;
      return 0;
    }
  };

  int fillSignal(ntupleType* ntuple , float customWeight ){
    if( signalHistoMap.find(ntuple) != signalHistoMap.end() ){
      return signalHistoMap[ntuple]->Fill(fillerFunc(ntuple),customWeight);
    }else{
      cout << "plot::fillSignal - ERROR: key not found, " << ntuple << endl;
      return 0;
    }
  };

  int fillData(ntupleType* ntuple ){
    //cout << "[[plot::fillData]]" << dataHist << endl;
    if( dataHist ){
      //cout << "result: " << fillerFunc(ntuple) << endl;
      return dataHist->Fill(fillerFunc(ntuple));
    }else{
      //cout << "dataHist: " << dataHist << endl;
      return 0;
    }
  };
    
  int fillData(ntupleType* ntuple , float customWeight ){
      if( dataHist )
          return dataHist->Fill(fillerFunc(ntuple),customWeight);
      else 
          return 0;
  };
  
  void setFillColor(ntupleType* ntuple , int color=1){
    if( histoMap.find(ntuple) != histoMap.end() ){
      histoMap[ntuple]->SetFillColor(color);
      histoMap[ntuple]->SetLineWidth(2);
    }else{
      cout << "plot::setFillColor - ERROR: couldn't find key, " << ntuple << endl;
    }
  };

  void setLineColor(ntupleType* ntuple , int color=1){
    if( signalHistoMap.find(ntuple) != signalHistoMap.end() ){
      signalHistoMap[ntuple]->SetLineColor(color);
      signalHistoMap[ntuple]->SetLineWidth(2);
    }else{
      cout << "plot::setFillColor - ERROR: couldn't find key, " << ntuple << endl;
    }
  };
  
  void buildStack(vector<ntupleType*> ntuples){
    if( stack ){
      for( int iSample = 0 ; iSample < ntuples.size() ; iSample++ ){
	stack->Add(histoMap[ntuples[iSample]]);
      }    
    }
  };

  void buildStack(){
    
    vector<ntupleType*> ntuples;
    for( typename map<ntupleType*,TH1D*>::iterator it = histoMap.begin() ;
	 it != histoMap.end() ; 
	 ++it){
      ntuples.push_back(it->first);
    }
    buildStack(ntuples);
  };

  void buildSum(){
    sum = NULL;
    for( typename map<ntupleType*,TH1D*>::iterator it = histoMap.begin() ;
	 it != histoMap.end() ; 
	 ++it){
      if( sum == NULL ){
	sum = new TH1D(*(it->second));
	sum->SetNameTitle(label+"_sum",label+"_sum");
      }else
	sum->Add(it->second);
    }
  };

  void Write(){

    for( typename map<ntupleType*,TH1D*>::iterator it = histoMap.begin() ;
	 it != histoMap.end() ; 
	 ++it ){
      if( it->second ) 
	it->second->Write();
    }

    for( typename map<ntupleType*,TH1D*>::iterator it = signalHistoMap.begin() ;
	 it != signalHistoMap.end() ; 
	 ++it ){
      if( it->second ) 
	it->second->Write();
    }
    if( dataHist )
      dataHist->Write();

  };

  void Write(ntupleType* ntuple){

    if( histoMap.find(ntuple) != histoMap.end() )
      histoMap[ntuple]->Write();

    if( signalHistoMap.find(ntuple) != signalHistoMap.end() )
      signalHistoMap[ntuple]->Write();

  };

  void WriteData(){
    dataHist->Write();
  };

  void Draw(TCanvas* can,
	    vector<ntupleType*>ntuples,
	    vector<ntupleType*>signalNtuples,
	    TString dir = "./" ){

    if( ! can ) return ;
    
    TPad* topPad = new TPad("topPad","topPad",0.,0.4,.99,.99);
    TPad* botPad = new TPad("botPad","botPad",0.,0.01,.99,.39);
    botPad->SetBottomMargin(0.25);
    botPad->SetTopMargin(0.02);
    topPad->SetTopMargin(0.06);
    topPad->SetBottomMargin(0.17);
    topPad->Draw();
    botPad->Draw();
    topPad->cd();
    
    double max = 0 ;
    if( histoMap.size() ){
      buildStack(ntuples);
      stack->Draw("histo");
      stack->GetXaxis()->SetTitle(xlabel);
      stack->GetXaxis()->SetNdivisions(505);
      stack->GetYaxis()->SetTitle("Events");
      buildSum();
      max = sum->GetMaximum(); 
    }
    
    for(int iSample = 0 ; iSample < signalNtuples.size() ; iSample++){
        TH1D* temp = signalHistoMap[signalNtuples[iSample]];
        if( temp ){
            //temp->Scale(sum->Integral()/temp->Integral());
            temp->Draw("histo,SAME");
            if( temp->GetMaximum() > max ) 
                max = temp->GetMaximum();
        }
    }

    if( dataHist ){
        dataHist->Draw("e1,SAME");
        if( dataHist->GetMaximum() > max ) 
            max = dataHist->GetMaximum();
    }
    
    stack->SetMaximum(max);
    stack->SetMinimum(0.1);

    stack->GetYaxis()->SetLabelFont(63);
    stack->GetYaxis()->SetLabelSize(14);
    stack->GetYaxis()->SetTitleFont(63);
    stack->GetYaxis()->SetTitleSize(20);
    stack->GetYaxis()->SetTitleOffset(1.6);

    stack->GetXaxis()->SetLabelFont(63);
    stack->GetXaxis()->SetLabelSize(14);
    stack->GetXaxis()->SetTitleFont(63);
    stack->GetXaxis()->SetTitleSize(20);
    stack->GetXaxis()->SetTitleOffset(1.7);

    writeExtraText = true;
    extraText="Supplementary";
    char lumiString[4];
    sprintf(lumiString,"%.1f fb^{-1}",lumi/1000.);
    lumi_13TeV = lumiString;
    CMS_lumi( can , 4 , 0 );
    can->Update();
    can->RedrawAxis();
    can->GetFrame()->Draw();
    
    botPad->cd();
    TH1D* ratio = new TH1D(*dataHist);
    ratio->SetNameTitle(sum->GetName()+TString("ratio"),sum->GetTitle());
    ratio->Divide(sum);
    //TGraphAsymmErrors* ratio = new TGraphAsymmErrors(sum,dataHist,"pois");
    ratio->SetMarkerStyle(8);
    ratio->GetYaxis()->SetTitle("Data/MC");
    ratio->GetYaxis()->SetRangeUser(0.2,2.0);
    ratio->GetXaxis()->SetRangeUser(lower,upper);
    ratio->GetXaxis()->SetTitle(xlabel);

    ratio->GetYaxis()->SetLabelFont(63);
    ratio->GetYaxis()->SetLabelSize(14);
    ratio->GetYaxis()->SetTitleFont(63);
    ratio->GetYaxis()->SetTitleSize(20);
    ratio->GetYaxis()->SetTitleOffset(1.6);

    ratio->GetXaxis()->SetLabelFont(63);
    ratio->GetXaxis()->SetLabelSize(14);
    ratio->GetXaxis()->SetTitleFont(63);
    ratio->GetXaxis()->SetTitleSize(20);
    ratio->GetXaxis()->SetTitleOffset(2.3);

    ratio->Draw("e1");
    TLine one(ratio->GetBinCenter(1)-ratio->GetBinWidth(1)/2.,1.,ratio->GetBinCenter(ratio->GetNbinsX())+ratio->GetBinWidth(ratio->GetNbinsX())/2.,1.);
    TLine avg(ratio->GetBinCenter(1)-ratio->GetBinWidth(1)/2.,dataHist->Integral()/sum->Integral(),ratio->GetBinCenter(ratio->GetNbinsX())+ratio->GetBinWidth(ratio->GetNbinsX())/2.,dataHist->Integral()/sum->Integral());
    avg.SetLineColor(2);
    avg.SetLineStyle(2);
    one.SetLineStyle(2);
    one.Draw();
    avg.Draw();

    //can->cd();
    topPad->cd();
    char SF[16];
    sprintf(SF,"data/MC=%1.1f",dataHist->Integral()/sum->Integral());
    //TText* scaleFactor = new TText(ratio->GetBinCenter(1)-ratio->GetBinWidth(1)/2.,2.1,SF);
    TText* scaleFactor = new TText(0.17,0.02,SF);
    scaleFactor->SetNDC();
    scaleFactor->SetTextFont(43);
    scaleFactor->SetTextSize(16);
    scaleFactor->Draw();

    //TCanvas* legCan = new TCanvas("legCan","legCan",500,500);
    TLegend* leg = new TLegend(0.8,.6,.9,.9);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    for(int iSample = 0 ; iSample < ntuples.size() ; iSample++){
      leg->AddEntry(histoMap[ntuples[iSample]],tagMap[ntuples[iSample]],"f");
    }
    for(int iSample = 0 ; iSample < signalNtuples.size() ; iSample++){
      leg->AddEntry(signalHistoMap[signalNtuples[iSample]],tagMap[signalNtuples[iSample]],"f");
    }
    if( dataHist ) 
      leg->AddEntry(dataHist,"data","p");
    leg->Draw();
    //legCan->SaveAs(dir+"/legend.png");

    can->cd();
    topPad->cd();
    gPad->SetLogy(false);
    can->SaveAs(dir+"/"+label+".png");
    can->SaveAs(dir+"/"+label+".eps");
    can->SaveAs(dir+"/"+label+".pdf");
    gPad->SetLogy(true);
    can->SaveAs(dir+"/"+label+"_LogY.png");
    can->SaveAs(dir+"/"+label+"_LogY.eps");
    can->SaveAs(dir+"/"+label+"_LogY.pdf");
  }

  void Print(){
    for( typename map<ntupleType*,TH1D*>::iterator it = histoMap.begin() ; it != histoMap.end() ; it++ ){
      cout << " -- " << it->second->GetName() << " -- " << endl;
      cout << "ntuple: " << it->first << " histo: " << it->second << endl;
      it->second->Print("all");
    }
  }
  
  TString label;
  TString xlabel;
  int nbins;
  double lower,upper;
  vector<double> binEdges;
  map<ntupleType*,TH1D*> histoMap;
  map<ntupleType*,TString> tagMap;
  map<ntupleType*,TH1D*> signalHistoMap;
  TH1D* dataHist{NULL};
  double (*fillerFunc)(ntupleType*);
  THStack* stack{NULL};
  TH1D* sum{NULL}; 
};
