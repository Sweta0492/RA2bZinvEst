from ROOT import *
gROOT.SetBatch(True)
from optparse import OptionParser
from math import sqrt

# - - - - ROOT STYLES  - - - - - -
gROOT.ProcessLine(".L ~/tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")

# - - - - - COMMAND LINE PARSER - - - - - 
parser = OptionParser()
parser.add_option("-s", "--saveTag", dest="saveTag", default="inclusive",
                  help="tag added to output file names")
parser.add_option("-e", "--encap", dest="endcap", default=False,action="store_true",
                  help="fit endcap photons")
parser.add_option("-b", "--binLabel", dest="binLabel", default="MHT_200",
                  help="histo bin label")
parser.add_option("-d", "--inputDir", dest="inputDir", default="",
                  help="path to the location of input root files")
(options, args) = parser.parse_args()

saveTag = options.saveTag
endcap = options.endcap
binLabel = options.binLabel
inputDir = options.inputDir
print "saveTag:",saveTag
print "endcap:",endcap
print "binLabel:",binLabel
print "inputDir:",inputDir
EEorEB=""
if endcap : 
    EEorEB="EE"
else:
    EEorEB="EB"

################################################
# - - - - - - - - - FITTER - - - - - - - - - - #
################################################
def runFit( promptHisto,nonPromptHisto,dataHisto,isEndcap,tag) : 

    IsoChrg = RooRealVar("iso","Iso_{chrg} [GeV]",0,10.)
    if isEndcap : 
        IsoChrg.setRange("SR",0.0,1.694);
    else :
        IsoChrg.setRange("SR",0.0,2.089);

    promptRooHist = RooDataHist("promptRooHist","promptaRooHist",RooArgList(IsoChrg),promptHisto)
    promptPdf = RooHistPdf("promptPdf","promptPdf",RooArgSet(IsoChrg),promptRooHist)

    nonpromptRooHist = RooDataHist("nonpromptRooHist","nonpromptaRooHist",RooArgList(IsoChrg),nonPromptHisto)
    nonpromptPdf = RooHistPdf("nonpromptPdf","nonpromptPdf",RooArgSet(IsoChrg),nonpromptRooHist)
    
    dataRooHist = RooDataHist("dataRooHist","dataRooHist",RooArgList(IsoChrg),dataHisto)

    frac = RooRealVar("frac","frac",.85,0.,1.)
    
    model = RooAddPdf("model","model",promptPdf,nonpromptPdf,frac)

    RooFit.gErrorIgnoreLevel=0
    model.fitTo(dataRooHist,RooFit.Verbose(kFALSE))

    can = TCanvas("can","can",500,500)
    
    fitPad = TPad("fit","fit",0.01,0.3,0.99,0.99)
    fitPad.SetTopMargin(.1)
    fitPad.Draw()
    ratioPad = TPad("ratio","ratio",0.01,0.01,0.99,0.29)
    ratioPad.SetBottomMargin(.07)
    ratioPad.SetTopMargin(0.)
    ratioPad.Draw()
    
    fitPad.cd()
    
    plot = IsoChrg.frame()
    dataRooHist.plotOn(plot)
    model.plotOn(plot)
    nonpromptPdf.plotOn(plot,RooFit.FillColor(2),RooFit.FillStyle(3144),RooFit.LineColor(2),RooFit.Normalization(1.0-frac.getVal()))
    plot.GetXaxis().SetNdivisions(505)
    plot.Draw()
    ratioPad.cd()

    modelHisto = TH1D(promptHisto)
    modelHisto.Scale(dataHisto.Integral()*frac.getVal()/promptHisto.Integral())
    modelHisto.Add(nonPromptHisto,float(dataHisto.Integral()*(1.-frac.getVal())/nonPromptHisto.Integral()))
    
    # for i in range(1,modelHisto.GetNbinsX()+1):
    #     print "modelHisto:",modelHisto.GetBinContent(i),"+/-",modelHisto.GetBinError(i)
    #     print "dataHisto:",dataHisto.GetBinContent(i),"+/-",dataHisto.GetBinError(i)
    #     print "r:",modelHisto.GetBinContent(i)/dataHisto.GetBinContent(i),"+/-",modelHisto.GetBinContent(i)/dataHisto.GetBinContent(i)*sqrt(modelHisto.GetBinError(i)*modelHisto.GetBinError(i)/modelHisto.GetBinContent(i)/modelHisto.GetBinContent(i)+dataHisto.GetBinError(i)*dataHisto.GetBinError(i)/dataHisto.GetBinContent(i)/dataHisto.GetBinContent(i))

    ratio = TGraphAsymmErrors(modelHisto,dataHisto,"pois")
    ratio.SetMarkerStyle(8)
    ratio.GetYaxis().SetRangeUser(0.,3.)
    ratio.GetXaxis().SetTitleOffset(1.5)
    ratio.GetXaxis().SetTitleSize(1.2)
    ratio.GetXaxis().SetLabelOffset(1.5)
    ratio.GetXaxis().SetNdivisions(505)
    ratio.GetYaxis().SetLabelSize(.1)
    ratio.GetYaxis().SetTitleSize(.15)
    ratio.GetYaxis().SetTitleOffset(.5)
    ratio.GetYaxis().SetTitle("Model/Data")
    ratio.Draw("Ap")

    promptIntegral   = promptPdf.createIntegral(RooArgSet(IsoChrg));
    promptSRIntegral = promptPdf.createIntegral(RooArgSet(IsoChrg),"SR");
    nonpromptIntegral = nonpromptPdf.createIntegral(RooArgSet(IsoChrg));
    nonpromptSRIntegral = nonpromptPdf.createIntegral(RooArgSet(IsoChrg),"SR");
    modelIntegral    = model.createIntegral(RooArgSet(IsoChrg));
    modelSRIntegral  = model.createIntegral(RooArgSet(IsoChrg),"SR");

    r_prompt = promptSRIntegral.getVal()/promptIntegral.getVal() 
    r_nonprompt = nonpromptSRIntegral.getVal()/nonpromptIntegral.getVal()

    print "REGION: {0} {1}".format(EEorEB,binLabel)
    purity=frac.getVal()/( (1.-frac.getVal()) * r_nonprompt/r_prompt + frac.getVal() )
    print "purity in SR: ",purity,"+/-",frac.getError()/frac.getVal()*purity

    ##############################################################
    # - - - - - - - - - - - make pretty plots - - - - - - - - -  #
    ##############################################################
    can.SetTopMargin(.1)

    cmsText = TPaveText(.17,.9,.47,.96,"NDC")
    cmsText.SetFillColor(0)
    cmsText.SetBorderSize(0)
    cmsText.AddText("CMS Preliminary")

    fitPad.cd()
    cmsText.Draw()
    
    lumiText = TPaveText(.72,.9,.97,.96,"NDC")
    lumiText.SetFillColor(0)
    lumiText.SetBorderSize(0)
    lumiText.AddText("35.9 fb^{-1} (13 TeV)")
    
    fitPad.cd()
    lumiText.Draw()

    can.SaveAs("../plots/purityFits/purityFits_{0}{1}_{2}.png".format(tag,EEorEB,saveTag))
    can.SaveAs("../plots/purityFits/purityFits_{0}{1}_{2}.eps".format(tag,EEorEB,saveTag))
    can.SaveAs("../plots/purityFits/purityFits_{0}{1}_{2}.pdf".format(tag,EEorEB,saveTag))

    fitPad.SetLogy()

    can.SaveAs("../plots/purityFits/purityFits_{0}{1}_{2}_LogY.png".format(tag,EEorEB,saveTag))
    can.SaveAs("../plots/purityFits/purityFits_{0}{1}_{2}_LogY.eps".format(tag,EEorEB,saveTag))
    can.SaveAs("../plots/purityFits/purityFits_{0}{1}_{2}_LogY.pdf".format(tag,EEorEB,saveTag))

def main() :
    ###########################################################
    # - - - - - - - - - - grab histo inputs - - - - - - - - - #
    ###########################################################
    inputFile = TFile("purityInputs.root","READ")
    prompt_hist = inputFile.Get("photonIsoChrgLowSieie_{0}_{1}_photonLoose_GJets".format(EEorEB,binLabel))
    assert(prompt_hist != None)
        
    ### ADD HACK BECAUSE THERE AREN'T ENOUGH EVENTS IN THE
    ### MC ALT MHT_600 BIN....
    if binLabel == "MHT_600" :
        nonPrompt_hist = inputFile.Get("photonIsoChrgHighSieie_{0}_{1}_photonLoose_Data".format(EEorEB,"MHT_350"))
        assert( nonPrompt_hist != None )
        nonPromptMC_hist = inputFile.Get("photonIsoChrgLowSieie_{0}_{1}_photonLoose_QCD".format(EEorEB,"MHT_350"))
        assert( nonPromptMC_hist != None )
        nonPromptMCALT_hist = inputFile.Get("photonIsoChrgHighSieie_{0}_{1}_photonLoose_QCD".format(EEorEB,"MHT_350"))
        assert( nonPromptMCALT_hist != None )
    else : 
        nonPrompt_hist = inputFile.Get("photonIsoChrgHighSieie_{0}_{1}_photonLoose_Data".format(EEorEB,binLabel))
        assert( nonPrompt_hist != None )
        nonPromptMC_hist = inputFile.Get("photonIsoChrgLowSieie_{0}_{1}_photonLoose_QCD".format(EEorEB,binLabel))
        assert( nonPromptMC_hist != None )
        nonPromptMCALT_hist = inputFile.Get("photonIsoChrgHighSieie_{0}_{1}_photonLoose_QCD".format(EEorEB,binLabel))
        assert( nonPromptMCALT_hist != None )
    data_hist = inputFile.Get("photonIsoChrgLowSieie_{0}_{1}_photonLoose_Data".format(EEorEB,binLabel))
    assert(data_hist!=None)

    ####################################################
     # - - - - - - - - - run fits - - - - - - - - - - - #
     ####################################################
    runFit(prompt_hist,nonPromptMC_hist,data_hist,endcap,"ALT_")
    runFit(prompt_hist,nonPrompt_hist,data_hist,endcap,"")
    runFit(prompt_hist,nonPromptMCALT_hist,data_hist,endcap,"MC_ALT_")

if __name__ == "__main__":
    main()
