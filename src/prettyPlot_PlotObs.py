from ROOT import *
from array import array
from sys import argv
import ROOT as r
r.gROOT.SetBatch(True)

r.gROOT.ProcessLine(".L ~/tdrstyle.C")
r.gROOT.ProcessLine("setTDRStyle()")
plot_dir="plotObs_plots"
input_file_name = "plotObs_photon_baseline.root"

input_file = r.TFile(input_file_name,"READ")    

def plot(plot_var = "MHT_photon_baseline_EE" ):
    
    samples=[["QCD_1000to1500",
              "QCD_1500to2000",
              "QCD_2000toInf",
              "QCD_200to300",
              "QCD_300to500",
              "QCD_500to700",
              "QCD_700to1000"],
             ["GJets0p4_100to200",
              "GJets0p4_200to400",
              "GJets0p4_400to600",
              "GJets0p4_600toInf"],
             ["TT_600to800",
              "TT_800to1200",
              "TT_1200to2500",
              "TT_2500toInf"],
             ["Other_WWTo1L1Nu2Q",
              "Other_WZTo1L3Nu",
              "Other_ZZTo2L2Q",
              "Other_TTWJetsToLNu",
              "Other_TTWJetsToQQ",
              "Other_TTGJets",
              "Other_TTZToLLNuNu",
              "Other_TTZToQQ",
              "Other_ST_s",
              "Other_ST_tW_antitop",
              "Other_ST_tW_top"]]
    
    data_samples=["SinglePhoton_2018A",
                  "SinglePhoton_2018B",
                  "SinglePhoton_2018C",
                  "SinglePhoton_2018D"]
    samples_labels = ["QCD","GJets","TT","Others"]
    samples_fill_color = [r.kGray,r.kGreen,r.kCyan,r.kRed+1]
    samples_line_color = [1,1,1,1]
    
    samples_histo=[]
    
    sum = None
    stack = r.THStack("stack","stack")
    
    for i,sample_names in enumerate(samples) :   
        for j,sample_name in enumerate(sample_names): 
            if len(samples_histo) <= i : 
                samples_histo.append(input_file.Get(plot_var+"_"+sample_name))
                if samples_histo[-1]==None :
                    print "looking for:",plot_var+"_"+sample_name
                    print input_file.ls(plot_var+"_"+sample_name)
                    input_file.ls()
                    assert(samples_histo[-1]!=None)
                elif samples_histo[-1].Integral() < 0.0001 :
                    print "oops.",plot_var+"_"+sample_name,"is empty"
                samples_histo[-1].SetLineColor(samples_line_color[i])
                samples_histo[-1].SetFillColor(samples_fill_color[i])
                samples_histo[-1].SetName(plot_var+"_"+samples_labels[i])
            else : 
                samples_histo[-1].Add(input_file.Get(plot_var+"_"+sample_name))

    data_histo=[]
    for i,s in enumerate(data_samples):
        if i == 0 : 
            data_histo.append(input_file.Get(plot_var+"_"+s))
            data_histo[-1].SetMarkerStyle(8)
            data_histo[-1].SetMarkerSize(.8)
            data_histo[-1].SetName(plot_var+"_Data")
            if data_histo[-1]==None :
                print "looking for:",plot_var+"_"+sample_name
                input_file.ls(plot_var+"*")
                input_file.ls("*"+sample_name)
                assert(data_histo[-1]!=None)
        else : 
            data_histo[-1].Add(input_file.Get(plot_var+"_"+s))

    for i,h in enumerate(samples_histo) : 
        stack.Add(h)
        if i == 0 : 
            sum = r.TH1D(h)
            sum.SetName(plot_var+"_sum")
        else : 
            sum.Add(h)

    can = r.TCanvas("can","can",500,500)
    topPad = r.TPad("topPad","topPad",0.,0.4,.99,.99);
    botPad = r.TPad("botPad","botPad",0.,0.01,.99,.39);
    botPad.SetBottomMargin(0.35);
    botPad.SetTopMargin(0.02);
    topPad.SetTopMargin(0.06);
    topPad.SetBottomMargin(0.17);
    topPad.Draw();
    botPad.Draw();
    topPad.cd();
    
    stack.Draw("histo")    
    data_histo[0].Draw("SAME,p")

    stack.SetMaximum(1.3*max(sum.GetMaximum(),data_histo[0].GetMaximum()))
    stack.SetMinimum(0.1)

    stack.GetYaxis().SetTitle("Events")
    stack.GetXaxis().SetTitle(data_histo[0].GetTitle())
    stack.GetYaxis().SetLabelFont(63);
    stack.GetYaxis().SetLabelSize(14);
    stack.GetYaxis().SetTitleFont(63);
    stack.GetYaxis().SetTitleSize(20);
    stack.GetYaxis().SetTitleOffset(1.6);

    stack.GetXaxis().SetLabelFont(63);
    stack.GetXaxis().SetLabelSize(14);
    stack.GetXaxis().SetTitleFont(63);
    stack.GetXaxis().SetTitleSize(20);
    stack.GetXaxis().SetTitleOffset(1.7);

    CMStext = r.TText(.17,.95,"CMS")
    CMStext.SetNDC()
    CMStext.SetTextFont(61)
    CMStext.SetTextSize(0.07)
    CMStext.Draw()
    
    SIMtext = r.TText(.28,.95,"preliminary")
    SIMtext.SetNDC()
    SIMtext.SetTextFont(52)
    SIMtext.SetTextSize(0.07)
    SIMtext.Draw()
    
    #LUMItext = r.TText(.65,.95,"13 TeV (21X/fb)")
    LUMItext = r.TText(.65,.95,"13 TeV (59.6X/fb)")
    LUMItext.SetNDC()
    LUMItext.SetTextFont(51)
    LUMItext.SetTextSize(0.07)
    LUMItext.Draw()
   
    SF = (1.1*(data_histo[0].Integral()/sum.Integral()))
    scaleFactor = r.TText(.17,.02,"data/MC = "+str(round(data_histo[0].Integral()/sum.Integral(),1)))
    scaleFactor.SetNDC()
    scaleFactor.SetTextFont(43)
    scaleFactor.SetTextSize(16)
    scaleFactor.Draw()

    leg = TLegend(.8,.6,.9,.9)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    for i,h in enumerate(samples_histo) :
	    leg.AddEntry(h,samples_labels[i],"f")
    
    leg.AddEntry(data_histo[0],"data","p")
    leg.Draw(); 


    botPad.cd()
    ratio = r.TH1D(data_histo[0])
    ratio.SetMarkerStyle(8)
    ratio.SetMarkerSize(.8)
    ratio.SetName(plot_var+"_ratio")
    ratio.Divide(sum)
    ratio.GetYaxis().SetRangeUser(0,2)
    ratio.GetYaxis().SetTitle("Obs/Exp")
    ratio.GetXaxis().SetTitle(data_histo[0].GetTitle())
    ratio.GetYaxis().SetLabelFont(63);
    ratio.GetYaxis().SetLabelSize(14);
    ratio.GetYaxis().SetTitleFont(63);
    ratio.GetYaxis().SetTitleSize(20);
    ratio.GetYaxis().SetNdivisions(505);
    ratio.GetYaxis().SetTitleOffset(1.6);

    ratio.GetXaxis().SetLabelFont(63);
    ratio.GetXaxis().SetLabelSize(14);
    ratio.GetXaxis().SetTitleFont(63);
    ratio.GetXaxis().SetTitleSize(20);
    ratio.GetXaxis().SetTitleOffset(2.3);

    ratio.Draw()
    
    one = r.TLine(ratio.GetBinCenter(1)-ratio.GetBinWidth(1)/2.,1.,ratio.GetBinCenter(ratio.GetNbinsX())+ratio.GetBinWidth(ratio.GetNbinsX())/2.,1.);
    avg = r.TLine(ratio.GetBinCenter(1)-ratio.GetBinWidth(1)/2.,data_histo[0].Integral()/sum.Integral(),ratio.GetBinCenter(ratio.GetNbinsX())+ratio.GetBinWidth(ratio.GetNbinsX())/2.,data_histo[0].Integral()/sum.Integral());
    avg.SetLineColor(2);
    avg.SetLineStyle(2);
    one.SetLineStyle(2);
    one.Draw();
    avg.Draw();


    can.SaveAs("../plots/"+plot_dir+"/"+plot_var+".png")
    can.SaveAs("../plots/"+plot_dir+"/"+plot_var+".pdf")
    topPad.SetLogy()
    can.SaveAs("../plots/"+plot_dir+"/"+plot_var+"_LogY.png")
    can.SaveAs("../plots/"+plot_dir+"/"+plot_var+"_LogY.pdf")

    output_file.cd()
    for h in samples_histo :
        r.TH1D(h).Write()
    data_histo[0].Write()

output_file = r.TFile("baselineInputs.root","RECREATE")

vars = []
list = input_file.GetListOfKeys()
next = r.TIter(list);
key = next()
while(key != None ) :
    obj = key.ReadObj();
    if obj.InheritsFrom("TH1") :
        name = r.TString(obj.GetName())
        if name.Contains("_QCD_300to500") : 
            #print name
            #print "integral:",obj.Integral()
            vars.append(name.ReplaceAll("_QCD_300to500","").Data())
    else :
        print obj.Print()
    key = next()

print vars
for var in vars : 
    plot(var)

    
output_file.Close()
input_file.Close()
