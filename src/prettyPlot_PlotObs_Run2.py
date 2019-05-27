import CMS_lumi
from ROOT import *
from array import array
from sys import argv
import ROOT as r
r.gROOT.SetBatch(True)

r.gROOT.ProcessLine(".L ~/tdrstyle.C")
r.gROOT.ProcessLine("setTDRStyle()")
plot_dir="plotObs_plots"
input_file_name = "Combined.root"

input_file = r.TFile(input_file_name,"READ")    

def plot(plot_var = "BTags_photon_baseline_EE" ):
#def plot(plot_var = "BTags_photon_baseline" ):
    
    data_samples=["Data"]
    '''samples=[["QCD"],
             ["GJets"],
             ["TT"],
             ["Others"]]

    samples_labels = ["QCD","GJets","TT","Others"]
    samples_fill_color = [r.kGray,r.kGreen,r.kCyan,r.kRed+1]
    samples_line_color = [1,1,1,1]'''
    
    samples=[["Others"],
             ["TT"],
             ["QCD"],
             ["GJets"]]

    samples_labels = ["Others","TT","QCD","GJets"]
    samples_fill_color = [r.kRed+1,r.kCyan,r.kGray,r.kGreen]
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
    botPad.SetTopMargin(0.0);
    topPad.SetTopMargin(0.06);
    topPad.SetBottomMargin(0.0);
    '''botPad.SetBottomMargin(0.35);
    botPad.SetTopMargin(0.02);
    topPad.SetTopMargin(0.06);
    topPad.SetBottomMargin(0.17);'''
    topPad.Draw();
    botPad.Draw();
    topPad.cd();
    
    stack.Draw("histo")    
    data_histo[0].Draw("SAME,pe")

    stack.SetMaximum(2*max(sum.GetMaximum(),data_histo[0].GetMaximum()))
    stack.SetMinimum(1)
    
    stack.GetYaxis().SetTitle("Events")
    stack.GetYaxis().SetLabelSize(0.055*1.15);
    stack.GetYaxis().SetTitleFont(42);
    stack.GetYaxis().SetTitleSize(0.06*1.15);
    stack.GetYaxis().SetTitleOffset(0.9);

    stack.GetXaxis().SetLabelSize(0.055*1.15);
    stack.GetXaxis().SetTitleFont(42);
    stack.GetXaxis().SetTitleSize(0.06*1.15);
    stack.GetXaxis().SetTitleOffset(0.9);

    '''stack.GetYaxis().SetTitle("Events")
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
    stack.GetXaxis().SetTitleOffset(1.7);'''

    CMS_lumi.extraText = "       Preliminary"
    CMS_lumi.writeExtraText = True
    CMS_lumi.channelText = ""
    CMS_lumi.CMS_lumi(can, 4, 0)

    '''CMStext = r.TText(.17,.95,"CMS")
    CMStext.SetNDC()
    CMStext.SetTextFont(61)
    CMStext.SetTextSize(0.07)
    CMStext.Draw()
    
    SIMtext = r.TText(.28,.95,"preliminary")
    SIMtext.SetNDC()
    SIMtext.SetTextFont(52)
    SIMtext.SetTextSize(0.07)
    SIMtext.Draw()
    
    LUMItext = r.TText(.65,.95,"13 TeV (137 X/fb)")
    LUMItext.SetNDC()
    LUMItext.SetTextFont(51)
    LUMItext.SetTextSize(0.07)
    LUMItext.Draw()'''
   
    SF = (1.1*(data_histo[0].Integral()/sum.Integral()))
    scaleFactor = r.TText(.17,.02,"Data/MC = "+str(round(data_histo[0].Integral()/sum.Integral(),1)))
    scaleFactor.SetNDC()
    scaleFactor.SetTextFont(43)
    scaleFactor.SetTextSize(16)
    #scaleFactor.Draw()

    leg = TLegend(.68,.71,.89,.93)    # for twiki
    leg.SetTextSize(0.03);
    leg.SetLineWidth(1)
    leg.SetFillStyle(1001)
    leg.SetFillColor(0)
    leg.SetLineColor(1)
    leg.SetBorderSize(1)
    leg.SetFillColor(0)
    
    leg.AddEntry(data_histo[0],"Data","pe")
    leg.AddEntry(samples_histo[3],"Signal \gamma+jets","f")
    leg.AddEntry(samples_histo[2],"QCD","f")
    leg.AddEntry(samples_histo[1],"t#bar{t}","f")
    leg.AddEntry(samples_histo[0],"Others","f") 
    leg.Draw(); 

    botPad.cd()
    ratio = r.TH1D(data_histo[0])

    ratio.GetYaxis().SetLabelSize(0.055*1.15);
    ratio.GetYaxis().SetTitleFont(42);
    ratio.GetYaxis().SetTitleSize(0.09*1.15);
    ratio.GetYaxis().SetTitleOffset(0.5);

    ratio.GetXaxis().SetLabelSize(0.055*1.15);
    ratio.GetXaxis().SetTitleFont(42);
    ratio.GetXaxis().SetTitleSize(0.09*1.15);
    ratio.GetXaxis().SetTitleOffset(0.8); 
    
    ratio.SetMarkerStyle(8)
    ratio.SetMarkerSize(.8)
    ratio.SetName(plot_var+"_ratio")
    ratio.Divide(sum)
    ratio.GetYaxis().SetRangeUser(0,1.9)
    ratio.GetYaxis().SetTitle("Data/MC")
    ratio.GetXaxis().SetTitle(data_histo[0].GetTitle())

    '''ratio.GetYaxis().SetLabelFont(63);
    ratio.GetYaxis().SetLabelSize(14);
    ratio.GetYaxis().SetTitleFont(63);
    ratio.GetYaxis().SetTitleSize(20);
    ratio.GetYaxis().SetNdivisions(505);
    ratio.GetYaxis().SetTitleOffset(1.6);

    ratio.GetXaxis().SetLabelFont(63);
    ratio.GetXaxis().SetLabelSize(14);
    ratio.GetXaxis().SetTitleFont(63);
    ratio.GetXaxis().SetTitleSize(20);
    ratio.GetXaxis().SetTitleOffset(2.3);'''

    ratio.Draw()
    
    one = r.TLine(ratio.GetBinCenter(1)-ratio.GetBinWidth(1)/2.,1.,ratio.GetBinCenter(ratio.GetNbinsX())+ratio.GetBinWidth(ratio.GetNbinsX())/2.,1.);
    avg = r.TLine(ratio.GetBinCenter(1)-ratio.GetBinWidth(1)/2.,data_histo[0].Integral()/sum.Integral(),ratio.GetBinCenter(ratio.GetNbinsX())+ratio.GetBinWidth(ratio.GetNbinsX())/2.,data_histo[0].Integral()/sum.Integral());
    avg.SetLineColor(2);
    avg.SetLineStyle(2);
    one.SetLineStyle(2);
    one.Draw();
    avg.Draw();


    topPad.SetLogy()
    can.SaveAs(plot_var+"_Run2_LogY.pdf")

    output_file.cd()
    for h in samples_histo :
        r.TH1D(h).Write()
    data_histo[0].Write()

output_file = r.TFile("baselineInputs_Run2.root","RECREATE")

vars = []
list = input_file.GetListOfKeys()
next = r.TIter(list);
key = next()
while(key != None ) :
    obj = key.ReadObj();
    if obj.InheritsFrom("TH1") :
        name = r.TString(obj.GetName())
        if name.Contains("_Others") : 
            #print name
            #print "integral:",obj.Integral()
            vars.append(name.ReplaceAll("_Others","").Data())
    else :
        print obj.Print()
    key = next()

print vars
for var in vars : 
    plot(var)

    
output_file.Close()
input_file.Close()
