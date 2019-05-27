import CMS_lumi
from ROOT import *
from array import array
from sys import argv
import ROOT as r
import RA2b
r.gROOT.SetBatch(True)

r.gROOT.ProcessLine(".L ~/tdrstyle.C")
r.gROOT.ProcessLine("setTDRStyle()")
def addDivisionsLow(ymin, ymax, axisTitleSize):
  # Putting lines and labels explaining search region definitions
  ls = 0.5

  # Njet separation lines
  tl_njet = r.TLine();
  tl_njet.SetLineStyle(2);
  tl_njet.DrawLine(10.+ls,ymin,10.+ls,10*ymax);
  tl_njet.DrawLine(20.+ls,ymin,20.+ls,10*ymax); 
  tl_njet.DrawLine(30.+ls,ymin,30.+ls,10*ymax);
  tl_njet.DrawLine(38.+ls,ymin,38.+ls,10*ymax);

  # Nb separation lines

def addDivisionsUp(ymin, ymax, axisTitleSize):
  # Putting lines and labels explaining search region definitions
  ymax2 = ymax/10
  ls = 0.5
  # Njet labels
  ttext_njet = r.TLatex();
  ttext_njet.SetTextFont(42);
  ttext_njet.SetTextSize(0.6*axisTitleSize);
  ttext_njet.SetTextAlign(22);
  ttext_njet.DrawLatex(6.2 , ymax/4.7 , "2 #leq N_{#scale[0.2]{ }jet} #leq 3");
  ttext_njet.DrawLatex(15.5 , ymax/4.7 , "4 #leq N_{#scale[0.2]{ }jet} #leq 5");
  ttext_njet.DrawLatex(25.5 , ymax/4.7 , "6 #leq N_{#scale[0.2]{ }jet} #leq 7");
  ttext_njet.DrawLatex(34.6 , ymax/4.7 , "8 #leq N_{#scale[0.2]{ }jet} #leq 9");
  ttext_njet.DrawLatex(42. , ymax/4.7 , "N_{#scale[0.2]{ }jet} #geq 10");

  # Njet separation lines
  tl_njet = r.TLine();
  tl_njet.SetLineStyle(2);
  tl_njet.DrawLine(10.+ls,ymin/10,10.+ls,0.4*ymax);
  tl_njet.DrawLine(20.+ls,ymin/10,20.+ls,0.4*ymax); 
  tl_njet.DrawLine(30.+ls,ymin/10,30.+ls,0.4*ymax);
  tl_njet.DrawLine(38.+ls,ymin/10,38.+ls,0.4*ymax);

  # Nb labels
  ttext_nb = r.TLatex();
  ttext_nb.SetTextFont(42);
  ttext_nb.SetTextSize(1.3*axisTitleSize);
  ttext_nb.SetTextAlign(22);
    
plot_dir="plotObs_plots"
#input_file_name = "Combined_forpaper_morecompact.root"
input_file_name = "Combined.root"

input_file_name_RzG = "RzGamma_Run2_Only_dRweight.root"
input_file = r.TFile(input_file_name,"READ")
input_file_RzG = r.TFile(input_file_name_RzG,"READ")
hist = input_file_RzG.Get('h')

def plot(plot_var = "AnalysisBins_BTag0_photon_baseline" ):

    '''samples=[["QCD"],
             ["TT"],
             ["Others"],
             ["GJets"]]
    samples_labels = ["QCD","TT","Others","GJets"]
    samples_fill_color = [r.kGray,r.kCyan,r.kRed+1,r.kGreen]'''

    samples=[["Others"],
             ["TT"],
             ["QCD"],
             ["GJets"]]

    samples_labels = ["Others","TT","QCD","GJets"]
    samples_fill_color = [r.kRed+1,r.kCyan,r.kGray,r.kGreen]
    samples_line_color = [1,1,1,1]

    data_samples=["Data"]
    
    '''samples=[["QCD"],
             ["GJets"]]

    samples_labels = ["QCD","GJets"]
    samples_fill_color = [r.kGray,r.kGreen]
    samples_line_color = [1,1]'''

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
            h.Scale(1.36353619186)
            sum = r.TH1D(h)
            sum.SetName(plot_var+"_sum")
        else : 
            h.Scale(1.36353619186)
            sum.Add(h)

    can = r.TCanvas("can","can",500,500)
    topPad = r.TPad("topPad","topPad",0.,0.4,.99,.99);
    botPad = r.TPad("botPad","botPad",0.,0.01,.99,.39);
    botPad.SetBottomMargin(0.35);
    botPad.SetTopMargin(0.0);
    topPad.SetTopMargin(0.06);
    topPad.SetBottomMargin(0.0);
    topPad.Draw();
    botPad.Draw();
    topPad.cd();
    
    stack.Draw("histo")    
    data_histo[0].Draw("SAME,p")

    # Upper panel
    p1b = 0.4  # pad1 vertical dimensions
    p1t = 0.99
    p2b = 0.01  # pad2 vertical dimensions
    p2t = 0.39
    p2mt = 0.11  # pad2 vertical margins
    p2mb = 0.005
    p2min = .5
    p2max = 5.0e6
    axisTitleSize = 0.055
    titleScale = (p1t-p1b)/(p2t-p2b)
    print "title scale factor = "+str(titleScale)
    axisTitleSize *= titleScale
    addDivisionsUp(p2min, p2max, axisTitleSize)

    stack.SetMaximum(8*max(sum.GetMaximum(),data_histo[0].GetMaximum()))
    stack.SetMinimum(1)

    #stack.GetYaxis().SetTitle("Events / bin   ")
    stack.GetYaxis().SetTitle("Events")


    #stack.GetYaxis().SetLabelFont(63);
    #stack.GetYaxis().SetLabelSize(14);
    #stack.GetYaxis().SetTitleFont(42);
    #stack.GetYaxis().SetTitleSize(20);
    #stack.GetYaxis().SetTitleOffset(1.6);

    #stack.GetXaxis().SetLabelFont(63);
    #stack.GetXaxis().SetLabelSize(14);
    #stack.GetXaxis().SetTitleFont(42);
    #stack.GetXaxis().SetTitleSize(20);
    #stack.GetXaxis().SetTitleOffset(1.7);

    stack.GetYaxis().SetLabelSize(0.055*1.15);
    stack.GetYaxis().SetTitleFont(42);
    stack.GetYaxis().SetTitleSize(0.06*1.15);
    stack.GetYaxis().SetTitleOffset(0.9);

    stack.GetXaxis().SetLabelSize(0.055*1.15);
    stack.GetXaxis().SetTitleFont(42);
    stack.GetXaxis().SetTitleSize(0.06*1.15);
    stack.GetXaxis().SetTitleOffset(0.9);

    CMStext = r.TText(.17,.95,"CMS")
    CMStext.SetNDC()
    CMStext.SetTextFont(61)
    CMStext.SetTextSize(0.07)
    #CMStext.Draw()
    
    #SIMtext = r.TText(.28,.95,"preliminary")
    SIMtext = r.TText(.28,.95,"")
    SIMtext.SetNDC()
    SIMtext.SetTextFont(52)
    SIMtext.SetTextSize(0.07)
    #SIMtext.Draw()
    
    #LUMItext = r.TText(.65,.95,"137 #fb^{-1} (13 Tev)") # old triveni
    #stamp = r.TLatex()
    #LUMItext = stamp.DrawLatex(.65,.95,"137 #fb^{-1} (13 TeV)")
    #LUMItext.SetNDC()
    #LUMItext.SetTextFont(51)
    #LUMItext.SetTextSize(0.07)
    #LUMItext.Draw()


    CMS_lumi.extraText = "  "
    CMS_lumi.writeExtraText = True
    CMS_lumi.channelText = ""
    CMS_lumi.CMS_lumi(can, 4, 0)
   
    SF = (1.1*(data_histo[0].Integral()/sum.Integral()))
    #scaleFactor = r.TText(.17,.02,"data/MC = "+str(round(data_histo[0].Integral()/sum.Integral(),1)))
    scaleFactor = r.TText(.17,.35,"  Simulation")
    scaleFactor.SetNDC()
    scaleFactor.SetTextFont(63)
    scaleFactor.SetTextSize(14)
    scaleFactor.Draw()

    #leg = TLegend(.7,.7,.9,.87)   # for paper
    leg = TLegend(.74,.66,.91,.88)    # for twiki
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
    ratio = r.TH1D(hist)
    ratio.SetMarkerStyle(8)
    ratio.SetMarkerSize(.8)
    ratio.SetName(plot_var+"_ratio")
    ratio.GetYaxis().SetRangeUser(0,1.39)
    ratio.GetYaxis().SetTitle("R_{Z->#nu#bar{#nu}/\gamma}^{sim}")
    ratio.GetXaxis().SetTitle("                                          N_{jet}, H_{T}^{miss}, H_{T} bin index")
    
    ratio.GetYaxis().SetLabelSize(0.055*1.15);
    ratio.GetYaxis().SetTitleFont(42);
    ratio.GetYaxis().SetTitleSize(0.09*1.15);
    ratio.GetYaxis().SetTitleOffset(0.6);

    ratio.GetXaxis().SetLabelSize(0.055*1.15);
    ratio.GetXaxis().SetTitleFont(42);
    ratio.GetXaxis().SetTitleSize(0.09*1.15);
    ratio.GetXaxis().SetTitleOffset(1.);
    
    #ratio.GetYaxis().SetLabelFont(63);
    #ratio.GetYaxis().SetLabelSize(14);
    #ratio.GetYaxis().SetTitleFont(42);
    #ratio.GetYaxis().SetTitleSize(20);
    #ratio.GetYaxis().SetNdivisions(505);
    #ratio.GetYaxis().SetTitleOffset(1.6);

    #ratio.GetXaxis().SetLabelFont(63);
    #ratio.GetXaxis().SetLabelSize(14);
    #ratio.GetXaxis().SetTitleFont(63);
    #ratio.GetXaxis().SetTitleSize(20);
    #ratio.GetXaxis().SetTitleOffset(2.3);

    ratio.Draw()
    # Lower panel
    p1b = 0.4  # pad1 vertical dimensions
    p1t = 0.99
    p1mt = 0.005  # pad1 vertical margins
    p1mb = 0.25
    p1min = 0.001
    p1max = 2
    axisTitleSize = 0.055
    addDivisionsLow(p1min, p1max, axisTitleSize)
    
    one = r.TLine(ratio.GetBinCenter(1)-ratio.GetBinWidth(1)/2.,1.,ratio.GetBinCenter(ratio.GetNbinsX())+ratio.GetBinWidth(ratio.GetNbinsX())/2.,1.);
    avg = r.TLine(ratio.GetBinCenter(1)-ratio.GetBinWidth(1)/2.,data_histo[0].Integral()/sum.Integral(),ratio.GetBinCenter(ratio.GetNbinsX())+ratio.GetBinWidth(ratio.GetNbinsX())/2.,data_histo[0].Integral()/sum.Integral());
    avg.SetLineColor(2);
    avg.SetLineStyle(2);
    one.SetLineStyle(2);
    #one.Draw();
    #avg.Draw();


    topPad.SetLogy()
    #can.SaveAs("Ngamma_RzG_Only_dRweight_for_Paper_.pdf")
    can.SaveAs("Ngamma_RzG_Only_dRweight_for_Twiki_.pdf")


    output_file.cd()
    for h in samples_histo :
        r.TH1D(h).Write()
    data_histo[0].Write()

output_file = r.TFile("baselineInputs_paper.root","RECREATE")

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
