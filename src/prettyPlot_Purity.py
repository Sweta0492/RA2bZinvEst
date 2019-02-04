from sys import argv
import ROOT as r
r.gROOT.SetBatch(True)

r.gROOT.ProcessLine(".L ~/tdrstyle.C")
r.gROOT.ProcessLine("setTDRStyle()")

plot_dir="plotPurityProperties_plots"

input_file_name = "plotPurityProperties_photonLoose.root"
#input_file_name = "plotObs_photon_baseline.root"
#input_file_name = "plotObs_photonLoose_baseline.root"
input_file = r.TFile(input_file_name,"READ")    

def plot(plot_var = "photonIsoChrgLowSieie_EB_photonLoose" ):
    
    samples=[["QCD_1000to1500",
              "QCD_1500to2000",
              "QCD_2000toInf",
              "QCD_200to300",
              "QCD_300to500",
              "QCD_500to700",
              "QCD_700to1000"],
             ["GJets_100to200",
              "GJets_100to200",
              "GJets_200to400",
              "GJets_400to600",
              "GJets_600toInf"]]
    
    data_samples=["SinglePhoton_2018A",
                  "SinglePhoton_2018B",
                  "SinglePhoton_2018C",
                  "SinglePhoton_2018D"]

    samples_labels = ["QCD","GJets"]
    samples_fill_color = [r.kGray,r.kGreen]
    samples_line_color = [1,1]
    
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
    topPad.SetBottomMargin(0.01);
    topPad.Draw();
    botPad.Draw();
    topPad.cd();
    
    stack.Draw("histo")    
    data_histo[0].Draw("SAME,p")

    stack.SetMaximum(1.3*max(sum.GetMaximum(),data_histo[0].GetMaximum()))
    stack.SetMinimum(0.1)

    stack.GetYaxis().SetTitle("Events")
    stack.GetYaxis().SetLabelFont(43);
    stack.GetYaxis().SetLabelSize(24);
    stack.GetYaxis().SetTitleFont(43);
    stack.GetYaxis().SetTitleSize(24);
    stack.GetYaxis().SetTitleOffset(1.5);

    stack.GetXaxis().SetLabelFont(43);
    stack.GetXaxis().SetLabelSize(24);
    stack.GetXaxis().SetTitleFont(43);
    stack.GetXaxis().SetTitleSize(24);
    stack.GetXaxis().SetTitleOffset(1.7);

    CMStext = r.TText(.17,.95,"CMS")
    CMStext.SetNDC()
    CMStext.SetTextFont(61)
    CMStext.SetTextSize(0.08)
    CMStext.Draw()
    
    SIMtext = r.TText(.28,.95,"preliminary")
    SIMtext.SetNDC()
    SIMtext.SetTextFont(52)
    SIMtext.SetTextSize(0.08)
    SIMtext.Draw()
    
    LUMItext = r.TText(.65,.95,"13 TeV (41.X/fb)")
    LUMItext.SetNDC()
    LUMItext.SetTextFont(51)
    LUMItext.SetTextSize(0.08)
    LUMItext.Draw()

    botPad.cd()
    ratio = r.TH1D(data_histo[0])
    ratio.SetMarkerStyle(8)
    ratio.SetName(plot_var+"_ratio")
    ratio.Divide(sum)
    ratio.GetYaxis().SetRangeUser(0,2)
    ratio.GetYaxis().SetTitle("Obs/Exp")
    ratio.GetXaxis().SetTitle(data_histo[0].GetTitle())
    ratio.GetYaxis().SetLabelFont(43);
    ratio.GetYaxis().SetLabelSize(24);
    ratio.GetYaxis().SetTitleFont(43);
    ratio.GetYaxis().SetTitleSize(24);
    ratio.GetYaxis().SetNdivisions(505);
    ratio.GetYaxis().SetTitleOffset(1.5);

    ratio.GetXaxis().SetLabelFont(43);
    ratio.GetXaxis().SetLabelSize(24);
    ratio.GetXaxis().SetTitleFont(43);
    ratio.GetXaxis().SetTitleSize(24);
    ratio.GetXaxis().SetTitleOffset(2.8);

    ratio.Draw()

    can.SaveAs("../plots/"+plot_dir+"/"+plot_var+".png")
    topPad.SetLogy()
    can.SaveAs("../plots/"+plot_dir+"/"+plot_var+"_LogY.png")

    output_file.cd()
    for h in samples_histo :
        r.TH1D(h).Write()
    data_histo[0].Write()

#output_file = r.TFile("baselineInputs.root","RECREATE")
output_file = r.TFile("purityInputs.root","RECREATE")

vars = []
list = input_file.GetListOfKeys()
next = r.TIter(list);
key = next()
while(key != None ) :
    obj = key.ReadObj();
    if obj.InheritsFrom("TH1") :
        name = r.TString(obj.GetName())
        if name.Contains("_QCD_300to500") : 
            vars.append(name.ReplaceAll("_QCD_300to500","").Data())
    else :
        print obj.Print()
    key = next()

print vars
for var in vars : 
    plot(var)

    
output_file.Close()
input_file.Close()
