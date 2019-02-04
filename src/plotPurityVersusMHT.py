from ROOT import *
from array import array

input_file = open('purity.log')

region=""
purity=0.0
purity_error=0.0

x=[275,325,475,800]
xErr=[25,25,125,200]
EB={"nom":[],"alt":[],"mcalt":[],"nomErr":[],"altErr":[],"mcaltErr":[]}
EE={"nom":[],"alt":[],"mcalt":[],"nomErr":[],"altErr":[],"mcaltErr":[]}

temp_purity=[]
temp_error=[]

for line in input_file:
    line = line[:-1]

    if "REGION" in line and region != line :
        print zip(temp_purity,temp_error)
        print "old:",region,"new:",line
        if len(temp_purity)==3:
            print "check"
            if "EB" in region :
                EB["alt"].append(temp_purity[0])
                EB["altErr"].append(temp_error[0])
                EB["nom"].append(temp_purity[1])
                EB["nomErr"].append(temp_error[1])
                EB["mcalt"].append(temp_purity[2])
                EB["mcaltErr"].append(temp_error[2])
            if "EE" in region :
                EE["alt"].append(temp_purity[0])
                EE["altErr"].append(temp_error[0])
                EE["nom"].append(temp_purity[1])
                EE["nomErr"].append(temp_error[1])
                EE["mcalt"].append(temp_purity[2])
                EE["mcaltErr"].append(temp_error[2])
            temp_purity=[]
            temp_error=[]
        region=line

    if "purity in SR:" in line :
        purity=float(line.split(":")[1].split("+/-")[0])
        purity_error=float(line.split(":")[1].split("+/-")[1])
        print region,purity,"+/-",purity_error
        temp_purity.append(purity)
        temp_error.append(purity_error)

if len(temp_purity)==3:
    print "check"
    print zip(temp_purity,temp_error)
    print region
    if "EB" in region :
        EB["alt"].append(temp_purity[0])
        EB["altErr"].append(temp_error[0])
        EB["nom"].append(temp_purity[1])
        EB["nomErr"].append(temp_error[1])
        EB["mcalt"].append(temp_purity[2])
        EB["mcaltErr"].append(temp_error[2])
    if "EE" in region :
        EE["alt"].append(temp_purity[0])
        EE["altErr"].append(temp_error[0])
        EE["nom"].append(temp_purity[1])
        EE["nomErr"].append(temp_error[1])
        EE["mcalt"].append(temp_purity[2])
        EE["mcaltErr"].append(temp_error[2])

for i in EB:
    print i,EB[i]
for i in EE:
    print i,EE[i]
EB["avg"] = []
EB["avgErr"] = []
EE["avg"] = []
EE["avgErr"] = []

for i in range(len(EE["nom"])):
    EB["avg"].append((EB["nom"][i]+EB["mcalt"][i]+EB["alt"][i])/3.)
    EE["avg"].append((EE["nom"][i]+EE["mcalt"][i]+EE["alt"][i])/3.)

for i in range(len(EE["nom"])):
    EB["avgErr"].append(0.)
    EE["avgErr"].append(0.)
    if( abs(EB["avg"][i] - EB["nom"][i]) > EB["avgErr"][i] ):  
        EB["avgErr"][i]=abs(EB["avg"][i] - EB["nom"][i])
    if( abs(EB["avg"][i] - EB["alt"][i]) > EB["avgErr"][i] ):  
        EB["avgErr"][i]=abs(EB["avg"][i] - EB["alt"][i])
    if( abs(EB["avg"][i] - EB["mcalt"][i]) > EB["avgErr"][i] ):  
        EB["avgErr"][i]=abs(EB["avg"][i] - EB["mcalt"][i])

    if( abs(EE["avg"][i] - EE["nom"][i]) > EE["avgErr"][i] ):  
        EE["avgErr"][i]=abs(EE["avg"][i] - EE["nom"][i])
    if( abs(EE["avg"][i] - EE["alt"][i]) > EE["avgErr"][i] ):  
        EE["avgErr"][i]=abs(EE["avg"][i] - EE["alt"][i])
    if( abs(EE["avg"][i] - EE["mcalt"][i]) > EE["avgErr"][i] ):  
        EE["avgErr"][i]=abs(EE["avg"][i] - EE["mcalt"][i])

print "x:",x
print "EB['nom']:",EB["nom"]
print "EB['nomErr']:",EB["nomErr"]
print "EB['mcalt']:",EB["mcalt"]
print "EB['mcaltErr']:",EB["mcaltErr"]
print "EB['alt']:",EB["alt"]
print "EB['altErr']:",EB["altErr"]

print "EE['nom']:",EE["nom"]
print "EE['nomErr']:",EE["nomErr"]
print "EE['mcalt']:",EE["mcalt"]
print "EE['mcaltErr']:",EE["mcaltErr"]
print "EE['alt']:",EE["alt"]
print "EE['altErr']:",EE["altErr"]

print "---------------------------------"
print "---------------------------------"
print "EB['avg']",EB["avg"]
print "EB['avgErr']",EB["avgErr"]
print "EE['avg']",EE["avg"]
print "EE['avgErr']",EE["avgErr"]

gROOT.ProcessLine(".L ~/tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")

can = TCanvas("can","can",500,500)
can.SetTopMargin(0.08)

avgEB = TGraphErrors(len(x),array("f",x),array("f",EB["avg"]),array("f",xErr),array("f",EB["avgErr"]))
avgEB.SetMarkerStyle(8)
avgEB.GetXaxis().SetTitle("H_{T}^{miss} [GeV]")
avgEB.GetYaxis().SetRangeUser(.5,1.)
avgEB.GetYaxis().SetTitle("Photon Purity")

avgEE = TGraphErrors(len(x),array("f",x),array("f",EE["avg"]),array("f",xErr),array("f",EE["avgErr"]))
avgEE.SetMarkerStyle(4)

# - - - - - - - -
nomEB = TGraphErrors(len(x),array("f",x),array("f",EB["nom"]),array("f",xErr),array("f",EB["nomErr"]))
nomEB.SetMarkerStyle(8)
nomEB.GetXaxis().SetTitle("H_{T}^{miss} [GeV]")
nomEB.GetYaxis().SetRangeUser(.5,1.)
nomEB.GetYaxis().SetTitle("Photon Purity")

nomEE = TGraphErrors(len(x),array("f",x),array("f",EE["nom"]),array("f",xErr),array("f",EE["nomErr"]))
nomEE.SetMarkerStyle(4)

# - - - - - - - - 

altEB = TGraphErrors(len(x),array("f",x),array("f",EB["alt"]),array("f",xErr),array("f",EB["altErr"]))
altEB.SetMarkerStyle(8)
altEB.SetMarkerColor(2)
altEB.SetLineColor(2)
altEB.GetXaxis().SetTitle("H_{T}^{miss} [GeV]")
altEB.GetYaxis().SetRangeUser(.5,1.)
altEB.GetYaxis().SetTitle("Photon Purity")

altEE = TGraphErrors(len(x),array("f",x),array("f",EE["alt"]),array("f",xErr),array("f",EE["altErr"]))
altEE.SetMarkerStyle(4)
altEE.SetMarkerColor(2)
altEE.SetLineColor(2)

# - - - - - - - - 

mcaltEB = TGraphErrors(len(x),array("f",x),array("f",EB["mcalt"]),array("f",xErr),array("f",EB["mcaltErr"]))
mcaltEB.SetMarkerStyle(8)
mcaltEB.SetMarkerColor(4)
mcaltEB.SetLineColor(4)
mcaltEB.GetXaxis().SetTitle("H_{T}^{miss} [GeV]")
mcaltEB.GetYaxis().SetRangeUser(.5,1.)
mcaltEB.GetYaxis().SetTitle("Photon Purity")

mcaltEE = TGraphErrors(len(x),array("f",x),array("f",EE["mcalt"]),array("f",xErr),array("f",EE["mcaltErr"]))
mcaltEE.SetMarkerStyle(4)
mcaltEE.SetMarkerColor(4)
mcaltEE.SetLineColor(4)

# - - - - - - - - 

leg = TLegend(.2,.2,.5,.5)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.AddEntry(nomEB,"EB","p")
leg.AddEntry(nomEE,"EE","p")
leg.AddEntry(altEB,"altEB","p")
leg.AddEntry(altEE,"altEE","p")
leg.AddEntry(mcaltEB,"mcaltEB","p")
leg.AddEntry(mcaltEE,"mcaltEE","p")

nomEB.Draw("Ap")
nomEE.Draw("p")
altEB.Draw("p")
altEE.Draw("p")
mcaltEB.Draw("p")
mcaltEE.Draw("p")
leg.Draw()

CMStext = TText(200,1.01,"CMS")
CMStext.SetTextFont(61)
CMStext.SetTextSize(0.045)
CMStext.Draw()

SIMtext = TText(285,1.01,"preliminary")
SIMtext.SetTextFont(52)
SIMtext.SetTextSize(0.045)
SIMtext.Draw()

LUMItext = TText(600,1.01,"13 TeV (41.5/fb)")
LUMItext.SetTextFont(51)
LUMItext.SetTextSize(0.045)
LUMItext.Draw()

can.SaveAs("../plots/purityResults/photonPurity_ALL_MHT.pdf")
can.SaveAs("../plots/purityResults/photonPurity_ALL_MHT.png")
can.SaveAs("../plots/purityResults/photonPurity_ALL_MHT.eps")

# - - - - - - - - 
canAvg = TCanvas("canAvg","canAvg",500,500)
canAvg.SetTopMargin(0.08)

legAvg = TLegend(.2,.2,.5,.5)
legAvg.SetBorderSize(0)
legAvg.SetFillColor(0)
legAvg.AddEntry(nomEB,"EB","p")
legAvg.AddEntry(nomEE,"EE","p")

avgEB.Draw("Ap")
avgEE.Draw("p")

CMStext.Draw()
SIMtext.Draw()
LUMItext.Draw()

canAvg.SaveAs("../plots/purityResults/photonPurity_MHT.pdf")
canAvg.SaveAs("../plots/purityResults/photonPurity_MHT.png")
canAvg.SaveAs("../plots/purityResults/photonPurity_MHT.eps")
