from multiprocessing import Process
import os

backgroundSamples=["QCD_200to300",
                   "QCD_300to500",
                   "QCD_500to700",
                   "QCD_700to1000",
                   "QCD_1000to1500",
                   "QCD_1500to2000",
                   "QCD_2000toInf",
                   "GJets_100to200",
                   "GJets_200to400",
                   "GJets_400to600",
                   "GJets_600toInf"
                    ]

dataSamples=["SinglePhoton_2018A",
             "SinglePhoton_2018B",
             "SinglePhoton_2018C",
             "SinglePhoton_2018D"]

def runPlotPurityProperties(bkg,data):
    print '../bin/plotPurityProperties 3 "{0}" "{1}"'.format(bkg,data)
    os.system('../bin/plotPurityProperties 3 "{0}" "{1}"'.format(bkg,data))

processes=[]
for sample in backgroundSamples : 
    p = Process(target=runPlotPurityProperties, args=(sample, "") )
    p.start()
    processes.append(p)

for sample in dataSamples : 
    p = Process(target=runPlotPurityProperties, args=("", sample) )
    p.start()
    processes.append(p)

for p in processes : 
    p.join()

#os.system("hadd -f plotPurityProperties_photonLoose.root plotPurityProperties_photonLoose_*.root")
#os.system("rm plotPurityProperties_photonLoose_*.root")
    


