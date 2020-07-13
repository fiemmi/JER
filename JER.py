import ROOT
import math

# =================================================================================== #

def iterativeFit(responseHisto) :
    str_histoname = str(responseHisto.GetName())
    str_fitheader = "** Fitting  **"
    str_totheader = str_histoname + str_fitheader
    print '=' * len(str_totheader)
    print("** Fitting {0} **".format(responseHisto.GetName()))
    print '=' * len(str_totheader)
    resp = ROOT.RooRealVar("response", "response", -1.0, 2.0) #response variable: clusteredJetPt/genJetPt
    l = ROOT.RooArgList(resp)
    responseDistribution = ROOT.RooDataHist("pt0to100_distribution", "pt0to100_distribution", l, responseHisto)
    mean = ROOT.RooRealVar("mean","Mean of Gaussian",responseHisto.GetMean(1),-2.0,2.0)
    sigma = ROOT.RooRealVar("sigma","Width of Gaussian",responseHisto.GetStdDev(1),-2.0,2.0)
    gauss = ROOT.RooGaussian("gauss","gauss(x,mean,sigma)",resp,mean,sigma)
    #perform a first fit as a reference for the iterative procedure
    result = gauss.fitTo(responseDistribution,ROOT.RooFit.PrintLevel(-1))
    previousMean = mean.getValV()
    previousSigma = sigma.getValV()
    str_previousMean = str(previousMean)
    str_previousSigma = str(previousSigma)
    str_header = "** Starting values for the iterative fit procedure: mean = ; sigma =  **"
    str_tot = str_previousMean + str_previousSigma + str_header
    print '=' * len(str_tot)
    print("** Starting values for the iterative fit procedure: mean = {0}; sigma = {1} **".format(previousMean, previousSigma))
    print '=' * len(str_tot)
    for idx, i in enumerate(range(20)) :
        mean = ROOT.RooRealVar("mean","Mean of Gaussian",previousMean,-2.0,2.0)
        sigma = ROOT.RooRealVar("sigma","Width of Gaussian",previousSigma,-2.0,2.0)
        gauss = ROOT.RooGaussian("gauss","gauss(x,mean,sigma)",resp,mean,sigma)
        result = gauss.fitTo(responseDistribution,ROOT.RooFit.PrintLevel(-10),ROOT.RooFit.Range(previousMean - 2.0*previousSigma, previousMean + 2.0*previousSigma))
        if abs(1 - sigma.getValV()/previousSigma) < 0.01 :
            previousMean = mean.getValV()
            previousSigma = sigma.getValV()
            print ("ITERATIVE FITTING STOPPED AT ITERATION {0}".format(idx))
            break
        else :
            previousMean = mean.getValV()
            previousSigma = sigma.getValV()
        print ("Iteration number {0}".format(idx))
        print("mean = {0}".format(previousMean))
        print("sigma = {0}".format(previousSigma))
        print("~~~~~ o ~~~~~~")
    #mean.Print()
    #sigma.Print()
    print '=' * 32
    print '=' * 30
    print("=======  FINAL RESULT ======")
    print("mean = {0} ".format(previousMean))
    print("sigma = {0}".format(previousSigma))
    print '=' * 30
    print '=' * 32
    respframe = resp.frame()
    responseDistribution.plotOn(respframe,ROOT.RooLinkedList())
    gauss.plotOn(respframe)
    can1 = ROOT.TCanvas("can1","Signal maximum likelihood fit", 900, 600);
    respframe.Draw()
    can1.Print("{0}.png".format(responseHisto.GetName()))
    
def matching(event) :
    deltaRmin = 10.0
    whichClusteredJet = 100
    for idx, CJet in enumerate(range(event.nClusteredJets)) :
        deltaEta = abs(event.ClusteredJetEta.at(CJet) - event.GenJetEta.at(GJet))
        deltaPhi = abs(event.ClusteredJetPhi.at(CJet) - event.GenJetPhi.at(GJet))
        if deltaPhi > 3.14159265358979323846 :
            deltaPhi = 2*3.14159265358979323846 - deltaPhi
        deltaR = math.sqrt(deltaEta**2 + deltaPhi**2)
        if deltaR < deltaRmin :
            deltaRmin = deltaR
            whichClusteredJet = idx
    return [deltaRmin, whichClusteredJet]
    
# =================================================================================== #

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.PROGRESS) #make fitting procedure silent up to progress messages; warning and errors are still shown

inputfile  = ROOT.TFile.Open("cluster_output.root")

ptResponseHisto_pt0to100 = ROOT.TH1F ("ptResponseHisto_pt0to100", "ptResponseHisto_pt0to100", 80, -2, 2)
ptResponseHisto_pt100to200 = ROOT.TH1F ("ptResponseHisto_pt100to200", "ptResponseHisto_pt100to200", 80, -2, 2)
ptResponseHisto_pt200to300 = ROOT.TH1F ("ptResponseHisto_pt200to300", "ptResponseHisto_pt200to300", 80, -2, 2)
ptResponseHisto_nPU0to10 = ROOT.TH1F ("ptResponseHisto_nPU0to10", "ptResponseHisto_nPU0to10", 80, -2, 2)
ptResponseHisto_nPU10to20 = ROOT.TH1F ("ptResponseHisto_nPU10to20", "ptResponseHisto_nPU10to20", 80, -2, 2)
ptResponseHisto_nPU20to30 = ROOT.TH1F ("ptResponseHisto_nPU20to30", "ptResponseHisto_nPU20to30", 80, -2, 2)
ptResponseHisto_nPU30to40 = ROOT.TH1F ("ptResponseHisto_nPU30to40", "ptResponseHisto_nPU30to40", 80, -2, 2)
etaResponseHisto_pt0to100 = ROOT.TH1F ("etaResponseHisto_pt0to100", "etaResponseHisto_pt0to100", 80, -2, 2)
etaResponseHisto_pt100to200 = ROOT.TH1F ("etaResponseHisto_pt100to200", "etaResponseHisto_pt100to200", 80, -2, 2)
etaResponseHisto_pt200to300 = ROOT.TH1F ("etaResponseHisto_pt200to300", "etaResponseHisto_pt200to300", 80, -2, 2)


#loop over events in inputfile
for event in inputfile.jets :
#fill histograms for pt and eta distributions
    for GJet in range (event.nGenJets) :
        #JER and eta resolution, 0 < pt < 100
        if event.GenJetPt.at(GJet) > 0 and event.GenJetPt.at(GJet) < 100 :
            matching_result = matching(event)
            if matching_result[0] < 0.4 :
                ptResponseHisto_pt0to100.Fill(event.ClusteredJetPt.at(matching_result[1])/event.GenJetPt.at(GJet))
                etaResponseHisto_pt0to100.Fill(event.ClusteredJetEta.at(matching_result[1])/event.GenJetEta.at(GJet))
        
        #JER and eta resolution, 100 < pt < 200
        if event.GenJetPt.at(GJet) > 100 and event.GenJetPt.at(GJet) < 200 :
            matching_result = matching(event)
            if matching_result[0] < 0.4 :
                ptResponseHisto_pt100to200.Fill(event.ClusteredJetPt.at(matching_result[1])/event.GenJetPt.at(GJet))
                etaResponseHisto_pt100to200.Fill(event.ClusteredJetEta.at(matching_result[1])/event.GenJetEta.at(GJet))
        
        #JER and eta resolution, 200 < pt < 300        
        if event.GenJetPt.at(GJet) > 200 and event.GenJetPt.at(GJet) < 300 :
            matching_result = matching(event)
            if matching_result[0] < 0.4 :
                ptResponseHisto_pt200to300.Fill(event.ClusteredJetPt.at(matching_result[1])/event.GenJetPt.at(GJet))
                etaResponseHisto_pt200to300.Fill(event.ClusteredJetEta.at(matching_result[1])/event.GenJetEta.at(GJet))
        
        #JER resolution, 0 < nPUint < 10
        if event.nPUint > 0 and event.nPUint < 10 :
            matching_result = matching(event)
            if matching_result[0] < 0.4 :
                ptResponseHisto_nPU0to10.Fill(event.ClusteredJetPt.at(matching_result[1])/event.GenJetPt.at(GJet))
        
        #JER resolution, 10 < nPUint < 20
        if event.nPUint > 10 and event.nPUint < 20 :
            matching_result = matching(event)
            if matching_result[0] < 0.4 :
                ptResponseHisto_nPU10to20.Fill(event.ClusteredJetPt.at(matching_result[1])/event.GenJetPt.at(GJet))
                
        #JER resolution, 20 < nPUint < 30
        if event.nPUint > 20 and event.nPUint < 30 :
            matching_result = matching(event)
            if matching_result[0] < 0.4 :
                ptResponseHisto_nPU20to30.Fill(event.ClusteredJetPt.at(matching_result[1])/event.GenJetPt.at(GJet))
                
        #JER resolution, 20 < nPUint < 30
        if event.nPUint > 30 and event.nPUint < 40 :
            matching_result = matching(event)
            if matching_result[0] < 0.4 :
                ptResponseHisto_nPU30to40.Fill(event.ClusteredJetPt.at(matching_result[1])/event.GenJetPt.at(GJet))
                
outputfile = ROOT.TFile( 'JER_output.root', 'RECREATE' )
ptResponseHisto_pt0to100.Write()
ptResponseHisto_pt100to200.Write()
ptResponseHisto_pt200to300.Write()
ptResponseHisto_nPU0to10.Write()
ptResponseHisto_nPU10to20.Write()
ptResponseHisto_nPU20to30.Write()
ptResponseHisto_nPU30to40.Write()
etaResponseHisto_pt0to100.Write()
etaResponseHisto_pt100to200.Write()
etaResponseHisto_pt200to300.Write()
outputfile.Close()

#iteratively fit the histograms to obtain JER
inputfile_bis = ROOT.TFile.Open("JER_output.root")
hpt_pt0to100 = inputfile_bis.Get("ptResponseHisto_pt0to100")
iterativeFit(hpt_pt0to100)
hpt_pt100to200 = inputfile_bis.Get("ptResponseHisto_pt100to200")
iterativeFit(hpt_pt100to200)
hpt_pt200to300 = inputfile_bis.Get("ptResponseHisto_pt200to300")
iterativeFit(hpt_pt200to300)
hpt_nPU0to10 = inputfile_bis.Get("ptResponseHisto_nPU0to10")
iterativeFit(hpt_nPU0to10)
hpt_nPU10to20 = inputfile_bis.Get("ptResponseHisto_nPU10to20")
iterativeFit(hpt_nPU10to20)
hpt_nPU20to30 = inputfile_bis.Get("ptResponseHisto_nPU20to30")
iterativeFit(hpt_nPU20to30)
hpt_nPU30to40 = inputfile_bis.Get("ptResponseHisto_nPU30to40")
iterativeFit(hpt_nPU30to40)
heta_pt0to100 = inputfile_bis.Get("etaResponseHisto_pt0to100")
iterativeFit(heta_pt0to100)
heta_pt100to200 = inputfile_bis.Get("etaResponseHisto_pt100to200")
iterativeFit(heta_pt100to200)
heta_pt200to300 = inputfile_bis.Get("etaResponseHisto_pt200to300")
iterativeFit(heta_pt200to300)
inputfile_bis.Close()
