#! /usr/bin/env python
import ROOT
import optparse
import json
import sys
import os
from TopLJets2015.TopAnalysis.storeTools import getEOSlslist
import MT2Calculator

"""
Take the ratio of two Breit-Wigner functions at fixed mass as a reweighting factor
"""
def weightTopWidth(tmassList,bwigner,targetWidth,origWidth=1.324):
    bwigner.FixParameter(2,origWidth)
    origNorm=bwigner.Integral(0,300)

    bwigner.FixParameter(2,targetWidth)
    targetNorm=bwigner.Integral(0,300)

    wgt=1.0
    for m in tmassList:
        bwigner.FixParameter(2,origWidth)
        origVal=bwigner.Eval(m)
        bwigner.FixParameter(2,targetWidth)
        targetVal=bwigner.Eval(m)
        wgt *= (targetVal/targetNorm) / (origVal/origNorm)
    return wgt

"""
Analysis loop
"""
def runTopWidthAnalysis(fileName,
                        outFileName,
                        widthList=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0],
                        systs=['','puup','pudn','btagup','btagdn','ltagup','ltagdn','jerup','jerdn','jesup','jesdn','lesup','lesdn','trigup','trigdn','selup','seldn','topptup','topptdn']):

    print '....analysing',fileName,'with output @',outFileName

    #check if this is data beforehand
    isData=False if 'MC13TeV' in fileName else True
    if isData:
        widthList=[1.0]
        systs=['']
    else:
        #generator level systematics for ttbar
        if 'TTJets' in fileName and 'powheg' in fileName:
            for i in xrange(1,120):  systs += ['gen%d'%i]

    smWidth=1.324 # all powheg samples have been generated with the width @ 172.5 GeV
    smMass=172.5
    if '166v5' in fileName : smMass=166.5
    if '169v5' in fileName : smMass=169.5
    if '175v5' in fileName : smMass=175.5
    if '178v5' in fileName : smMass=178.5

    #define the relativistic Breit-Wigner function
    bwigner=ROOT.TF1('bwigner',
                     '[0]*([1]*[2]*sqrt([1]*[1]*([1]*[1]+[2]*[2]))/sqrt([1]*[1]+sqrt([1]*[1]*([1]*[1]+[2]*[2]))))/(TMath::Power(x*x-[1]*[1],2)+TMath::Power([1]*[2],2))',
                     0,300)
    bwigner.SetParName(0,"N")
    bwigner.SetParameter(0,1.0)
    bwigner.SetParName(1,"m_{t}")
    bwigner.FixParameter(1,smMass)
    bwigner.SetParName(2,"#Gamma_{t}")
    bwigner.FixParameter(2,smWidth)

    #book histograms
    observablesH={}

    #MC truth control histograms
    for w in widthList:
        var='tmass_%3.1fw'%w
        observablesH[var]=ROOT.TH1F(var,';Top quark mass [GeV];Top quarks',100,150,200)
        for assig in ['cor','wro']:
            var='%sgenmlbvsmtop_%3.1fw'%(assig,w)
            observablesH[var]=ROOT.TH2F(var,';Mass(lepton,jet) [GeV];Top mass;l+j pairs',30,0,300,100,150,200)

    #RECO level histograms
    for j in ['EE','MM','EM']:
        var=j+'_evcount'
        observablesH[var]=ROOT.TH1F(var,';Category;Events',2,0,2)
        observablesH[var].GetXaxis().SetBinLabel(1,'=1b')
        observablesH[var].GetXaxis().SetBinLabel(2,'#geq2b')
        for b in ['1b','2b']:

            var=j+b+'_ptlb'
            observablesH[var]=ROOT.TH1F(var,';p_{T}(lepton,jet) [GeV];l+j pairs',30,0,300)
            var=j+b+'_met'
            observablesH[var]=ROOT.TH1F(var,';Missing transverse energy [GeV];Events',30,0,300)
            var=j+b+'_mll'
            observablesH[var]=ROOT.TH1F(var,';Dilepton invariant mass [GeV];Events',30,0,300)
            var=j+b+'_njets'
            observablesH[var]=ROOT.TH1F(var,';Jet multiplicity;Events',6,2,8)
            var=j+b+'_njetsnonboost'
            observablesH[var]=ROOT.TH1F(var,';Jet multiplicity;Events',6,2,8)
            var=j+b+'_njetsboost'
            observablesH[var]=ROOT.TH1F(var,';Jet multiplicity;Events',6,2,8)

            for s in systs:
                for i in ['lowpt','highpt']:
                    for w in widthList:
                        var=s+i+j+b+'_mlb_%3.1fw'%w
                        observablesH[var]=ROOT.TH1F(var,';Mass(lepton,jet) [GeV];l+j pairs',30,0,300)
                        var=s+i+j+b+'_incmlb_%3.1fw'%w
                        observablesH[var]=ROOT.TH1F(var,';Mass(lepton,jet) (Inclusive) [GeV];l+j pairs',30,0,300)
                        #var=s+i+j+b+'_sncmlb_%3.1fw'%w
                        #observablesH[var]=ROOT.TH1F(var,';Mass(lepton,jet) (Semi-Inclusive) [GeV];l+j pairs',30,0,300)
                        #var=s+i+j+b+'_mdrmlb_%3.1fw'%w
                        #observablesH[var]=ROOT.TH1F(var,';Mass(lepton,jet) (Minimum #Delta R) [GeV];l+j pairs',30,0,300)
                        #var=s+i+j+b+'_minmlb_%3.1fw'%w
                        #observablesH[var]=ROOT.TH1F(var,';Mass(lepton,jet) (Minimum) [GeV];l+j pairs',30,0,300)
                        var=s+i+j+b+'_mt2mlb_%3.1fw'%w
                        observablesH[var]=ROOT.TH1F(var,';Mass(lepton,jet) (M_{T2} Method) [GeV];l+j pairs',30,0,300)

                        if w!=1.0 or len(s)>0 : continue
                        var=i+j+b+'_pairing'
                        observablesH[var]=ROOT.TH1F(var,';Pairing;l+j pairs',2,0,2)
                        observablesH[var].GetXaxis().SetBinLabel(1,'correct')
                        observablesH[var].GetXaxis().SetBinLabel(2,'wrong')

    for var in observablesH:
        observablesH[var].SetDirectory(0)
        observablesH[var].Sumw2()

    #open file
    puNormSF=1.0
    if isData:
        fIn=ROOT.TFile.Open(fileName)
        puCorrH=fIn.Get('puwgtctr')
        nonWgt=puCorrH.GetBinContent(1)
        wgt=puCorrH.GetBinContent(2)
        if wgt>0 : puNormSF=nonWgt/wgt
        fIn.Close()

    tree=ROOT.TChain('twev')
    tree.AddFile(fileName)

    #loop over events in the tree and fill histos
    totalEntries=tree.GetEntries()
    for i in xrange(0,totalEntries):

        tree.GetEntry(i)

        if i%100==0 : sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )

        #focus on dilepton analysis for the moment
        if abs(tree.cat)==11 or abs(tree.cat)==13 : continue
        if abs(tree.cat)==11*11 : evcat='EE'
        if abs(tree.cat)==11*13 : evcat='EM'
        if abs(tree.cat)==13*13 : evcat='MM'

        #base event weight
        baseEvWeight=puNormSF*tree.weight[0]

        #determine weighting factors for the width
        tops={}
        tmassList=[]
        for it in xrange(0,tree.nt):
            if it>1 : break
            tid=tree.t_id[it]
            tops[ tid ] = ROOT.TLorentzVector()
            tops[ tid ].SetPtEtaPhiM(tree.t_pt[it],tree.t_eta[it],tree.t_phi[it],tree.t_m[it])
            tmassList.append( tops[tid].M() )
        widthWeight={}
        for w in widthList:
            # eac edit
            widthCorrFactor= 1 if (isData or "MC13TeV_TTJets_dilpowheg" not in filename) else 4
            widthWeight[w]=weightTopWidth(tmassList,bwigner,w/widthCorrFactor*smWidth,smWidth)

            #paranoid check
            var='tmass_%3.1fw'%w
            for mtop in tmassList:
                observablesH[var].Fill(mtop,baseEvWeight*widthWeight[w])


        #preselect the b-jets (central b-tag, b-tag up, b-tag dn, l-tag up, l-tag dn, jer up, jer dn, jes up, jes dn)
        bjets=( [], [], [], [], [], [], [], [], [] )
        otherjets=( [], [], [], [], [], [], [], [], [] )

        #define b-tag variations (bit to use, do heavy flavour shift with None=use nominal)
        btagVars=[ (0,None), (1,True), (2,True), (1,False), (2,False) ]
        for ij in xrange(0,tree.nj):

            jp4=ROOT.TLorentzVector()
            jp4.SetPtEtaPhiM(tree.j_pt[ij],tree.j_eta[ij],tree.j_phi[ij],tree.j_m[ij])

            for ibvar in xrange(0,len(btagVars)):

                #reset b-tag bit to 0 if flavour is not the one to vary
                ibit, shiftHeavyFlav = btagVars[ibvar]
                if shiftHeavyFlav is not None:
                    hadFlav=abs(tree.gj_hadflav[ij])
                    if shiftHeavyFlav :
                        if hadFlav!=4 and hadFlav!=5 : ibit=0
                    else:
                        if hadFlav==4 or hadFlav==5 : ibit=0

                #get b-tag decision
                btagVal=((tree.j_btag[ij] >> ibit) & 0x1)

                if btagVal > 0:
                    bjets[ibvar].append( (ij,jp4) )

                    #use standard b-tag decision for JES/JER variations
                    if shiftHeavyFlav is not None : continue

                    jres=ROOT.TMath.Abs(1-tree.j_jer[ij])
                    bjets[5].append( (ij,jp4*(1+jres)) )
                    bjets[6].append( (ij,jp4*(1-jres)) )

                    jscale=tree.j_jes[ij]
                    bjets[7].append( (ij,jp4*(1+jscale)) )
                    bjets[8].append( (ij,jp4*(1-jscale)) )
                elif btagVal == 0 :
                    otherjets[ibvar].append( (ij,jp4) )

                    #use standard b-tag decision for JES/JER variations
                    if shiftHeavyFlav is not None : continue

                    jres=ROOT.TMath.Abs(1-tree.j_jer[ij])
                    otherjets[5].append( (ij,jp4*(1+jres)) )
                    otherjets[6].append( (ij,jp4*(1-jres)) )

                    jscale=tree.j_jes[ij]
                    otherjets[7].append( (ij,jp4*(1+jscale)) )
                    otherjets[8].append( (ij,jp4*(1-jscale)) )

        #build the dilepton
        dilepton=ROOT.TLorentzVector()
        for il in xrange(0,2):
            stdlp4=ROOT.TLorentzVector()
            stdlp4.SetPtEtaPhiM(tree.l_pt[il],tree.l_eta[il],tree.l_phi[il],tree.l_m[il])
            dilepton+=stdlp4

        #global control histos
        nbtags=len(bjets[0])
        if nbtags>2 : nbtags=2
        btagcat='1b' if nbtags==1 else '2b'
        if nbtags==0 : btagcat='0b'

        if nbtags>0:
            var=evcat+btagcat+'_mll'
            observablesH[var].Fill(dilepton.M(),baseEvWeight)

        #remove Z/quarkonia candidates
        if abs(tree.cat)==11*11 or abs(tree.cat)==13*13:
            if ROOT.TMath.Abs(dilepton.M()-91)<15 : continue
            if dilepton.M()<20: continue

        if nbtags>0:
            var=evcat+"_evcount"
            observablesH[var].Fill(nbtags-1,baseEvWeight)
            var=evcat+btagcat+'_met'
            observablesH[var].Fill(tree.met_pt,baseEvWeight)
            var=evcat+btagcat+'_njets'
            observablesH[var].Fill(tree.nj,baseEvWeight)

            #check if event has a boosted pair
            hasBoostedPair=False
            for il in xrange(0,2):
                lp4=ROOT.TLorentzVector()
                lp4.SetPtEtaPhiM(tree.l_pt[il],tree.l_eta[il],tree.l_phi[il],tree.l_m[il])
                for ib in xrange(0,nbtags):
                    ij,jp4 = bjets[0][ib]
                    ptlb=(lp4+jp4).Pt()
                    if ptlb>100 : hasBoostedPair=True
            var=evcat+btagcat+'_njetsboost' if hasBoostedPair else evcat+btagcat+'_njetsnonboost'
            observablesH[var].Fill(tree.nj,baseEvWeight)

        # setup mlb calculations
        mlbTypes  = ["inc","mt2"]#,"snc","min","mdr"]
        ptCatList = ["highpt", "lowpt"]
        bTagCats  = ["1b", "2b"]
        evCatList = ["EE", "EM", "MM"]
        mlbMap = {}
        for s,ptC,evC,btC,mC,w in [(s,a,b,c,d,e) for s in systs
                for a in ptCatList
                for b in evCatList
                for c in bTagCats
                for d in mlbTypes
                for e in widthList]:
            initVals = [(float('inf'),0)]
            if mC == "inc" or mC == "snc" : initVals = []
            mlbMap["%s%s%s%s_%smlb_%3.1fw"%(s,ptC,evC,btC,mC,w)] = initVals

        metForMT2 = ROOT.TLorentzVector()
        metForMT2.SetPtEtaPhiM(tree.met_pt, 0, tree.met_phi, 0)


        #pair with the leptons
        for il in xrange(0,2):

            stdlp4=ROOT.TLorentzVector()
            stdlp4.SetPtEtaPhiM(tree.l_pt[il],tree.l_eta[il],tree.l_phi[il],tree.l_m[il])
            lscale=tree.l_les[il]

            for s in systs:

                #event weight
                evWeight=baseEvWeight

                #base lepton kinematics
                lp4=ROOT.TLorentzVector(stdlp4)

                #experimental uncertainties
                ijhyp=0
                if s=='btagup' : ijhyp=1
                if s=='btagdn' : ijhyp=2
                if s=='ltagup' : ijhyp=3
                if s=='ltagdn' : ijhyp=4
                if s=='jerup'  : ijhyp=5
                if s=='jerdn'  : ijhyp=6
                if s=='jesup'  : ijhyp=7
                if s=='jesdn'  : ijhyp=8
                if s=='lesup'  : lp4 *= (1.0+lscale)
                if s=='lesdn'  : lp4 *= (1.0-lscale)
                if s=='puup'   : evWeight=puNormSF*tree.weight[1]
                if s=='pudn'   : evWeight=puNormSF*tree.weight[2]
                if s=='trigup':  evWeight=puNormSF*tree.weight[3]
                if s=='trigdn':  evWeight=puNormSF*tree.weight[4]
                if s=='selup':  evWeight=puNormSF*tree.weight[5]
                if s=='seldn':  evWeight=puNormSF*tree.weight[6]
                if s=='topptup':  evWeight=puNormSF*tree.weight[7]
                if s=='topptdn':  evWeight=puNormSF*tree.weight[8]
                if 'gen' in s : evWeight=puNormSF*tree.weight[8+int(s[3:])]

                #btag hypothesis
                nbtags=len(bjets[ijhyp])
                if nbtags<1 : continue
                if nbtags>2 : nbtags=2
                btagcat='1b' if nbtags==1 else '2b'
                for ib in xrange(0,nbtags):

                    ij,jp4 = bjets[ijhyp][ib]

                    #MC truth for this pair
                    pairFullyMatchedAtGen = (tree.gl_id[il]!=0 and abs(tree.gj_flav[ij])==5)
                    assignmentType,tmass,genmlb=1,0.0,0.0
                    if pairFullyMatchedAtGen and tree.nt>0:

                        #MC truth  kinematkcs
                        glp4=ROOT.TLorentzVector()
                        glp4.SetPtEtaPhiM(tree.gl_pt[il],tree.gl_eta[il],tree.gl_phi[il],tree.gl_m[il])
                        gjp4=ROOT.TLorentzVector()
                        gjp4.SetPtEtaPhiM(tree.gj_pt[ij],tree.gj_eta[ij],tree.gj_phi[ij],tree.gj_m[ij])
                        genmlb=(glp4+gjp4).M()

                        #correctness of the assignment can be checked by op. charge
                        if tree.gl_id[il]*tree.gj_flav[ij]<0 : assignmentType=0

                        #top mass (parton level)
                        try:
                            if tree.gl_id[il]<0 : tmass=tops[6].M()
                            else                : tmass=tops[-6].M()
                        except:
                            pass

                    #save MC truth distribution
                    if s=='' and pairFullyMatchedAtGen:
                        for w in widthList:
                            var='cor' if assignmentType==0 else 'wro'
                            var+='genmlbvsmtop_%3.1fw'%w
                            observablesH[var].Fill(genmlb,tmass,baseEvWeight*widthWeight[w])

                    #RECO kinematics of the l,b system
                    mlb=(lp4+jp4).M()
                    ptlb=(lp4+jp4).Pt()
                    ptCat='lowpt' if ptlb<100 else 'highpt'
                    dRlb=lp4.DeltaR(jp4)

                    # calculate mt2 only in the last loop
                    mt2=float('inf')
                    mt2ptCat=''
                    if il == 1 and ib == nbtags-1:
                        # get other lepton
                        tlp=ROOT.TLorentzVector()
                        tlp.SetPtEtaPhiM(tree.l_pt[0],tree.l_eta[0],tree.l_phi[0],tree.l_m[0])
                        if s=="lesup" :
                            tlp*=(1+tree.l_les[0])
                        if s=="lesdown" :
                            tlp*=(1-tree.l_les[0])

                        # get other b-jet (or next highest csv)
                        tb=None
                        if nbtags==1 :
                            _,tb = otherjets[ijhyp][0]
                        elif nbtags==2 :
                            _,tb = bjets[ijhyp][0]

                        # calculate mt2 for the highest HT
                        tlpjp4IsHighestHT = ((tlp+jp4).Pt() + (lp4+tb).Pt()) >= ((lp4+jp4).Pt() + (tlp+tb).Pt())

                        try :
                            if tlpjp4IsHighestHT :
                                mt2 = MT2Calculator.calcMt2(tlp+jp4,lp4+tb,metForMT2)
                                mt2ptCat= 'highpt' if ((tlp+jp4).Pt() > 100 or (lp4+tb).Pt() > 100) else 'lowpt'
                            else :
                                mt2 = MT2Calculator.calcMt2(tlp+tb,lp4+jp4,metForMT2)
                                mt2ptCat= 'highpt' if ((lp4+jp4).Pt() > 100 or (tlp+tb).Pt() > 100) else 'lowpt'
                        except ZeroDivisionError:
                            mt2=float('inf')
                            print '\t\t - event has incalculable mt2, skipping... \n'



                    #fill histos
                    if s=='':
                        var=evcat+btagcat+'_ptlb'
                        observablesH[var].Fill(ptlb,evWeight)

                    # keep the old way of storing mlb
                    for w in widthList:
                        mlbWt = evWeight * widthWeight[w]

                        var=s+ptCat+evcat+btagcat+'_mlb_%3.1fw'%w
                        observablesH[var].Fill(mlb,mlbWt)

                        # fill inclusive mlb
                        var=s+ptCat+evcat+btagcat+'_incmlb_%3.1fw'%w
                        mlbMap[var] += [(mlb,mlbWt)]

                        # fill semi-inclusive mlb (to sort later)
                        #var=s+ptCat+evcat+btagcat+'_sncmlb_%3.1fw'%w
                        #mlbMap[var] += [(mlb,mlbWt)]+[(dRlb,0)]

                        # fill min mlb
                        #var=s+ptCat+evcat+btagcat+'_minmlb_%3.1fw'%w
                        #if mlb < mlbMap[var][0] :
                        #    mlbMap[var] = [(mlb,mlbWt)]

                        # fill mdr mlb
                        #var=s+ptCat+evcat+btagcat+'_mdrmlb_%3.1fw'%w
                        #if mlb < mlbMap[var][0] and dRlb<=1 :
                        #    mlbMap[var] = [(mlb,mlbWt)]

                        # fill mt2 mlb
                        if il == 1 and ib == nbtags-1 and (mt2ptCat=='highpt' or mt2ptCat=='lowpt'):
                            var=s+mt2ptCat+evcat+btagcat+'_mt2mlb_%3.1fw'%w
                            if mt2 < mlbMap[var][0] :
                                mlbMap[var] = [(mt2,mlbWt)]

                        #only for standard width and syst variations
                        if w!=1.0 or len(s)>0 : continue
                        var=s+ptCat+evcat+btagcat+'_pairing'
                        observablesH[var].Fill(assignmentType,evWeight*widthWeight[w])

        # fill all relevant histos
        for histoName in mlbMap :
            # make sure there's one entry for appropriate mlb
            #if len(mlbMap[histoName]) > 1:
            #    if "min" in histoName or "mdr" in histoName or "mt2" in histoName :
            #        print "WARNING: storing more than one mlb for %s"%histoName
            # perform snc reduction
            if "snc" in histoName:
                if len(mlbMap[histoName]) % 2 == 1 :
                    print "WARNING: incorrect structure for %s"%histoName
                minDR2=[(float('inf'),0,float('inf')),(float('inf'),0,float('inf'))]

                # get pair with lowest DR separation
                for i in range(0,len(mlbMap[histoName]),2) :
                    mlb,wt = mlbMap[histoName][i]
                    mdr,_  = mlbMap[histoName][i+1]

                    if mdr < minDR2[0][2] :
                        minDR2[1] = minDR2[0]
                        minDR2[0] = (mlb,wt,mdr)
                    elif mdr < minDR2[1][2] :
                        minDR2[1] = (mlb,wt,mdr)

                mlbMap[histoName] = [(mlb,wt) for mlb,wt,_ in minDR2]

            # fill
            for mlb,wt in mlbMap[histoName] :
                if mlb == float('inf') : continue
                observablesH[histoName].Fill(mlb,wt)

    #save results
    fOut=ROOT.TFile.Open(outFileName,'RECREATE')
    for var in observablesH: observablesH[var].Write()
    fOut.Close()


"""
Wrapper for when the analysis is run in parallel
"""
def runTopWidthAnalysisPacked(args):
    try:
        fileNames,outFileName=args
        runTopWidthAnalysis(fileNames,outFileName)
    except : # ReferenceError:
        print 50*'<'
        print "  Problem with", name, "continuing without"
        print 50*'<'
        return False

"""
Create analysis tasks
"""
def createAnalysisTasks(opt):

    onlyList=opt.only.split('v')

    ## Local directory
    file_list=[]
    if os.path.isdir(opt.input):
        for file_path in os.listdir(opt.input):
            if file_path.endswith('.root'):
                file_list.append(os.path.join(opt.input,file_path))
    elif opt.input.startswith('/store/'):
        file_list = getEOSlslist(opt.input)
    elif '.root' in opt.input:
        file_list.append(opt.input)

    #list of files to analyse
    tasklist=[]
    for filename in file_list:
        baseFileName=os.path.basename(filename)
        tag,ext=os.path.splitext(baseFileName)
        if len(onlyList)>0:
            processThis=False
            for filtTag in onlyList:
                if filtTag in tag:
                    processThis=True
            if not processThis : continue
        tasklist.append((filename,'%s/%s'%(opt.output,baseFileName)))

    #loop over tasks
    if opt.queue=='local':
        if opt.jobs>1:
            print ' Submitting jobs in %d threads' % opt.jobs
            import multiprocessing as MP
            pool = MP.Pool(opt.jobs)
            pool.map(runTopWidthAnalysisPacked,tasklist)
        else:
            for fileName,outFileName in tasklist:
                runTopWidthAnalysis(fileName,outFileName)
    else:
        cmsswBase=os.environ['CMSSW_BASE']
        for fileName,_ in tasklist:
            localRun='python %s/src/TopLJets2015/TopAnalysis/scripts/runTopWidthAnalysis.py -i %s -o %s -q local'%(cmsswBase,fileName,opt.output)
            cmd='bsub -q %s %s/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh \"%s\"' % (opt.queue,cmsswBase,localRun)
            print cmd
            os.system(cmd)

"""
steer
"""
def main():
	usage = 'usage: %prog [options]'
	parser = optparse.OptionParser(usage)
	parser.add_option('-i', '--input',
                          dest='input',
                          default='/afs/cern.ch/user/p/psilva/work/TopWidth',
                          help='input directory with the files [default: %default]')
	parser.add_option('--jobs',
                          dest='jobs',
                          default=1,
                          type=int,
                          help='# of jobs to process in parallel the trees [default: %default]')
	parser.add_option('--only',
                          dest='only',
                          default='',
                          type='string',
                          help='csv list of tags to process')
	parser.add_option('-o', '--output',
                          dest='output',
                          default='analysis',
                          help='Output directory [default: %default]')
	parser.add_option('-q', '--queue',
                          dest='queue',
                          default='local',
                          help='Batch queue to use [default: %default]')
	(opt, args) = parser.parse_args()

        ROOT.FWLiteEnabler.enable()
	os.system('mkdir -p %s' % opt.output)

        createAnalysisTasks(opt)


if __name__ == "__main__":
	sys.exit(main())
