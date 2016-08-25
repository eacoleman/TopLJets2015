#!/bin/bash

WHAT=$1; 
ERA=$2
if [ "$#" -ne 2 ]; then 
    echo "steerTOPWidthAnalysis.sh <SEL/MERGESEL/PLOTSEL/WWWSEL/ANA/MERGE/BKG/PLOT/WWW> <ERA>";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGESEL     - merge the output of the jobs";
    echo "        PLOTSEL      - runs the plotter tool on the selection";
    echo "        WWWSEL       - moves the plots to an afs-web-based area";
    echo "        ANA          - analyze the selected events";
    echo "        MERGE        - merge the output of the analysis jobs";
    echo "        BKG          - estimate DY scale factor from data";
    echo "        PLOT         - runs the plotter tool on the analysis outputs";
    echo "        WWW          - moves the analysis plots to an afs-web-based area";
    echo " "
    echo "        ERA          - era2015/era2016";
    exit 1; 
fi

export LSB_JOB_REPORT_MAIL=N

queue=2nw
githash=8db9ad6
lumi=12870
lumiUnc=0.062
eosdir=/store/cmst3/user/psilva/LJets2016/${githash}
summaryeosdir=/store/cmst3/group/top/summer2016/TopWidth_${ERA}_ichepv2
case $ERA in
    era2015)
	githash=8c1e7c9;
	lumi=2267.84
	lumiUnc=0.027
	eosdir=/store/cmst3/user/psilva/LJets2015/${githash}
	summaryeosdir=/store/cmst3/group/top/summer2016/TopWidth_${ERA}_notrig 
        #/store/cmst3/group/top/summer2016/TopWidth_${ERA}
	;;
esac

outdir=/afs/cern.ch/work/e/ecoleman/public/TopWidth/TopWidth_${ERA}_widthx4
wwwdir=~/www/TopWidth_${ERA}_widthx4


RED='\e[31m'
NC='\e[0m'
case $WHAT in
    SEL )
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${outdir} --era ${ERA} -m TOP-16-019::RunTop16019 --ch 0;
	;;
    MERGESEL )
	./scripts/mergeOutputs.py ${outdir} True;	
	;;
    PLOTSEL )
	python scripts/plotter.py -i ${outdir} --puNormSF puwgtctr  -j data/${ERA}/samples.json -l ${lumi} --saveLog --mcUnc ${lumiUnc};	
	;;
    WWWSEL )
	mkdir -p ${wwwdir}/sel
	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel
	cp test/index.php ${wwwdir}/sel
	;;
    ANA )
	python scripts/runTopWidthAnalysis.py -i ${summaryeosdir} -o ${outdir}/analysis -q ${queue};
	;;
    MERGE )
	./scripts/mergeOutputs.py ${outdir}/analysis;
	;;
    BKG )
	python scripts/plotter.py -i ${outdir}/analysis  -j data/${ERA}/samples.json  -l ${lumi} --onlyData --only mll -o dy_plotter.root;        
	python scripts/runDYRinRout.py --in ${outdir}/analysis/plots/dy_plotter.root --categs 1b,2b --out ${outdir}/analysis/plots/;
	;;
    PLOT )
	#python scripts/plotter.py -i ${outdir}/analysis  -j data/${ERA}/samples.json      -l ${lumi} --mcUnc ${lumiUnc} --only count --saveTeX -o count_plotter.root --procSF DY:${outdir}/analysis/plots/.dyscalefactors.pck; 
        #python scripts/plotter.py -i ${outdir}/analysis  -j data/${ERA}/samples.json      -l ${lumi} --mcUnc ${lumiUnc} --onlyData --procSF DY:${outdir}/analysis/plots/.dyscalefactors.pck;
	#python scripts/plotter.py -i ${outdir}/analysis  -j data/${ERA}/syst_samples.json -l ${lumi} --mcUnc ${lumiUnc} --silent -o syst_plotter.root;        
	#combined plots
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py ptlb
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py ptlb EE1b,MM1b,EM1b
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py ptlb EE2b,MM2b,EM2b
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_1.0w lowptEE1b,lowptMM1b,lowptEM1b
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_1.0w lowptEE2b,lowptMM2b,lowptEM2b
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_1.0w highptEE1b,highptMM1b,highptEM1b
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_1.0w highptEE2b,highptMM2b,highptEM2b
        ;;
    WWW )
        mkdir -p ${wwwdir}/ana
        cp ${outdir}/analysis/plots/*.{png,pdf} ${wwwdir}/ana        
        cp test/index.php ${wwwdir}/ana
	mkdir -p ${wwwdir}/comb
        cp plots/*.{root,png,pdf} ${wwwdir}/comb 
        cp test/index.php ${wwwdir}/comb
	;;
esac
