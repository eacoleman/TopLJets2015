#!/bin/bash

WHAT=$1; 
ERA=$2;
if [ "$#" -ne 2 ]; then 
    echo "steerTOPWidthAnalysis.sh <SEL/MERGESEL/PLOTSEL/WWWSEL/ANA/MERGE/PLOT/WWW>";
    echo "        WORKSPACE       - create a workspace with the TopHypoTest model";
    echo "        SCAN            - scans the likelihood";
    echo "        CLs             - calculates CLs using higgs_combine";
    echo "        TOYS            - plots toys and saves locally";
    echo "        MERGE_DATACARDS - merge all datacards for one LFS/WID into one WID datacard";
    echo "        QUANTILES       - plots the quantiles chart for all distributions, by lepton final state";
    echo "        WWW             - moves the analysis plots to an afs-web-based area";
    echo " "
    echo "        ERA          - era2015/era2016";
    exit 1; 
fi

queue=2nd
outdir=/afs/cern.ch/work/e/ecoleman/public/TopWidth/TopWidth_${ERA}_old/
cardsdir=${outdir}/datacards
wwwdir=~/www/TopWidth_${ERA}/
CMSSW_7_4_7dir=~/CMSSW_7_4_7/src/
CMSSW_7_6_3dir=~/CMSSW_8_0_11/src/

unblind=true
nPseudo=1000

lfs=(EE EM MM)
wid=(0p2w 0p4w 0p6w 0p8w 1p0w 1p2w 1p4w 1p6w 1p8w 2p0w 2p2w 2p4w 2p6w 2p8w 3p0w 3p5w 4p0w)
#dists=(mdrmlb minmlb incmlb sncmlb mt2mlb)
dists=(incmlb sncmlb mt2mlb)
cat=(1b 2b)
lbCat=(highpt lowpt)

nuisances=(jes jer pu btag les trig sel toppt Mtop ttPartonShower MEmuR MEmuF MEtot Herwig amcnloFxFx)

RED='\e[31m'
NC='\e[0m'

# Helpers: getting comma-separated lists of variables
distStr=""
lfsStr=""
widStr=""
for dist in ${dists[*]}
do
    if [[ "${distStr}" == "" ]];
    then
        distStr="${dist}"
    else
        distStr="${distStr},${dist}"
    fi
done

for tlfs in ${lfs[*]} ; do
    if [[ "${lfsStr}" == "" ]];
    then
        lfsStr="${tlfs}"
    else
        lfsStr="${lfsStr},${tlfs}"
    fi
done

for twid in ${wid[*]} ; do
    if [[ "${widStr}" == "" ]];
    then
        widStr="${twid}"
    else
        widStr="${widStr},${twid}"
    fi
done
       


case $WHAT in
    SHAPES ) # to get the shapes file / datacards for the full analysis
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/
        export PYTHONPATH=$PYTHONPATH:/usr/lib64/python2.6/site-packages/
        python test/TopWidthAnalysis/createShapesFromPlotter.py \
                -s tbart,tW \
                --dists ${distStr} \
                -o ${outdir}/datacards/ \
                -i ${outdir}/analysis/plots/plotter.root \
                --systInput ${outdir}/analysis/plots/syst_plotter.root \
                --nomorph
    ;;
    MORPH ) # to get the morphs for the full analysis
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/
        export PYTHONPATH=$PYTHONPATH:/usr/lib64/python2.6/site-packages/
        python test/TopWidthAnalysis/createShapesFromPlotter.py \
                -s tbart,tW \
                --dists ${distStr} \
                -o ${outdir}/datacards/ \
                -i ${outdir}/analysis/plots/plotter.root \
                --systInput ${outdir}/analysis/plots/syst_plotter.root \
                --noshapes 
    ;;
    MORPHALL ) # get shapes/morphs replacing data with a given width template 
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/
        export PYTHONPATH=$PYTHONPATH:/usr/lib64/python2.6/site-packages/

        # for each width, generate a truth dataset and produce the proper morph datacards for it
        # --> output in datacards_<wid>/
        for twid in ${wid[*]}; do
           rm -rf ${outdir}/datacards_${twid}
           nohup python test/TopWidthAnalysis/createShapesFromPlotter.py \
                   -s tbart,tW \
                   --dists ${distStr} \
                   -o ${outdir}/datacards_${twid}/ \
                   -i ${outdir}/analysis/plots/plotter.root \
                   --systInput ${outdir}/analysis/plots/syst_plotter.root \
                   --truth ${twid} &
        done
    ;;
    2DBIAS ) # get the bias plot for the 2D scan
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`

        # Merge datacards
        for twid in ${wid[*]}; do
        for dist in ${dist[*]}; do
            cmd="python ${CMSSW_7_4_7dir}/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "

            for tlbCat in ${lbCat[*]} ; do
            for tlfs in ${lfs[*]} ; do
            for tCat in ${cat[*]} ; do
                cardname="${tlbCat}${tlfs}${tCat}=${outdir}/datacards_${twid}"
                cardname="${cardname}/datacard_widfit__${tlbCat}${tlfs}${tCat}_${dist}.dat"
                cmd="${cmd} ${cardname} "
            done
            done
            done
        done
        done
            cmd="${cmd} > ${outdir}/datacards_${twid}/datacard_widfit__${twid}_${dist}.dat"
            eval `${cmd}`

        # Get workspaces 
        for twid in ${wid[*]}; do 
        for dist in ${dist[*]}; do
            cd ${outdir}/datacards_${twid}
            text2workspace.py ${outdir}/datacards_${twid}/datacard_widfit__${twid}_${dist}.dat \
                -o ${outdir}/datacards_${twid}/workspace_biascheck_${twid}_${dist}.root
        done
        done
       
        # Run combine
        for twid in ${wid[*]}; do 
        for dist in ${dist[*]}; do
            cd ${outdir}/datacards_${twid}
            combine ${outdir}/datacards_${twid}/workspace_biascheck_${twid}_${dist}.root \
                -M MultiDimFit --expectSignal=1 --poi r --poi alpha --poi beta \
                --setPhysicsModelParameterLimits alpha=0,9:beta=0,2 \
                --algo=grid --points=1000
        done
        done
    ;;
    MERGE_DATACARDS ) # get all the datacards you could possibly desire for the analysis
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do

        # for a given width, merge all
        allcmd="python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "
        for tlbCat in ${lbCat[*]} ; do

            # for a given width and lbcat, merge all
            lbccmd="python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "
            for tlfs in ${lfs[*]} ; do

                # for a given width, lfs, and lbcat, merge all
                lfscmd="python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "
                for tCat in ${cat[*]} ; do
                    cardname="${tlbCat}${tlfs}${tCat}=${cardsdir}"
                    cardname="${cardname}/datacard__${twid}_${tlbCat}${tlfs}${tCat}_${dist}.dat"
                    allcmd="${allcmd} ${cardname} "
                    lbccmd="${lbccmd} ${cardname} "
                    lfscmd="${lfscmd} ${cardname} "
                done

                lfscmd="${lfscmd} > ${cardsdir}/datacard__${twid}_${tlbCat}${tlfs}_${dist}.dat"
            done

            lbccmd="${lbccmd} > ${cardsdir}/datacard__${twid}_${tlbCat}_${dist}.dat"
        done

        allcmd="${allcmd} > ${cardsdir}/datacard__${twid}_${dist}.dat"
            
        eval $allcmd
        eval $lbccmd
        eval $lfscmd

        done
        done
    ;;
    WORKSPACE ) # generate combine workspaces
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        mkdir ${outdir}
        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do
            
            # All datacards
            echo "Creating workspace for ${twid}${dist}" 
            text2workspace.py ${cardsdir}/datacard__${twid}_${dist}.dat -P \
                HiggsAnalysis.CombinedLimit.TopHypoTest:twoHypothesisTest \
                -m 172.5 --PO verbose --PO altSignal=${twid} --PO muFloating \
                -o ${outdir}/${twid}_${dist}.root 
        done
        done
    ;;
    SCAN ) # perform likelihood scans on Asimov datasets 
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${outdir}
        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do

            # All datacards
            cmd="combine ${outdir}/${twid}_${dist}.root -M MultiDimFit" 
            cmd="${cmd} -m 172.5 -P x --floatOtherPOI=1 --algo=grid --points=200"
            cmd="${cmd} --expectSignal=1 --setPhysicsModelParameters x=0,r=1"
            cmd="${cmd} -n x0_scan_Asimov_${twid}_${dist}"
            if [[ ${unblind} == false ]] ; then 
                echo "Analysis is blinded"
                cmd="${cmd} -t -1"
            fi
            
            bsub -q ${queue} \
                ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                "${outdir}/" "${cmd}" 
        done
        done
    ;;
    CLs ) # get CLs statistics from combine
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${outdir}
        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do

            # All datacards
            echo "Making CLs for ${twid} ${dist}"
            cmd="combine ${outdir}/${twid}_${dist}.root -M HybridNew --seed 8192 --saveHybridResult" 
            cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 6"
            cmd="${cmd} --clsAcc 0 --fullBToys  --generateExt=1 --generateNuis=0"
            cmd="${cmd} --expectedFromGrid 0.5 -n x_pre-fit_exp__${twid}_${dist}"
            #cmd="${cmd} &> ${outdir}/x_pre-fit_exp__${twid}_${dist}.log"

            bsub -q ${queue} ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh "${outdir}/" "${cmd}" 
        done
        done
    ;;
    TOYS ) # get toys distributions from the pseudoexperiments
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${outdir}
        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do
            
            rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
            rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx(\"x_pre-fit_exp__${twid}_${dist}.qvals.root"
            rootcmds="${rootcmds}\",172.5,1,\"x\",1000,\"\",\"${twid}\",\"${dist}\",${unblind})"
            # All datacards
            root -l -q -b \
                higgsCombinex_pre-fit_exp__${twid}_${dist}.HybridNew.mH172.5.8192.quant0.500.root \
                ${rootcmds}
        done
        done
    ;;
    QUANTILES ) # plot quantiles distributions of all toys, get CLsPlot
        cd ${outdir}
        for dist in ${dists[*]} ; do

        # TODO Add unblind command

            # All datacards
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getQuantilesPlot.py \
                -i ${outdir}/ -o ${outdir}/ \
                --wid ${widStr} \
                --dist ${dist}
            
            # Get Separation plots / LaTeX tables as well
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py\
                -i ${outdir}/ -o ${outdir}/ \
                --wid ${widStr} \
                --dist ${dist}
        done
        
        # Get Separation plots for all dists 
        python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py \
            -i ${outdir}/ -o ${outdir}/ \
            --dist ${distStr} \
            --doAll
    ;;
    NUISANCES ) # run once for each nuisance, doing a full analysis to understand effects of systematics 
        for dist in ${dists[*]} ; do
            systList=""
            i=0
        for syst in ${nuisances[*]} ; do
            if [[ "${syst}" != "${nuisances: 0 : 1}" ]]; then
               systList="${systList}," 
            fi
            systList="${systList}${syst}"

            echo "Frozen systematics are now: ${systList[@]}"

            mkdir ${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}
            cd ${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}
            
            # run likelihood scan
            cmd="echo 'Starting job ${dist} ${syst}'"
            for twid in ${wid[*]}; do
                cmd="${cmd}; combine ${outdir}/${twid}_${dist}.root -M MultiDimFit" 
                cmd="${cmd} -m 172.5 -P x --floatOtherPOI=1 --algo=grid --points=200"
                cmd="${cmd} --expectSignal=1 --setPhysicsModelParameters x=0,r=1"
                cmd="${cmd} -n x0_scan_Asimov_${twid}_${dist}"
                cmd="${cmd} --freezeNuisances ${systList}"
                if [[ ${unblind} == false ]] ; then
                  cmd="${cmd} -t -1"   
                fi
            done

            # run CLs
            for twid in ${wid[*]}; do
                cmd="${cmd}; combine ${outdir}/${twid}_${dist}.root -M HybridNew --seed 8192 --saveHybridResult" 
                cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 6"
                cmd="${cmd} --clsAcc 0 --fullBToys  --generateExt=1 --generateNuis=0"
                cmd="${cmd} --expectedFromGrid 0.5 -n x_pre-fit_exp_${twid}_${dist}"
                cmd="${cmd} --freezeNuisances ${systList}"
            done

            # run toy macro 
            for twid in ${wid[*]}; do
                rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
                rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx(\"x_pre-fit_exp_${twid}_${dist}.qvals.root"
                rootcmds="${rootcmds}\",172.5,1,\"x\",1000,\"\",\"${twid}\",\"${dist}\",${unblind})"

                # All datacards
                cmd="${cmd}; root -l -q -b"
                cmd="${cmd} ${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}"
                cmd="${cmd}/higgsCombinex_pre-fit_exp_${twid}_${dist}.HybridNew.mH172.5.8192.quant0.500.root"
                cmd="${cmd} ${rootcmds}"
            done

            # launch job
            bsub -q ${queue} \
                ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                "${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}/" "${cmd}" 

            let "i += 1" 

        done
        done
    ;;
    WWW )
        mkdir -p ${wwwdir}/ana_${ERA}
        cp ${outdir}/analysis/plots/*.{png,pdf} ${wwwdir}/ana_${ERA} 
        cp ${outdir}/*.{png,pdf} ${wwwdir}/ana_${ERA} 
        cp test/index.php ${wwwdir}/ana_${ERA}
    ;;
esac
