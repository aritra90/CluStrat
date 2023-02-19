#!/bin/bash
#SBATCH -A pdrineas
#SBATCH -n 128
#SBATCH -t 1-00:00:00

# set name of job
#SBATCH --job-name=metaclustrat
#SBATCH -o metaclustrat.out
#SBATCH -e metaclustrat.err

# mail alert at (b)eginning, (e)nd and (a)bortion of execution
# SBATCH --mail-type=ALL
# send mail to the following address
# SBATCH --mail-user=mcburch@purdue.edu

# use submission environment
#SBATCH --export=ALL
PATHTOPLINK="/depot/pdrineas/data/gwas_software/plink2.0/plink2"
PATHTOGCTA="/depot/pdrineas/data/GCTA/gcta_v194/gcta64"
# DATA="/scratch/bell/mcburch/MOSAIC_SIM/case_control_glmmsim_5k_5k/merged_simulated_data"
# DATA="/scratch/bell/mcburch/MOSAIC_SIM/case_control_glmmsim_2k_8k/merged_simulated_data"
OLDIFS=$IFS

PVAL="1e-5"
ITERS=100

models=( "BN" "PSD" "TGP" )
mharr=( 10 50 100 )
clustarr=( 4 5 6 )

# mharr=( 10 )
# clustarr=( 5 )
for m in "${models[@]}"
do
  DATA="/depot/pdrineas/data/DIMs/CluStrat/META_CluStrat/${m}_simdata/simdata"
  for l in "${!mharr[@]}"
  do
    for j in "${!clustarr[@]}"
    do

      rm -rf cluster_iids_* pheno_cluster_iids_* metal_input causal_* spurious_*
      for i in $(eval echo "{1..$ITERS}")
      do

        ##################################### SIMULATE DATA #####################################
        rm -rf ${m}_simdata/ metal_input

        echo "simulating data..."
        python data_simulate.py --model ${m}
        awk '{print $2}' /depot/pdrineas/data/DIMs/CluStrat/META_CluStrat/${m}_simdata/simdata.bim | shuf -n 100 > /depot/pdrineas/data/DIMs/CluStrat/META_CluStrat/${m}_simdata/causal_variants.txt
        $PATHTOGCTA --bfile /depot/pdrineas/data/DIMs/CluStrat/META_CluStrat/${m}_simdata/simdata \
                    --simu-causal-loci /depot/pdrineas/data/DIMs/CluStrat/META_CluStrat/${m}_simdata/causal_variants.txt \
                    --simu-hsq 0.5 --simu-k 0.2 --simu-cc 200 856 --out /depot/pdrineas/data/DIMs/CluStrat/META_CluStrat/${m}_simdata/simdata > /dev/null 2>&1


        echo " "
        echo "running pca..."
        "$PATHTOPLINK" --bfile $DATA --pca 2 --silent

        ##################################### META-CLUSTRAT #####################################
        echo " "
        echo "running meta-clustrat..."
        python3 META_CluStrat.py --fn $DATA --maxclust ${clustarr[j]} --mhdistk ${mharr[l]} #> /dev/null 2>&1
        # python3 pca_and_project.py "clustrat"
        # python3 pca_and_project.py "random"

        array=($(ls cluster_iids_*))

        # grab allele frequencies for smaller clusters
        $PATHTOPLINK --bfile ${DATA} --freq --silent

        # plink lmm steps
        echo " "
        echo "running plink lmm per cluster..."
        for i in "${array[@]}"
        do
          echo $i
          # extract genotype data for cluster
          $PATHTOPLINK \
                --bfile ${DATA} \
                --keep $i \
                --make-bed \
                --out ${i}_data \
                --silent
          # get covariates
          testing=$(wc -l $i)
          read -a strarr <<< "$testing"
          if (($strarr < 51));
          then
            $PATHTOPLINK \
                  --bfile ${i}_data \
                  --pca 2 \
                  --read-freq plink2.afreq \
                  --silent
          else
            $PATHTOPLINK \
                  --bfile ${i}_data \
                  --pca 2 \
                  --silent
          fi
          # run plink lmm
          $PATHTOPLINK \
                --bfile ${i}_data \
                --glm 'cols=chrom,pos,ref,alt,test,se,p,beta' \
                --covar plink2.eigenvec \
                --pheno pheno_${i} \
                --pheno-name 'status' \
                --out ${i}_plink_lmm \
                --silent

          # filter association output and add frequencies
          awk '($7 == "ADD" || $7 == "TEST") {print}' ${i}_plink_lmm.status.glm.logistic  > temp # USE THIS FOR BINARY TRAIT
          # awk '($7 == "ADD" || $7 == "TEST") {print}' ${i}_plink_lmm.status.glm.linear  > temp # USE THIS FOR QUANT TRAIT
          awk '{print $5}' plink2.afreq > temp2
          paste temp temp2 > ${i}_plink_lmm.status.glm.logistic.ADDONLYANDFREQ
          rm -rf temp temp2

          echo " "
          echo "writing to metal script..."
          # grabbing just the file name
          IFS='/'
          read -a strarr <<< "$i"
          clusterfilnm=${strarr[${#strarr[@]} - 1]}_plink_lmm.status.glm.logistic.ADDONLYANDFREQ
          processtext=/depot/pdrineas/data/DIMs/CluStrat/META_CluStrat/${clusterfilnm}
          IFS=$OLDIFS

          # create input for METAL
          echo "MARKER   ID" >> metal_input
          echo "ALLELE   REF ALT" >> metal_input
          echo "FREQ     ALT_FREQS" >> metal_input
          echo "EFFECT   BETA" >> metal_input
          echo "STDERR   SE" >> metal_input
          echo "PVAL     P" >> metal_input
          metal_weight=$(wc -l $i | tr ' ' '\n' | head -1)
          echo "DEFAULT  ${metal_weight}" >> metal_input
          # echo "COLUMNCOUNTING     LENIENT" >> metal_input
          echo " " >> metal_input
          echo "PROCESS " $processtext >> metal_input
          echo " " >> metal_input

        done

        echo "OUTFILE   /depot/pdrineas/data/DIMs/CluStrat/META_CluStrat/METAANALYSIS .TBL" >> metal_input
        echo " " >> metal_input
        echo "ANALYZE" >> metal_input
        echo " "
        echo "running metal..."

        rm -rf /depot/pdrineas/data/DIMs/CluStrat/META_CluStrat/METAANALYSIS1.TBL*
        /depot/pdrineas/data/METAL/generic-metal/metal /depot/pdrineas/data/DIMs/CluStrat/META_CluStrat/metal_input > /dev/null 2>&1

        # sort -k 6 -g METAANALYSIS1.TBL > temp
        # mv temp METAANALYSIS1.TBL
        python3 metal_results_clustrat.py METAANALYSIS1.TBL "metal" ${PVAL} ${m}



        ##################################### OTHER METHODS #####################################

        ##################################### PLINK2 #####################################
        # run plink lmm
        echo " "
        echo "running plink2..."
        awk '{print $1, $2, $3}' ${DATA}.phen > pheno_all_iids
        awk 'BEGIN{print "FID IID status"}1' pheno_all_iids > temp; mv temp pheno_all_iids
        $PATHTOPLINK \
              --bfile ${DATA} \
              --pca 2 \
              --silent
        $PATHTOPLINK \
              --bfile ${DATA} \
              --glm 'cols=chrom,pos,ref,alt,test,se,p,beta' \
              --covar plink2.eigenvec \
              --pheno pheno_all_iids \
              --pheno-name 'status' \
              --out PLINK2_GLM \
              --silent

        # quant trait plink
        awk '($7 == "ADD" || $7 == "TEST") {print}' PLINK2_GLM.status.glm.logistic >  temp # | sort -k 10 -g > temp
        # awk '($7 == "ADD" || $7 == "TEST") {print}' PLINK2_GLM.status.glm.linear >  temp # | sort -k 10 -g > temp
        mv temp PLINK2_GLM.status.glm.logistic
        python3 metal_results_clustrat.py PLINK2_GLM.status.glm.logistic "plink2" ${PVAL} ${m}

        # echo "number of significant SNPs in common:"
        # comm -12 <(sort sig_snps_metal) <(sort sig_snps_plink2) | wc -l

        ##################################### ARMITAGE CHI-SQ / GEMMA / EMMAX / EIGENSTRAT #####################################
        echo " "
        python3 compare_others.py --fn ${DATA} --pval ${PVAL} --model ${m}

        rm -rf cluster_iids_* pheno_cluster_iids_*

      done
      python3 compare_plots.py ${clustarr[j]} ${mharr[l]} ${m}

    done
  done

  rm -rf ${m}_compare_spurious/ ${m}_compare_causal/
  mkdir ${m}_compare_spurious/
  mkdir ${m}_compare_causal/
  mv compare_others_${m}_mhdist*spurious.png ${m}_compare_spurious/
  mv compare_others_${m}_mhdist*causal.png ${m}_compare_causal/

done






echo "DONE" | mail -s "META_CluStrat.sh finished" mcburch@purdue.edu

# # end of script
