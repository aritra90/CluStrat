#!/bin/sh -l
source /depot/pdrineas/data/gwas_scripts/software_paths.conf

######---------------------------------------------------------
######    Last updated: 2017 Jul 17
######---------------------------------------------------------


######---------------------------------------------------------
######    Input arguments

function usage()
{
    echo "Usage: $0 [--geno1 <float>][--geno2 <float>][--ibrc <float>]"
}

PARAMS=""

# Default values
geno1=0.05
geno2=0.02
ibrc=0.2

float_re='^[+-]?[0-9]+([.][0-9]+)?$'
# int_re='^[0-9]+$'

while (( "$#" )); do
  case "$1" in
    --geno1)
        #echo "$2"
        geno1=$2
        if ! [[ $geno1 =~ $float_re ]]
        then
            usage
            echo "Not a float!!"
            exit 1
        fi
        if (( $(echo "$geno1 > 1" |bc -l) ||  $(echo "$geno1 < 0" |bc -l)  )); then
            usage
            echo "Enter a float between 0 and 1"
            exit 1
        fi
        shift 2
        ;;
    --geno2)
      #echo "$2"
      geno2=$2
      if ! [[ $geno2 =~ $float_re ]]
      then
          usage
          echo "Not a float!!"
          exit 1
      fi
      if (( $(echo "$geno2 > 1" |bc -l) ||  $(echo "$geno2 < 0" |bc -l)  )); then
          usage
          echo "Enter a float between 0 and 1"
          exit 1
      fi
      shift 2
      ;;
    --ibrc)
        #echo "$2"
        ibrc=$2
        if ! [[ $ibrc =~ $float_re ]]
        then
            usage
            echo "Not a float!!"
            exit 1
        fi
        if (( $(echo "$ibrc > 1" |bc -l) ||  $(echo "$ibrc < 0" |bc -l)  )); then
            usage
            echo "Enter a float between 0 and 1"
            exit 1
        fi
        shift 2
        ;;
    --) # end argument parsing
      shift
      break
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done
# set positional arguments in their proper place
eval set -- "$PARAMS"
echo "Value of SNP call rate: "$geno1
echo "Value of individual call rate: "$geno2
echo "Value of interbreeding coefficent: "$ibrc

dataset=`readlink -e $1.bed | sed 's/\.bed//'`
status=$2
memory=$3
######---------------------------------------------------------



######---------------------------------------------------------
######    Keep a log
date >> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
echo 'Running dataset_qc.sh with the following arguments:' \
					>> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
for i; do
	echo -e '\t'$i >> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
done
echo -e '\n' >> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
######---------------------------------------------------------



######---------------------------------------------------------
######     Start analysis  ---  Part 1 - Individual QC
######---------------------------------------------------------

mkdir ${dataset##*/}_qc
cd ${dataset##*/}_qc
echo '\documentclass{article}

\usepackage[utf8]{inputenc}
\usepackage[a4paper, margin=1in]{geometry}
\usepackage{booktabs}
\usepackage{authblk}
\usepackage{graphicx}
\usepackage{float}
\usepackage{listings}
\usepackage{tabularx}
\usepackage{array}
\usepackage{underscore}
\usepackage{adjustbox}
\newcolumntype{L}{>{$}l<{$}}

\begin{document}

' > report.tex

echo '\documentclass{beamer}

\usepackage[utf8]{inputenc}
\usepackage{booktabs}
\usepackage{graphicx}
\usepackage{float}
\usepackage{listings}
\usepackage{tabularx}
\usepackage{array}
\usepackage{underscore}
\usepackage{adjustbox}
\newcolumntype{L}{>{$}l<{$}}

\usetheme{default}

\begin{document}
' > slides.tex


echo -e '\n'
echo -e '🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔'
echo -e '🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔'
echo -e '🙒🙜 🙒🙜 🙒🙜                                                             🙞🙔 🙞🙔 🙞🙔'
echo -e '🙒🙜 🙒🙜                                                                   🙞🙔 🙞🙔'
echo -e '🙒🙜                                                                         🙞🙔'



###################
###          Get crude sample numbers
###################
echo `awk '{if($6==2) print}' ${dataset}.fam | \
		wc -l`' Cases ( '`awk '{if($5==1 && $6==2) print}' ${dataset}.fam | \
		wc -l`' Males and '`awk '{if($5==2 && $6==2) print}' ${dataset}.fam | \
		wc -l`' Females )
'`awk '{if($6==1) print}' ${dataset}.fam | \
		wc -l`' Controls ( '`awk '{if($5==1 && $6==1) print}' ${dataset}.fam | \
		wc -l`' Males and '`awk '{if($5==2 && $6==1) print}' ${dataset}.fam | \
		wc -l`' Females )
'`awk '{if($6==0) print}' ${dataset}.fam | \
		wc -l`' Unknown phenotype ( '`awk '{if($5==1 && $6==0) print}' ${dataset}.fam | \
		wc -l`' Males and '`awk '{if($5==2 && $6==0) print}' ${dataset}.fam | \
		wc -l`' Females )

'`awk '{if($5==1) print}' ${dataset}.fam | \
		wc -l`' Males ( '`awk '{if($5==1 && $6==2) print}' ${dataset}.fam | \
		wc -l`' Cases and '`awk '{if($5==1 && $6==1) print}' ${dataset}.fam | \
		wc -l`' Controls )
'`awk '{if($5==2) print}' ${dataset}.fam | \
		wc -l`' Females ( '`awk '{if($5==2 && $6==2) print}' ${dataset}.fam | \
		wc -l`' Cases and '`awk '{if($5==2 && $6==1) print}' ${dataset}.fam | \
		wc -l`' Controls )
'`awk '{if($5==0) print}' ${dataset}.fam | \
		wc -l`' Unknown sex ( '`awk '{if($5==0 && $6==2) print}' ${dataset}.fam | \
		wc -l`' Cases and '`awk '{if($5==0 && $6==1) print}' ${dataset}.fam | \
		wc -l`' Controls )'			 > status_summary.txt

awk '{if($5==0) print}' ${dataset}.fam > no_sex.ind
awk '{if($6==0) print}' ${dataset}.fam > no_pheno.ind









###################
###          First filter SNPs with call rate < 0.95 (geno1)
###################

echo -e '\nStep ❶'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --bfile ${dataset} \
				--allow-no-sex \
				--missing \
				--out ${dataset##*/}_miss \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

############################################################################
####################### MYSON CMD LINE CHANGES #############################
############################################################################

awk -v var="$geno1" '{if($5>var) print $2, $0}' ${dataset##*/}_miss.lmiss \
									> ${dataset##*/}_over005.snps
awk -v var="$geno2" '{if($6>var) print $2, $0}' ${dataset##*/}_miss.imiss \
									> ${dataset##*/}_over005.ind

############################################################################
####################### MYSON CMD LINE CHANGES #############################
############################################################################

echo -e '\nStep ❷'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --bfile ${dataset} \
				--cluster missing \
				--out ${dataset##*/}_miss_clustermiss \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

gzip ${dataset##*/}_miss_clustermiss.cluster3.missing
gzip ${dataset##*/}_miss_clustermiss.mdist.missing

echo -e '\n\n'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
echo -e '                            Starting plotting round'
python ${SCRIPTPATH}heatmap_from_matrix.py ${dataset##*/}_miss_clustermiss.mdist.missing.gz

python ${SCRIPTPATH}histogram.py ${dataset##*/}_miss.lmiss \
									-col=4 \
									-header \
									-vline=0.05 \
									-binsize=0.01
python ${SCRIPTPATH}histogram.py ${dataset##*/}_miss.lmiss \
									-col=4 \
									-header \
									-vline=0.05 \
									-binsize=0.01 \
									-log
python ${SCRIPTPATH}histogram.py ${dataset##*/}_miss.imiss \
									-col=5 \
									-header \
									-vline=0.02 \
									-binsize=0.01
python ${SCRIPTPATH}histogram.py ${dataset##*/}_miss.imiss \
									-col=5 \
									-header \
									-vline=0.02 \
									-binsize=0.01 \
									-log
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'





###################
###          Filter Individuals with call rate < 0.98, inbreeding coefficients > 0.2
###                and sexcheck, if there are sex chromosome data
###################

# >&2 echo "MARKER 1"

if [ $(awk '{if($1==23) print}' ${dataset}.bim | wc -l) -eq 0 ]; then
	echo -e '\nStep ❸'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --bfile ${dataset} \
					--allow-no-sex \
					--missing \
					--exclude ${dataset##*/}_over005.snps \
					--het \
					--out ${dataset##*/}_over005_mhs \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
else
	echo -e '\nStep ❸'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --bfile ${dataset} \
					--allow-no-sex \
					--missing \
					--exclude ${dataset##*/}_over005.snps \
					--het \
					--check-sex \
					--out ${dataset##*/}_over005_mhs \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	python ${SCRIPTPATH}histogram.py ${dataset##*/}_over005_mhs.sexcheck \
										-col=5 \
										-header \
										-binsize=0.01 \
										-log
fi

echo -e '\n\n'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
echo -e '                            Starting plotting round'
python ${SCRIPTPATH}histogram.py ${dataset##*/}_over005_mhs.het \
									-col=5 \
									-header \
									-vline=0.2,-0.2 \
									-binsize=0.01 \
									-log

python ${SCRIPTPATH}histogram.py ${dataset##*/}_over005_mhs.imiss \
									-col=5 \
									-header \
									-vline=0.02 \
									-binsize=0.01 \
									-log
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

############################################################################
####################### MYSON CMD LINE CHANGES #############################
############################################################################


awk -v var="$geno2" '{if($6>var) print}' ${dataset##*/}_over005_mhs.imiss \
				| grep -v 'F_MISS'	> ${dataset##*/}_over002.ind
awk -v var="$ibrc" '{if($6>var || $6<-var) print}' ${dataset##*/}_over005_mhs.het \
				| grep -v 'O(HOM)'	> ${dataset##*/}_inbrcoeff.ind

############################################################################
####################### MYSON CMD LINE CHANGES #############################
############################################################################


if [ $(awk '{if($1==23) print}' ${dataset}.bim | wc -l) -eq 0 ]; then
	cat ${dataset##*/}_over002.ind ${dataset##*/}_inbrcoeff.ind > qc_remove.ind
else
	grep PROBLEM ${dataset##*/}_over005_mhs.sexcheck \
				| grep -v 'PEDSEX' | sort -nk6,6	> sexcheck_problematic.ind
	awk '{if($6>0.25 && $6<0.75) print}' ${dataset##*/}_over005_mhs.sexcheck \
													> sexcheck_remove.ind
	cat ${dataset##*/}_over002.ind ${dataset##*/}_inbrcoeff.ind \
								sexcheck_remove.ind > qc_remove.ind
fi

echo -e '\nStep ❹'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --bfile ${dataset} \
				--allow-no-sex \
				--remove qc_remove.ind \
				--make-bed \
				--out ${dataset##*/}_qcind \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'





###################
###          Get crude after QC sample numbers
###################
echo `awk '{if($6==2) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Cases ( '`awk '{if($5==1 && $6==2) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Males and '`awk '{if($5==2 && $6==2) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Females )
'`awk '{if($6==1) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Controls ( '`awk '{if($5==1 && $6==1) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Males and '`awk '{if($5==2 && $6==1) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Females )
'`awk '{if($6==0) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Unknown phenotype ( '`awk '{if($5==1 && $6==0) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Males and '`awk '{if($5==2 && $6==0) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Females )

'`awk '{if($5==1) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Males ( '`awk '{if($5==1 && $6==2) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Cases and '`awk '{if($5==1 && $6==1) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Controls )
'`awk '{if($5==2) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Females ( '`awk '{if($5==2 && $6==2) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Cases and '`awk '{if($5==2 && $6==1) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Controls )
'`awk '{if($5==0) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Unknown sex ( '`awk '{if($5==0 && $6==2) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Cases and '`awk '{if($5==0 && $6==1) print}' ${dataset##*/}_qcind.fam | \
		wc -l`' Controls )' 		> status_summary_afterqc.txt

echo '' >> status_summary_afterqc.txt
echo `wc -l  ${dataset##*/}_over002.ind \
	| awk '{print $1}'`' Individuals with over 2% missingness - removed' \
									>> status_summary_afterqc.txt
echo `wc -l  ${dataset##*/}_inbrcoeff.ind \
	| awk '{print $1}'`' Individuals with over 20% absolute inbreeding coefficient - removed' \
									>> status_summary_afterqc.txt
if [ $(awk '{if($1==23) print}' ${dataset}.bim | wc -l) -eq 0 ]; then
	echo 'No sex chromosome data to check for discrepancies between genomic and phenotypic sex'
else
	echo `wc -l sexcheck_problematic.ind \
	| awk '{print $1}'`' Individuals with genomic and phenotypic sex discrepancies - recheck' \
									>> status_summary_afterqc.txt
fi











######---------------------------------------------------------
######     Start analysis  ---  Part 2 - SNP QC
######---------------------------------------------------------



###################
###          Remove SNPs with call rate < 0.98 (geno2)
###################


echo -e '\nStep ❺'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --bfile ${dataset##*/}_qcind \
				--allow-no-sex \
				--missing \
				--test-missing \
				--test-mishap \
				--cluster missing \
				--hardy \
				--freq \
				--out ${dataset##*/}_qcind_snp \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

echo -e '\n\n'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
echo -e '                            Starting plotting round'
gzip -f ${dataset##*/}_qcind_snp.cluster3.missing
gzip -f ${dataset##*/}_qcind_snp.mdist.missing
python ${SCRIPTPATH}heatmap_from_matrix.py ${dataset##*/}_qcind_snp.mdist.missing.gz


############################################################################
####################### MYSON CMD LINE CHANGES #############################
############################################################################

awk -v var="$geno2" '{if($5>var) print $2, $0}' ${dataset##*/}_qcind_snp.lmiss \
									> ${dataset##*/}_qcind_over002.snps

############################################################################
####################### MYSON CMD LINE CHANGES #############################
############################################################################

python ${SCRIPTPATH}histogram.py ${dataset##*/}_qcind_snp.lmiss \
									-col=4 \
									-header \
									-vline=0.02 \
									-binsize=0.01
python ${SCRIPTPATH}histogram.py ${dataset##*/}_qcind_snp.lmiss \
									-col=4 \
									-header \
									-vline=0.02 \
									-binsize=0.01 \
									-log
python ${SCRIPTPATH}histogram.py ${dataset##*/}_qcind_snp.imiss \
									-col=5 \
									-header \
									-vline=0.02 \
									-binsize=0.01
python ${SCRIPTPATH}histogram.py ${dataset##*/}_qcind_snp.imiss \
									-col=5 \
									-header \
									-vline=0.02 \
									-binsize=0.01 \
									-log


if [ "${status}" == "mixed" ] ; then
	###   Check missingness and differential missingness between cases and controls
	awk '{print $0,$3 - $4}' ${dataset##*/}_qcind_snp.missing \
										> ${dataset##*/}_qcind_snp_casecontrol.diff

    ############################################################################
    ####################### MYSON CMD LINE CHANGES #############################
    ############################################################################


	awk -v var="$geno2" '{if($6 > var || $6 < -var) print}' ${dataset##*/}_qcind_snp_casecontrol.diff \
										> casecontrol_over002.snps


    ############################################################################
    ####################### MYSON CMD LINE CHANGES #############################
    ############################################################################

	python ${SCRIPTPATH}histogram.py ${dataset##*/}_qcind_snp_casecontrol.diff \
										-col=5 \
										-header \
										-vline=0.02,-0.02 \
										-binsize=0.01
	python ${SCRIPTPATH}histogram.py ${dataset##*/}_qcind_snp_casecontrol.diff \
										-col=5 \
										-header \
										-vline=0.02,-0.02 \
										-binsize=0.01 \
										-log
fi

###   Use HWE exact test to filter
###   Check for how many cases and controls and filter accordingly

if [ "$(awk '{if($6==1) print}' ${dataset}.fam | wc -l)" -lt 10 ] ; then
	awk '{if($9<0.000001 && $3=="AFF") print}' ${dataset##*/}_qcind_snp.hwe \
									> ${dataset##*/}_qcind_hardy_over1e6.snps
	cp ${dataset##*/}_qcind_hardy_over1e6.snps ${dataset##*/}_qcind_hardy.snps
elif [ "$(awk '{if($6==2) print}' ${dataset}.fam | wc -l)" -lt 10 ] ; then
	awk '{if($9<0.000001 && $3=="UNAFF") print}' ${dataset##*/}_qcind_snp.hwe \
									> ${dataset##*/}_qcind_hardy_over1e6.snps
	cp ${dataset##*/}_qcind_hardy_over1e6.snps ${dataset##*/}_qcind_hardy.snps
else
	cat <(awk '{if($9<0.0000000001 && $3=="AFF") print}' ${dataset##*/}_qcind_snp.hwe) \
		<(awk '{if($9<0.000001 && $3=="UNAFF") print}' ${dataset##*/}_qcind_snp.hwe) \
									> ${dataset##*/}_qcind_hardyover1e61e10.snps
	cp ${dataset##*/}_qcind_hardyover1e61e10.snps ${dataset##*/}_qcind_hardy.snps
fi

awk '{if($9!=NA) print $1,$2,$3,$4,$5,$6,$7,$8,-log($9)/log(10)}' ${dataset##*/}_qcind_snp.hwe \
									> ${dataset##*/}_qcind_hardy_log.hwe

python ${SCRIPTPATH}histogram.py ${dataset##*/}_qcind_hardy_log.hwe \
									-col=8 \
									-header \
									-binsize=0.1 \
									-vline=6,10 \
									-log

###   Check frequencies
awk '{if($5 == 0) print $2,$0}' ${dataset##*/}_qcind_snp.frq \
									> monomorphic.snps

python ${SCRIPTPATH}histogram.py ${dataset##*/}_qcind_snp.frq \
									-col=4 \
									-header \
									-binsize=0.01
python ${SCRIPTPATH}histogram.py ${dataset##*/}_qcind_snp.frq \
									-col=4 \
									-header \
									-binsize=0.01 \
									-log
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'



###   Remove unplaced SNPs
awk '{if($1==0) print $2,$0}' ${dataset##*/}_qcind.bim \
									> unplaced.snps




###   Gather all SNPs to be removed in a single file
if [ "${status}" == "mixed" ] ; then
	cat casecontrol_over002.snps ${dataset##*/}_qcind_over002.snps \
		${dataset##*/}_qcind_hardy.snps unplaced.snps \
										> qc_remove.snps
else
	cat ${dataset##*/}_qcind_over002.snps ${dataset##*/}_qcind_hardy.snps unplaced.snps \
										> qc_remove.snps
fi

# >&2 echo "MARKER 2"

###################
###          Filter dataset to get the final clean file
###################
echo -e '\nStep ❻'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --bfile ${dataset##*/}_qcind \
				--allow-no-sex \
				--exclude qc_remove.snps \
				--make-bed \
				--out ${dataset##*/}_qcind_qcsnp \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

echo -e '\nStep ❼'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --bfile ${dataset##*/}_qcind_qcsnp \
				--allow-no-sex \
				--cluster missing \
				--out ${dataset##*/}_qcind_qcsnp_clustermiss \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
gzip -f ${dataset##*/}_qcind_qcsnp_clustermiss.mdist.missing
python ${SCRIPTPATH}heatmap_from_matrix.py ${dataset##*/}_qcind_qcsnp_clustermiss.mdist.missing.gz




###################
###          Get crude after QC SNP numbers
###################
echo `wc -l ${dataset}.bim \
	| awk '{print $1}'`' SNPs in the original dataset' \
									> snp_summary_after_qc.txt
echo `wc -l ${dataset##*/}_qcind_over002.snps \
	| awk '{print $1}'`' SNPs with over 2% missingness - Removed' \
									>> snp_summary_after_qc.txt
if [ "${status}" == "mixed" ] ; then
	echo `wc -l casecontrol_over002.snps \
		| awk '{print $1}'`' SNPs with over 2% difference between cases and controls - Removed' \
										>> snp_summary_after_qc.txt
fi
echo `wc -l ${dataset##*/}_qcind_hardy.snps \
	| awk '{print $1}'`' SNPs failing HWE threshold - Removed' \
									>> snp_summary_after_qc.txt
echo `wc -l unplaced.snps | awk '{print $1}'`' SNPs with no position - Removed' \
									>> snp_summary_after_qc.txt
echo `wc -l monomorphic.snps | awk '{print $1}'`' monomorphic SNPs - Not removed' \
									>> snp_summary_after_qc.txt
echo `wc -l qc_remove.snps | awk '{print $1}'`' SNPs removed in total' \
									>> snp_summary_after_qc.txt
echo `wc -l ${dataset##*/}_qcind_qcsnp.bim \
	| awk '{print $1}'`' SNPs after removing those that failed' \
									>> snp_summary_after_qc.txt





###################
###          Gzip files for less space
###################
gzip -f *.hwe &
gzip -f *.lmiss &
gzip -f *.frq &
gzip -f *.diff &
gzip -f *.missing &
gzip -f *.hap
wait


echo '
\newpage
\section{QC of the dataset \detokenize{'${dataset##*/}'}}
\subsection{Brief introduction}

In this section of the report, lies the description of the QC steps for this particular dataset. For this part, we analyse all datasets separately, focusing on individual-level QC and SNP-level QC. \newline \newline
The steps are outlined as follows: \\
 \\
Step 1: Initial filtering for SNPs with call rate $<0.95$ \newline
 \\
Step 2: Remove samples with:
\begin{itemize}
	\item Call rate $<0.98$
	\item Inbreeding coefficients $>|0.2|$
	\item Reported sex and genomic sex deviation
\end{itemize}

Step 3: Remove SNPs with:
\begin{itemize}
	\item Call rate $<0.98$
	\item Call rate difference between cases and controls $>0.02$
	\item Hardy Weinberg Equilibrium p-value
	\begin{itemize}
		\item $<10^{-6}$ if only cases
		\item $<10^{-6}$ if only controls
		\item $<10^{-6}$ for controls and $<10^{-10}$ for cases if cases and controls
	\end{itemize}
\end{itemize}
\vspace{1em}

\subsection{Individual QC stats}

\begin{table}[H]
	\caption{End results of the QC}
	\centering
	\begin{tabular}{lll}\toprule
		Step&Samples&Action \\ \midrule
		Before QC&'`wc -l ${dataset}.fam | awk '{print $1}'`'& None \\
		Missingness $> 0.02$ &'`grep missingness status_summary_afterqc.txt \
									| awk '{print $1}'`'&Removed \\
		$|$Inbreeding Coeff$| > 0.2$ &'`grep inbreeding status_summary_afterqc.txt \
									| awk '{print $1}'`'&Removed \\
		Sex Discrepancies&'`grep discrepancies status_summary_afterqc.txt \
						| awk '{print $1}'`'&Re-check \\
		After QC&'`wc -l ${dataset##*/}_qcind_qcsnp.fam | awk '{print $1}'`'& None \\
		\bottomrule
	\end{tabular}
\end{table}

%\lstinputlisting[breaklines,
%				firstline=1,
%				frame=single,
%				title={Sample summary in the dataset
%						as delivered in the original file}
%						]{'`readlink -e status_summary.txt`'}

%\lstinputlisting[breaklines,
%				firstline=1,
%				frame=single,
%				title={Sample summary in the dataset after Individual QC}
%						]{'`readlink -e status_summary_afterqc.txt`'}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_miss.imiss_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of initial individual missingness.}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_over005_mhs.imiss_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of individual missingness after removing low-quality SNPs}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_qcind_snp.imiss_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of individual missingness after individual QC}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_over005_mhs.het_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of individual inbreeding coefficents before individual QC}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_over005_mhs.sexcheck_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of individual X chromosome heterozygosity before individual QC}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e heatmap-plot-${dataset##*/}_miss_clustermiss.mdist.missing.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Heatmap of initial missingness clustering of individuals}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e heatmap-plot-${dataset##*/}_qcind_snp.mdist.missing.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Heatmap of missingness clustering of individuals after individual QC}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e heatmap-plot-${dataset##*/}_qcind_qcsnp_clustermiss.mdist.missing.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Heatmap of missingness clustering of individuals after Individual and SNP standard QC}
\end{figure}



\newpage
\subsection{SNP QC stats}

\begin{table}[H]
	\caption{End results of the QC}
	\centering
	\begin{tabular}{lll}\toprule
		Step&SNPs&Action \\ \midrule
		Before QC&'`grep original snp_summary_after_qc.txt | awk '{print $1}'`'& None \\
		Missingness $>0.02$&'`grep missingness snp_summary_after_qc.txt \
								| awk '{print $1}'`'&Removed \\' >> report.tex
if [ "${status}" == "mixed" ] ; then
	echo '		Differential missingness $>0.02$&'`grep difference snp_summary_after_qc.txt \
								| awk '{print $1}'`'&Removed \\' >> report.tex
fi
echo '		Failing Hardy-Weinberg Equilibrium&'`grep HWE snp_summary_after_qc.txt \
													| awk '{print $1}'`'&Removed \\
		Unplaced&'`grep position snp_summary_after_qc.txt | awk '{print $1}'`'&Removed \\
		Removed in total&'`grep total snp_summary_after_qc.txt | awk '{print $1}'`'&Removed \\
		Monomorphic&'`grep monomorphic snp_summary_after_qc.txt | awk '{print $1}'`'&None \\
		After QC&'`grep removing snp_summary_after_qc.txt | awk '{print $1}'`'&None \\
		\bottomrule
	\end{tabular}
\end{table}

%\lstinputlisting[breaklines,
%				firstline=1,
%				frame=single,
%				title={SNP summary in the dataset after SNP QC}
%						]{'`readlink -e snp_summary_after_qc.txt`'}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_miss.lmiss_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of initial SNP missingness}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_qcind_snp.lmiss_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of SNP missingness after individual QC}
\end{figure}
' >> report.tex

if [ "${status}" == "mixed" ] ; then
echo '
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_qcind_snp_casecontrol.diff_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of differential missingness between cases and controls}
\end{figure}
' >> report.tex
fi

echo '
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_qcind_hardy_log.hwe_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of the HWE test p-values}
\end{figure}
' >> report.tex



echo '

\section{QC of the dataset \detokenize{'${dataset##*/}'}}

\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'}}
%\subsection{Individual QC stats}
%\lstinputlisting[breaklines,
%				firstline=1,
%				frame=single,
%				title={Sample summary in the dataset
%						as delivered in the original file}
%						]{'`readlink -e status_summary.txt`'}
\begin{table}[H]
	\caption{End results of the QC}
	\centering
	\begin{tabular}{lll}\toprule
		Step&Samples&Action \\ \midrule
		Before QC&'`wc -l ${dataset}.fam | awk '{print $1}'`'& None \\
		Missingness $> 0.2$ &'`grep missingness status_summary_afterqc.txt \
									| awk '{print $1}'`'&Removed \\
		Inbreeding Coeff $> 0.2$ &'`grep inbreeding status_summary_afterqc.txt \
									| awk '{print $1}'`'&Removed \\
		Sex Discrepancie&'`grep discrepancies status_summary_afterqc.txt \
						| awk '{print $1}'`'&Re-check \\
		Before QC&'`wc -l ${dataset##*/}_qcind_qcsnp.fam | awk '{print $1}'`'& None \\
		\bottomrule
	\end{tabular}
\end{table}
\end{frame}

%\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'}}
%\lstinputlisting[breaklines,
%				firstline=1,
%				frame=single,
%				title={Sample summary in the dataset after Individual QC}
%						]{'`readlink -e status_summary_afterqc.txt`'}
%\end{frame}

\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_miss.imiss_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of initial individual missingness}
\end{figure}
\end{frame}

\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_over005_mhs.imiss_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of individual missingness after removing low-quality SNPs}
\end{figure}
\end{frame}

\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_qcind_snp.imiss_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of individual missingness after individual QC}
\end{figure}
\end{frame}

\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_over005_mhs.het_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of individual inbreeding coefficents before individual QC}
\end{figure}
\end{frame}

\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_over005_mhs.sexcheck_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of individual X chromosome heterozygosity before individual QC}
\end{figure}
\end{frame}

\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e heatmap-plot-${dataset##*/}_miss_clustermiss.mdist.missing.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Heatmap of initial missingness clustering}
\end{figure}
\end{frame}

\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e heatmap-plot-${dataset##*/}_qcind_snp.mdist.missing.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Heatmap of missingness clustering of individuals after individual QC}
\end{figure}
\end{frame}

\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e heatmap-plot-${dataset##*/}_qcind_qcsnp_clustermiss.mdist.missing.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Heatmap of missingness clustering of individuals after Individual and SNP standard QC}
\end{figure}
\end{frame}




\subsection{SNP QC stats}

\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'}}
%\lstinputlisting[breaklines,
%				firstline=1,
%				frame=single,
%				title={SNP summary in the dataset after SNP QC}
%						]{'`readlink -e snp_summary_after_qc.txt`'}
\begin{table}[H]
	\caption{End results of the QC}
	\centering
	\begin{tabular}{lll}\toprule
		Step&SNPs&Action \\ \midrule
		Before QC&'`grep original snp_summary_after_qc.txt | awk '{print $1}'`'& None \\
		Missingness $>0.02$&'`grep missingness snp_summary_after_qc.txt \
								| awk '{print $1}'`'&Removed \\' >> slides.tex
if [ "${status}" == "mixed" ] ; then
	echo '		Differential missingness $>0.02$&'`grep difference snp_summary_after_qc.txt \
								| awk '{print $1}'`'&Removed \\' >> slides.tex
fi
echo '		Failing Hardy-Weinberg Equilibrium&'`grep HWE snp_summary_after_qc.txt \
													| awk '{print $1}'`'&Removed \\
		Unplaced&'`grep position snp_summary_after_qc.txt | awk '{print $1}'`'&Removed \\
		Removed in total&'`grep total snp_summary_after_qc.txt | awk '{print $1}'`'&Removed \\
		Monomorphic&'`grep monomorphic snp_summary_after_qc.txt | awk '{print $1}'`'&None \\
		After QC&'`grep removing snp_summary_after_qc.txt | awk '{print $1}'`'&None \\
		\bottomrule
	\end{tabular}
\end{table}
\end{frame}

\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_miss.lmiss_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of initial SNP missingness}
\end{figure}
\end{frame}

\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_qcind_snp.lmiss_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of SNP missingness after individual QC}
\end{figure}
\end{frame}
' >> slides.tex

if [ "${status}" == "mixed" ] ; then
echo '
\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_qcind_snp_casecontrol.diff_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of differential missingness between cases and controls}
\end{figure}
\end{frame}
' >> slides.tex
fi

echo '
\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-${dataset##*/}_qcind_hardy_log.hwe_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of the HWE test p-values}
\end{figure}
\end{frame}
' >> slides.tex




###################
###          Run association without filtering
###################

if [ "${status}" == "mixed" ] ; then
	echo -e '\nStep ❽'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --bfile ${dataset} \
					--allow-no-sex \
					--assoc \
					--out ${dataset##*/}_assoc \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

	python ${SCRIPTPATH}manhattan_plot.py ${dataset##*/}_assoc.assoc 8,0,2,1 0.05 1
	python ${SCRIPTPATH}qq-plot.py ${dataset##*/}_assoc.assoc 8 1 \
							$(awk '{if($6==2) print}' ${dataset}.fam | wc -l),$(awk '{if($6==1) print}' ${dataset}.fam | wc -l)

	######    Get top results ( p <= 1e-3 )
	cat <(head -1 ${dataset##*/}_assoc.assoc) \
		<(sort -gk9 ${dataset##*/}_assoc.assoc | grep -v NA \
							| awk '{if($9<=1e-3) print}') \
			> ${dataset##*/}_assoc_top.assoc

	echo -e '\nStep ❽'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --annotate ${dataset##*/}_assoc_top.assoc \
								ranges=${REFPATH}plink/glist-hg19 \
					--border 20 \
					--out ${dataset##*/}_assoc_top \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	gzip -f ${dataset##*/}_assoc.assoc

	echo '
	\newpage
	\subsection{Preliminary associations of the dataset}

	After all the quality control procedures, we run two simple associations, one on the dataset, as it was before QC, and one on the dataset after individual and SNP QC.
	\textbf{These associations are not the final ones, because we need to remove IBD and PCA outliers.}\\

	\begin{figure}[H]
		\centering
		\includegraphics[width=\textwidth,keepaspectratio]{{'`readlink -e manhattan-plot-${dataset##*/}_assoc.assoc.png | sed 's/\.png//' `'}.png}
		\caption{Manhattan plot of the association results before Individual and SNP standard QC}
	\end{figure}

	\begin{figure}[H]
		\centering
		\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e qq-plot-${dataset##*/}_assoc.assoc.png | sed 's/\.png//' `'}.png}
		\caption{QQ-plot of the association results before Individual and SNP standard QC}
	\end{figure}

	\begin{table}[H]
		\caption{}
		\centering
		\begin{scriptsize}
		\begin{tabularx}{\textwidth}{lLLLLLLl}\toprule
			'`awk '{print $2,$1,$3,$4,$5,$6,$9,$10,$11}' ${dataset##*/}_assoc_top.annot \
					| awk -F"|" '{print $1}' | head -20 \
					| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
			\bottomrule
		\end{tabularx}
		\end{scriptsize}
	\end{table}
	' >> report.tex

	echo '
	\subsection{Preliminary associations of the dataset}
	\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'} - Crude Association}
		\begin{figure}[H]
			\centering
			\includegraphics[width=\textwidth,keepaspectratio]{{'`readlink -e manhattan-plot-${dataset##*/}_assoc.assoc.png | sed 's/\.png//' `'}.png}
			\caption{Manhattan plot of the association results before
						Individual and SNP standard QC}
		\end{figure}
	\end{frame}

	\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'} - Crude Association}
		\begin{figure}[H]
			\centering
			\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e qq-plot-${dataset##*/}_assoc.assoc.png | sed 's/\.png//' `'}.png}
			\caption{QQ-plot of the association results before
						Individual and SNP standard QC}
		\end{figure}
	\end{frame}

	\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'} - Crude Association}
	\begin{table}[H]
		\caption{}
		\centering
		\begin{tiny}
		\begin{tabularx}{\textwidth}{lLLLLLLl}\toprule
			'`awk '{print $2,$1,$3,$4,$5,$6,$9,$10,$11}' ${dataset##*/}_assoc_top.annot \
					| awk -F"|" '{print $1}' | head -20 \
					| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
			\bottomrule
		\end{tabularx}
		\end{tiny}
	\end{table}
	\end{frame}
	' >> slides.tex
fi


###################
###          Run association on clean data
###################
if [ "${status}" == "mixed" ] ; then
	echo -e '\nStep ❾'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --bfile ${dataset##*/}_qcind_qcsnp \
					--allow-no-sex \
					--assoc \
					--out ${dataset##*/}_qcind_qcsnp_assoc \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

	python ${SCRIPTPATH}manhattan_plot.py ${dataset##*/}_qcind_qcsnp_assoc.assoc 8,0,2,1 0.05 1
	python ${SCRIPTPATH}qq-plot.py ${dataset##*/}_qcind_qcsnp_assoc.assoc 8 1 \
							$(awk '{if($6==2) print}' ${dataset}.fam | wc -l),$(awk '{if($6==1) print}' ${dataset}.fam | wc -l)

	######    Get top results ( p <= 1e-3 )
	cat <(head -1 ${dataset##*/}_qcind_qcsnp_assoc.assoc) \
		<(sort -gk9 ${dataset##*/}_qcind_qcsnp_assoc.assoc | grep -v NA \
							| awk '{if($9<=1e-3) print}') \
			> ${dataset##*/}_qcind_qcsnp_assoc_top.assoc

	echo -e '\nStep ❾'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --annotate ${dataset##*/}_qcind_qcsnp_assoc_top.assoc \
								ranges=${REFPATH}plink/glist-hg19 \
					--border 20 \
					--out ${dataset##*/}_qcind_qcsnp_assoc_top \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	gzip -f ${dataset##*/}_qcind_qcsnp_assoc.assoc

	echo '
	\newpage
	\begin{figure}[H]
		\centering
		\includegraphics[width=\textwidth,keepaspectratio]{{'`readlink -e manhattan-plot-${dataset##*/}_qcind_qcsnp_assoc.assoc.png | sed 's/\.png//' `'}.png}
		\caption{Manhattan plot of the association results after Individual and SNP standard QC}
	\end{figure}

	\begin{figure}[H]
		\centering
		\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e qq-plot-${dataset##*/}_qcind_qcsnp_assoc.assoc.png | sed 's/\.png//' `'}.png}
		\caption{QQ-plot of the association results after Individual and SNP standard QC}
	\end{figure}

	\begin{table}[H]
		\caption{}
		\centering
		\begin{scriptsize}
		\begin{tabularx}{\textwidth}{lLLLLLLl}\toprule
			'`awk '{print $2,$1,$3,$4,$5,$6,$9,$10,$11}' ${dataset##*/}_qcind_qcsnp_assoc_top.annot \
					| awk -F"|" '{print $1}' | head -20 \
					| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
			\bottomrule
		\end{tabularx}
		\end{scriptsize}
	\end{table}
	' >> report.tex

	echo '
	\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'} - Clean Association}
	\begin{figure}[H]
		\centering
		\includegraphics[width=\textwidth,keepaspectratio]{{'`readlink -e manhattan-plot-${dataset##*/}_qcind_qcsnp_assoc.assoc.png | sed 's/\.png//' `'}.png}
		\caption{Manhattan plot of the association results after Individual and SNP standard QC}
	\end{figure}
	\end{frame}

	\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'} - Clean Association}
	\begin{figure}[H]
		\centering
		\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e qq-plot-${dataset##*/}_qcind_qcsnp_assoc.assoc.png | sed 's/\.png//' `'}.png}
		\caption{QQ-plot of the association results after Individual and SNP standard QC}
	\end{figure}
	\end{frame}

	\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'} - Clean Association}
	\begin{table}[H]
		\caption{}
		\centering
		\begin{tiny}
		\begin{tabularx}{\textwidth}{lLLLLLLl}\toprule
			'`awk '{print $2,$1,$3,$4,$5,$6,$9,$10,$11}' ${dataset##*/}_qcind_qcsnp_assoc_top.annot \
					| awk -F"|" '{print $1}' | head -20 \
					| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
			\bottomrule
		\end{tabularx}
		\end{tiny}
	\end{table}
	\end{frame}
	' >> slides.tex

fi



###################
###          Make a final table with what has been left
###################

echo '


\subsection{Final QC table of the dataset \detokenize{'${dataset##*/}'}}

These are the results of the quality control procedure. Listed are the number of SNPs and samples before and after QC.\\

\begin{table}[H]
	\caption{End results of the QC}
	\centering
	\begin{tabular}{lllll}\toprule
		Step&Samples&Cases&Controls&SNPs \\ \midrule
		Before QC&'`wc -l ${dataset}.fam | awk '{print $1}'`'&'`awk '{if($6==2) print}' ${dataset}.fam | wc -l`'&'`awk '{if($6==1) print}' ${dataset}.fam | wc -l`'&'`wc -l ${dataset}.bim | awk '{print $1}' `' \\
		After QC&'`wc -l ${dataset##*/}_qcind_qcsnp.fam | awk '{print $1}'`'&'`awk '{if($6==2) print}' ${dataset##*/}_qcind_qcsnp.fam | wc -l`'&'`awk '{if($6==1) print}' ${dataset##*/}_qcind_qcsnp.fam | wc -l`'&'`wc -l ${dataset##*/}_qcind_qcsnp.bim | awk '{print $1}' `' \\
		\bottomrule
	\end{tabular}
\end{table}

Next to follow is the detection of outliers in the IBD and the PCA analysis. Removal of those samples will account for confouding and stratification issues that may inflate the association results.
' >> report.tex

echo '
\begin{frame}{QC of the dataset \detokenize{'${dataset##*/}'}}
\begin{table}[H]
	\caption{End results of the QC}
	\centering
	\begin{tabular}{lllll}\toprule
		Step&Samples&Cases&Controls&SNPs \\ \midrule
		Before QC&'`wc -l ${dataset}.fam | awk '{print $1}'`'&'`awk '{if($6==2) print}' ${dataset}.fam | wc -l`'&'`awk '{if($6==1) print}' ${dataset}.fam | wc -l`'&'`wc -l ${dataset}.bim | awk '{print $1}' `' \\
		After QC&'`wc -l ${dataset##*/}_qcind_qcsnp.fam | awk '{print $1}'`'&'`awk '{if($6==2) print}' ${dataset##*/}_qcind_qcsnp.fam | wc -l`'&'`awk '{if($6==1) print}' ${dataset##*/}_qcind_qcsnp.fam | wc -l`'&'`wc -l ${dataset##*/}_qcind_qcsnp.bim | awk '{print $1}' `' \\
		\bottomrule
	\end{tabular}
\end{table}
\end{frame}
' >> slides.tex

echo -e '\n\n'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
echo -e '               Writing the report and the slides on LaTeX ✄ ✏\n'

echo '\end{document}' >> report.tex
pdflatex -interaction=batchmode report.tex
rm -rf report.aux report.log

echo '\end{document}' >> slides.tex
pdflatex -interaction=batchmode slides.tex
rm -rf slides.aux slides.log slides.nav slides.out slides.snm slides.toc
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

echo -e '🙐🙝                                                                         🙟🙖'
echo -e '🙐🙝 🙐🙝                                                                   🙟🙖 🙟🙖'
echo -e '🙐🙝 🙐🙝 🙐🙝                                                             🙟🙖 🙟🙖 🙟🙖'
echo -e '🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖'
echo -e '🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙐🙝 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖 🙟🙖'
echo -e '\n'
cd ..
