#!/bin/sh -l
source /depot/pdrineas/data/gwas_scripts/software_paths.conf

######---------------------------------------------------------
######    Last updated: 2017 Jul 17
######---------------------------------------------------------



######---------------------------------------------------------
######    Input arguments
dataset=`readlink -e $1.bed | sed 's/\.bed//'`
fixdup=$2
fixxy=$3
fixindel=$4
memory=$5
######---------------------------------------------------------



######---------------------------------------------------------
######    Keep a log
date >> /depot/pdrineas/data/gwas_scripts/scripthistory_updatedatasets.log
echo 'Running illumina_psycharray.sh with the following arguments:' \
					>> /depot/pdrineas/data/gwas_scripts/scripthistory_updatedatasets.log
for i; do
	echo -e '\t'$i >> /depot/pdrineas/data/gwas_scripts/scripthistory_updatedatasets.log
done
echo -e '\n' >> /depot/pdrineas/data/gwas_scripts/scripthistory_updatedatasets.log
######---------------------------------------------------------



######---------------------------------------------------------
######     Start analysis
######---------------------------------------------------------

######    Start by creating a new folder with the updated file
mkdir ${dataset##*/}_array
cd ${dataset##*/}_array

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

######    Update names to rs nonmeclature, based on the manufacturer official files
unzip ${REFPATH}Beadchips/PsychArray/v1-2/infinium-psycharray-24-v1-2-a1-b144-rsids.zip -d ./
unzip ${REFPATH}Beadchips/PsychArray/v1-2/infinium-psycharray-24-v1-2-a1-annotated.zip -d ./
	### This step allows snps without name (="." in the annotation file) to remain unchanged
awk '{if($2==".") $2=$1; print}' InfiniumPsychArray-24v1-2_A1_b144_rsids.txt \
									> InfiniumPsychArray-24v1-2_A1_b144_rsids.txt.fixed
rm -rf InfiniumPsychArray-24v1-2_A1_b144_rsids.txt
	### This step creates an annotation file to be used for the further associations
paste <(sort -gk1,1 InfiniumPsychArray-24v1-2_A1_b144_rsids.txt.fixed) \
		<(sort -gk1,1 InfiniumPsychArray-24v1-2_A1.annotated.txt) \
		| cut -f2,4- | sed 's/	/ /g' | gzip -c > snp_annotation.txt.gz
rm -rf InfiniumPsychArray-24v1-2_A1.annotated.txt



echo -e '\nStep ❶'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --bfile ${dataset} \
				--allow-no-sex \
				--update-name InfiniumPsychArray-24v1-2_A1_b144_rsids.txt.fixed \
				--make-bed \
				--out ${dataset##*/}_rs \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

gzip InfiniumPsychArray-24v1-2_A1_b144_rsids.txt.fixed

mv ${dataset##*/}_rs.fam ${dataset##*/}_rs.fam0
awk '{print $2,$2,$3,$4,$5,$6}' ${dataset##*/}_rs.fam0 > ${dataset##*/}_rs.fam
mv ${dataset##*/}_rs.bim ${dataset##*/}_rs.bim0

sed 's/psy_rs/rs/g' ${dataset##*/}_rs.bim0 | sed 's/newrs/rs/g' \
	| sed 's/exm-rs/rs/g' | sed 's/snv-rs/rs/g' | sed 's/_/-/g' \
									> ${dataset##*/}_rs.bim.upd


### Keep some stats of what changed after all these name-updates
python ${SCRIPTPATH}bim-analyser.py ${dataset}.bim \
									> bim_before_update.txt
python ${SCRIPTPATH}bim-analyser.py ${dataset##*/}_rs.bim0 \
									> bim_after_update1.txt
python ${SCRIPTPATH}bim-analyser.py ${dataset##*/}_rs.bim.upd \
									> bim_after_update2.txt

echo `grep -w 'I\|D' ${dataset##*/}_rs.bim.upd | wc -l`" indels in the dataset" \
									>> bim_after_update2.txt
echo `awk '{if($1>22) print}' ${dataset##*/}_rs.bim.upd | \
		wc -l`" non-autosomal SNPs in the dataset" \
									>> bim_after_update2.txt


######    After updates, discover the name-duplicates and merge them
python ${SCRIPTPATH}bim_dup_marker.py ${dataset##*/}_rs.bim.upd
mv ${dataset##*/}_rs.bim.upd.markdup ${dataset##*/}_rs.bim

nextdatasetname=${dataset##*/}_rs








if [ "${fixdup}" == "fixdup" ] ; then
	####    Check concordance of duplicate SNPs
	grep DUP ${nextdatasetname}.bim | awk '{print $2}' \
										> dups.snps
	sed 's/_DUP.//g' dups.snps | sort | uniq \
										> orig.snps
	cat dups.snps orig.snps 			> check_dup.snps
	
	echo `wc -l dups.snps | awk '{print $1}'`" duplicate markers in the dataset" \
										>> bim_after_update2.txt
	
	if [ "$(wc -l check_dup.snps | awk '{print $1}')" -gt 0 ] ; then 
		echo -e '\nStep ❷'
		echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
		${PLINK2PATH} --bfile ${nextdatasetname} \
						--allow-no-sex \
						--extract check_dup.snps \
						--cluster missing \
						--out ${nextdatasetname}_checkdup \
						--memory ${memory}
		echo -e '╚════════════════════════════════════════════════════════════════════════════╝'
		echo -e ''
		gzip ${nextdatasetname}_checkdup.mdist.missing
		python ${SCRIPTPATH}heatmap_from_matrix.py ${nextdatasetname}_checkdup.mdist.missing.gz
		
		echo -e '\nStep ❸'
		echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
		${PLINK2PATH} --bfile ${nextdatasetname} \
						--allow-no-sex \
						--extract check_dup.snps \
						--recode vcf-iid bgz \
						--out ${nextdatasetname}_checkdupvcf \
						--memory ${memory}
		echo -e '╚════════════════════════════════════════════════════════════════════════════╝'
		echo -e ''
		
		python ${SCRIPTPATH}vcf_flag_fix.py ${nextdatasetname}_checkdupvcf.vcf.gz > update.snps
		grep -v 'needs_flip\|different_alleles\|missing' update.snps \
										> update_chrpos.snps
		grep 'needs_flip' update.snps > flip.snps
		grep 'different_alleles' update.snps | awk '{print $1}' > exclude.snps
		awk '{print $1,$2}' update_chrpos.snps | grep -v -w -f exclude.snps > update_name.snps
		grep -wf exclude.snps ${nextdatasetname}.bim \
								| awk '{print $2,$2"_"$1"_"$4"_"$5"_"$6"_3AL"}' \
														>> update_name.snps
		sed -i 's/_FIXINFILE/_FIXDUP/g' update_name.snps
		awk '{print $1,$3}' update_chrpos.snps | grep -v -w -f exclude.snps > update_chr.snps
		awk '{print $1,$4}' update_chrpos.snps | grep -v -w -f exclude.snps > update_pos.snps
		
		echo -e '\nStep ❹'
		echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
		${PLINK2PATH} --bfile ${nextdatasetname} \
						--allow-no-sex \
						--extract dups.snps \
						--update-chr update_chr.snps \
						--update-map update_pos.snps \
						--flip flip.snps \
						--out ${nextdatasetname}_dups1 \
						--make-bed \
						--memory ${memory}
		echo -e '╚════════════════════════════════════════════════════════════════════════════╝'
		echo -e ''
		
		echo -e '\nStep ❺'
		echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
		${PLINK2PATH} --bfile ${nextdatasetname}_dups1 \
						--allow-no-sex \
						--update-name update_name.snps \
						--out ${nextdatasetname}_dups \
						--make-bed \
						--memory ${memory}
		echo -e '╚════════════════════════════════════════════════════════════════════════════╝'
		echo -e ''
		echo -e '\nStep ❻'
		echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
		${PLINK2PATH} --bfile ${nextdatasetname} \
						--allow-no-sex \
						--exclude dups.snps \
						--update-chr update_chr.snps \
						--update-map update_pos.snps \
						--flip flip.snps \
						--make-bed \
						--out ${nextdatasetname}_nodup1 \
						--memory ${memory}
		echo -e '╚════════════════════════════════════════════════════════════════════════════╝'
		echo -e ''
		
		echo -e '\nStep ❼'
		echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
		${PLINK2PATH} --bfile ${nextdatasetname}_nodup1 \
						--allow-no-sex \
						--update-name update_name.snps \
						--make-bed \
						--out ${nextdatasetname}_nodup \
						--memory ${memory}
		echo -e '╚════════════════════════════════════════════════════════════════════════════╝'
		echo -e ''
		
		echo -e '\nStep ❽'
		echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
		${PLINK2PATH} --bfile ${nextdatasetname}_nodup \
						--allow-no-sex \
						--bmerge ${nextdatasetname}_dups \
						--out ${nextdatasetname}_mrgdup \
						--make-bed \
						--memory ${memory}
		echo -e '╚════════════════════════════════════════════════════════════════════════════╝'
		echo -e ''
		
		echo -e '\nStep ❾'
		echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
		${PLINK2PATH} --bfile ${nextdatasetname}_nodup \
						--allow-no-sex \
						--bmerge ${nextdatasetname}_dups \
						--merge-mode 6 \
						--out ${nextdatasetname}_checkmrgdup \
						--memory ${memory}
		echo -e '╚════════════════════════════════════════════════════════════════════════════╝'
		echo -e ''
		
		nextdatasetname=${nextdatasetname}_mrgdup
	else
	nextdatasetname=${dataset##*/}_rs
	fi
fi








######    Explore and fix PAR regions and polymorphisms
if [ "${fixxy}" == "fixxy" ] ; then
	echo -e '\nStep ❿'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --bfile ${nextdatasetname} \
					--allow-no-sex \
					--merge-x no-fail \
					--make-bed \
					--out ${nextdatasetname}_noxy \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	echo -e '\nStep ❶❶'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --bfile ${nextdatasetname}_noxy \
					--allow-no-sex \
					--split-x hg19 no-fail \
					--make-bed \
					--out ${nextdatasetname}_splitxy \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	echo -e '\nStep ❶❷'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --bfile ${nextdatasetname}_splitxy \
					--allow-no-sex \
					--freq \
					--out ${nextdatasetname}_frq \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	awk '{print $3}' ${nextdatasetname}_frq.hh | sort | uniq -c | sort -k1 \
										> hh.snps
	awk '{if(($1==23)||($1==24)) print }' ${nextdatasetname}_splitxy.bim \
										> xy.bim
	(while read line; do 
			paste <(grep -w `echo $line | awk '{print $2}'` xy.bim) \
					<(echo $line | awk '{print $1}') ; done < hh.snps ) \
					| sort -k1,1 -gk4,4 > hh.snplist
	nextdatasetname=${nextdatasetname}_splitxy
fi








if [ "${fixindel}" == "fixindel" ] ; then
	####    Fix INDEL notation, if available
	zgrep -v '#' ${REFPATH}hg19/human_9606_b147_GRCh37p13/edits/All_20160601_indels.vcf.gz \
		| awk '{if(length($4)>1) print $3,"I D",$4,$5 ; \
				else if(length($5)>1) print $3,"D I",$4,$5}' \
		| grep -v ',' 					> indels.update
	
	echo -e '\nStep ❶❸'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --bfile ${nextdatasetname} \
					--allow-no-sex \
					--update-alleles indels.update \
					--make-bed \
					--out ${nextdatasetname}_updindel \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	gzip indels.update
	
	nextdatasetname=${nextdatasetname}_updindel
fi






######    When all is said and done, rename the last files to the original name,
######        but in this annotated folder
echo "Renaming files from "${nextdatasetname}" to "${dataset##*/}
mv ${nextdatasetname}.bed ${dataset##*/}.bed
mv ${nextdatasetname}.bim ${dataset##*/}.bim
mv ${nextdatasetname}.fam ${dataset##*/}.fam






echo -e '\n\n'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
echo -e '               Writing the report and the slides on LaTeX ✄ ✏\n'

echo '
\newpage
\section{Inital overview of \detokenize{'${dataset##*/}'}}

This section describes the original overview of the SNP composition of the dataset. This dataset was produced by genotyping the samples using the Illumina Infinium PsychArray. The first frame shows the SNP composition per name and per feature as produced by the original Illumina manifest (bpm file), the second frame shows the SNP composition per name and feature, after annotating with the updated Illumina annotation files (infinium-psycharray-24-v1-1-a1-b144-rsids).\newline \newline
'$( awk '{if($6==2) print $6}' ${dataset##*/}.fam | wc -l )' cases in the dataset.\newline
'$( awk '{if($6==1) print $6}' ${dataset##*/}.fam | wc -l )' controls in the dataset.\newline

\lstinputlisting[breaklines, 
				firstline=3, 
				frame=single, 
				title={Composition of the SNP content in the dataset 
						as delivered in the original file}
						]{'`readlink -e bim_before_update.txt`'}

\lstinputlisting[breaklines, 
				firstline=3, 
				frame=single, 
				title={Composition of the SNP content in the dataset 
				after using Illumina annotation and extracting the rs names}
				]{'`readlink -e bim_after_update2.txt`'}
' >> report.tex

echo '
\section{Inital overview of \detokenize{'${dataset##*/}'}}

\begin{frame}{Inital overview of \detokenize{'${dataset##*/}'} }
	\vspace{2em}
	\large{Initial annotation of the dataset}
\end{frame}

\begin{frame}{Inital overview of \detokenize{'${dataset##*/}'}}
	\lstinputlisting[breaklines, 
					firstline=3, 
					basicstyle=\tiny, 
					frame=single, 
					title={Composition of the SNP content in the dataset 
							as delivered in the original file}
							]{'`readlink -e bim_before_update.txt`'}
\end{frame}

\begin{frame}{Inital overview of \detokenize{'${dataset##*/}'}}
	\lstinputlisting[breaklines, 
					firstline=3, 
					basicstyle=\tiny, 
					frame=single, 
					title={Composition of the SNP content in the dataset 
					after using Illumina annotation and extracting the rs names}
					]{'`readlink -e bim_after_update2.txt`'}
\end{frame}
' >> slides.tex


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
