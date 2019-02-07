#!/bin/sh -l
source /depot/pdrineas/data/gwas_scripts/software_paths.conf

######---------------------------------------------------------
######    Last updated: 2017 Jul 17
######---------------------------------------------------------



######---------------------------------------------------------
######    Input arguments
data1=`readlink -e ${1}.bed | sed 's/\.bed//'`
data2=`readlink -e ${2}.bed | sed 's/\.bed//'`
name1=${3}
name2=${4}
samp1=`readlink -e ${5}`
samp2=`readlink -e ${6}`
outl3=`readlink -e ${7}`
memory=${8}
######---------------------------------------------------------



######---------------------------------------------------------
######    Keep a log
date >> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
echo 'Running merge_qc.sh with the following arguments:' \
					>> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
for i; do
	echo -e '\t'$i >> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
done
echo -e '\n' >> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
######---------------------------------------------------------

######---------------------------------------------------------
######     Start analysis
######---------------------------------------------------------

mkdir ${name1}_${name2}_merge
cd ${name1}_${name2}_merge
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



##### This checks for overlapping SNPs, problematic SNPs and SNPs that need flipping
python ${SCRIPTPATH}check-bims.py ${data1}.bim ${data2}.bim | tee bims_overlap.log

#####       Prepare the first dataset
echo -e '\n\nStep ❶'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --allow-no-sex \
				--bfile ${data1} \
				--extract ${data1##*/}.${data2##*/}-common.snps \
				--exclude ${data1##*/}.${data2##*/}-problem.snps \
				--make-bed \
				--keep ${samp1} \
				--out ${name1}_${name2} \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
#####       Prepare the second dataset
echo -e '\nStep ❷'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --allow-no-sex \
				--bfile ${data2} \
				--extract ${data1##*/}.${data2##*/}-common.snps \
				--exclude ${data1##*/}.${data2##*/}-problem.snps \
				--make-bed \
				--flip ${data1##*/}.${data2##*/}-needflip.snps \
				--keep ${samp2} \
				--out ${name2}_${name1} \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
#####       Merge
echo -e '\nStep ❸'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --allow-no-sex \
				--bfile ${name1}_${name2} \
				--bmerge ${name2}_${name1} \
				--make-bed \
				--out ${name1}${name2} \
				--remove ${outl3} \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

#####       Find the ambiguous SNPs
python ${SCRIPTPATH}bim-analyser.py ${name1}${name2}.bim -ambiguous > ${name1}${name2}.bim.analysis.log
echo -e '\n'


#########     BATCH EFFECT TEST, SPLIT IN CASES AND CONTROLS
echo -e '\nStep ❹'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --allow-no-sex \
				--bfile ${name1}${name2} \
				--filter-cases \
				--make-bed \
				--out ${name1}${name2}_cases \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

echo -e '\nStep ❺ '
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --allow-no-sex \
				--bfile ${name1}${name2} \
				--filter-controls \
				--make-bed \
				--out ${name1}${name2}_controls \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

#####       And for each subset, test association, HWE and missingness by batch

for pheno in cases controls; do
	echo -e '\nStep ❻ '
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --allow-no-sex \
					--bfile ${name1}${name2}_${pheno} \
					--make-pheno ${data1}.fam '*' \
					--assoc \
					--hardy \
					--test-missing \
					--out batch_${pheno}_assoc \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

	awk '{print $0,$3 - $4}' batch_${pheno}_assoc.missing \
										> batch_${pheno}_assoc.missing.diff

	awk '{if((($4=="A")&&($7=="T")) \
			||(($4=="T")&&($7=="A")) \
			||(($4=="G")&&($7=="C")) \
			||(($4=="C")&&($7=="G"))) print}' batch_${pheno}_assoc.assoc \
									> batch_${pheno}_assoc.assoc.ambiguous
	
	grep -F -x -v -f batch_${pheno}_assoc.assoc.ambiguous batch_${pheno}_assoc.assoc \
									> batch_${pheno}_assoc.assoc.nonambiguous
	
	awk '{if($5<0.00001) print $2,$0}' batch_${pheno}_assoc.missing.diff \
									> batch_${pheno}_assoc.missing.diff.remove
	
	awk '{if($9<0.00001) print $2,$0}' batch_${pheno}_assoc.assoc.nonambiguous \
									> batch_${pheno}_assoc.assoc.nonambiguous.remove
	
	awk '{if((0.6>$5&&$5>0.4) \
			|| (0.6>$6&&$6>0.4))print $2}' batch_${pheno}_assoc.assoc.ambiguous \
									> batch_${pheno}_assoc.assoc.ambiguous.remove
	
	awk '{if($9<1e-20) print $2,$0}' batch_${pheno}_assoc.assoc.ambiguous \
									> batch_${pheno}_assoc.assoc.ambiguous.flip
	
	awk '{if($9!=NA) print $1,-log($9)/log(10)}' batch_${pheno}_assoc.assoc \
						| grep -v inf > batch_${pheno}_assoc_log.assoc
	awk '{if($5!=NA) print $1,-log($5)/log(10)}' batch_${pheno}_assoc.missing \
						| grep -v inf > batch_${pheno}_assoc_log.missing
	
	python ${SCRIPTPATH}histogram.py batch_${pheno}_assoc_log.assoc \
									-col=1 \
									-header \
									-binsize=0.1 \
									-vline=5 \
									-log
	python ${SCRIPTPATH}histogram.py batch_${pheno}_assoc_log.missing \
									-col=1 \
									-header \
									-binsize=0.1 \
									-vline=5 \
									-log
done

cat ${data1##*/}.${data2##*/}-needflip.snps \
		<(awk '{print $1}' batch_*_assoc.assoc.ambiguous.flip) \
									> all_for_flip.snps

cat batch_*.remove > to_remove.snps

echo -e '\nStep ❼ '
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --allow-no-sex \
				--bfile ${data2} \
				--extract ${data1##*/}.${data2##*/}-common.snps \
				--exclude ${data1##*/}.${data2##*/}-problem.snps \
				--make-bed \
				--flip all_for_flip.snps \
				--keep ${samp2} \
				--out ${name2}_${name1}_ambflip \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
#####       Merge
echo -e '\nStep ❽ '
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --allow-no-sex \
				--bfile ${name1}_${name2} \
				--bmerge ${name2}_${name1}_ambflip  \
				--make-bed \
				--out ${name1}${name2}_ambflip \
				--remove ${outl3} \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
echo -e '\nStep ❾ '
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --allow-no-sex \
				--bfile ${name1}${name2}_ambflip \
				--exclude to_remove.snps \
				--out ${name1}${name2}_batcheff \
				--make-bed \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

######      Remove ambiguous
echo -e '\nStep ❿ '
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --allow-no-sex \
				--bfile ${name1}${name2}_batcheff \
				--exclude ${name1}${name2}.bim.ambiguous \
				--make-bed \
				--out ${name1}${name2}_batcheff_noamb \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'


###################
###          Gzip files for less space
###################
gzip -f *.hwe &
gzip -f *.assoc &
gzip -f *.nonambiguous &
gzip -f *.diff &
gzip -f *.remove &
gzip -f *.snps &
gzip -f *.missing &
wait



echo '
\newpage
\section{Merging overview}

In this section we have an overview of the datasets merged.

\subsection{Datasets merged}

\begin{table}[H]
	\caption{Overview of datasets merged}
	\centering
	\begin{tabular}{lllllll}\toprule
		Dataset&Cases&Controls&SNPs&Overlap&Problematic&Flipped \\ \midrule
		'${name1}'&'`awk '{if($6==2) print}' ${data1}.fam | wc -l`'&'`awk '{if($6==1) print}' ${data1}.fam | wc -l`'&'`wc -l ${data1}.bim | awk '{print $1}'`'&'`zcat ${data1##*/}.${data2##*/}-common.snps.gz | wc -l`'&'`zcat ${data1##*/}.${data2##*/}-problem.snps.gz | wc -l`'&0 \\
		'${name2}'&'`awk '{if($6==2) print}' ${data2}.fam | wc -l`'&'`awk '{if($6==1) print}' ${data2}.fam | wc -l`'&'`wc -l ${data2}.bim | awk '{print $1}'`'&'`zcat ${data1##*/}.${data2##*/}-common.snps.gz | wc -l`'&'`zcat ${data1##*/}.${data2##*/}-problem.snps.gz | wc -l`'&'`zcat ${data1##*/}.${data2##*/}-needflip.snps.gz | wc -l`'\\
		\bottomrule
	\end{tabular}
\end{table}

\subsection{Batch effect}

For the batch effects, we split the datasets into control-only and case-only subsets and then perform association using the dataset membership as the phenotype. For this part, we also perform a differential missingness test. The p-value cut-off for these tests is $10^{-5}$.

\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-batch_cases_assoc_log.assoc_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of the association test p-values in cases}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-batch_cases_assoc_log.missing_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of the differential missingess test p-values in cases}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-batch_controls_assoc_log.assoc_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of the association test p-values in controls}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-batch_controls_assoc_log.missing_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of the differential missingess test p-values in controls}
\end{figure}

\begin{table}[H]
	\caption{SNPs removed from the batch effect tests}
	\centering
	\begin{tabular}{lllll}\toprule
		Subset&SNPs \\ \midrule
		Failing batch assoc in cases&'`zcat batch_cases_assoc.assoc.nonambiguous.remove.gz | wc -l`'\\
		Failing missingness in cases&'`zcat batch_cases_assoc.missing.diff.remove.gz | wc -l`'\\
		Failing batch assoc in controls&'`zcat batch_controls_assoc.assoc.nonambiguous.remove.gz | wc -l`'\\
		Failing missingness in controls&'`zcat batch_controls_assoc.missing.diff.remove.gz | wc -l`'\\
		Ambiguous&'`wc -l ${name1}${name2}.bim.ambiguous | awk '{print $1}'`' \\
		\bottomrule
	\end{tabular}
\end{table}

' >> report.tex


echo '
\section{Merging overview}

\begin{frame}
\begin{table}[H]
	\caption{Overview of datasets merged}
	\centering
	\begin{tabular}{lllll}\toprule
		Dataset&Cases&Controls&SNPs&Overlap \\ \midrule
		'${name1}'&'`awk '{if($6==2) print}' ${data1}.fam | wc -l`'&'`awk '{if($6==1) print}' ${data1}.fam | wc -l`'&'`wc -l ${data1}.bim | awk '{print $1}'`'&'`zcat ${data1##*/}.${data2##*/}-common.snps | wc -l`'\\
		'${name2}'&'`awk '{if($6==2) print}' ${data2}.fam | wc -l`'&'`awk '{if($6==1) print}' ${data2}.fam | wc -l`'&'`wc -l ${data2}.bim | awk '{print $1}'`'&'`zcat ${data1##*/}.${data2##*/}-common.snps| wc -l`'\\
		\bottomrule
	\end{tabular}
\end{table}
\end{frame}

\begin{frame}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-batch_cases_assoc_log.assoc_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of the association test p-values in cases}
\end{figure}
\end{frame}

\begin{frame}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-batch_cases_assoc_log.missing_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of the differential missingess test p-values in cases}
\end{figure}
\end{frame}

\begin{frame}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-batch_controls_assoc_log.assoc_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of the association test p-values in controls}
\end{figure}
\end{frame}

\begin{frame}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-batch_controls_assoc_log.missing_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of the differential missingess test p-values in controls}
\end{figure}
\end{frame}

\begin{frame}
\begin{table}[H]
	\caption{SNPs removed from the batch effect tests}
	\centering
	\begin{tabular}{lllll}\toprule
		Subset&SNPs \\ \midrule
		Failing batch assoc in cases&'`zcat batch_cases_assoc.assoc.nonambiguous.remove.gz | wc -l`'\\
		Failing missingness in cases&'`zcat batch_cases_assoc.missing.diff.remove.gz | wc -l`'\\
		Failing batch assoc in controls&'`zcat batch_controls_assoc.assoc.nonambiguous.remove.gz | wc -l`'\\
		Failing missingness in controls&'`zcat batch_controls_assoc.missing.diff.remove.gz | wc -l`'\\
		Ambiguous&'`wc -l ${name1}${name2}.bim.ambiguous | awk '{print $1}'`' \\
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
