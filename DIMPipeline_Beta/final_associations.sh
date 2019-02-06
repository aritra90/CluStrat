#!/bin/sh -l
source /depot/pdrineas/data/gwas_scripts/software_paths.conf

######---------------------------------------------------------
######    Last updated: 2017 Jul 17
######---------------------------------------------------------



######---------------------------------------------------------
######    Input arguments
dataset=`readlink -e $1.bed | sed 's/\.bed//'`
covarpath=`readlink -e $2`
covars=$3
mode=$4
outliers=`readlink -e $5`
threads=$6
memory=$7

pvalThresh=0.05


float_re='^[+-]?[0-9]+([.][0-9]+)?$'
int_re='^[0-9]+$'

while (( "$#" )); do
  case "$1" in
    --pvalThresh)
        #echo "$2"
        pvalThresh=$2
        if ! [[ $pvalThresh =~ $float_re ]]
        then
            usage
            echo "Not a float!!"
            exit 1
        fi
        if (( $(echo "$pvalThresh > 1" |bc -l) ||  $(echo "$pvalThresh < 0" |bc -l)  )); then
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
######---------------------------------------------------------



######---------------------------------------------------------
######    Keep a log
date >> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
echo 'Running final_associations.sh with the following arguments:' \
					>> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
for i; do
	echo -e '\t'$i >> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
done
echo -e '\n' >> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
######---------------------------------------------------------



######---------------------------------------------------------
######     Start analysis
######---------------------------------------------------------

mkdir ${dataset##*/}_assoc_$(echo ${outliers##*/} | cut -d'.' -f1)
cd ${dataset##*/}_assoc_$(echo ${outliers##*/} | cut -d'.' -f1)
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

\begin{document}' > report.tex
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

\begin{document}' > slides.tex



echo -e '\n'
echo -e '🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔'
echo -e '🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔'
echo -e '🙒🙜 🙒🙜 🙒🙜                                                             🙞🙔 🙞🙔 🙞🙔'
echo -e '🙒🙜 🙒🙜                                                                   🙞🙔 🙞🙔'
echo -e '🙒🙜                                                                         🙞🙔'




###############################################
####    First a simple association using PLINK
###############################################

if [ "$mode" == "simple" ] || [ "$mode" == "simplelogistic" ] || [ "$mode" == "full" ] ; then
	echo -e '\nStep ❶'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --bfile ${dataset} \
					--allow-no-sex \
					--assoc \
					--remove ${outliers} \
					--out ${dataset##*/}_assoc \
					--threads ${threads} \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	
	######    Get top results ( p <= 1e-3 )
	cat <(head -1 ${dataset##*/}_assoc.assoc) \
		<(sort -gk9 ${dataset##*/}_assoc.assoc | grep -v NA \
							| awk -v var="$pvalThresh" '{if($9<=var) print}') \
			> ${dataset##*/}_assoc_top.assoc
	
	echo -e '\nStep ❶'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --annotate ${dataset##*/}_assoc_top.assoc \
								ranges=${REFPATH}plink/glist-hg19 \
					--border 20 \
					--out ${dataset##*/}_assoc_top \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	
	echo -e '\nStep ❶'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --clump ${dataset##*/}_assoc.assoc \
					--bfile ${dataset} \
					--remove ${outliers} \
					--clump-range ${REFPATH}plink/glist-hg19 \
					--clump-range-border 20 \
					--out ${dataset##*/}_assoc_clump \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	
	######    Plot and zip
	python ${SCRIPTPATH}manhattan_plot.py ${dataset##*/}_assoc.assoc 8,0,2,1 0.05 1
	python ${SCRIPTPATH}qq-plot.py ${dataset##*/}_assoc.assoc 8 1 \
							$(awk '{if($6==2) print}' ${dataset}.fam | wc -l),$(awk '{if($6==1) print}' ${dataset}.fam | wc -l)
	gzip -f ${dataset##*/}_assoc.assoc
	
	echo '
	\newpage
	\section{Associations of \detokenize{'${dataset##*/}'}}
	
	In this section we describe the results of the associations as performed on the clean data.
	
	\subsection{Simple Association}
	Associations without using regression on PCs.
	\begin{figure}[H]
		\centering
		\includegraphics[width=\textwidth,keepaspectratio]{{'`readlink -e \
						manhattan-plot-${dataset##*/}_assoc.assoc.png | sed 's/\.png//' `'}.png}
		\caption{Manhattan plot without PC correction}
	\end{figure}

	\begin{figure}[H]
		\centering
		\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e \
								qq-plot-${dataset##*/}_assoc.assoc.png | sed 's/\.png//' `'}.png}
		\caption{QQ-plot without PC correction}
	\end{figure}
	
	\begin{table}[H]
		\caption{Top annotated results without PC correction}
		\centering
		\begin{scriptsize}
		\begin{tabularx}{\textwidth}{lLLLLLLl}\toprule
			'`awk '{print $2,$1,$3,$4,$5,$6,$9,$10,$11}' ${dataset##*/}_assoc_top.annot \
					| head -20 \
					| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
			\bottomrule
		\end{tabularx}
		\end{scriptsize}
	\end{table}
	
	\begin{table}[H]
		\caption{Top annotated clumped results without PC correction}
		\centering
		\begin{scriptsize}
		\begin{tabularx}{\textwidth}{lLLLlLl}\toprule
			'`awk '{print $2,$1,$3,$4,$5,$6,$7}' ${dataset##*/}_assoc_clump.clumped.ranges \
					| head -20 \
					| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
			\bottomrule
		\end{tabularx}
		\end{scriptsize}
	\end{table}
	
	' >> report.tex
	
	echo '
	\section{Associations of \detokenize{'${dataset##*/}'}}
	
	\begin{frame}{Association of \detokenize{'${dataset##*/}'}}
	\begin{figure}[H]
		\centering
		\includegraphics[width=\textwidth,keepaspectratio]{{'`readlink -e \
						manhattan-plot-${dataset##*/}_assoc.assoc.png | sed 's/\.png//' `'}.png}
		\caption{Manhattan plot without PC correction}
	\end{figure}
	\end{frame}
	
	\begin{frame}{Association of \detokenize{'${dataset##*/}'}}
	\begin{figure}[H]
		\centering
		\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e \
								qq-plot-${dataset##*/}_assoc.assoc.png | sed 's/\.png//' `'}.png}
		\caption{QQ-plot without PC correction}
	\end{figure}
	\end{frame}
	
	\begin{frame}{Association of \detokenize{'${dataset##*/}'}}
	\begin{table}[H]
		\caption{Top annotated results without PC correction}
		\centering
		\begin{tiny}
		\begin{tabularx}{\textwidth}{lLLLLLLl}\toprule
			'`awk '{print $2,$1,$3,$4,$5,$6,$9,$10,$11}' ${dataset##*/}_assoc_top.annot \
					| head -20 \
					| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
			\bottomrule
		\end{tabularx}
		\end{tiny}
	\end{table}
	\end{frame}
	
	\begin{frame}{Association of \detokenize{'${dataset##*/}'}}
	\begin{table}[H]
		\caption{Top annotated clumped results without PC correction}
		\centering
		\begin{scriptsize}
		\begin{tabularx}{\textwidth}{lLLLlLl}\toprule
			'`awk '{print $2,$1,$3,$4,$5,$6,$7}' ${dataset##*/}_assoc_clump.clumped.ranges \
					| head -20 \
					| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
			\bottomrule
		\end{tabularx}
		\end{scriptsize}
	\end{table}
	\end{frame}
	' >> slides.tex
fi


###############################################
####    Logistic Regression using PLINK
###############################################

if [ "$mode" == "logistic" ] || [ "$mode" == "simplelogistic" ] || [ "$mode" == "full" ] ; then
	
	######    Association using the first PC    ######
	echo -e '\nStep ❷'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --bfile ${dataset} \
					--allow-no-sex \
					--logistic \
					--remove ${outliers} \
					--covar ${covarpath} \
					--out ${dataset##*/}_logistic_1 \
					--covar-number 1 \
					--hide-covar \
					--threads ${threads} \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	
	######    Get top results ( p <= 1e-3 )
	cat <(head -1 ${dataset##*/}_logistic_1.assoc.logistic) \
		<(sort -gk9 ${dataset##*/}_logistic_1.assoc.logistic | grep -v NA \
							| awk '{if($9<=1e-3) print}') \
			> ${dataset##*/}_logistic_1_top.assoc.logistic
	######
	
	echo -e '\nStep ❷'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --annotate ${dataset##*/}_logistic_1_top.assoc.logistic \
								ranges=${REFPATH}plink/glist-hg19 \
					--border 20 \
					--out ${dataset##*/}_logistic_1_top \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	
	echo -e '\nStep ❷'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --clump ${dataset##*/}_logistic_1.assoc.logistic \
					--bfile ${dataset} \
					--remove ${outliers} \
					--clump-range ${REFPATH}plink/glist-hg19 \
					--clump-range-border 20 \
					--out ${dataset##*/}_logistic_1_clump \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	
	######    Plot and zip
	python ${SCRIPTPATH}manhattan_plot.py ${dataset##*/}_logistic_1.assoc.logistic 8,0,2,1 0.05 1
	python ${SCRIPTPATH}qq-plot.py ${dataset##*/}_logistic_1.assoc.logistic 8 1 \
							$(awk '{if($6==2) print}' ${dataset}.fam | wc -l),$(awk '{if($6==1) print}' ${dataset}.fam | wc -l)
	gzip -f ${dataset##*/}_logistic_1.assoc.logistic
	######
	
	echo '
	\newpage
	\subsection{Associations using logistic regression}
	In this section we present the results after running logistic regression on the PCs generated by the PCA. We add each PC in every run.
	\begin{figure}[H]
		\centering
		\includegraphics[width=\textwidth,keepaspectratio]{{'`readlink -e \
						manhattan-plot-${dataset##*/}_logistic_1.assoc.logistic.png \
															| sed 's/\.png//' `'}.png}
		\caption{Manhattan plot correcting for the 1st PC}
	\end{figure}

	\begin{figure}[H]
		\centering
		\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e \
				qq-plot-${dataset##*/}_logistic_1.assoc.logistic.png | sed 's/\.png//' `'}.png}
		\caption{QQ-plot correcting for the 1st PC}
	\end{figure}
	
	\begin{table}[H]
		\caption{Top annotated results correcting for the 1st PC}
		\centering
		\begin{scriptsize}
		\begin{tabularx}{\textwidth}{lLLLLLLl}\toprule
			'`awk '{print $2,$1,$3,$4,$7,$9,$10}' ${dataset##*/}_logistic_1_top.annot \
					| head -20 \
					| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
			\bottomrule
		\end{tabularx}
		\end{scriptsize}
	\end{table}
	
	\begin{table}[H]
		\caption{Top annotated clumped results correcting for the 1st PC}
		\centering
		\begin{scriptsize}
		\begin{tabularx}{\textwidth}{lLLLlLl}\toprule
			'`awk '{print $2,$1,$3,$4,$5,$6,$7}' ${dataset##*/}_logistic_1_clump.clumped.ranges \
					| head -20 \
					| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
			\bottomrule
		\end{tabularx}
		\end{scriptsize}
	\end{table}
	
	' >> report.tex
	
	echo '
	\begin{frame}{Association of \detokenize{'${dataset##*/}'}}
	\begin{figure}[H]
		\centering
		\includegraphics[width=\textwidth,keepaspectratio]{{'`readlink -e \
						manhattan-plot-${dataset##*/}_logistic_1.assoc.logistic.png \
															| sed 's/\.png//' `'}.png}
		\caption{Manhattan plot correcting for the 1st PC}
	\end{figure}
	\end{frame}
	
	\begin{frame}{Association of \detokenize{'${dataset##*/}'}}
	\begin{figure}[H]
		\centering
		\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e \
				qq-plot-${dataset##*/}_logistic_1.assoc.logistic.png | sed 's/\.png//' `'}.png}
		\caption{QQ-plot correcting for the 1st PC}
	\end{figure}
	\end{frame}
	
	\begin{frame}{Association of \detokenize{'${dataset##*/}'}}
	\begin{table}[H]
		\caption{Top annotated results correcting for the 1st PC}
		\centering
		\begin{tiny}
		\begin{tabularx}{\textwidth}{lLLLLLLl}\toprule
			'`awk '{print $2,$1,$3,$4,$7,$9,$10}' ${dataset##*/}_logistic_1_top.annot \
					| head -20 \
					| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
			\bottomrule
		\end{tabularx}
		\end{tiny}
	\end{table}
	\end{frame}
	' >> slides.tex
	
	######    Association using more PCs, if they are to be used    ######
	if [ "${covars}" -gt 1 ]; then
	
	
		for covar in `seq 2 ${covars}`; do
			echo -e '\nStep ❷'
			echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
			${PLINK2PATH} --bfile ${dataset} \
							--allow-no-sex \
							--remove ${outliers} \
							--logistic \
							--covar ${covarpath} \
							--out ${dataset##*/}_logistic_1-${covar} \
							--covar-number 1-${covar} \
							--hide-covar \
							--threads ${threads} \
							--memory ${memory}
			echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
		
			######    Get top results ( p <= 1e-3 )
			cat <(head -1 ${dataset##*/}_logistic_1-${covar}.assoc.logistic) \
				<(sort -gk9 ${dataset##*/}_logistic_1-${covar}.assoc.logistic | grep -v NA \
									| awk '{if($9<=1e-3) print}') \
					> ${dataset##*/}_logistic_1-${covar}_top.assoc.logistic
			
			echo -e '\nStep ❷'
			echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
			${PLINK2PATH} --annotate ${dataset##*/}_logistic_1-${covar}_top.assoc.logistic \
										ranges=${REFPATH}plink/glist-hg19 \
							--border 20 \
							--out ${dataset##*/}_logistic_1-${covar}_top \
							--memory ${memory}
			echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
			echo -e '\nStep ❷'
			echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
			${PLINK2PATH} --clump ${dataset##*/}_logistic_1-${covar}.assoc.logistic \
							--bfile ${dataset} \
							--remove ${outliers} \
							--clump-range ${REFPATH}plink/glist-hg19 \
							--clump-range-border 20 \
							--out ${dataset##*/}_logistic_1-${covar}_clump \
							--memory ${memory}
			echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
			
			
			######    Plot and zip
			python ${SCRIPTPATH}manhattan_plot.py \
							${dataset##*/}_logistic_1-${covar}.assoc.logistic 8,0,2,1 0.05 1
			python ${SCRIPTPATH}qq-plot.py ${dataset##*/}_logistic_1-${covar}.assoc.logistic 8 1 \
							$(awk '{if($6==2) print}' ${dataset}.fam | wc -l),$(awk '{if($6==1) print}' ${dataset}.fam | wc -l)
			gzip -f ${dataset##*/}_logistic_1-${covar}.assoc.logistic
			
			echo '
			\newpage
			\begin{figure}[H]
				\centering
				\includegraphics[width=\textwidth,keepaspectratio]{{'`readlink -e \
						manhattan-plot-${dataset##*/}_logistic_1-${covar}.assoc.logistic.png \
																| sed 's/\.png//' `'}.png}
				\caption{Manhattan plot correcting for PCs 1-'${covar}'}
			\end{figure}

			\begin{figure}[H]
				\centering
				\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e \
						qq-plot-${dataset##*/}_logistic_1-${covar}.assoc.logistic.png \
																| sed 's/\.png//' `'}.png}
				\caption{QQ-plot correcting for PCs 1-'${covar}'}
			\end{figure}
			
			\begin{table}[H]
				\caption{Top annotated results correcting for PCs 1-'${covar}'}
				\centering
				\begin{scriptsize}
				\begin{tabularx}{\textwidth}{lLLLLLLl}\toprule
					'`awk '{print $2,$1,$3,$4,$7,$9,$10}' \
							${dataset##*/}_logistic_1-${covar}_top.annot \
							| head -20 \
							| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
					\bottomrule
				\end{tabularx}
				\end{scriptsize}
			\end{table}
			
			\begin{table}[H]
				\caption{Top annotated clumped results correcting for PCs 1-'${covar}'}
				\centering
				\begin{scriptsize}
				\begin{tabularx}{\textwidth}{lLLLlLl}\toprule
					'`awk '{print $2,$1,$3,$4,$5,$6,$7}' ${dataset##*/}_logistic_1-${covar}_clump.clumped.ranges \
							| head -20 \
							| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
					\bottomrule
				\end{tabularx}
				\end{scriptsize}
			\end{table}
			
			' >> report.tex
			
			echo '
			\begin{frame}{Association of \detokenize{'${dataset##*/}'}}
			\begin{figure}[H]
				\centering
				\includegraphics[width=\textwidth,keepaspectratio]{{'`readlink -e \
					manhattan-plot-${dataset##*/}_logistic_1-${covar}.assoc.logistic.png \
																| sed 's/\.png//' `'}.png}
				\caption{Manhattan plot correcting for PCs 1-'${covar}'}
			\end{figure}
			\end{frame}
			
			\begin{frame}{Association of \detokenize{'${dataset##*/}'}}
			\begin{figure}[H]
				\centering
				\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e \
					qq-plot-${dataset##*/}_logistic_1-${covar}.assoc.logistic.png \
																| sed 's/\.png//' `'}.png}
				\caption{QQ-plot correcting for PCs 1-'${covar}'}
			\end{figure}
			\end{frame}
			
			\begin{frame}{Association of \detokenize{'${dataset##*/}'}}
			\begin{table}[H]
				\caption{Top annotated results correcting for PCs 1-'${covar}'}
				\centering
				\begin{tiny}
				\begin{tabularx}{\textwidth}{lLLLLLLl}\toprule
					'`awk '{print $2,$1,$3,$4,$7,$9,$10}' \
							${dataset##*/}_logistic_1-${covar}_top.annot \
							| head -20 \
							| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
					\bottomrule
				\end{tabularx}
				\end{tiny}
			\end{table}
			\end{frame}
			
			' >> slides.tex
		done
		
		
	fi
fi

###  Using Eigensoft's proposed covariates
if [ "$mode" == "eig" ] || [ "$mode" == "full" ] ; then
	
	covar=$(grep Control_Case $(dirname ${covarpath})/pca.pcalog | sort -gk2,2 \
			| awk '{if($2<0.001) print}'| awk -F"_" '{print $2}' |xargs | tr ' ' ',')
	
	if [ -n "$covar" ]; then
		echo -e '\nStep ❷'
		echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
		${PLINK2PATH} --bfile ${dataset} \
						--allow-no-sex \
						--remove ${outliers} \
						--logistic \
						--covar ${covarpath} \
						--out ${dataset##*/}_eig_logistic_${covar} \
						--covar-number ${covar} \
						--hide-covar \
						--memory ${memory}
		echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
		
		######    Get top results ( p <= 1e-3 )
		cat <(head -1 ${dataset##*/}_eig_logistic_${covar}.assoc.logistic) \
			<(sort -gk9 ${dataset##*/}_eig_logistic_${covar}.assoc.logistic | grep -v NA \
								| awk '{if($9<=1e-3) print}') \
				> ${dataset##*/}_eig_logistic_${covar}_top.assoc.logistic
		
		echo -e '\nStep ❷'
		echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
		${PLINK2PATH} --annotate ${dataset##*/}_eig_logistic_${covar}_top.assoc.logistic \
									ranges=${REFPATH}plink/glist-hg19 \
						--border 20 \
						--out ${dataset##*/}_eig_logistic_${covar}_top \
						--memory ${memory}
		echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
		mv ${dataset##*/}_eig_logistic_${covar}.assoc.logistic \
			${dataset##*/}_eig_logistic_$(echo ${covar} | sed 's/,/-/g' ).assoc.logistic
		echo -e '\nStep ❷'
		echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
		${PLINK2PATH} --clump ${dataset##*/}_eig_logistic_$(echo ${covar} | sed 's/,/-/g' ).assoc.logistic \
						--bfile ${dataset} \
						--remove ${outliers} \
						--clump-range ${REFPATH}plink/glist-hg19 \
						--clump-range-border 20 \
						--out ${dataset##*/}_eig_logistic_${covar}_clump \
						--memory ${memory}
		echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
		mv ${dataset##*/}_eig_logistic_$(echo ${covar} | sed 's/,/-/g' ).assoc.logistic \
					${dataset##*/}_eig_logistic_${covar}.assoc.logistic
		
		######    Plot and zip
		python ${SCRIPTPATH}manhattan_plot.py \
						${dataset##*/}_eig_logistic_${covar}.assoc.logistic 8,0,2,1 0.05 1
		python ${SCRIPTPATH}qq-plot.py ${dataset##*/}_eig_logistic_${covar}.assoc.logistic 8 1 \
						$(awk '{if($6==2) print}' ${dataset}.fam | wc -l),$(awk '{if($6==1) print}' ${dataset}.fam | wc -l)
		gzip -f ${dataset##*/}_eig_logistic_${covar}.assoc.logistic
		
		echo '
		\newpage
		\subsection{By EIGENSOFTs proposed covariates}
		In this section we use the covariates that have been proposed for EIGENSOFT to correct for. These covariates have been identified by an ANOVA on the dataset, and have a p-value of less than 0.001. \\
		
		Covariates (PCs produced by EIGENSOFT): '`echo ${covar}`' (in order of significance) \\
		
		\begin{figure}[H]
			\centering
			\includegraphics[width=\textwidth,keepaspectratio]{{'`readlink -e \
						manhattan-plot-${dataset##*/}_eig_logistic_${covar}.assoc.logistic.png | sed 's/\.png//' `'}.png}
			\caption{Manhattan plot using EIGENSOFTs proposed covariates}
		\end{figure}

		\begin{figure}[H]
			\centering
			\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e \
								qq-plot-${dataset##*/}_eig_logistic_${covar}.assoc.logistic.png | sed 's/\.png//' `'}.png}
			\caption{QQ-plot using EIGENSOFTs proposed covariates}
		\end{figure}
	
		\begin{table}[H]
			\caption{Top annotated results using EIGENSOFTs proposed covariates}
			\centering
			\begin{scriptsize}
			\begin{tabularx}{\textwidth}{lLLLLLLl}\toprule
				'`awk '{print $2,$1,$3,$4,$7,$9,$10}' \
						${dataset##*/}_eig_logistic_${covar}_top.annot \
						| head -20 \
						| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
				\bottomrule
			\end{tabularx}
			\end{scriptsize}
		\end{table}
		
		\begin{table}[H]
			\caption{Top annotated clumped results using EIGENSOFTs proposed covariates}
			\centering
			\begin{scriptsize}
			\begin{tabularx}{\textwidth}{lLLLlLl}\toprule
				'`awk '{print $2,$1,$3,$4,$5,$6,$7}' ${dataset##*/}_eig_logistic_${covar}_clump.clumped.ranges \
						| head -20 \
						| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
				\bottomrule
			\end{tabularx}
			\end{scriptsize}
		\end{table}
		' >> report.tex
	fi
fi

#######################################################
####    Also do an MLM association analysis using GCTA
#######################################################

if [ "$mode" == "mlma" ] || [ "$mode" == "fll" ] ; then

######### Removed this since it's essentially the same as the next part
###	${GCTAPATH} --bfile ${dataset} \
###				--autosome \
###				--maf 0.01 \
###				--remove ${outliers} \
###				--make-grm \
###				--out ${dataset##*/}_grm \
###				--thread-num ${threads}
###			
###	${GCTAPATH} --mlma \
###				--autosome \
###				--bfile ${dataset} \
###				--remove ${outliers} \
###				--grm ${dataset##*/}_grm \
###				--pheno ${dataset}.fam \
###				--mpheno 4 \
###				--out ${dataset##*/}_gcta_mlma \
###				--thread-num ${threads}
###				
###				
###	######    Plot and zip
###	python ${SCRIPTPATH}manhattan_plot.py \
###					${dataset##*/}_gcta_mlma.mlma 8,0,2,1 0.05 1
###	python ${SCRIPTPATH}qq-plot.py ${dataset##*/}_gcta_mlma.mlma 8 1 \
###							$(awk '{if($6==2) print}' ${dataset}.fam | wc -l),$(awk '{if($6==1) print}' ${dataset}.fam | wc -l)
###	gzip -f ${dataset##*/}_gcta_mlma.mlma
	
	
	
	
	####
	#### MLM based association analysis including the candidate SNP (MLMi)
	####
	
	${GCTAPATH} --mlma \
				--autosome \
				--bfile ${dataset} \
				--remove ${outliers} \
				--pheno ${dataset}.fam \
				--mpheno 4 \
				--out ${dataset##*/}_gcta_mlmi \
				--thread-num ${threads}
	
	sed -i -e '1!b;s/Chr/CHR/' \
		-e '1!b;s/bp/BP/' \
		-e '1!b;s/p/P/' \
		-e '1!b;s/se/SE/' \
		-e '1!b;s/b/STAT/' ${dataset##*/}_gcta_mlmi.mlma
	
	######    Get top results ( p <= 1e-3 )
	cat <(head -1 ${dataset##*/}_gcta_mlmi.mlma ) \
		<(sort -gk9 ${dataset##*/}_gcta_mlmi.mlma | grep -v nan \
							| awk '{if($9<=1e-3) print}') \
			> ${dataset##*/}_gcta_mlmi_top.mlma
	
	echo -e '\nStep ❸'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --annotate ${dataset##*/}_gcta_mlmi_top.mlma \
								ranges=${REFPATH}plink/glist-hg19 \
					--border 20 \
					--out ${dataset##*/}_gcta_mlmi_top \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	echo -e '\nStep ❸'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --clump ${dataset##*/}_gcta_mlmi.mlma \
					--bfile ${dataset} \
					--remove ${outliers} \
					--clump-range ${REFPATH}plink/glist-hg19 \
					--clump-range-border 20 \
					--out ${dataset##*/}_gcta_mlmi_clump \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	
	
	######    Plot and zip
	python ${SCRIPTPATH}manhattan_plot.py \
					${dataset##*/}_gcta_mlmi.mlma 8,0,2,1 0.05 1
	python ${SCRIPTPATH}qq-plot.py ${dataset##*/}_gcta_mlmi.mlma 8 1 \
							$(awk '{if($6==2) print}' ${dataset}.fam | wc -l),$(awk '{if($6==1) print}' ${dataset}.fam | wc -l)
	gzip -f ${dataset##*/}_gcta_mlmi.mlma
	
	
	
	
	
	####
	#### MLM leaving-one-chromosome-out (LOCO) analysis
	####
	
	${GCTAPATH} --mlma-loco \
				--autosome \
				--bfile ${dataset} \
				--remove ${outliers} \
				--pheno ${dataset}.fam \
				--mpheno 4 \
				--out ${dataset#*/}_gcta_mlmaloco \
				--thread-num ${threads}
	
	sed -i -e '1!b;s/Chr/CHR/' \
		-e '1!b;s/bp/BP/' \
		-e '1!b;s/p/P/' \
		-e '1!b;s/se/SE/' \
		-e '1!b;s/b/STAT/' ${dataset##*/}_gcta_mlmaloco.loco.mlma
	
	######    Get top results ( p <= 1e-3 )
	cat <(head -1 ${dataset##*/}_gcta_mlmaloco.loco.mlma ) \
		<(sort -gk9 ${dataset##*/}_gcta_mlmaloco.loco.mlma | grep -v nan \
							| awk '{if($9<=1e-3) print}') \
			> ${dataset##*/}_gcta_mlmaloco_top.loco.mlma
	
	echo -e '\nStep ❸'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --annotate ${dataset##*/}_gcta_mlmaloco_top.loco.mlma \
								ranges=${REFPATH}plink/glist-hg19 \
					--border 20 \
					--out ${dataset##*/}_gcta_mlmaloco_top \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	echo -e '\nStep ❸'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --clump ${dataset##*/}_gcta_mlmaloco.loco.mlma \
					--bfile ${dataset} \
					--remove ${outliers} \
					--clump-range ${REFPATH}plink/glist-hg19 \
					--clump-range-border 20 \
					--out ${dataset##*/}_gcta_mlmaloco_clump \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	
	
	######    Plot and zip
	python ${SCRIPTPATH}manhattan_plot.py \
					${dataset##*/}_gcta_mlmaloco.loco.mlma 8,0,2,1 0.05 1
	python ${SCRIPTPATH}qq-plot.py ${dataset##*/}_gcta_mlmaloco.loco.mlma 8 1 \
							$(awk '{if($6==2) print}' ${dataset}.fam | wc -l),$(awk '{if($6==1) print}' ${dataset}.fam | wc -l)
	gzip -f ${dataset##*/}_gcta_mlmaloco.loco.mlma
	
	
	echo '
	\newpage
	\subsection{Mixed Linear Models Association}
	In this section we desribed the MLM association as utilized by GCTA. For the analysis to have maximum accuracy, a total number of more than 3000 samples is required.
%%%	\begin{figure}[H]
%%%		\centering
%%%		\includegraphics[width=\textwidth,keepaspectratio]{{'`readlink -e \
%%%					manhattan-plot-${dataset##*/}_gcta_mlma.mlma.png | sed 's/\.png//' `'}.png}
%%%		\caption{Manhattan plot using MLM and a GRM}
%%%	\end{figure}

%%%	\begin{figure}[H]
%%%		\centering
%%%		\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e \
%%%							qq-plot-${dataset##*/}_gcta_mlma.mlma.png | sed 's/\.png//' `'}.png}
%%%		\caption{QQ-plot using MLM and a GRM}
%%%	\end{figure}
	
	\begin{figure}[H]
		\centering
		\includegraphics[width=\textwidth,keepaspectratio]{{'`readlink -e \
					manhattan-plot-${dataset##*/}_gcta_mlmi.mlma.png | sed 's/\.png//' `'}.png}
		\caption{Manhattan plot using MLMi}
	\end{figure}

	\begin{figure}[H]
		\centering
		\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e \
							qq-plot-${dataset##*/}_gcta_mlmi.mlma.png | sed 's/\.png//' `'}.png}
		\caption{QQ-plot using MLMi}
	\end{figure}
	
	\begin{table}[H]
		\caption{Top annotated results using MLMi}
		\centering
		\begin{scriptsize}
		\begin{tabularx}{\textwidth}{lLLLLLLl}\toprule
			'`awk '{print $2,$1,$3,$4,$7,$9,$10}' \
					${dataset##*/}_gcta_mlmi_top.annot \
					| head -20 \
					| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
			\bottomrule
		\end{tabularx}
		\end{scriptsize}
	\end{table}
	
	\begin{table}[H]
		\caption{Top annotated clumped using MLMi}
		\centering
		\begin{scriptsize}
		\begin{tabularx}{\textwidth}{lLLLlLl}\toprule
			'`awk '{print $2,$1,$3,$4,$5,$6,$7}' ${dataset##*/}_gcta_mlmi_clump.clumped.ranges \
					| head -20 \
					| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
			\bottomrule
		\end{tabularx}
		\end{scriptsize}
	\end{table}
	
	\newpage
	
	\begin{figure}[H]
		\centering
		\includegraphics[width=\textwidth,keepaspectratio]{{'`readlink -e \
			manhattan-plot-${dataset##*/}_gcta_mlmaloco.loco.mlma.png | sed 's/\.png//' `'}.png}
		\caption{Manhattan plot using MLMe-LOCO}
	\end{figure}

	\begin{figure}[H]
		\centering
		\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e \
					qq-plot-${dataset##*/}_gcta_mlmaloco.loco.mlma.png | sed 's/\.png//' `'}.png}
		\caption{QQ-plot using MLMe-LOCO}
	\end{figure}
	
	\begin{table}[H]
		\caption{Top annotated results using MLMi}
		\centering
		\begin{scriptsize}
		\begin{tabularx}{\textwidth}{lLLLLLLl}\toprule
			'`awk '{print $2,$1,$3,$4,$7,$9,$10}' \
					${dataset##*/}_gcta_mlmaloco_top.annot \
					| head -20 \
					| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
			\bottomrule
		\end{tabularx}
		\end{scriptsize}
	\end{table}
	
	\begin{table}[H]
		\caption{Top annotated clumped using MLMi}
		\centering
		\begin{scriptsize}
		\begin{tabularx}{\textwidth}{lLLLlLl}\toprule
			'`awk '{print $2,$1,$3,$4,$5,$6,$7}' ${dataset##*/}_gcta_mlmaloco_clump.clumped.ranges \
					| head -20 \
					| sed 's/ /\&/g' | awk '{print $0,"\\\\\\\\"}'`'
			\bottomrule
		\end{tabularx}
		\end{scriptsize}
	\end{table}
	
	' >> report.tex
	
	echo '
%	\begin{frame}{Association of \detokenize{'${dataset##*/}'}}
%	\begin{figure}[H]
%		\centering
%		\includegraphics[width=\textwidth,keepaspectratio]{{'`readlink -e \
%					manhattan-plot-${dataset##*/}_gcta_mlma.mlma.png | sed 's/\.png//' `'}.png}
%		\caption{Manhattan plot using MLM and a GRM}
%	\end{figure}
%	\end{frame}
	
%	\begin{frame}{Association of \detokenize{'${dataset##*/}'}}
%	\begin{figure}[H]
%		\centering
%		\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e \
%							qq-plot-${dataset##*/}_gcta_mlma.mlma.png | sed 's/\.png//' `'}.png}
%		\caption{QQ-plot using MLM and a GRM}
%	\end{figure}
%	\end{frame}
	
	\begin{frame}{Association of \detokenize{'${dataset##*/}'}}
	\begin{figure}[H]
		\centering
		\includegraphics[width=\textwidth,keepaspectratio]{{'`readlink -e \
					manhattan-plot-${dataset##*/}_gcta_mlmi.mlma.png | sed 's/\.png//' `'}.png}
		\caption{Manhattan plot using MLMi}
	\end{figure}
	\end{frame}
	
	\begin{frame}{Association of \detokenize{'${dataset##*/}'}}
	\begin{figure}[H]
		\centering
		\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e \
							qq-plot-${dataset##*/}_gcta_mlmi.mlma.png | sed 's/\.png//' `'}.png}
		\caption{QQ-plot using MLMi}
	\end{figure}
	\end{frame}
	
	\begin{frame}{Association of \detokenize{'${dataset##*/}'}}
	\begin{figure}[H]
		\centering
		\includegraphics[width=\textwidth,keepaspectratio]{{'`readlink -e \
			manhattan-plot-${dataset##*/}_gcta_mlmaloco.loco.mlma.png | sed 's/\.png//' `'}.png}
		\caption{Manhattan plot using MLMe-LOCO}
	\end{figure}
	\end{frame}
	
	\begin{frame}{Association of \detokenize{'${dataset##*/}'}}
	\begin{figure}[H]
		\centering
		\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e \
					qq-plot-${dataset##*/}_gcta_mlmaloco.loco.mlma.png | sed 's/\.png//' `'}.png}
		\caption{QQ-plot using MLMe-LOCO}
	\end{figure}
	\end{frame}
	
	' >> slides.tex
fi


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
