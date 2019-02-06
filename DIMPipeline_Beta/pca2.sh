#!/bin/sh -l
source /depot/pdrineas/data/gwas_scripts/software_paths.conf

######---------------------------------------------------------
######    Last updated: 2017 Jul 17
######---------------------------------------------------------



######---------------------------------------------------------
######    Input arguments

function usage()
{
    echo "Usage: $0 [--maf <float>][--hwe <float>][--snpmissing <float>][--remainingSNPS <integer>][--mhcchrom <integer>][--mhcstart <integer>][--mhcend <integer>][--chr8invchrom <integer>][--chr8invstart <integer>][--chr8invstop <integer>][--chr17invchrom <integer>][--chr17invstart <integer>][--chr17invstop <integer>]"
}

PARAMS=""

# Default values - TAKE SPECIAL CARE IF YOU CHANGE THESE
maf=0.05
hwe=0.001
snpmissing=0.02

mhcchrom=6
mhcstart=25000000
mhcend=35000000

chr8invchrom=8
chr8invstart=7000000
chr8invstop=13000000

chr17invchrom=17
chr17invstart=40900000
chr17invstop=44900000

remainingSNPS=150000

float_re='^[+-]?[0-9]+([.][0-9]+)?$'
int_re='^[0-9]+$'

while (( "$#" )); do
  case "$1" in
    --maf)
        #echo "$2"
        maf=$2
        if ! [[ $maf =~ $float_re ]]
        then
            usage
            echo "Not a float!!"
            exit 1
        fi
        if (( $(echo "$maf > 1" |bc -l) ||  $(echo "$maf < 0" |bc -l)  )); then
            usage
            echo "Enter a float between 0 and 1"
            exit 1
        fi
        shift 2
        ;;
    --hwe)
      #echo "$2"
      hwe=$2
      if ! [[ $hwe =~ $float_re ]]
      then
          usage
          echo "Not a float!!"
          exit 1
      fi
      if (( $(echo "$hwe > 1" |bc -l) ||  $(echo "$hwe < 0" |bc -l)  )); then
          usage
          echo "Enter a float between 0 and 1"
          exit 1
      fi
      shift 2
      ;;
    --snpmissing)
        #echo "$2"
        snpmissing=$2
        if ! [[ $snpmissing =~ $float_re ]]
        then
            usage
            echo "Not a float!!"
            exit 1
        fi
        if (( $(echo "$snpmissing > 1" |bc -l) ||  $(echo "$snpmissing < 0" |bc -l)  )); then
            usage
            echo "Enter a float between 0 and 1"
            exit 1
        fi
        shift 2
        ;;
	--remainingSNPS)
        #echo "$2"
        remainingSNPS=$2
        if ! [[ $remainingSNPS =~ $int_re ]]
        then
            usage
            echo "Not a positive integer!!"
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
echo "Parameter values"
echo "---------------------------------------------------------------------------"
echo "maf: "$maf
echo "hwe: "$hwe
echo "snpmissing: "$snpmissing
echo "mhcchrom: "$mhcchrom
echo "mhcstart: "$mhcstart
echo "mhcend: "$mhcend
echo "chr8invchrom: "$chr8invchrom
echo "chr8invstart: "$chr8invstart
echo "chr8invstop: "$chr8invstop
echo "chr17invchrom: "$chr17invchrom
echo "chr17invstart: "$chr17invstart
echo "chr17invstop: "$chr17invstop
echo "remainingSNPS: "$remainingSNPS
echo "---------------------------------------------------------------------------"



dataset=`readlink -e $1.bed | sed 's/\.bed//'`
outliers=`readlink -e $2`
projectfile=$3
threads=$4
memory=$5
######---------------------------------------------------------



######---------------------------------------------------------
######    Keep a log
date >> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
echo 'Running pca.sh with the following arguments:' \
					>> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
for i; do
	echo -e '\t'$i >> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
done
echo -e '\n' >> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
######---------------------------------------------------------


######---------------------------------------------------------
######     Start analysis
######---------------------------------------------------------

mkdir ${dataset##*/}_pca_$(echo ${outliers##*/} | cut -d'.' -f1)
cd ${dataset##*/}_pca_$(echo ${outliers##*/} | cut -d'.' -f1)
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

commandvar="awk '{if((\$1==$mhcchrom && \$4>$mhcstart && \$4<$mhcend) || (\$1==$chr8invchrom && \$4>$chr8invstart && \$4<$chr8invstop) || (\$1==$chr17invchrom && \$4>$chr17invstart && \$4<$chr17invstop)) print \$2, \$1, \$4}' ${dataset}.bim > mhc817.snps"
eval $commandvar
awk '{if($1>22 || $1==0 ) print}' ${dataset}.bim > nonautosomal.snps
cat mhc817.snps nonautosomal.snps > mhc817nonauto.snps

echo -e '\n'
echo -e '🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔'
echo -e '🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙒🙜 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔 🙞🙔'
echo -e '🙒🙜 🙒🙜 🙒🙜                                                             🙞🙔 🙞🙔 🙞🙔'
echo -e '🙒🙜 🙒🙜                                                                   🙞🙔 🙞🙔'
echo -e '🙒🙜                                                                         🙞🙔'

if [ "${projectfile}" != "no" ] ; then
	awk '{print $1,$2}' ${projectfile} > project.ind
	awk '{print $1,$2,3}' ${projectfile} > project.pheno
	cp ${dataset}.fam projected.pheno1
	while read line ; do
		sed -i "/$line/d" projected.pheno1
	done < project.ind
	awk '{print $1,$2,$6}' projected.pheno1 > projected.pheno
	cat project.pheno projected.pheno > all.pheno
else
	awk '{print $1,$2,$6}' ${dataset}.fam > all.pheno
fi

echo -e '\nStep ❶'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --bfile ${dataset} \
				--geno ${snpmissing} \
				--maf ${maf} \
				--hwe ${hwe} \
				--exclude mhc817nonauto.snps \
				--indep-pairwise 200 100 0.2 \
				--remove ${outliers} \
				--out pca1 \
				--snps-only \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

############   If there are still more than 150K SNPs in the dataset,
############   then prune again using the same settings.


############################################################################
####################### MYSON CMD LINE CHANGES #############################
############################################################################


if [ "$(wc -l pca1.prune.in | awk '{print $1}')" -gt $remainingSNPS ] ; then
	echo -e '\nStep ❷'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --bfile ${dataset} \
					--extract pca1.prune.in \
					--indep-pairwise 200 100 0.2 \
					--remove ${outliers} \
					--out pca2 \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	echo -e '\nStep ❷'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --bfile ${dataset} \
					--pheno all.pheno \
					--extract pca2.prune.in \
					--make-bed \
					--remove ${outliers} \
					--out pca \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'


	############################################################################
	####################### MYSON CMD LINE CHANGES #############################
	############################################################################


else
	echo -e '\nStep ❷'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --bfile ${dataset} \
					--pheno all.pheno \
					--extract pca1.prune.in \
					--make-bed \
					--out pca \
					--remove ${outliers} \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	python /depot/pdrineas/data/gwas_pipelines/pca_bims.py pca.bim
	cp pca.bim pca.bim0
	cp pca2.bim pca.bim
fi

echo 'genotypename: pca.bed
snpname:      pca.bim
indivname:    pca.fam
evecoutname:  pca.evec
evaloutname:  pca.eval
numoutevec:   20
numoutlieriter: 0
numthreads: '${threads} > pca.parfile
if [ "$projectfile" != "no" ] ; then
echo '3' > projectpop.txt
echo 'poplistname: projectpop.txt' >> pca.parfile
fi
echo 'familynames: NO
altnormstyle: NO' >> pca.parfile

echo -e '\nStep ❸'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${EIGPATH}bin/smartpca -p pca.parfile | tee pca.pcalog
echo -e '╚════════════════════════════════════════════════════════════════════════════╝'


echo -e '\n\n'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
echo -e '                            Starting plotting round'
python ${SCRIPTPATH}pca-2d-plot_v2.py pca.evec 1:2
python ${SCRIPTPATH}pca-2d-plot_v2.py pca.evec 2:3
python ${SCRIPTPATH}pca-2d-plot_v2.py pca.evec 3:4
python ${SCRIPTPATH}pca-2d-plot_v2.py pca.evec 4:5
python ${SCRIPTPATH}pca-2d-plot_v2.py pca.evec 5:6
python ${SCRIPTPATH}pca-2d-plot_v2.py pca.evec 7:8
python ${SCRIPTPATH}pca-2d-plot_v2.py pca.evec 9:10
python ${SCRIPTPATH}pca-3d-plot_v2.py pca.evec 1:2:3
python ${SCRIPTPATH}pca-3d-plot_v2.py pca.evec 3:4:5
python ${SCRIPTPATH}pca-3d-plot_v2.py pca.evec 6:7:8
python ${SCRIPTPATH}pca-1d-plot_v2.py pca.evec 1:2:3
python ${SCRIPTPATH}pca-1d-plot_v2.py pca.evec 4:5:6
python ${SCRIPTPATH}pca-1d-plot_v2.py pca.evec 7:8:9
python ${SCRIPTPATH}pca-1d-plot_v2.py pca.evec 1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20
echo -e '╚════════════════════════════════════════════════════════════════════════════╝'
awk '{print $1,$0}' pca.evec > pca.covar

echo -e '\n\n'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
echo -e '               Writing the report and the slides on LaTeX ✄ ✏\n'


echo '
\newpage
\section{PCA of \detokenize{'${dataset##*/}'}}

In this section we describe the results of the Principal Component Analysis. The PCA step is run after the IBD analysis, removing the IBD individuals that could lead to confounding. The PCA step is also rerun a second time to remove population outliers from the PCA. For the PCA, we first run an LD-pruning step.

\begin{table}[H]
	\caption{Overview of samples in the analysis}
	\centering
	\begin{tabular}{lllll}\toprule
		Dataset&Cases&Controls&Total&SNPs \\ \midrule
		Input dataset&'`awk '{if($6==2) print}' ${dataset}.fam | wc -l`'&'`awk '{if($6==1) print}' ${dataset}.fam | wc -l`'&'`wc -l ${dataset}.fam | awk '{print $1}'`'&'`wc -l ${dataset}.bim | awk '{print $1}'`'\\
		PCA-ready dataset&'`awk '{if($6==2) print}' pca.fam | wc -l`'&'`awk '{if($6==1) print}' pca.fam | wc -l`'&'`wc -l pca.fam | awk '{print $1}'`'&'`wc -l pca.bim | awk '{print $1}'`' \\
		\bottomrule
	\end{tabular}
\end{table}

\begin{table}[H]
	\caption{Eigenvalues of the eigenvectors generated}
	\centering
	\begin{tabular}{l|llllllllll}\toprule
		Eigenvectors&1&2&3&4&5&6&7&8&9&10 \\ \midrule
		Eigenvalues&'`head -1 pca.evec | xargs | cut -d" " -f2-11 | tr " " "&" `' \\
		\bottomrule
	\end{tabular}
\end{table}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e pca-plot-pca-1:2-2d.pdf| sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the first two PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e pca-plot-pca-3:4-2d.pdf| sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the third and fourth PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e pca-plot-pca-5:6-2d.pdf| sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the fifth and sixth PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e pca-plot-pca-7:8-2d.pdf| sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the seventh and eighth PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e pca-plot-pca-9:10-2d.pdf| sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the ninth and tenth PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e pca-plot-pca-1:2:3-3d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the first three PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e pca-plot-pca-3:4:5-3d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the next three PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e pca-plot-pca-1:2:3-1d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the first three PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e pca-plot-pca-4:5:6-1d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the next three PCs generated by the PCA}
\end{figure}

\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e pca-plot-pca-7:8:9-1d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the next three PCs generated by the PCA}
\end{figure}
' >> report.tex

echo '
\section{PCA of \detokenize{'${dataset##*/}'}}

\begin{frame}{PCA of \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e pca-plot-pca-1:2-2d.pdf| sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the first two PCs generated by the PCA}
\end{figure}
\end{frame}

\begin{frame}{PCA of \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e pca-plot-pca-3:4-2d.pdf| sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the third and fourth PCs generated by the PCA}
\end{figure}
\end{frame}

\begin{frame}{PCA of \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.75\textwidth,keepaspectratio]{{'`readlink -e pca-plot-pca-1:2:3-3d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the first three PCs generated by the PCA}
\end{figure}
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
