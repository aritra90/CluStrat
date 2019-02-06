#!/bin/sh -l
source /depot/pdrineas/data/gwas_scripts/software_paths.conf

######---------------------------------------------------------
######    Last updated: 2017 Jul 17
######---------------------------------------------------------



######---------------------------------------------------------
######    Input arguments

function usage()
{
    echo "Usage: $0 [--maf <float>][--hwe <float>][--snpmissing <float>][--remainingSNPS <integer>][--relatednessThresh <float>][--mhcchrom <integer>][--mhcstart <integer>][--mhcend <integer>][--chr8invchrom <integer>][--chr8invstart <integer>][--chr8invstop <integer>][--chr17invchrom <integer>][--chr17invstart <integer>][--chr17invstop <integer>]"
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
relatednessThresh=0.1

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
	--relatednessThresh)
        #echo "$2"
        relatednessThresh=$2
        if ! [[ $relatednessThresh =~ $float_re ]]
        then
            usage
            echo "Not a float!!"
            exit 1
        fi
        if (( $(echo "$relatednessThresh > 1" |bc -l) ||  $(echo "$relatednessThresh < 0" |bc -l)  )); then
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
echo "relatednessThresh: "$relatednessThresh
echo "---------------------------------------------------------------------------"


dataset=`readlink -e $1.bed | sed 's/\.bed//'`
outliers=`readlink -e $2`
threads=$3
memory=$4
######---------------------------------------------------------



######---------------------------------------------------------
######    Keep a log
date >> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
echo 'Running ibd_check.sh with the following arguments:' \
					>> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
for i; do
	echo -e '\t'$i >> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
done
echo -e '\n' >> /depot/pdrineas/data/gwas_scripts/scripthistory_gwas.log
######---------------------------------------------------------


######---------------------------------------------------------
######     Start analysis
######---------------------------------------------------------

mkdir ${dataset##*/}_ibd_$(echo ${outliers##*/} | cut -d'.' -f1)
cd ${dataset##*/}_ibd_$(echo ${outliers##*/} | cut -d'.' -f1)
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

# >&2 echo "MARKER 1"

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

echo -e '\n\nStep ❶'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --bfile ${dataset} \
				--geno ${snpmissing} \
				--maf ${maf} \
				--hwe ${hwe} \
				--exclude mhc817nonauto.snps \
				--indep-pairwise 200 100 0.2 \
				--remove ${outliers} \
				--out ibd1 \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'


############################################################################
####################### MYSON CMD LINE CHANGES #############################
############################################################################

# >&2 echo "MARKER 2"

############   If there are still more than 150K SNPs in the dataset,
############   then prune again using the same settings.
if [ "$(wc -l ibd1.prune.in | awk '{print $1}')" -gt $remainingSNPS ] ; then
	echo -e '\nStep ❷'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --bfile ${dataset} \
					--extract ibd1.prune.in \
					--indep-pairwise 200 100 0.2 \
					--remove ${outliers} \
					--out ibd2 \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
	echo -e '\nStep ❷'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --bfile ${dataset} \
					--extract ibd2.prune.in \
					--make-bed \
					--remove ${outliers} \
					--out ibd \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

############################################################################
####################### MYSON CMD LINE CHANGES #############################
############################################################################

else
	echo -e '\nStep ❷'
	echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
	${PLINK2PATH} --bfile ${dataset} \
					--extract ibd1.prune.in \
					--make-bed \
					--remove ${outliers} \
					--out ibd \
					--memory ${memory}
	echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'
fi

# >&2 echo "MARKER 3"

echo -e '\nStep ❸'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --bfile ibd \
				--nonfounders \
				--freq \
				--remove ${outliers} \
				--out ibdfreq \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

# >&2 echo "MARKER 4"

echo -e '\nStep ❹'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
${PLINK2PATH} --bfile ibd \
				--genome full unbounded \
				--min 0.05 \
				--remove ${outliers} \
				--threads ${threads} \
				--allow-no-sex \
				--out ibdgenome \
				--memory ${memory}
echo -e '╚════════════════════════════════════════════════════════════════════════════╝\n'

############   Just for check, keep a list of couples related >0.1
############...Normal threshold is 0.1875


############################################################################
####################### MYSON CMD LINE CHANGES #############################
############################################################################

# >&2 echo "MARKER 5"

awk -v var="$relatednessThresh" '{if($10>var) print}' ibdgenome.genome > ibdcheck01.rel
awk '{if($10>0.1875) print}' ibdcheck01.rel > ibdcheck02.rel
awk '{print $3,$4}' ibdcheck02.rel > ibdout.ind

awk '{if($10>0.95) print}' ibdcheck01.rel > ibdcheck_dup.rel        #### See the duplicates
awk '{if(($10>0.4)&&($10<0.6)) print}' ibdcheck01.rel > ibdcheck_50.rel    ####See the siblings, parents
awk '{if($8>0.9) print}' ibdcheck_50.rel > ibdcheck_parents.rel
awk '{if($8<0.9) print}' ibdcheck_50.rel > ibdcheck_siblings.rel

############################################################################
####################### MYSON CMD LINE CHANGES #############################
############################################################################


cat <(awk '{print $1,$2}' ibdcheck01.rel) <(awk '{print $3,$4}' ibdcheck01.rel) | sort | uniq -c | sort -gk1 | awk '{print $2,$3,$1}' > most_pairings01_double.ind
awk '{if($3>1) print}' most_pairings01_double.ind > top01.ind
cat <(awk '{print $1,$2}' ibdcheck02.rel) <(awk '{print $3,$4}' ibdcheck02.rel) | sort | uniq -c | sort -gk1 | awk '{print $2,$3,$1}' > most_pairings02_double.ind
awk '{print $1,$2}' ibdcheck02.rel | sort | uniq -c | sort -gk1 | awk '{print $2,$3,$1}' > most_pairings02.ind
awk '{print $1,$2}' ibdcheck01.rel | sort | uniq -c | sort -gk1 | awk '{print $2,$3,$1}' > most_pairings01.ind

echo `grep -v 'HOMHOM' ibdcheck02.rel | wc -l `' pairs with over 18.75% relatedness' \
	> ibd_check_summary.txt
echo `grep -v 'HOMHOM' ibdcheck_dup.rel |  wc -l `' duplicate pairs' >> ibd_check_summary.txt
echo `grep -v 'HOMHOM' ibdcheck_parents.rel |  wc -l `' parent pairs' >> ibd_check_summary.txt
echo `grep -v 'HOMHOM' ibdcheck_siblings.rel |  wc -l `' sibling pairs' >> ibd_check_summary.txt

echo -e '\n\n'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
echo -e '                            Starting plotting round'
python ${SCRIPTPATH}histogram.py ibdgenome.genome -col=9 -header -binsize=0.01 -log
python ${SCRIPTPATH}ibd_scatter_plot.py ibdcheck01.rel 6:7
python ${SCRIPTPATH}ibd_scatter_plot.py ibdcheck01.rel 7:8
echo -e '╚════════════════════════════════════════════════════════════════════════════╝'


gzip -f ibdgenome.genome

echo -e '\n\n'
echo -e '╔════════════════════════════════════════════════════════════════════════════╗'
echo -e '               Writing the report and the slides on LaTeX ✄ ✏\n'


echo '
\newpage
\section{IBD analysis of \detokenize{'${dataset##*/}'}}

In this section we describe the IBD analysis. The IBD analysis preceeds the Principal Component Analysis, and is used to detect high-IBD pairs. The cut-off for the pairs is a pi-hat of $0.1875$. We remove one sample of each pair, focusing on removing the least cases, and the samples that aggregate the most IBD pairings. For the IBD analysis, we first run an LD-pruning step.\\

\begin{table}[H]
	\caption{Overview of samples in the analysis}
	\centering
	\begin{tabular}{lllll}\toprule
		Dataset&Cases&Controls&Total&SNPs \\ \midrule
		Input dataset&'`awk '{if($6==2) print}' ${dataset}.fam | wc -l`'&'`awk '{if($6==1) print}' ${dataset}.fam | wc -l`'&'`wc -l ${dataset}.fam | awk '{print $1}'`'&'`wc -l ${dataset}.bim | awk '{print $1}'`'\\
		IBD-ready dataset&'`awk '{if($6==2) print}' ibd.fam | wc -l`'&'`awk '{if($6==1) print}' ibd.fam | wc -l`'&'`wc -l ibd.fam | awk '{print $1}'`'&'`wc -l ibd.bim | awk '{print $1}'`' \\
		\bottomrule
	\end{tabular}
\end{table}

We show a histogram of the pair frequencies of all the IBD rates.\\

\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-ibdgenome.genome_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of the pairs generated by the IBD analysis}
\end{figure}

The next plots show the Z0,Z1 and Z2 scatter plots. The color represents the pihat $(=Z1+1/2*Z2)$

\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e ibd-plot-ibdcheck01-6:7-2d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the Z probability values generated by the IBD analysis}
\end{figure}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e ibd-plot-ibdcheck01-7:8-2d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the Z probability values generated by the IBD analysis}
\end{figure}

\lstinputlisting[breaklines,
				firstline=1,
				frame=single,
				title={Sample summary in the dataset after Individual QC}
						]{'`readlink -e ibd_check_summary.txt`'}

' >> report.tex

echo '
\section{IBD analysis of \detokenize{'${dataset##*/}'}}

\begin{frame}{IBD check of \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e hist-plot-ibdgenome.genome_log.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Histogram of the pairs generated by the IBD analysis}
\end{figure}
\end{frame}

\begin{frame}{IBD check of \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e ibd-plot-ibdcheck01-6:7-2d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the Z probability values generated by the IBD analysis}
\end{figure}
\end{frame}

\begin{frame}{IBD check of \detokenize{'${dataset##*/}'}}
\begin{figure}[H]
	\centering
	\includegraphics[width=0.5\textwidth,keepaspectratio]{{'`readlink -e ibd-plot-ibdcheck01-7:8-2d.pdf | sed 's/\.pdf//' `'}.pdf}
	\caption{Scatterplot of the Z probability values generated by the IBD analysis}
\end{figure}
\end{frame}

\begin{frame}{IBD check of \detokenize{'${dataset##*/}'}}
\lstinputlisting[breaklines,
				firstline=1,
				frame=single,
				title={Sample summary in the dataset after Individual QC}
						]{'`readlink -e ibd_check_summary.txt`'}
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
