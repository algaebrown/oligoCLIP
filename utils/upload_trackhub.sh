module load maketrackhubs
cd $1
folder=$2 #CITS or COV
name=$3
rep=$4

maketrackhubs \
--hub $name \
--genome hg38 \
$rep/bw/$folder/*.bw \