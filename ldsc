

setwd("~/Desktop/Applied Genomics/BMEC/db_db_brain_integrated/DEG_celltypes")

###reading the DEG files
tmp = list.files(pattern="*_markers.txt")

tmp

#read and separate the files
Allfiles = lapply(tmp,read.delim,sep="\t")

filenames=list.files(pattern="*_markers.txt", full.names=TRUE)

filenames

# get rid of unnessasary information
names(Allfiles) = gsub("_db_markers.txt", "", filenames)

names(Allfiles)

Allfiles2 <- Allfiles

# leaving only the celltype and gene names
library(tidyr)
A=lapply(Allfiles2,function(x)x[,1])
str(A)

# melting all the data and putting them together
All_Genes=melt(A)
head(All_Genes)
All_Genes$setID=All_Genes$L1
All_Genes$geneName=All_Genes$value
All_Genes$value=NULL
All_Genes$L1=NULL

# get rid of the ./ before celltype
All_Genes$setID = gsub("./", "", All_Genes$setID)

Rah_Gene=read.csv("ENSG_GeneSymbol_GeneName_HUMAN.txt",sep=",")
head(Rah_Gene)

# how many genes are found in human
length(intersect(All_Genes$geneName,Rah_Gene$Gene.name))

# create a column that contains the human gene name
All_Genes$geneID= Rah_Gene$Human.gene.stable.ID[match(paste(All_Genes$geneName),paste(Rah_Gene$Gene.name))]
head(All_Genes)

# identify the blank ones and get rid of them
All_Genes[All_Genes==""]<-NA
B.B=All_Genes %>% drop_na()
B.B$geneName=NULL
head(B.B)
dim(B.B)

#creating a txt for each celltype
B.B.LIST=unstack(B.B, (B.B$geneID) ~ B.B$setID)

 

AllGenes=B.B.LIST

 

for (i in 1:length(AllGenes)) { write.table(AllGenes[i], file=paste0(names(AllGenes)[i], ".ENTID.txt"),sep="\t",row.names=F, col.names = F, quote=FALSE)
}



############################################
############################################
############################################
############################################
############################################
#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -m bea
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=240gb
#PBS -V


cd /rds/general/user/zw4419/home/LDSC/ldsc
#singularity shell *.sif


singularity exec -B $TMPDIR:/tmp *.sif bash -c "source activate ldsc && /rds/general/user/zw4419/home/LDSC/A_ALL.sh"



############################################
############################################
############################################
############################################
############################################
############################################

cd /rds/general/user/zw4419/home/LDSC

cd ldsc




for i in *.ENTID.txt;do

    basefile_name=$(basename $i .txt)

    for((j=1; j<=22; j=j+1));do

    ./make_annot.py \

    --gene-set-file $i \

    --gene-coord-file /rds/general/user/zw4419/home/LDSC/REF/ENSG_coord.txt \

    --windowsize 100000 \

    --bimfile /rds/general/user/zw4419/home/LDSC/REF/1000G_EUR_Phase3_plink/1000G.EUR.QC.$j.bim \

    --annot-file $basefile_name.$j.annot.gz




    mv $basefile_name.$j.annot.gz /rds/general/user/zw4419/home/LDSC/Annot_ALL

    done

done




for i in *.ENTID.txt;do

    basefile_name=$(basename $i .txt)

    for((j=1; j<=22; j=j+1));do

    ./ldsc.py \

    --l2 \

    --bfile /rds/general/user/zw4419/home/LDSC/REF/1000G_EUR_Phase3_plink/1000G.EUR.QC.$j \

    --ld-wind-cm 1 \

    --thin-annot \

    --annot /rds/general/user/zw4419/home/LDSC/Annot_ALL/$basefile_name.$j.annot.gz \

    --out $basefile_name.$j \

    --print-snps /rds/general/user/zw4419/home/LDSC/REF/hapmap3_snps/hm.$j.snp




    mv $basefile_name.$j.l2.M $basefile_name.$j.l2.M_5_50 $basefile_name.$j.l2.ldscore.gz -t /rds/general/user/zw4419/home/LDSC/Annot_ALL




        done

done


###################################
###################################
###################################
###################################


qsub -I -l select=01:ncpus=8:mem=96gb -l walltime=08:00:00


ssh zw4419@login-a.hpc.ic.ac.uk
export TMUX_TMPDIR=$EPHEMERAL
singularity shell *.sif
source activate ldsc



