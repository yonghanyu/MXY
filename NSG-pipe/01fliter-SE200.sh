#!/bin/sh

# usage : sh 01fliter-SE200.sh infotab outpath
####################	     *** PARA SET ***		###################

outpath=$2
infotab=$1

adaptor_SE="AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAANNNNNNNNNNCAACTCCTTGGCTCACAGAACGACATGGCTACGATCCGACTT"
infofile="/share/FGI2017B/users/yuhuiyun/data/flow_pe50/xinyao_STR/20180718CG_CL100077652_SE200/pydata/STR144.v1.info.tab"
snpinfo="/ldfssz1/FGI/CCS/libowen/01_PROJECT/MultiPCR/20181022.workflow_se200/shell/snp.near25.info.v2.txt"
adapterRM_para="--trimns --trimqualities --minlength 50 "
#qsub_para="-S /bin/bash -cwd -l vf=1g -l num_proc=5 -q fgi.q -P FGI -binding linear:5" 	 #fgi
qsub_para="-S /bin/bash -cwd -l vf=6g -l num_proc=5 -q fgi.q -P NIPPT -binding linear:5"  #cngb

#################### 		MAIN SHELL 		####################
########## MAKE DIR & check
#for sam in `cat $1`
while read sam cell lane index storedir
do
mkdir -p $outpath/$sam/01clean

########## MAKE ALIGN.SH
#infq_1=$storedir/${cell}_${date}/$lane/${cell}_${date}${lane}_$index.fq.gz
fqgz_1=($(echo $index | awk -vRS="," '{split("'$lane'",ln,",");for(i in ln){print "'$storedir'/'$cell'/"ln[i]"/'$cell'_"ln[i]"_"$1".fq.gz"}}')) 

echo \
"
echo --- clean start \`date\` --- && 
zcat ${fqgz_1[@]} > $outpath/$sam/01clean/${sam}_1.merge.fq && 
gzip $outpath/$sam/01clean/${sam}_1.merge.fq && 
/share/FGI2017B/users/yuhuiyun/soft/adapterremoval-2.1.7/build/AdapterRemoval --file1 $outpath/$sam/01clean/${sam}_1.merge.fq.gz --adapter1 $adaptor_SE $adapterRM_para --basename $outpath/$sam/01clean/$sam &&
mv $outpath/$sam/01clean/$sam.truncated $outpath/$sam/01clean/${sam}_R1.clean.fq &&
echo --- clean done \`date\` --- && 
echo ---------------   CLEAN COMPLETED   --------------- && 

echo ---uniq start \`date\` --- && 
awk '{if(NR%4 == 1){print \">\" substr(\$0, 2)}}{if(NR%4 == 2){print}}' $outpath/$sam/01clean/${sam}_R1.clean.fq > $outpath/$sam/01clean/${sam}_R1.clean.fasta &&
cat $outpath/$sam/01clean/${sam}_R1.clean.fasta | sort | uniq -c | awk '{if(\$1 > 5){print \$1\"\t\"\$2}}' | sort -nr -k1 > $outpath/$sam/01clean/${sam}_R1_d5.reads && 
awk '{if(NR%2 == 0){print}}' $outpath/$sam/01clean/${sam}_R1.clean.fasta |  sort | uniq -c | sort -nr -k1 > $outpath/$sam/01clean/${sam}_R1.sort.reads && 
echo --- uniq end \`date\` --- && 

echo ---------------   uniq COMPLETED   --------------- && 

echo ---------------  genotype start \`date\`  --------------- &&
#python /share/FGI2017B/users/yuhuiyun/data/flow_pe50/xinyao_STR/20180718CG_CL100077652_SE200/pydata/STRlen.near30.string.mismach.rc.multiprocess.py -i $infofile -f $outpath/$sam/01clean/${sam}_R1_d5.reads -o $outpath/$sam/01clean/${sam}_R1_d5.reads.mult15.tab -p 5 > $outpath/$sam/01clean/${sam}_R1_d5.reads.step1.log &&
#python /share/FGI2017B/users/yuhuiyun/data/flow_pe50/xinyao_STR/20180718CG_CL100077652_SE200/pydata/STR_genotype.py -i $outpath/$sam/01clean/${sam}_R1_d5.reads.mult15.tab -o $outpath/$sam/01clean/${sam}_R1_d5.reads.mult15.gt.tab > $outpath/$sam/01clean/${sam}_R1_d5.reads.mult15.gt.tab.log &&
python /ldfssz1/FGI/CCS/libowen/01_PROJECT/MultiPCR/20181022.workflow_se200/shell/find_SNP.v2.py -i $snpinfo -f $outpath/$sam/01clean/${sam}_R1_d5.reads -o $outpath/$sam/01clean/${sam}_R1_d5.reads.snp.gt.tab -p 5 > $outpath/$sam/01clean/${sam}_R1_d5.reads.snp.gt.log &&
python /share/FGI2017B/users/yuhuiyun/data/flow_pe50/xinyao_STR/20180718CG_CL100077652_SE200/pydata/STRlen.near30.string.mismach.rc.multiprocess.outseq.py -i $infofile -f $outpath/$sam/01clean/${sam}_R1_d5.reads -o $outpath/$sam/01clean/${sam}_R1_d5.reads.rc.outseq.out -p 5 > $outpath/$sam/01clean/${sam}_R1_d5.reads.outseq.log &&
python /share/FGI2017B/users/yuhuiyun/data/flow_pe50/xinyao_STR/20180718CG_CL100077652_SE200/pydata/STR_genotype.v1.py -i $outpath/$sam/01clean/${sam}_R1_d5.reads.rc.outseq.out -o $outpath/$sam/01clean/${sam}_R1_d5.reads.rc.outseq.out.gt > $outpath/$sam/01clean/${sam}_R1_d5.reads.rc.outseq.out.gt.log &&
echo ---------------  genotype COMPLETED \`date\`  --------------- && echo

" > $outpath/$sam/$sam.sh

########## QSUB
qsub $qsub_para -o $outpath/$sam/$sam.log -e $outpath/$sam/$sam.err $outpath/$sam/$sam.sh
#sh $outpath/$sam/$sam.sh >>$outpath/$sam/$sam.log 2>>$outpath/$sam/$sam.err &

done < $infotab
