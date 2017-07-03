#!/bin/bash
## Author: Kun Lu, drlukun@swu.edu.cn
## Version: 1.0,  Date: 2017-2-27
## The final result will be Ath.best.primers in Step 11 and Ath.all.primers in Step 12.
## Prepare coding/mRNA sequences and genomic sequence in the ~/qPrimerDB/example fold, such as Ath.fa and Ath.genome.fa.
## Please install e-PCR, BLAT, Primer3, BBmap, faSomeRecords, MFEprimer-2.0 and MPprimer before analysis.

cd ~/qPrimerDB/example

# The first section used for shortening gene name. You'd better don't use symbol "_" in the gene name. In the following analysis, We use the model plant ‘Arabidopsis thaliana’ as an example to demonstrate primer design pipeline.

perl -pe 's/\_/./g' Ath.fa | awk -F " " '/>/{$0=$1}1' > Ath.shortname.fa;

## Step 1. split the cDNA sequences into small fragment with 300-bp windows and 50-bp step.

~/qPrimerDB/bbmap/shred.sh in=Ath.shortname.fa out=Ath.shortname.frag.fa length=300 overlap=250 minlength=200;

## Step 2. Create primer3 input file for large fasta sequences; split each transcriptome into smaller fasta files, each with 100000 splited fasta sequences

perl ~/qPrimerDB/MPprimer/CreateMPprimerInput.py -i Ath.shortname.frag.fa -o Ath.shortname.p3in;

## Step 3. Design primers for each p3in file
primer3_core < Ath.shortname.p3in > Ath.shortname.p3out;

## Step 4. Format output file to e-PCR file using perl script extr_primer3new.pl, and Ath.shortname.uniepcr will be used for e-PCR program.
perl ~/qPrimerDB/scripts/extr_primer3new.pl -p Ath.shortname.p3out -l Ath.shortname.list -s Ath.shortname.blast -e Ath.shortname.epcr;
awk '{print $0"\t"$2$3}' Ath.shortname.epcr | awk  '!a[$5]++' | awk '{print $1"\t"$2"\t"$3"\t"$4}' > Ath.shortname.uniepcr;

## Step 5. Check the specificity of all the primers, and filter out non-specific primers based on e-PCR results
e-PCR -w9 -f1 -n1 -g1 -m100 Ath.shortname.uniepcr D=50-1500 Ath.shortname.fa N=1 G=1 T=3 O=Ath.shortname.STS;
awk -F "_" '{print $1"\t"$2"\t"$3"\t"$4}' Ath.shortname.STS | awk '{if($8!=0 || $9!=0) print $2"_"$3}'|sort -u > Ath.uniSTS.remove.id.tmp;
awk -F "_" '{print $1"\t"$2"\t"$3"\t"$4}' Ath.shortname.STS | awk '{if($1==$2 && $8==0 && $9==0) print $2"_"$3}' > Ath.uniSTS.id.tmp;
awk 'NR==FNR{a[$1]=1}NR!=FNR{if(a[$1]!=1) print $0}'  Ath.uniSTS.remove.id.tmp Ath.uniSTS.id.tmp > Ath.uniSTS.id;

##extract results for MFEprimer analysis
grep -F -f Ath.uniSTS.id Ath.shortname.list > Ath.uniSTS.unilist;
grep -F -f Ath.uniSTS.id Ath.shortname.blast | awk -F ">" '{print $2}' > Ath.uniSTS.uniblast.id;
faSomeRecords Ath.shortname.blast Ath.uniSTS.uniblast.id Ath.uniSTS.uniblast;

## Step 6. Prepare primer files for each gene, and move them into individual fold.
mkdir Ath;
grep "^>" Ath.shortname.fa | awk -F '>' '{print $2}' > Ath.shortname.name;
for file1 in $(cat Ath.shortname.name); 
do
grep $file1 Ath.shortname.blast | awk -F ">" '{print$2}' > $file1.blast.in;
faSomeRecords Ath.shortname.blast $file1.blast.in Ath.$file1.blast.fa;
rm -rf $file1.blast.in;
mv Ath.*.blast.fa ./Ath;
done

cd ./Ath;
for file2 in $(ls *.blast.fa);
do
perl ~/qPrimerDB/scripts/fastaDeal.pl -cuts 2 $file2;
done;
rm -rf *.blast.fa;


## Step 7. Index original cDNA fasts file. Then, check the specificity of primers using MFPrimer

~/qPrimerDB/MFEprimer/IndexDb.sh Ath.shortname.fa;

## Step8. Check the specificity of primers using MFPrimer

cd ./Ath;
mkdir primers;
for file3 in $(cat ../Ath.shortname.name); 
do 
cd ./Ath.$file3.blast.fa.cut; 
for file4 in $(ls -1);
do
~/qPrimerDB/MFEprimer/MFEprimer.py --size_stop=1500 --size_start=50 --tab --ppc=10 -i $file4 -d ~/qPrimerDB/example/Ath.shortname.fa -o $file4.txt;
done; 
cat *.txt > $file3.MFP;
mv $file3.MFP ~/qPrimerDB/example/Ath/primers; 
rm -rf *.txt; 
cd ../; 
done; 

cd ~/qPrimerDB/example/Ath/primers;
cat *.MFP > Ath.primers;
mv Ath.primers ../../;
rm -f *.MFP


## Step9 format the result file to obtain gene-specific primers. l1: PPC<10; l2: FpDg <-9 && RpDG <-9; l3: FpDg <-11 && RpDG <-11
# All candidate primers
cd ../../;
awk '$1=="AmpID" {print a} {a=$0}' Ath.primers | awk '$1 =="1" {print $0}' > Ath.primers.l1;
awk '{if($1=="AmpID" || $10 < -9) print $0}' Ath.primers | awk '{if($1=="AmpID" || $11 < -9)  print $0}' > Ath.primers.l2;
awk '{if($1=="AmpID" || $10 < -11) print $0}' Ath.primers | awk '{if($1=="AmpID" || $11 < -11)  print $0}' > Ath.primers.l3;

# obtain the unique l1 primers 
awk '{print $0"\t"$10+$11}' Ath.primers.l1 | sort -k 4,4 -k 15,15n -k 12,12nr | awk  '!a[$4]++' > Ath.primers.unil1;
cut -f 4 Ath.primers.unil1 > Ath.primers.unil1.id;
# Uniq the l2 primers
awk '$1=="AmpID" {print a} {a=$0}' Ath.primers.l2 | awk '$1 =="1" {print $0}' | awk '{print $0"\t"$10+$11}'  | sort -k 4,4 -k 15,15n -k 12,12nr | awk  '!a[$4]++' > Ath.primers.unil2; 
cut -f 4 Ath.primers.unil2  > Ath.primers.unil2.id;
# Uniq the l3 primers
awk '$1=="AmpID" {print a} {a=$0}' Ath.primers.l3 | awk '$1 =="1" {print $0}' | awk '{print $0"\t"$10+$11}'  | sort -k 4,4 -k 15,15n -k 12,12nr | awk  '!a[$4]++' > Ath.primers.unil3; 
cut -f 4 Ath.primers.unil3  > Ath.primers.unil3.id;


###Step10 Extract uniq primers and primer parameters
# Extract l1 uniq primers and primer parameters
awk '{print $2}' Ath.primers.unil1 | sed 's/^L_//g' > Ath.primers.unil1.primer.id;
grep -F -f Ath.primers.unil1.primer.id Ath.shortname.epcr | cut -f 1,2,3 | sort -k 1,1 > Ath.primers.unil1.epcr.sorted;
sort -k 4,4 Ath.primers.unil1 | cut -f 2-15 | sed 's/^L_//g' | sort -k 1,1 > Ath.primers.unil1.sorted;
join -1 1 -2 1 -a 1 Ath.primers.unil1.sorted Ath.primers.unil1.epcr.sorted | perl -pe 's/ /\t/g' > Ath.primers.unil1.primers.tmp;
awk '{print $3"\tL_"$1"\t"$2"\t"$15"\t"$16"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' Ath.primers.unil1.primers.tmp > Ath.unil1.primers;
# Extract l2 uniq primer and primer parameters
awk '{print $2}' Ath.primers.unil2 | sed 's/^L_//g' > Ath.primers.unil2.primer.id;
grep -F -f Ath.primers.unil2.primer.id Ath.shortname.epcr | cut -f 1,2,3 | sort -k 1,1 > Ath.primers.unil2.epcr.sorted; 
sort -k 4,4 Ath.primers.unil2 | cut -f 2-15 | sed 's/^L_//g' | sort -k 1,1 > Ath.primers.unil2.sorted; 
join  -1 1 -2 1 -a 1 -o auto  Ath.primers.unil2.sorted Ath.primers.unil2.epcr.sorted | perl -pe 's/ /\t/g' > Ath.primers.unil2.primers.tmp; 
awk '{print $3"\tL_"$1"\t"$2"\t"$15"\t"$16"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' Ath.primers.unil2.primers.tmp > Ath.unil2.primers;
# Extract l3 uniq primer and primer parameters
awk '{print $2}' Ath.primers.unil3 | sed 's/^L_//g' > Ath.primers.unil3.primer.id;  
grep -F -f Ath.primers.unil3.primer.id Ath.shortname.epcr | cut -f 1,2,3 | sort -k 1,1 > Ath.primers.unil3.epcr.sorted; 
sort -k 4,4 Ath.primers.unil3 | cut -f 2-15 | sed 's/^L_//g' | sort -k 1,1 > Ath.primers.unil3.sorted; 
join -1 1 -2 1 -a 1 -o auto  Ath.primers.unil3.sorted Ath.primers.unil3.epcr.sorted | perl -pe 's/ /\t/g' >  Ath.primers.unil3.primers.tmp; 
awk '{print $3"\tL_"$1"\t"$2"\t"$15"\t"$16"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' Ath.primers.unil3.primers.tmp > Ath.unil3.primers;

## All candidate primers
awk '$1=="AmpID" {print a} {a=$0}' Ath.primers | awk '$1 =="1" {print $0}'  > Ath.primers.l1.all;
awk '$1=="AmpID" {print a} {a=$0}' Ath.primers.l2 | awk '$1 =="1" {print $0}'  > Ath.primers.l2.all; 
awk '$1=="AmpID" {print a} {a=$0}' Ath.primers.l3 | awk '$1 =="1" {print $0}' > Ath.primers.l3.all; 

##extract l1, l2 and l3 primers and primer parameters.
awk '{print $2}' Ath.primers.l1.all | sed 's/^L_//g' > Ath.primers.l1.all.primer.id;
grep -F -f Ath.primers.l1.all.primer.id Ath.shortname.epcr | cut -f 1,2,3 | sort -k 1,1 > Ath.primers.l1.all.epcr.sorted;
sort -k 4,4 Ath.primers.l1.all | cut -f 2-15 | sed 's/^L_//g' | sort -k 1,1 > Ath.primers.l1.all.sorted; 
join -1 1 -2 1 -a 1 -o auto  Ath.primers.l1.all.sorted Ath.primers.l1.all.epcr.sorted | perl -pe 's/ /\t/g' >  Ath.primers.l1.primers.tmp; 
awk '{print $3"\tL_"$1"\t"$2"\t"$14"\t"$15"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' Ath.primers.l1.primers.tmp > Ath.l1.all.primers;
awk '{print $2}' Ath.primers.l2.all | sed 's/^L_//g' > Ath.primers.l2.all.primer.id;  
grep -F -f Ath.primers.l2.all.primer.id Ath.shortname.epcr | cut -f 1,2,3 | sort -k 1,1 > Ath.primers.l2.all.epcr.sorted; 
sort -k 4,4 Ath.primers.l2.all | cut -f 2-15 | sed 's/^L_//g'| sort -k 1,1 > Ath.primers.l2.all.sorted; 
join -1 1 -2 1 -a 1 -o auto  Ath.primers.l2.all.sorted Ath.primers.l2.all.epcr.sorted | perl -pe 's/ /\t/g' >  Ath.primers.l2.all.primers.tmp; 
awk '{print $3"\tL_"$1"\t"$2"\t"$14"\t"$15"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' Ath.primers.l2.all.primers.tmp > Ath.l2.all.primers;
awk '{print $2}' Ath.primers.l3.all | sed 's/^L_//g' > Ath.primers.l3.all.primer.id;  
grep -F -f Ath.primers.l3.all.primer.id Ath.shortname.epcr | cut -f 1,2,3 | sort -k 1,1 > Ath.primers.l3.all.epcr.sorted; 
sort -k 4,4 Ath.primers.l3.all | cut -f 2-15 | sed 's/^L_//g'| sort -k 1,1 > Ath.primers.l3.all.sorted; 
join -1 1 -2 1 -a 1 -o auto  Ath.primers.l3.all.sorted Ath.primers.l3.all.epcr.sorted | perl -pe 's/ /\t/g' >  Ath.primers.l3.all.primers.tmp; 
awk '{print $3"\tL_"$1"\t"$2"\t"$14"\t"$15"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' Ath.primers.l3.all.primers.tmp > Ath.l3.all.primers; 


###Step11 Obtain the best gene-specific primers, which will be put into pPrimerDB.
mkdir best.primers
mv Ath.unil1.primers ./best.primers
mv Ath.unil2.primers ./best.primers
mv Ath.unil3.primers ./best.primers
cd ./best.primers

mv Ath.unil1.primers Ath.l1;
awk '{print $1}' Ath.l1 > Ath.l1.id;
awk 'NR==FNR{a[$2]=1}NR!=FNR{if(a[$1]!=1) print $0}' Ath.l1.id Ath.unil2.primers > Ath.l2; 
awk '{print $1}' Ath.l2 > Ath.l2.id; cat Ath.l1.id Ath.l2.id > Ath.l12.id; 
awk 'NR==FNR{a[$2]=1}NR!=FNR{if(a[$1]!=1) print $0}' Ath.l12.id Ath.unil3.primers > Ath.l3; 
rename "s/.l1/%%l1/" *  
rename "s/.l2/%%l2/" * 
rename "s/.l3/%%l3/" * 
awk '{print FILENAME"\t"$0}' Ath%%l1 > Ath.l1.tmp ; 
awk '{print FILENAME"\t"$0}' Ath%%l2 > Ath.l2.tmp; 
awk '{print FILENAME"\t"$0}' Ath%%l3 > Ath.l3.tmp;  
cat Ath.l1.tmp Ath.l2.tmp Ath.l3.tmp | awk -F "%%" '{print $1"\t"$2}' > Ath.4db.tmp;

# Obtain amplification parameters using blat against genome sequence.
cut -f 3,17 Ath.4db.tmp | awk '{print ">"$1"\n"$2}' > Ath.4db.best.tmp; 
blat ~/qPrimerDB/example/Ath.genome.fa Ath.4db.best.tmp -noHead Ath.4db.best.blat.tmp; 
sort -k 10,10 -k 1,1nr Ath.4db.best.blat.tmp | awk '!a[$10]++' > Ath.4db.best.blat.uniq.tmp;
sort -k3,3 Ath.4db.tmp | awk '!a[$3]++' >  Ath.4db.best.sorted.tmp;
join -1 3 -2 10 -a 1  Ath.4db.best.sorted.tmp  Ath.4db.best.blat.uniq.tmp | perl -pe 's/ /\t/g' > Ath.4db.best.records.tmp;
awk '{print $0"\tbest"}' Ath.4db.best.records.tmp > Ath.best.4db;
awk -F "\t" 'BEGIN {print "GeneID\tpLevel\tFpID\tRpID\tFprimer\tRprimer\tPPC\tAmpSize\tAmpGC\tFpTm\tRpTm	\tFpDg\tRpDg\tNumExonCorss"} {print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$34}'  Ath.best.4db > Ath.best.primers


###Step12 Obtain all the candidate gene-specific primers.
cd ..
mkdir all.primers
mv Ath.l1.all.primers ./all.primers
mv Ath.l2.all.primers ./all.primers
mv Ath.l3.all.primers ./all.primers

cd ./all.primers
rename "s/.l1.all.primers/%%l1/" *  
rename "s/.l2.all.primers/%%l2/" * 
rename "s/.l3.all.primers/%%l3/" * 
awk '{print FILENAME"\t"$0}' Ath%%l1 > Ath.l1.all.tmp ; 
awk '{print FILENAME"\t"$0}' Ath%%l2 > Ath.l2.all.tmp; 
awk '{print FILENAME"\t"$0}' Ath%%l3 > Ath.l3.all.tmp;  
cat Ath.l1.all.tmp Ath.l2.all.tmp Ath.l3.all.tmp | awk -F "%%" '{print $1"\t"$2}' > Ath.4db.all.tmp;
sort -k 3,3 -k 2,2 Ath.4db.all.tmp | awk '!a[$4]++' > Ath.4db.all.uniq.tmp;

# Obtain amplification parameters using blat against genome sequence.
cut -f 4,17 Ath.4db.all.uniq.tmp | sed 's/^L_//g' | awk '{print ">"$1"\n"$2}' > Ath.4db.all.uniq.4blat.tmp; 
blat ~/qPrimerDB/example/Ath.genome.fa Ath.4db.all.uniq.4blat.tmp -noHead Ath.4db.all.blat.tmp; 
sort -k 10,10 -k 1,1nr Ath.4db.all.blat.tmp | awk '!a[$10]++' > Ath.4db.all.blat.uniq.tmp;
sort -k4,4 Ath.4db.all.uniq.tmp > Ath.4db.all.uniq.sorted.tmp;
cut -f 4 Ath.4db.all.uniq.sorted.tmp | sed 's/^L_//g' > Ath.4db.all.uniq.sorted.tmp.id;
paste  Ath.4db.all.uniq.sorted.tmp.id Ath.4db.all.uniq.sorted.tmp > Ath.4db.all.uniq.sorted.id.tmp
join -1 1 -2 10 Ath.4db.all.uniq.sorted.id.tmp  Ath.4db.all.blat.uniq.tmp | perl -pe 's/ /\t/g' > Ath.4db.all.records.tmp;
awk '{print $0"\tall"}' Ath.4db.all.records.tmp > Ath.all.4db;
awk -F "\t" '{print $1"\t"$4"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$35}'  Ath.all.qPDB > Ath.all.qPDB.download;
awk -F "\t" 'BEGIN {print "primerID\tGeneID\tpLevel\tFpID\tRpID\tFprimer\tRprimer\tPPC\tAmpSize\tAmpGC\tFpTm\tRpTm	\tFpDg\tRpDg\tNumExonCorss"} {print $1"\t"$4"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$35}' Ath.all.4db > Ath.all.primers