#!/usr/bin/bash
#usage bash get.table.snp.forge2.sh rsid
#e.g.  bash get.table.snp.forge2.sh rs189631248

perl /data/breezece/eforge.safe_copy_nonthreaded.table.only.version2.pl -dmps $1 -noplot -noproxy --data erc2-DHS --label $1
perl /data/breezece/eforge.safe_copy_nonthreaded.table.only.version2.pl -dmps $1 -noplot -noproxy --data erc2-H3-all --label $1
perl /data/breezece/eforge.safe_copy_nonthreaded.table.only.version2.pl -dmps $1 -noplot -noproxy --data erc --label $1
perl /data/breezece/eforge.safe_copy_nonthreaded.table.only.version2.pl -dmps $1 -noplot -noproxy --data blueprint --label $1
perl /data/breezece/eforge.safe_copy_nonthreaded.table.only.version2.pl -dmps $1 -noplot -noproxy -cato 1 -eqtl 1 -eqtlgen 1 -abc 1 -tf 1 -gene 1 --label $1
perl /data/breezece/eforge.safe_copy_nonthreaded.table.only.version2.pl -dmps $1 -noplot -noproxy --data encode --label $1

#ls -tr 0x*/$1*.gz | xargs gunzip

tail -n +1 0x*/$1.GWAS*|sed 's/==>/<b>/g'|sed 's/<==/<\/b>/g'|sed 's/DHS\t/\ \ DNase I hotspot\ \ \t/g' > summary.table.$1

#include code to separate tables and add headers

ls -tr 0x*/$1.GWAS*| awk -F"/" '{print "rm -rf "$1}'|bash

#split file on empty lines

mkdir /lscratch/$SLURM_JOBID/$1
cd /lscratch/$SLURM_JOBID/$1

awk -v RS= '{print > ("splitset-" NR ".txt")}' ../summary.table.$1

cp ../tableToHtml.pl .

perl tableToHtml.pl < splitset-1.txt > splitset-1.txt.html
perl tableToHtml.pl < splitset-2.txt > splitset-2.txt.html
perl tableToHtml.pl < splitset-3.txt > splitset-3.txt.html
perl tableToHtml.pl < splitset-4.txt > splitset-4.txt.html
perl tableToHtml.pl < splitset-5.txt > splitset-5.txt.html
perl tableToHtml.pl < splitset-6.txt > splitset-6.txt.html
perl tableToHtml.pl < splitset-7.txt > splitset-7.txt.html
perl tableToHtml.pl < splitset-8.txt > splitset-8.txt.html
perl tableToHtml.pl < splitset-9.txt > splitset-9.txt.html
perl tableToHtml.pl < splitset-10.txt > splitset-10.txt.html
perl tableToHtml.pl < splitset-11.txt > splitset-11.txt.html
perl tableToHtml.pl < splitset-12.txt > splitset-12.txt.html
perl tableToHtml.pl < splitset-13.txt > splitset-13.txt.html
perl tableToHtml.pl < splitset-14.txt > splitset-14.txt.html

cp ../minibr.txt .

echo "<b>Summary and number of annotations for $1:</b>" > table.count
wc -l splitset-10.txt.html |awk '{print "RefSeq closest gene data:\t"$1-14}' >> table.count
wc -l splitset-11.txt.html |awk '{print "UCSC genome browser data:\t1"}' >> table.count
wc -l splitset-12.txt.html |awk '{print "CADD v1.6 annotations:\t"$1-14}' >> table.count
wc -l splitset-7.txt.html |awk '{print "GTEx cis-eQTLs:\t"$1-14}' >> table.count
wc -l splitset-13.txt.html |awk '{print "eQTLGen blood cis-eQTLs:\t"$1-14}' >> table.count
wc -l splitset-8.txt.html |awk '{print "ABC contacts:\t"$1-14}' >> table.count
wc -l splitset-9.txt.html |awk '{print "FORGE2-TF motifs:\t"$1-14}' >> table.count
wc -l splitset-6.txt.html |awk '{print "CATO score:\t"$1-14}' >> table.count
wc -l splitset-1.txt.html |awk '{print "FORGE2 consolidated roadmap DNase I hotspots (erc2-DHS):\t"$1-14}' >> table.count
wc -l splitset-2.txt.html |awk '{print "FORGE2 consolidated roadmap H3 histone marks (erc2-H3-all):\t"$1-14}' >> table.count
wc -l splitset-3.txt.html |awk '{print "FORGE2 unconsolidated roadmap DNase I hotspots (erc):\t"$1-14}' >> table.count
wc -l splitset-4.txt.html |awk '{print "FORGE2 blueprint DNase I hotspots (blueprint):\t"$1-14}' >> table.count
wc -l splitset-14.txt.html |awk '{print "FORGE2 ENCODE DNase I hotspots (encode):\t"$1-14}' >> table.count
wc -l splitset-5.txt.html |awk '{print "FORGE2 consolidated roadmap chromatin states (erc2-15state):\t"$1-14}' >> table.count 

perl ../tableToHtml2.pl < table.count > table.count.html

cat ../miniabout.txt ../img.txt table.count.html minibr.txt minibr.txt minibr.txt splitset-10.txt.html minibr.txt minibr.txt minibr.txt splitset-11.txt.html minibr.txt minibr.txt minibr.txt splitset-12.txt.html minibr.txt minibr.txt minibr.txt splitset-7.txt.html minibr.txt minibr.txt minibr.txt splitset-13.txt.html minibr.txt minibr.txt minibr.txt splitset-8.txt.html minibr.txt minibr.txt minibr.txt splitset-9.txt.html minibr.txt minibr.txt minibr.txt splitset-6.txt.html minibr.txt minibr.txt minibr.txt splitset-1.txt.html minibr.txt minibr.txt minibr.txt splitset-2.txt.html minibr.txt minibr.txt minibr.txt splitset-3.txt.html minibr.txt minibr.txt minibr.txt splitset-4.txt.html minibr.txt minibr.txt minibr.txt splitset-14.txt.html minibr.txt minibr.txt minibr.txt splitset-5.txt.html minibr.txt minibr.txt minibr.txt | sed 's/Feta /Fetal /g' | fgrep -v "mouse" | sed 's/0x[a-zA-Z0-9]*\//<a href=\"https:\/\/forge2.altiusinstitute.org\"\/>FORGE2<\/a> /g'|sed 's/0x[a-zA-Z0-9]*\///g' |sed 's/\.GWAS\././g'|sed 's/.chart.tsv/:/g' > summary.table.$1.html

#cat ../img.txt table.count.html minibr.txt minibr.txt minibr.txt splitset-10.txt.html minibr.txt minibr.txt minibr.txt splitset-11.txt.html minibr.txt minibr.txt minibr.txt splitset-12.txt.html minibr.txt minibr.txt minibr.txt splitset-7.txt.html minibr.txt minibr.txt minibr.txt splitset-13.txt.html minibr.txt minibr.txt minibr.txt splitset-8.txt.html minibr.txt minibr.txt minibr.txt splitset-9.txt.html minibr.txt minibr.txt minibr.txt splitset-6.txt.html minibr.txt minibr.txt minibr.txt splitset-1.txt.html minibr.txt minibr.txt minibr.txt splitset-2.txt.html minibr.txt minibr.txt minibr.txt splitset-3.txt.html minibr.txt minibr.txt minibr.txt splitset-4.txt.html minibr.txt minibr.txt minibr.txt splitset-14.txt.html minibr.txt minibr.txt minibr.txt splitset-5.txt.html minibr.txt minibr.txt minibr.txt | sed 's/Feta /Fetal /g' | fgrep -v "mouse" | sed 's/0x[a-zA-Z0-9]*\//<a href=\"https:\/\/forge2.altiusinstitute.org\"\/>FORGE2<\/a> /g'|sed 's/0x[a-zA-Z0-9]*\///g' |sed 's/\.GWAS\././g'|sed 's/.chart.tsv/:/g' > summary.table.$1.html

folderset=`echo "$1"|awk '{ print substr($1,1,6)}'`

#cp summary.table.$1.html /var/www/forge2/production/src/html/files/

#cp summary.table.$1.html /data/breezece/web.rs.dirs/$folderset

cp summary.table.$1.html ../done.htmls

awk '{print $(NF)}' table.count > $1.test.txt

fgrep -vf /lscratch/$SLURM_JOBID/snps.to.avoid $1.test.txt > $1.test.txt1 #these SNPs break the tool bc of issues re forge2 hash construction, so we want to avoid that

perl -pe 's/\S+/abs($&) > 1 ? 1 : $&/ge' < $1.test.txt1 > $1.test.txt2

score=`python /lscratch/$SLURM_JOBID/scores.py $1.test.txt2 |sed 's/:<\/b>//g' | paste - - |awk '{print $NF}'|sed 's/\.0//g'`

#score=`grep -P "$1\t" /lscratch/$SLURM_JOBID/complete.forgedb.done.htmls.tar.gz.scores2.txt |awk '{print $NF}'|tr -d "\n"`

sed -i "s@$1:@$1 (FORGEdb score=<a href="https://forge2.altiusinstitute.org/files/FORGEdb_scores.html">$score</a>):@g" ../done.htmls/summary.table.$1.html

rm -rf /lscratch/$SLURM_JOBID/$1
#rm 0x*/*$1*GWAS*
#error checking code below:

if [ `grep "erc2-15state" /lscratch/$SLURM_JOBID/done.htmls/summary.table.$1.html |grep -v 127|awk '{print $(NF-1)+1}'` ] || [ `grep "RefSeq closest gene data:" /lscratch/$SLURM_JOBID/done.htmls/summary.table.$1.html |awk '$(NF-1)!=1'|awk '{print $(NF-1)+1}'` ] || [ `grep "UCSC genome browser data:" /lscratch/$SLURM_JOBID/done.htmls/summary.table.$1.html |awk '$(NF-1)!=1'|awk '{print $(NF-1)+1}'` ] || [ `grep "CADD v1.6 annotations:" /lscratch/$SLURM_JOBID/done.htmls/summary.table.$1.html |awk '$(NF-1)!=3'|awk '{print $(NF-1)+1}'` ]
then
    cd /lscratch/$SLURM_JOBID
    TMPDIR=/lscratch/$SLURM_JOBID; bash get.table.snp.forge2.2.sh $1

fi

bash /lscratch/$SLURM_JOBID/ad.hoc.sh ../done.htmls/summary.table.$1.html

#if [ `grep "RefSeq closest gene data:" /lscratch/$SLURM_JOBID/done.htmls/summary.table.$1.html |awk '$(NF-1)!=1'|awk '{print $(NF-1)+1}'` ]
#then
#    cd /lscratch/$SLURM_JOBID
#    TMPDIR=/lscratch/$SLURM_JOBID; bash get.table.snp.forge2.2.sh $1
#fi

#if [ `grep "UCSC genome browser data:" /lscratch/$SLURM_JOBID/done.htmls/summary.table.$1.html |awk '$(NF-1)!=1'|awk '{print $(NF-1)+1}'` ]
#then
#    cd /lscratch/$SLURM_JOBID
#    TMPDIR=/lscratch/$SLURM_JOBID; bash get.table.snp.forge2.2.sh $1
#fi


#if [ `grep "CADD v1.6 annotations:" /lscratch/$SLURM_JOBID/done.htmls/summary.table.$1.html |awk '$(NF-1)!=3'|awk '{print $(NF-1)+1}'` ]
#then
#    cd /lscratch/$SLURM_JOBID
#    TMPDIR=/lscratch/$SLURM_JOBID; bash get.table.snp.forge2.2.sh $1
#fi

cd /lscratch/$SLURM_JOBID

#cp summary.table.$1.html ..

#rm summary.table.$1.html


#echo "https://forge2.altiusinstitute.org/files/summary.table."$1".html"

#cat img.txt minibr.txt splitset-8.txt.html minibr.txt minibr.txt splitset-9.txt.html minibr.txt splitset-10.txt.html minibr.txt splitset-7.txt.html minibr.txt splitset-1.txt.html minibr.txt splitset-2.txt.html minibr.txt splitset-3.txt.html minibr.txt splitset-4.txt.html minibr.txt splitset-5.txt.html minibr.txt splitset-6.txt.html minibr.txt | sed 's/Feta /Fetal /g' | sed 's/0x[a-zA-Z0-9]*\//<a href=\"https:\/\/forge2.altiusinstitute.org\"\/>FORGE2<\/a> /g'|sed 's/0x[a-zA-Z0-9]*\///g' > summary.table.$1.html


#rm table.count*

#rm minibr.txt tableToHtml.pl

#rm splitset-*txt
#rm splitset-*html



#rm -rf $1 ../summary.table.$1

#perl tableToHtml.pl < summary.table.$1 > summary.table.$1.html

#echo "CATO coordinates are in hg19" >> summary.table.$1

#zcat 0x*/$1*.gz |tail -15 >> summary.table.$1

