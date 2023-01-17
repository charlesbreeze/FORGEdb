# FORGEdb
Tutorial: generating 37M FORGEdb pages 

This tutorial lists the different steps to generate 37 million webpages using FORGEdb standalone. The code presented here is optimised for the NIH Biowulf cluster, but can be amended to work on other large compute clusters. Note that different required FORGEdb databases can be downloaded from https://forgedb.altiusinstitute.org/?download. 

PREPARATORY STEPS:

The "file.commands.regenerate.forgedb.htmls.sh.using.local.memory.sh" is a file containing all forgedb standalone commands for each SNP, one per line. There are 37964763 SNPs

	[breezece@biowulf breezece]$ wc file.commands.regenerate.forgedb.htmls.sh.using.local.memory.sh &
	[3] 68788
	[breezece@biowulf breezece]$ 37964763 file.commands.regenerate.forgedb.htmls.sh.using.local.memory.sh

To split into 312 nodes:

	[breezece@biowulf breezece]$ bc 
	37964763/312
	121681.93269230769230769230
	quit

We then split the file:

	[breezece@biowulf breezece]$ split -121682 file.commands.regenerate.forgedb.htmls.sh.using.local.memory.sh test.localmem2/for.job &
	[3] 69391
	breezece@biowulf breezece]$ 
	[3]+  Done                    split -121682 file.commands.regenerate.forgedb.htmls.sh.using.local.memory.sh test.localmem2/for.job

List all the files made:

	[breezece@biowulf breezece]$ ls test.localmem2/ > ls.localmem
	[breezece@biowulf breezece]$ wc ls.localmem
	312 ls.localmem

Modify myjobscript template to include the real files from each directory:

	[breezece@biowulf breezece]$ cat ls.localmem |awk '{print "cat myjobscript|sed \"s/test.sh/"$0"/g\" > test.localmem2/myjobscript."$0}'|bash


RUN with (from directory /data/breezece/test.localmem2):

	sbatch --cpus-per-task=32 --mem=120g --gres=lscratch:700 --constraint=x2650 --partition norm --time=9-03:08:42 myjobscript.for.jobaa

	(to run all at once the commands are): 

	ls myjob* |grep -v head|awk '{print "sbatch --cpus-per-task=32 --mem=120g --gres=lscratch:700 --constraint=x2650 --partition norm --time=9-03:08:42 "$0}' > submits.jobs.sbatch.sh
	
	[test.localmem2]$ nohup bash submits.jobs.sbatch.sh & 

this runs:

	
	[breezece@biowulf breezece]$ cat test.localmem7/myjobscript.for.jobaa
	
	#!/bin/bash
	
	mkdir /lscratch/$SLURM_JOBID/done.htmls
	
	cp for.jobaa /lscratch/$SLURM_JOBID
	
	cd /lscratch/$SLURM_JOBID
	
	cp /data/breezece/tableToHtml2.pl /lscratch/$SLURM_JOBID
	cp /data/breezece/miniabout.txt /lscratch/$SLURM_JOBID
	cp /data/breezece/img.txt /lscratch/$SLURM_JOBID
	cp /data/breezece/tableToHtml.pl /lscratch/$SLURM_JOBID
	cp /data/breezece/minibr.txt /lscratch/$SLURM_JOBID
	cp /data/breezece/get.table.snp.forge2.2.sh /lscratch/$SLURM_JOBID
	
	cp /data/breezece/abc.db /lscratch/$SLURM_JOBID 
	cp /data/breezece/cato.db /lscratch/$SLURM_JOBID 
	cp /data/breezece/closest.refseq.gene.hg19.db /lscratch/$SLURM_JOBID 
	cp /data/breezece/eqtlgen.db /lscratch/$SLURM_JOBID 
	cp /data/breezece/forge_2.0.db /lscratch/$SLURM_JOBID 
	cp /data/breezece/forge2.tf.fimo.jaspar.1e-5.taipale.1e-5.taipaledimer.1e-5.uniprobe.1e-5.xfac.1e-5.db /lscratch/$SLURM_JOBID 
	cp /data/breezece/gtex.eqtl.db /lscratch/$SLURM_JOBID
	cp /data/breezece/dmp_GWAS_bins /lscratch/$SLURM_JOBID
	cp /data/breezece/dmp_GWAS_params /lscratch/$SLURM_JOBID
	
	cp /data/breezece/snps.to.avoid /lscratch/$SLURM_JOBID
	cp /data/breezece/scores.py /lscratch/$SLURM_JOBID
	
	export SLURM_JOBID
	
	# for.jobaa should run, say, 100 commands
	#cat for.jobaa |/usr/local/apps/parallel/20220422/bin/parallel -j $SLURM_CPUS_PER_TASK
	cat for.jobaa |/usr/local/apps/parallel/20220422/bin/parallel -j 16
	
	# after all the parallel jobs have finished
	
	tar -cvz -f $SLURM_JOBID.done.htmls.tar.gz /lscratch/$SLURM_JOBID/done.htmls
	
	cp $SLURM_JOBID.done.htmls.tar.gz /data/breezece/done.htmls

and for.jobaa

	[breezece@biowulf breezece]$ head test.localmem7/for.jobaa
	TMPDIR=/lscratch/$SLURM_JOBID; bash get.table.snp.forge2.2.sh rs10000242
	TMPDIR=/lscratch/$SLURM_JOBID; bash get.table.snp.forge2.2.sh rs10000384

After jobs are complete, the full list of 37M webpages (in compressed format) can be viewed at the "done.htmls" directory.
 
