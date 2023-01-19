#!/usr/bin/env perl

=head1 NAME

eforge.pl - Experimentally derived Functional element Overlap analysis of ReGions from EWAS.

=head1 SYNOPSIS

eforge.pl options (-f file) (-dmp dmplist)

=head1 DESCRIPTION

Analyse a set of DMPs for their overlap with DNase 1 hotspots compared to matched background DMPs. 
Identifies enrichment in DHS by tissue and plots graphs and table to display. Arbitrarily a minumum of 5* DMPs is required.  
Note that if no DMPs are given the script will run on A DEFAULT EWAS* as an example output.

Several outputs are made.

A straight base R graphics pdf chart of the data.

A polychart (https://github.com/Polychart/polychart2) interactive javascript graphic using rCharts (http://ramnathv.github.io/rCharts/).

A dimple (http://dimplejs.org) d3 interactive graphic using rCharts.

A table using the Datatables (https://datatables.net) plug-in for the jQuery Javascript library, again accessed through rCharts.

In each of the graphics the colouring should be consistent. Blue (p value > 0.05), light red or pink (0.05 => p value > 0.01), red or dark red (p value <= 0.01 ) for the 95% and 99% cIs. 
Or whatever other thresholds are specified. 

eForge functions, plotting options and stats are provided by eForge::eForge, eForge::ePlot and eForge::eStats modules.

=head1 OPTIONS

=over

=item B<--dataset TAG>

Set of functional data to look for enrichment. Either ENCODE data ('encode'), unconsolidated Roadmap
Epigenome data ('erc'), consolidated Roadmap Epigenome data ('erc2'), or Blueprint data ('blueprint').
erc by default.

Use --dataset ? to get a list of available datasets on your local install.

=item B<--array TAG>

Array (FKA background) is set at default to 450k array ('450k'), the Illumina Infinium HumanMethylation450 BeadChip.

For the time being, it is suficient for DMPs to be on the 450k array. Probes within 1kb of each other
will undergo filtering.

Use --array ? to get a list of available backgrounds on your local install.

=item B<--label STRING>

Supply a label that you want to use for the plotting titles, and filenames.

=item B<--f FILENAME>

Supply the name of a file containing a list of DMPs. 
Format must be given by the -format flag. 
If not supplied the analysis is performed either on dmps provided as probeids (cg or ch probes) in a
comma separated list through the dmps option or on a set of data from a default ewas study, namely a
set of monocyte tDMPs from Jaffe AE and Irizarry RA, Genome Biol 2014.

Note that at least 5 DMPs are required at a minimum by default.

=item B<--dmps probe_id,probe_id...>

Can provide the dmps as probeids in a comma separated list.

=item B<--min_dmps INT>

Specify the minimum number of DMPs to be allowed. Default is 5 now we are using binomial test.

=item B<--thresh FLOAT,FLOAT>

Alter the default binomial p value thresholds. Give a comma separate list of three e.g. 0.05,0.01 for the defaults

=item B<--format STRING>

If f is specified, specify the file format as follow:

probeid = list of dmps as probeids each on a separate line. Optionally can add other fields after the probeid which are ignored,
unless the pvalue filter is specified, in which case eForge assumes that the second field is the minus log10 pvalue

bed  = File given is a bed file of locations (chr\tbeg\tend).  bed format should be 0 based and the chromosome should be given as chrN.
However we will also accept chomosomes as just N (ensembl) and 1-based format where beg and end are the same*.

tabix = File contains DMPs in tabix format.

=item B<--filter FLOAT>

Set a filter on the DMPs based on the -log10 pvalue.  This works for files in the probeid' format.
Give a value as the lower threshold and only DMPs with -log10 pvalues >= to the threshold will be
analysed. Default is no filtering.

=item B<--save_stats>

Output annotation stats for the original and the random picks.

=item B<--reps>

The number of background matching sets to pick and analyse. Default 1000.

=item B<--proxy TAG>

Apply filter for DMPs in proximity (within 1 kb of another test DMP). With proximity filter specified,
eForge will report DMPs removed due to proximity with another DMP in the list and will randomly pick
one of the probes among the set of probes that are in proximity (within 1 kb of each other).

At the moment, this is a dummy flag as only one proximity filter is available for each array. It
will become useful if the database and code support more than one. At the moment to turn off
proximity filtering, simply specify -noproxy

=item B<--cato>                                                                                                                              
Show CATO score for individual SNP.

=item B<--eqtl>                                                                                                                                
Show GTEx eQTL for individual SNP.                                                                                                             

=item B<--abc>                                                                                                                                 
Show Activity-By-Contact (ABC) data and target gene for individual SNP.                                                                                    
=item B<--gene>                                                                                                                                   Show RefSeq closest gene for individual SNP. Also generates UCSC genome browser links.
      
=item B<--tf>                                                                                                                                     Show FORGE2-TF motif data from jaspar/taipale/taipaledimer/uniprobe/xfac at FIMO p<1e-5 for individual SNP.

=item B<--eqtlgen>                                                                                                                                Show eQTLGen eQTL for individual SNP.                   

=item B<--noproxy>

Turn off proximity filtering.

=item B<--depletion>

Analyse for depletion pattern instead of the default enrichment analysis. Use when dealing with
datasets suspected not to overlap with DHS (or the relevant functional assay). Specifying depletion
will be indicated on the label (the text "Depletion Analysis" will be added to the file label).

=item B<--noplot>

Just make the data file, don't plot.

=item B<--help|-h|-?>

Print a brief help message and exits.

=item B<--man|-m>

Print this perldoc and exit.

=back

=head1 LICENCE AND COPYRIGHT

eforge.pl Functional analysis of EWAS DMPs

Copyright (C) [2014-2015] EMBL - European Bioinformatics Institute and University College London

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; version 2 dated June, 1991 or at your option
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available in the source tree;
if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

=head1 CONTACT

Charles Breeze, C<< <c.breeze at ucl.ac.uk> >>

Javier Herrero, C<< <javier.herrero at ucl.ac.uk> >>

=head1 ACKNOWLEDGEMENTS

This software is based on the FORGE tool developed by Ian Dunham at the EMBL-EBI

Javier Herrero <javier.herrero@ucl.ac.uk>

=cut

use strict;
use 5.010;
use warnings;
use DBD::SQLite;
use DBI; #database link to sqlite database
use Sort::Naturally;
use Cwd;
use Getopt::Long; #check this module
use File::Basename;
use Config::IniFiles;
use Pod::Usage;
use Scalar::Util qw(looks_like_number);
# use eForge::eStats;
use eForge::ePlot;
use eForge::eForge;
use Forge2::Stats2;
use lib '/var/www/forge2/production/src/bin';
use Data::UUID;
use Statistics::Multtest qw(BH);
use Data::Dumper;

my $cwd = getcwd;

#my $dbname = "forge_2.0.2.db";
my $dbname = "forge_2.0.db";

my $array; # Default value
my $array_label;
my $format = 'probeid'; # Input format
my $label = 'Unnamed'; # Label for plots
my $reps = 1;
# set binomial p values, multiple test correction is used
my $thresh; # string for command line option
my $t_marginal = 0.05; # default marginal p-value threshold
my $t_strict = 0.01; # default strict p-value threshold

my $min_num_probes = 1; # the minimum number of probes allowed for test. Set to 5 as we have binomial p

my ($dataset, $filename, $save_probe_annotation_stats, $noplot,
    $help, $man, $proxy, $noproxy, $depletion, $filter, $out_dir, $probe_list,
    $web, $autoopen, $cato, $eqtl, $eqtlgen, $abc, $tf, $gene);

GetOptions (
    'dataset=s'  => \$dataset,
    'save_stats|bkgrd' => \$save_probe_annotation_stats,
    'array|bkgd=s' => \$array,
    'label=s'    => \$label,
    'f=s'        => \$filename,
    'format=s'   => \$format,
    'probes|dmps=s@' => \$probe_list,
    'min_num_probes|min_dmps=i' => \$min_num_probes,
    'noplot'     => \$noplot,
    'reps=i'     => \$reps,
    'thresh=s'   => \$thresh,
    'proxy=s'    => \$proxy,
    'noproxy'    => \$noproxy,
    'cato'       => \$cato,
    'gene'       => \$gene,
    'eqtl'       => \$eqtl,
    'eqtlgen'    => \$eqtlgen,
    'abc'        => \$abc,
    'tf'        => \$tf,
    'depletion'  => \$depletion,
    'filter=f'   => \$filter,
    'out_dir=s'  => \$out_dir,
    'web=s'      => \$web,
    'autoopen'   => \$autoopen,
    'help|h|?'   => \$help,
    'man|m'      => \$man,

);


pod2usage(1) if ($help);
pod2usage(-verbose => 2) if ($man);

if (!$out_dir) {
    my $ug = new Data::UUID;
    $out_dir = $ug->to_hexstring($ug->create());
}
mkdir $out_dir;


# Define the thresholds to use.
if ($thresh) {
    ($t_marginal, $t_strict) = parse_pvalue_thresholds($thresh);
}


## ============================================================================
## Connect to the DB
## ============================================================================
# This reads the config file and sets up the $datadir variable
#commented on 30-Aug-2022:
#my $dirname = dirname(__FILE__);
#my $cfg = Config::IniFiles->new( -file => "$dirname/eforge.ini" );
#my $datadir = $cfg->val('Files', 'datadir');

my $datadir = "/lscratch/$ENV{SLURM_JOBID}";

##warn "[".scalar(localtime())."] DBD::SQLite plugin version: $DBD::SQLite::sqlite_version\n";

unless (-s "$datadir/$dbname") {
    die "Database $dbname not found or empty [$datadir/$dbname]";
}
my $dsn = "dbi:SQLite:dbname=$datadir/$dbname";
my $dbh = DBI->connect($dsn, "", "") or die $DBI::errstr;



## ============================================================================


## ============================================================================
## Check the dataset against the info on the DB
## ============================================================================
my $all_datasets = get_all_datasets($dbh);
if (!defined($all_datasets)) {
    die "Empty database: no dataset loaded!\n";
} elsif (!defined($dataset)) {
    $dataset = $all_datasets->[0]->{tag};
    #print "Using default dataset: [$dataset] ".$all_datasets->[0]->{name}."\n";
} elsif ($dataset eq "?") {
    print "Available datasets:\n - [".join("\n - [", map {$_->{tag}."] ".$_->{name}} @$all_datasets)."\n";
    exit();
} elsif (!grep {$_ eq $dataset} map {$_->{tag}} @$all_datasets) {
    die "Dataset $dataset unknown\nAvailable datasets:\n - [".join("\n - [", map {$_->{tag}."] ".$_->{name}} @$all_datasets)."\n";
}
## ============================================================================


## ============================================================================
## Check the array name (A.K.A. background) against DB
## ============================================================================
my $all_arrays = get_all_arrays($dbh);
if (!defined($all_arrays)) {
    die "Empty database: no background loaded!\n";
} elsif (!defined($array)) {
    $array = $all_arrays->[0]->{tag};
    print "Using default background: [$array] ".$all_arrays->[0]->{name}."\n";
    $array_label = $all_arrays->[0]->{name};
} elsif ($array eq "?") {
    print "Available arrays:\n - [".join("\n - [", map {$_->{tag}."] ".$_->{name}} @$all_arrays)."\n";
    exit();
} elsif (!grep {$_ eq $array} map {$_->{tag}} @$all_arrays) {
    die "Array $array unknown\nAvailable arrays:\n - [".join("\n - [", map {$_->{tag}."] ".$_->{name}} @$all_arrays)."\n";
} else {
    foreach my $this_array (@$all_arrays) {
        if ($this_array->{tag} eq $array) {
            $array_label = $this_array->{name};
            last;
        }
    }
}
## ============================================================================


## ============================================================================
## Check the proxy_filter (A.K.A. filter) against DB
## ============================================================================
# Set proximity filter
if (defined $noproxy) {
    $proxy = undef;
} else {
    my $all_proxy_filters = get_all_proximity_filters($dbh);
    if ($all_proxy_filters->{$array}) {
        $proxy = $all_proxy_filters->{$array};
    }
}
## ============================================================================


## ============================================================================
## Append main options (depletion on/off; array; dataset) to $label
## ============================================================================
if (defined $depletion) {
    $label = "$label.depletion";
}
(my $lab = $label) =~ s/\s/_/g; # Avoid whitespaces on the label
$lab = "$lab.$array.$dataset";
## ============================================================================


## ============================================================================
## Read and process the input DMPs
## ============================================================================
##warn "[".scalar(localtime())."] Processing input...\n";
# This will read the probes from the file if provided, from the probe list otherwise or use the
# example data set as a last resort.
my $dmps = get_input_probes($filename, $probe_list);
my $original_dmps = [@$dmps];
my $num_of_input_dmps = scalar(@$dmps);

my ($proximity_excluded);
if (0) {
    # eForge proximity filter disabled in favor of LD filter
    # Apply the proximity filter if requested
    if (defined $proxy) {
        ##warn "Applying LD filter... ($proxy)\n";
        ($proximity_excluded, $dmps) = proximity_filter($dbh, $array, $dmps);
        while (my ($excluded_dmp, $dmp) = each %$proximity_excluded) {          
         ##warn "$excluded_dmp excluded for $proxy LD filter with $dmp\n";
        }
    } else {
        ##warn "Not applying LD filter...\n";
    }
} 

if (defined $noproxy) {
    ##warn "Not applying LD filter...\n";
} else {
    # LD method from earlier Forge2.pm/forge2.db
    my $LD_threshold = "r08";
    ##warn "Applying LD filter... ($LD_threshold)\n";
    ($proximity_excluded, $dmps) = ld_filter($dmps, $LD_threshold, $dbh);
    # Note:  sub ld_filter already prints the excluded SNPs
}

#else {
#    # LD method from earlier Forge2.pm/forge2.db
#    my $LD_threshold = "r08";
#    warn "Applying LD filter... ($LD_threshold)\n";
#    ($proximity_excluded, $dmps) = ld_filter($dmps, $LD_threshold, $dbh);
#    # Note:  sub ld_filter already prints the excluded SNPs
#}

# $annotated_probes is an arrayref with probe_id, sum, bit, gene_group, cgi_group for each input probe
my $annotated_probes = get_probe_annotations_and_overlap_for_dataset($dbh, $dataset, $array, $dmps);
my $existing_probes = {map {$_->[0] => 1} @$annotated_probes};
$dmps = [keys %$existing_probes];

## Detect and remove the missing probes.
my $num_missing_probes = find_missing_probes($original_dmps, $existing_probes, $proximity_excluded);

# Print summary of filtering and checks:
my $msg = "For $label, $num_of_input_dmps SNPs provided, ". scalar @$dmps.
        " retained: $num_missing_probes were not found";
if (defined $proxy) {
    $msg .= " and " . scalar(keys %$proximity_excluded) . " excluded using $proxy filter";
}
##warn $msg, ".\n";

# Check we have enough DMPs left
my $num_of_valid_probes = scalar @$dmps;
if ($num_of_valid_probes < $min_num_probes) {
    die "Fewer than $min_num_probes DMPs. Analysis not run\n";
}
## ============================================================================


# get the cell list array and the hash that connects the cells and tissues
# $samples is a hash whose keys are the $cells (short name for the cell type/lines) and value is
# another hash with 'tissue', 'datatype', 'file' and 'acc' keys.
# IMPORTANT: $cells contains the list of cells in the order defined in the DB. This is critical
# to correctly assign each bit to the right sample.
my ($cells, $samples) = get_samples_from_dataset($dbh, $dataset);

#print Dumper $cells;
#exit -1;

# unpack the bitstrings and store the overlaps by cell.
# $overlaps is a complex hash like:
# $overlaps->{'DMPS'}->{$probe_id}->{'SUM'} (total number of overlaps of this probe with features in this dataset)
# $overlaps->{'DMPS'}->{$probe_id}->{'PARAMS'} (gene and CGI annotations for this probe)
# $overlaps->{'CELLS'}->{$cell}->{'COUNT'} (number of input DMPs that overlap with the signal on this sample)
# $overlaps->{'CELLS'}->{$cell}->{'DMPS'} (list of input DMPs that overlap with the signal on this sample)
my $overlaps = process_overlaps($annotated_probes, $cells, $dataset);

# generate stats on the background selection
if (defined $save_probe_annotation_stats) {
    save_probe_annotation_stats($overlaps, $out_dir, $lab, "test");
}



# only pick background dmps matching dmps that had bitstrings originally.
#reference to hash key 'DMPS' is due to use of eforge.pm module from eForge tool
#(in subroutines process_overlaps, etc)


# Identify the MAF, TSS distance and GC content relationship, and then make random picks
##warn "[".scalar(localtime())."] Loading the $array background...\n";
my $random_picks = get_random_matching_picks($overlaps, $array, $datadir, $reps);

#print Dumper $random_picks;
#these lines below can be uncommented to get random picks
##print Dumper @$random_picks;
##exit -1;

########check below lines:
 
# for bkgrd set need to get distribution of counts instead
# make a hash of data -> cell -> bkgrd-Set -> overlap counts
my %overlaps_per_cell; #this hash is going to store the overlaps for the random picks, per cell

# Get the bits for the background sets and process
my $total_num_probes_in_random_picks;

##warn "[".scalar(localtime())."] Running the analysis with $num_of_valid_probes DMPs...\n";
my $count = 0;

#
# Intialize overlaps for cells, regardless of whether overlaps are discovered in random picks, or not.
# This may help prevent downstream statistical tests from failing by divide-by-zero errors.
# 8 Feb 2018, APR
#
foreach my $cell (@$cells) {
    push @{$overlaps_per_cell{$cell}}, 0;
}

foreach my $this_random_pick (@$random_picks) {
    ##warn "[".scalar(localtime())."] Repetition $count out of ".$reps."\n" if (++$count%100 == 0);
    $annotated_probes = get_probe_annotations_and_overlap_for_dataset($dbh, $dataset, $array, $this_random_pick);

    $total_num_probes_in_random_picks += scalar @$annotated_probes;

    unless (scalar @$annotated_probes == $num_of_valid_probes) {
        ##warn "Random pick #$count only has " . scalar @$annotated_probes . " probes compared to $num_of_valid_probes in the input set.\n";
    }

    my $this_pick_overlaps = process_overlaps($annotated_probes, $cells, $dataset);

    # accumulate the overlap counts by cell
    foreach my $cell (keys %{$this_pick_overlaps->{'CELLS'}}) {
        push @{$overlaps_per_cell{$cell}}, $this_pick_overlaps->{'CELLS'}->{$cell}->{'COUNT'}; 
        #print Dumper "cell: [ $cell ] overlaps: $this_pick_overlaps->{'CELLS'}->{$cell}->{'COUNT'}";
    }

    if (defined $save_probe_annotation_stats) {
        save_probe_annotation_stats($this_pick_overlaps, $out_dir, $lab, $count);
    }
}


#my $result_chr=fetch_all_chr_pos($dbh_2,$dmps);
#print Dumper $result_chr;

$dbh->disconnect();

##warn "[".scalar(localtime())."] All repetitions done.\n";

##warn "[".scalar(localtime())."] Calculating p-values...\n";
#Having got the test overlaps and the bkgd overlaps now calculate p values and output 
#the table to be read into R for plotting.


if (!$web) {
    open(BACKGROUND, "| gzip -9 > $out_dir/background.tsv.gz") or die "Cannot open background.tsv";
}



my @results;
my @pvalues;
###ncmp is a function from Sort::Naturally
foreach my $cell (sort {ncmp($$samples{$a}{'tissue'},$$samples{$b}{'tissue'}) || ncmp($a,$b)} @$cells){
    # above line sorts by the tissues alphabetically (from $samples hash values)

    # ultimately want a data frame of names(results)<-c("Zscore", "Cell", "Tissue", "File", "DMPs")
    if (!$web) {
        print BACKGROUND join("\t", @{$overlaps_per_cell{$cell}}), "\n";
    }
    my $teststat = ($overlaps->{'CELLS'}->{$cell}->{'COUNT'} or 0); #number of overlaps for the test DMPs

    # binomial pvalue, probability of success is derived from the background overlaps over the tests for this cell
    # $backdmps is the total number of background dmps analysed
    # $tests is the number of overlaps found over all the background tests
    my $total_num_overlaps_in_random_picks = 1;
    foreach (@{$overlaps_per_cell{$cell}}) {
        $total_num_overlaps_in_random_picks += $_;
    }
    my $p = sprintf("%.6f", $total_num_overlaps_in_random_picks / $total_num_probes_in_random_picks);

    # binomial probability for $teststat or more hits out of $dmpcount dmps
    # sum the binomial for each k out of n above $teststat
    my $pbinom;
    if (defined $depletion) {
        foreach my $k (0 .. $teststat) {
           # $pbinom += binomial($k, $num_of_valid_probes, $p);
        }
    } else {
        foreach my $k ($teststat .. $num_of_valid_probes) {
           # $pbinom += binomial($k, $num_of_valid_probes, $p);
        }
    }
#    if ($pbinom >1) {
    if (1){
        $pbinom=1;
    }
    # Store the p-values in natural scale (i.e. before log transformation) for FDR correction
    push(@pvalues, $pbinom);
    $pbinom = sprintf("%.2e", $pbinom);

    # Z score calculation (note: this is here only for legacy reasons. Z-scores assume normal distribution)
    
    if (! @{$overlaps_per_cell{$cell}}) {
      print Dumper "undefined overlaps_per_cell for cell $cell";
    }
    
    my $zscore = zscore($teststat, $overlaps_per_cell{$cell});

    my $dmp_string = "";
    $dmp_string = join(",", @{$overlaps->{'CELLS'}->{$cell}->{'DMPS'}})
            if defined $overlaps->{'CELLS'}->{$cell}->{'DMPS'};
    # This gives the list of overlapping DMPs for use in the tooltips. If there are a lot of them this can be a little useless
    my ($shortcell, undef) = split('\|', $cell); # undo the concatenation from earlier to deal with identical cell names.
    push(@results, [$shortcell, $$samples{$cell}{'tissue'}, $$samples{$cell}{'datatype'}, $dmp_string]);
#     push(@results, [$$samples{$cell}{'datatype'}, $shortcell, $$samples{$cell}{'tissue'}, $dmp_string}]);
    
#    push(@results, [$shortcell, $$samples{$cell}{'tissue'}, $$samples{$cell}{'datatype'}, $$samples{$cell}{'file'}, $dmp_string, $$samples{$cell}{'acc'}]);
}
if (!$web) {
    close(BACKGROUND);
}

## ============================================================================
## Correct the p-values for multiple testing using the Benjamini-Yekutieli FDR control method
## ============================================================================
my $qvalues = BH(\@pvalues);
$qvalues = [map {sprintf("%.2e", $_)} @$qvalues];
## ============================================================================


## ============================================================================
## Write the results to a tab-separated file
## ============================================================================
my $results_filename = "$lab.chart.tsv.gz";
#open(TSV, "| gzip -9 > $out_dir/$results_filename") or die "Cannot open $out_dir/$results_filename: $!";
open(TSV, " > $out_dir/$results_filename") or die "Cannot open $out_dir/$results_filename: $!";
#print TSV join("\t", "Cell", "Tissue", "Datatype", "File", "Variant", "Accession"), "\n";
print TSV join("\t", "<b>Cell</b>", "<b>Tissue</b>","<b>Datatype</b>","<b>Variant rsid</b>"), "\n";
for (my $i = 0; $i < @results; $i++) {

    foreach my $probe_id (keys %{$overlaps->{'DMPS'}}) {

	if ( grep( /^$probe_id/, @{$results[$i]} ) ) {
	    #if ( grep( /klsjfkldjsklfadjkfjsdkfjsdkjfd/, @{$results[$i]} ) ) {
	    #}else{
		print TSV join("\t", @{$results[$i]}), "\n";
	    #}
	}

    }
    
 
}

###get cato score
my $dbname2 = "cato.db";
my $dsn_2 = "dbi:SQLite:dbname=$datadir/$dbname2";
my $dbh_2 = DBI->connect($dsn_2, "", "") or die $DBI::errstr;


if (defined $cato) {
    my $result_chr=fetch_all_chr_pos($dbh_2,$dmps);
    print TSV "\n<b> <a href=\"https://forge2.altiusinstitute.org/files/cato.html\">CATO</a> (contextual analysis of TF occupancy) score (coordinates in hg19, for more details, see <a href=\"http://www.mauranolab.org/CATO/\">http://www.mauranolab.org/CATO/</a>):\n</b>";
    print TSV "<b>chr</b>\t<b>start</b>\t<b>end</b>\t<b>rsid</b>\t<b>CATO score</b>\t<b>strand</b>\t<b>TF motif</b>\t<b>score1</b>\t<b>allele1</b>\t<b>allele2</b>\t<b>score2</b>\t<b>biosample</b>\n";
    #$result_chr =~ y/\n//d;
    print TSV join("\t", @$result_chr) . "\n";    
    #print TSV Dumper $result_chr;
}

$dbh_2->disconnect();

###get eqtls

my $dbname3 = "gtex.eqtl.db";
my $dsn_3 = "dbi:SQLite:dbname=$datadir/$dbname3";
my $dbh_3 = DBI->connect($dsn_3, "", "") or die $DBI::errstr;

if (defined $eqtl) {
    print TSV "\n<b>GTEx cis-eQTL information (for more details, see <a href=\"https://gtexportal.org/home/\">https://gtexportal.org/home/</a>):\n</b>";
    print TSV "<b>gene symbol</b>\t<b>gene ID</b>\t<b>variant ID</b>\t<b>gene ID2</b>\t<b>TSS distance</b>\t<b>reference allele</b>\t<b>alternative allele</b>\t<b>MA samples</b>\t<b>MA count</b>\t<b>MAF</b>\t<b>p-value nominal</b>\t<b>slope</b>\t<b>slope se</b>\t<b>p-value nominal threshold</b>\t<b>min p-value nominal</b>\t<b>p-value beta</b>\t<b>biosample</b>\t<b>chr (hg38)</b>\t<b>variant position (hg38) </b>\t<b>num alternative alleles per site</b>\t<b>rs ID dbSNP151 GRCh38p7</b>\t<b>variant ID b37</b>\n";
    my $result_chr2=fetch_all_chr_pos_eqtls($dbh_3,$dmps);
#    print TSV join("\t", @$result_chr2) . "\n";
#    for (my $i = 0; $i < @result_chr2; $i++) {
#	print TSV join("\t", @{$result_chr2[$i]}), "\n";
#    }
}

sub fetch_all_chr_pos_eqtls {
    my ($dbh_3, $probe_ids) = @_;
    my $sth = $dbh_3->prepare("SELECT * FROM eqtl WHERE rs_id_dbSNP151_GRCh38p7 = ?");

    my $locations = [];

    foreach my $probe_id (@$probe_ids) {
	$sth->execute($probe_id);
	my $result_chr = $sth->fetchall_arrayref();
	my $location;
	foreach my $row (@{$result_chr}) {
	    push(@$locations,$$row[21],"\t",$$row[0],"\t",$$row[1],"\t",$$row[2],"\t",$$row[3],"\t",$$row[16],"\t",$$row[17],"\t",$$row[4],"\t",$$row[5],"\t",$$row[6],"\t",$$row[7],"\t",$$row[8],"\t",$$row[9],"\t",$$row[10],"\t",$$row[11],"\t",$$row[12],"\t",$$row[13],"\t",$$row[14],"\t",$$row[15],"\t",$$row[18],"\t",$$row[19],"\t",$$row[20],"\n");

#	    print TSV join($$row[0],"\t",$$row[1],"\t",$$row[2],"\t",$$row[3],"\t",$$row[4],"\t",$$row[5],"\t",$$row[6],"\t",$$row[7],"\t",$$row[8],"\t",$$row[9],"\t",$$row[10],"\t",$$row[11],"\t",$$row[12],"\t",$$row[13],"\t",$$row[14],"\t",$$row[15],"\t",$$row[16],"\t",$$row[17],"\t",$$row[18],"\t",$$row[19],"\n");
		
#	    foreach(@$locations)
#	    {
#		#s/txt/txt\n/g;
#	    }
	    
#	    push(@$locations, $$row[0],$$row[1],$$row[2],$$row[3],$$row[4],$$row[5],$$row[6],$$row[7],$$row[8],$$row[9],$$row[10],$$row[11],$$row[12],$$row[13],$$row[14],$$row[15],$$row[16],$$row[17],$$row[18],$$row[19]);

	}
	    print TSV @$locations,"\n" ;
    }
    $sth->finish();
    return $locations;
}

$dbh_3->disconnect();


###get abc

my $dbname4 = "abc.db";
my $dsn_4 = "dbi:SQLite:dbname=$datadir/$dbname4";
my $dbh_4 = DBI->connect($dsn_4, "", "") or die $DBI::errstr;

if (defined $abc) {
    print TSV "<b>Activity-By-Contact (ABC) information (hg19, for more details, see <a href=\"https://www.engreitzlab.org/resources\">https://www.engreitzlab.org/resources</a>):\n</b>";
    print TSV "<b>SNP chromosome</b>\t<b>SNP start</b>\t<b>SNP end</b>\t<b>rsid</b>\t<b>bait chromosome</b>\t<b>bait start</b>\t<b>bait end</b>\t<b>name</b>\t<b>class</b>\t<b>activity base</b>\t<b>TargetGene</b>\t<b>TargetGeneTSS</b>\t<b>TargetGeneExpression</b>\t<b>TargetGenePromoterActivityQuantile</b>\t<b>TargetGeneIsExpressed</b>\t<b>distance</b>\t<b>isSelfPromoter</b>\t<b>Hi-C contact</b>\t<b>powerlaw contact</b>\t<b>powerlaw contact reference</b>\t<b>Hi-C contact pl scaled</b>\t<b>Hi-C pseudocount</b>\t<b>Hi-C contact pl scaled adj</b>\t<b>ABC Score Numerator</b>\t<b>ABC Score</b>\t<b>powerlaw Score Numerator</b>\t<b>powerlaw Score</b>\t<b>CellType</b>\t<b>base overlap</b>\n";
    my @result_chr4=fetch_all_chr_pos_abc($dbh_4,$dmps);

    #print TSV join("\t", @$result_chr4) . "\n";
#    for (my $i = 0; $i < @result_chr4; $i++) {
#	print TSV join("\t", @{$result_chr4[$i]}), "\n";
#    }
}

sub fetch_all_chr_pos_abc {
    my ($dbh_4, $probe_ids) = @_;
    my $sth = $dbh_4->prepare("SELECT * FROM abc WHERE rsid = ?");

    my $locations = [];

    foreach my $probe_id (@$probe_ids) {
	$sth->execute($probe_id);
	my $result_chr = $sth->fetchall_arrayref();
	my $location;
	foreach my $row (@{$result_chr}) {
	    push(@$locations, $$row[0],"\t",$$row[1],"\t",$$row[2],"\t",$$row[3],"\t",$$row[4],"\t",$$row[5],"\t",$$row[6],"\t",$$row[7],"\t",$$row[8],"\t",$$row[9],"\t",$$row[10],"\t",$$row[11],"\t",$$row[12],"\t",$$row[13],"\t",$$row[14],"\t",$$row[15],"\t",$$row[16],"\t",$$row[17],"\t",$$row[18],"\t",$$row[19],"\t",$$row[20],"\t",$$row[21],"\t",$$row[22],"\t",$$row[23],"\t",$$row[24],"\t",$$row[25],"\t",$$row[26],"\t",$$row[27],"\t",$$row[28],"\n");
	    #push(@$locations, $$row[0],$$row[1],$$row[2],$$row[3],$$row[4],$$row[5],$$row[6],$$row[7],$$row[8],$$row[9],$$row[10],$$row[11],$$row[12],$$row[13],$$row[14],$$row[15],$$row[16],$$row[17],$$row[18],$$row[19],$$row[20],$$row[21],$$row[22],$$row[23],$$row[24],$$row[25],$$row[26],$$row[27],$$row[28]);
#	    print TSV join("\t", @$locations) . "\n"; 
#	    for (my $i = 0; $i < @$locations; $i++) {
#		print TSV join("\t", @{$locations[$i]}), "\n";
#	    }

	}
	    print TSV @$locations, "\n";
    }
    $sth->finish();
    return $locations;
}

$dbh_4->disconnect();

###get tf

my $dbname5 = "forge2.tf.fimo.jaspar.1e-5.taipale.1e-5.taipaledimer.1e-5.uniprobe.1e-5.xfac.1e-5.db";
my $dsn_5 = "dbi:SQLite:dbname=$datadir/$dbname5";
my $dbh_5 = DBI->connect($dsn_5, "", "") or die $DBI::errstr;

if (defined $tf) {
    print TSV "<b>FORGE2-TF motif information (hg19, for more details, see <a href=\"https://forge2-tf.altiusinstitute.org/\">FORGE2-TF</a>, <a href=\"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6853678/\">eFORGE-TF</a> and the <a href=\"https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02560-3\">FORGE2 main paper</a>):\n</b>";
    print TSV "<b>chr</b>\t<b>start</b>\t<b>end</b>\t<b>rsid</b>\t<b>chr2</b>\t<b>start2</b>\t<b>end2</b>\t<b>TF motif name</b>\t<b>p-value</b>\t<b>strand</b>\t<b>genomic sequence</b>\n";
    my @result_chr5=fetch_all_chr_pos_tf($dbh_5,$dmps);

    #    print TSV join("\t", @$result_chr5) . "\n";
    for (my $i = 0; $i < @result_chr5; $i++) {
	#print TSV join("\t", @{$result_chr5[$i]}), "\n";
    }

}

sub fetch_all_chr_pos_tf {
    my ($dbh_5, $probe_ids) = @_;
    my $sth = $dbh_5->prepare("SELECT * FROM tf WHERE rsid = ?");

    my $locations = [];

    foreach my $probe_id (@$probe_ids) {
	$sth->execute($probe_id);
	my $result_chr = $sth->fetchall_arrayref();
	my $location;
	foreach my $row (@{$result_chr}) {
	    push(@$locations, $$row[0],"\t",$$row[1],"\t",$$row[2],"\t",$$row[3],"\t",$$row[4],"\t",$$row[5],"\t",$$row[6],"\t",$$row[7],"\t",$$row[8],"\t",$$row[9],"\t",$$row[10],"\n");


	}
	    print TSV @$locations, "\n";
    }
    $sth->finish();
   # return $locations;
}

$dbh_5->disconnect();


###get closest gene

my $dbname6 = "closest.refseq.gene.hg19.db";
my $dsn_6 = "dbi:SQLite:dbname=$datadir/$dbname6";
my $dbh_6 = DBI->connect($dsn_6, "", "") or die $DBI::errstr;

if (defined $gene) {
    print TSV "<b>RefSeq closest gene information (hg19, for more details, see <a href=\"https://www.ncbi.nlm.nih.gov/refseq/about/\">NCBI RefSeq</a>):\n</b>";
    print TSV "<b>chr</b>\t<b>start</b>\t<b>end</b>\t<b>rsid</b>\t<b>chr2</b>\t<b>start2</b>\t<b>end2</b>\t<b>gene name</b>\n";
    my @result_chr6=fetch_all_chr_pos_gene($dbh_6,$dmps);

    #    print TSV join("\t", @$result_chr6) . "\n";
    for (my $i = 0; $i < @result_chr6; $i++) {
	#print TSV join("\t", @{$result_chr6[$i]}), "\n";
    }

}

sub fetch_all_chr_pos_gene {
    my ($dbh_6, $probe_ids) = @_;
    my $sth = $dbh_6->prepare("SELECT * FROM gene WHERE rsid = ?");

    my $locations = [];

    foreach my $probe_id (@$probe_ids) {
	$sth->execute($probe_id);
	my $result_chr = $sth->fetchall_arrayref();
	my $location;
	my $counter = 1;
	my $ucsc = "";
	my $cadd = "";
	my $cadd_out = "";
	my $cadd_chr = "";
	foreach my $row (@{$result_chr}) {
	    #only print first row:
	    
	    if ($counter==1){
	    push(@$locations, $$row[0],"\t",$$row[1],"\t",$$row[2],"\t",$$row[3],"\t",$$row[4],"\t",$$row[5],"\t",$$row[6],"\t",$$row[7],"\n");
	    $counter = 2;
	    $ucsc = join("", "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=",$$row[0],"%3A",$$row[2],"-",$$row[2]);

	    $cadd_chr = $$row[0];
	    $cadd_chr =~ tr/chr//d;
#	    $cadd = join("", "tabix /net/seq/data/projects/Forge2/f2.db/cadd.whole_genome_SNVs_hg19.tsv.gz ",$cadd_chr,":",$$row[2],"-",$$row[2]);
	    $cadd = join("", "tabix /data/breezece/cadd.whole_genome_SNVs_hg19.tsv.gz ",$cadd_chr,":",$$row[2],"-",$$row[2]);
	    $cadd_out = `$cadd`;
	    }
	    
	}
	print TSV @$locations, "\n";
	print TSV "<b>UCSC genome browser link (hg19):</b>\n";
	print TSV $ucsc,"\n";

	print TSV "\n<b><a href=\"https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00835-9\">CADD v1.6</a> annotations at SNP coordinates in <a href=\"https://cadd.gs.washington.edu\">hg19</a> (note that none or several variants may be shown at this position):</b>\n";
	print TSV "<b>chromosome</b>\t<b>position</b>\t<b>reference_allele</b>\t<b>alternative_allele</b>\t<b>raw_CADD_score</b>\t<b>PHRED_score</b>\n";
	print TSV $cadd_out,"\n";
	#print TSV $cadd,"\n";
    }
    $sth->finish();
    # return $locations;
}

$dbh_6->disconnect();

###get eqtlgen eqtls

my $dbname7 = "eqtlgen.db";
my $dsn_7 = "dbi:SQLite:dbname=$datadir/$dbname7";
my $dbh_7 = DBI->connect($dsn_7, "", "") or die $DBI::errstr;

if (defined $eqtlgen) {
    print TSV "\n<b><a href=\"https://www.nature.com/articles/s41588-021-00913-z\">eQTLGen</a> blood FDR 0.05 cis-eQTL information (hg19, for more details, see <a href=\"https://www.eqtlgen.org/\">https://www.eqtlgen.org</a>):\n</b>";
    print TSV "<b>Gene Symbol</b>\t<b>p-value</b>\t<b>rsid</b>\t<b>SNP Chromosome</b>\t<b>SNP Position</b>\t<b>Assessed Allele</b>\t<b>Other Allele</b>\t<b>Z-score</b>\t<b>Gene ID</b>\t<b>Gene Chromosome</b>\t<b>Gene Position </b>\t<b>Number of Cohorts</b>\t<b>Number of Samples</b>\t<b>FDR</b>\t<b>Bonferroni p-value</b>\n";
    my $result_chr2=fetch_all_chr_pos_eqtls2($dbh_7,$dmps);
}

sub fetch_all_chr_pos_eqtls2 {
    my ($dbh_7, $probe_ids) = @_;
    my $sth = $dbh_7->prepare("SELECT * FROM eqtlgen WHERE SNP = ?");

    my $locations = [];

    foreach my $probe_id (@$probe_ids) {
	$sth->execute($probe_id);
	my $result_chr = $sth->fetchall_arrayref();
	my $location;
	foreach my $row (@{$result_chr}) {
	    push(@$locations,$$row[8],"\t",$$row[0],"\t",$$row[1],"\t",$$row[2],"\t",$$row[3],"\t",$$row[4],"\t",$$row[5],"\t",$$row[6],"\t",$$row[7],"\t",$$row[9],"\t",$$row[10],"\t",$$row[11],"\t",$$row[12],"\t",$$row[13],"\t",$$row[14],"\n");

	}
	print TSV @$locations,"\n" ;
    }
    $sth->finish();
    return $locations;
}

$dbh_7->disconnect();


close(TSV);
## ============================================================================


## ============================================================================
## Generate plots
## ============================================================================
##warn "[".scalar(localtime())."] Generating plots...\n";
unless (defined $noplot){
    #Plotting and table routines
    ##warn "[".scalar(localtime())."] Charting against dataset [".$dataset."]...\n"; 
    if ($dataset eq "erc2-H3-all") {
      ChartAllH3Marks($results_filename, $lab, $out_dir, $samples, $cells, $label, $t_marginal, $t_strict, $dataset); # more complex pdf plot
    }
    elsif ($dataset eq "erc2-chromatin15state-all") {
      ChartAllChromatin15StateMarks($results_filename, $lab, $out_dir, $samples, $cells, $label, $t_marginal, $t_strict, $dataset); # more complex pdf plot
    }
    else {
      Chart($results_filename, $lab, $out_dir, $samples, $cells, $label, $t_marginal, $t_strict, $dataset); # basic pdf plot
    }
    dChart($results_filename, $lab, $out_dir, $dataset, $label, $t_marginal, $t_strict, $web); # rCharts Dimple chart
    table($results_filename, $lab, $out_dir, $web); # Datatables chart
  }
## ============================================================================

##warn "[".scalar(localtime())."] Done.\n";

if ($autoopen) {
    system("open $out_dir/$lab.table.html");
    system("open $out_dir/$lab.dchart.html");
    system("open $out_dir/$lab.chart.pdf");
}


####################################################################################################
####################################################################################################
##
##  Sub-functions
##
####################################################################################################
####################################################################################################


=head2 parse_pvalue_thresholds

 Arg[1]         : string $thresholds
 Returns        : arrayref of marginal and strict thresholds (floats)
 Example        : ($t_marginal, $t_strict) = parse_pvalue_thesholds("0.05,0.01");
 Description    : This function returns the both marginal and strict p-value thresholds as read from
                  the command line option. The input string should contain both numbers separated by
                  a comma.
 Exceptions     : Dies if $thresholds is empty, does not contain numbers or are not defined between
                  0 and 1 and/or the marginal threshold is not larger or equal to the strict one.

=cut

sub parse_pvalue_thresholds {
    my ($thresh) = @_;
    my ($t_marginal, $t_strict);

    if (!$thresh) {
        die "Cannot read p-value thresholds from an empty string\n";
    }

    ($t_marginal, $t_strict) = split(",", $thresh);
    unless (looks_like_number($t_marginal) && looks_like_number($t_strict)){
        die "You must specify numerical p-value thresholds in a comma separated list\n";
    }
    unless ((1 >= $t_marginal) && ($t_marginal >= $t_strict) && ($t_strict >= 0)) {
        die "The p-value thresholds should be 1 >= T.marginal >= T.strict >= 0\n";
    }
    return ($t_marginal, $t_strict);
}


=head2 get_input_probes

 Arg[1]         : string $filename
 Arg[2]         : arrayref $probe_list
 Returns        : arrayref of probe IDs (string)
 Example        : $dmps = get_input_probes("input.txt", undef);
 Example        : $dmps = get_input_probes(undef, ["cg13430807", "cg10480329,cg06297318,cg19301114"]);
 Example        : $dmps = get_input_probes(undef, undef);
 Description    : This function returns the list of input probe IDs. This can come from either
                  $filename if defined or from $probe_list otherwise. Each element in $probe_list is a
                  string which contains one or more probe IDs separated by commas (see Examples).
                  Falls back to the default data set from Jaffe and Irizarry.
                  The set of probe IDs is checked to remove redundant entries.
 Exceptions     : Dies if the file is not found or cannot be opened for whatever reason.

=cut

sub get_input_probes {
    my ($filename, $probe_list) = @_;
    my $probes;

    if (defined $filename) {
        my $fh;
        if ($filename =~ /\.gz$/) {
            open($fh, "gunzip -c $filename |") or die "cannot open file $filename : $!";
        } elsif ($filename =~ /\.bz2$/) {
            open($fh, "bunzip2 -c $filename |") or die "cannot open file $filename : $!";
        } else {
            open($fh, "$filename") or die "cannot open file $filename : $!";
        }
        $probes = process_file($fh, $format, $dbh, $array, $filter);

    } elsif ($probe_list and @$probe_list) {
        @$probes = split(/,/, join(',', @$probe_list));

    } else{
        # If no options are given it will run on the default set of DMPs
        warn "No SNP input given, so running on default set of SNPs (a set of PR interval-associated SNPs from the GWAS catalogue).";
        @$probes = qw(rs11773845 rs17744182 rs17026156 rs7433306 rs7135659 rs3733017 rs1524976 rs3103778 rs1994318 rs16926523 rs10447419 rs7604827 rs12595668 rs746265 rs11732231 rs11773845 rs267567 rs3891585 rs6801957 rs3922844 rs6763048 rs1895585  rs10865355 rs6798015 rs3922844 rs6599222 rs7312625 rs7692808 rs3807989 rs11897119 rs6800541 rs11708996 rs4944092 rs11047543 rs251253 rs1896312);
    }

    # Remove redundancy in the input
    my %probes_hash;
    foreach my $probe (@$probes) {
        $probes_hash{$probe}++;
    }

    while (my ($probe, $num) = each %probes_hash) {
        if ($num > 1) {
            warn "$probe is present $num times in the input. Analysing only once.\n";
        }
    }

    @$probes = keys %probes_hash;

    return($probes);
}


=head2 find_missing_probes

 Arg[1]         : arrayref of strings $original_probe_ids
 Arg[2]         : hashref $existing_probe_ids (keys are probe_ids, values are ignored)
 Arg[3]         : hashref $excluded_probe_ids (keys are probe_ids, values are ignored)
 Returns        : int $num_missing_probes
 Example        : my $num_missing_probes = find_missing_probes(['cg001', 'cg002', 'cg003', 'cg004'],
                    {'cg001' => 1, 'cg003 => 1}, {'cg002' => 'cg001');
 Description    : Detects and prints the list of missing probes if any.
 Exceptions     :

=cut

sub find_missing_probes {
    my ($original_probes, $existing_probes_hash, $excluded_probes_hash) = @_;
    my $num_missing_probes = 0;

    my $missing_probes = [];
    foreach my $probe_id (@$original_probes) {
        unless ($existing_probes_hash->{$probe_id} or $excluded_probes_hash->{$probe_id}) {
            push @$missing_probes, $probe_id;
        }
    }
    $num_missing_probes = scalar @$missing_probes;

    if ($num_missing_probes > 0) {
        warn "The following $num_missing_probes DMPs have not been analysed because they were not found on the selected array\n";
        warn join("\n", @$missing_probes) . "\n";
    }

    return $num_missing_probes;
}



    


sub fetch_all_chr_pos {
    my ($dbh_2, $probe_ids) = @_;

    # Get array ID
#    my $array_id_sth = $dbh->prepare("SELECT array_id FROM array WHERE array_tag = ?");
#    $array_id_sth->execute($array_tag);
#    my $array_id_result = $array_id_sth->fetchrow_arrayref();
#    my $array_id = @{$array_id_result}[0];
#    $array_id_sth->finish();

    
# Get probe mapping ID
#    my $probe_mapping_id_sth = $dbh->prepare("SELECT probe_mapping_id FROM probe_mapping_info WHERE array_id = ?");
#    $probe_mapping_id_sth->execute($array_id);
#    my $probe_mapping_id_result = $probe_mapping_id_sth->fetchrow_arrayref();
#    my $probe_mapping_id = @{$probe_mapping_id_result}[0];
#    $probe_mapping_id_sth->finish();

#        my $sth = $dbh->prepare("SELECT chr_name, position FROM probe_mapping WHERE probe_id = ?");
    my $sth = $dbh_2->prepare("SELECT chr,start,end,probe,catoscore,strand,motif,score1,allele1,allele2,score2,tissuesample FROM cato WHERE probe = ?");
    
    my $locations = [];

    foreach my $probe_id (@$probe_ids) {
	$sth->execute($probe_id);
	my $result_chr = $sth->fetchall_arrayref();
	my $location;
	foreach my $row (@{$result_chr}) {
	    push(@$locations, $$row[0],$$row[1],$$row[2],$$row[3],$$row[4],$$row[5],$$row[6],$$row[7],$$row[8],$$row[9],$$row[10],$$row[11]);
#	my ($chr, $pos) = @$locations;
	}
    }
    $sth->finish();
    return $locations;
}


