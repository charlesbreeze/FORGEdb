#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw( basename dirname );

my $useFirstFieldColoring = 0;
my $alternateRowColoring = 1;
if (@ARGV > 0) {
    if ($ARGV[0] eq "firstcol") {
        $useFirstFieldColoring = 1;
    } else {
        $alternateRowColoring = int( $ARGV[0] );
    }
}
my $showPNG = 1;
my $header1 = 0; # flag to change class of first row to "header"
if (@ARGV > 1) {
    if ($ARGV[1] eq "linkonly") {
        $showPNG = 0;
    } elsif ($ARGV[1] eq "header") {
        $header1 = 1;
    }
}

#print "<html>\n";
#print "<body>\n";

my $message = <<'END_MESSAGE';

<html>
  <head>
    <title>FORGEdb</title>
    <link rel="stylesheet" type="text/css" href="blog0.css" />
  </head>
    <body>

END_MESSAGE

	 
print $message;

my $message2 = <<'END_MESSAGE';                                                                                                                                                                                                                                                              
    </body>                                                                                                                                      </html>
                                                                                                                                              
END_MESSAGE

print $message2;

print '<table class="center" border="1" cellspacing="0" cellpadding="2" style="font-family:Helvetica" style="text-align:center">' . "\n";
my $count = 0;
my $prevFirstCol = "";
while (<STDIN>) {
    chomp;
    #my $trclass= (++$count % 2 == 0) ? "row_even" : "row_odd";
    my @things = split /\t/, $_;
    my $trclass;
    if ($useFirstFieldColoring) {
        my $firstcol = $things[0];
        if ($prevFirstCol ne $firstcol) {
            ++$count;
        }
        $trclass= ($count % 2 == 0) ? "row_even" : "row_odd";
        $prevFirstCol = $firstcol;
    } else {
        $trclass= (int((++$count - 1)/$alternateRowColoring) % 2 == 0) ? "row_even" : "row_odd";
    }
    if ($header1 and ($count < 2)) {
        $trclass = "header";
    }
    print "<tr class=$trclass>";
    #my @things = split; 
    for (my $i = 0; $i < @things; ++$i) {

	my $colspan= "";
	my $align="style=\"text-align:center\"";
	if (scalar(@things)<2){ 
	    $colspan= "colspan=10";
	    $align="style=\"text-align:left\"";
	}
	    
        my $tdclass= "";
        my $localpath = $things[$i]; # save this before fixpath changes it
        $things[$i] = fixpath( $things[$i] );
        if ($things[$i] =~ /^[0-9,\.Ee+-]+$/ or $things[$i] =~ /\>[0-9,\.Ee+-]+\</) {
            $tdclass= "class=numeric";

            if ($things[$i] =~ /\./) {
                my $precision = 10000;
                #$things[$i] = int( 0.5 + ($things[$i]*$precision) ) / $precision;
                #TODO  Preserve significant digits!
            } elsif ($things[$i] =~ /^([\d]+)$/) {
                $things[$i] = encomma($things[$i]);
            }
        } elsif (($things[$i] =~ /^http.*\/([^\/]+)$/) or ($things[$i] =~ /^\/\~.*\/([^\/]+)$/)) {
            my $shortname = $1;
            my $orig = $things[$i];
            if ($showPNG and (($orig =~ /\.png$/i) or ($orig =~ /\.svg$/i))) { 
                my $filename = basename( $localpath );
                my $localdirname = dirname( $localpath );
                my $webdirname = dirname( $orig );
                my $localthumbnail = "$localdirname/th_$filename";
                if (-s $localthumbnail) {
                    # show thumbnail as a link to the image
                    my $webthumbnail = "$webdirname/th_$filename";
                    $things[$i] = "<a href=\"$orig\"><img src=\"$webthumbnail\"/></a>";
                } else {
                    #warn "$localthumbnail not found\n";
                    # just show the image
                    $things[$i] = "<img src=\"$orig\"/>";
                }
            } else {
                if ($shortname =~ /greatStart.php/) {
                    $shortname = "GREAT"; # used to be a very long name
                }
                $things[$i] = "<a href=\"$orig\">$shortname</a>";
            }
        }
        print " <td $tdclass $colspan $align> " . $things[$i] . " </td>";
    }
    print "</tr>\n";
}
print "</table>\n";
#print "</body>\n";
#print "</html>\n";

sub fixpath {
    my $orig = shift;
    my $path = $orig;
    if ($orig =~ /^\/home\/([^\/]+)\/public_html\/(.*)$/) {
        my ($user,$relativePath) = ($1,$2);
        $path = "/~$user/$relativePath";
    } elsif ($orig =~ /.*\/([^\/]+)_enhancer\/(.*)$/) {
        my ($user,$relativePath) = ($1,$2);
        $path = "http://www.uwencode.org/proj/${user}_enhancer/$relativePath";
    } elsif ($orig =~ /.net.lebowski.vol1.work.proj.(.*)$/) {
        my $relativePath = $1;
        $path = "http://www.uwencode.org/proj/$relativePath";
    }
    return $path;
}
#/net/lebowski/vol1/work/proj/ehaugen/2012oct17

#<a href="http://www.uwencode.org/proj/ehaugen_enhancer/2011Sep27/clusters.html">Clustered motifs.</a>
#<a href="http://www.uwencode.org/proj/ehaugen_enhancer/2011Sep27/summary.fvi_OnlyPromoters.factorscores2.html"> in promoters (Gencode V7 TSS +/- 2.5kb)</a>

sub encomma {
    my $orig = shift;
    my $remain = "$orig";
    my $result = "";
    while (length($remain) > 3) {
        my $chunk = substr( $remain, length($remain) - 3 );
        $remain = substr( $remain, 0, length($remain) - 3 );
        if (length($result) > 0) {
            $result = ",$result";
        }
        $result = "$chunk$result";
        #warn "$remain | $result\n";
    }
    if ((length($result) > 0) and (length($remain)>0)) {
        $result = ",$result";
    }
    $result = "$remain$result";
    #warn "encomma($orig)=($result)\n";
    return $result;
}

