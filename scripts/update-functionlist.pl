#! /usr/bin/perl

# ---------------------------------------------------------------------------
# Run this script in the top-level Runuran directory.
# The script updates list of distributions in
#    man/Runuran.distributions.Rd
#    inst/doc/src/tab-distributions.tex
# ---------------------------------------------------------------------------

use strict;
use File::Find;
use Getopt::Std;

# --- Constants -------------------------------------------------------------

my $ur_Vignette = "./inst/doc/src/tab-generators.tex";
my $ud_Vignette = "./inst/doc/src/tab-distributions.tex";
my $generators_Rd_file = "./man/Runuran.special.generators.Rd";
my $distributions_Rd_file = "./man/Runuran.distributions.Rd";

my $man_dir = "./man";

# --- Usage -----------------------------------------------------------------

sub usage {
    my $progname = $0;
    $progname =~ s#^.*/##g;
        
    print STDERR <<EOM;

usage: $progname -u

  -u ... update list of 'ur' calls for distributions in vignette
         and package help file 'Runuran-package.Rd'

EOM

    exit -1;
}

# --- Read command line options ---------------------------------------------

my %opts;
getopts('u', \%opts) or usage();
my $update = $opts{'u'};

usage unless $update;

# --- Read file 'DESCRIPTION' -----------------------------------------------

open DESC, "DESCRIPTION" 
    or die "You must run this script in the top-level Runuran directory";
my $description;
while (<DESC>) {
    $description .= $_;
}
close DESC; 

# check name of package
die "You must run this script in the top-level Runuran directory"
    unless $description =~ /^\s*Package:\s+Runuran\s*\n/;

# --- Get list of .Rd files -------------------------------------------------

opendir MAN, $man_dir
    or die "Cannot open directory $man_dir: $!";

my @ud_files = grep { /^ud.*\.Rd/ && -f "$man_dir/$_" } readdir MAN;
rewinddir MAN;

my @ur_files = grep { /^ur.*\.Rd/ && -f "$man_dir/$_" } readdir MAN;
rewinddir MAN;

my @new_files = grep { /^.*\.new\.Rd/ && -f "$man_dir/$_" } readdir MAN;
closedir MAN;

# --- Get list of distributions ---------------------------------------------

my %udcont;
my %uddiscr;

my $udmax = "Function";
foreach my $ud (sort @ud_files) {
    my $udcall = $ud;
    $udcall =~ s/\.Rd\s*$//;
    $udmax = $udcall if length($udcall) > length($udmax);
    open UD, "$man_dir/$ud" 
	or die "Cannot open file '$man_dir/$ud' for reading: $!";
    while (<UD>) { 
	chomp;
	next unless $_ =~/^\s*\[Distribution]\s+\-\-\s+(.*)$/;
	my $tmp = $1;
	$tmp =~ /^(.*)\s*\.\s*\%\%\s*(Continuous|Discrete)\s*$/
	    or die "Format error in '$_'";
	my $distr = $1;
	my $type = $2;
	if ($type eq "Continuous") {
	    $udcont{$udcall} = $distr;
	}
	if ($type eq "Discrete") {
	    $uddiscr{$udcall} = $distr;
	}
	last;
    }
}

my $n_udcont  = scalar keys %udcont;
my $n_uddiscr = scalar keys %uddiscr;

print "# Distributions\n";
print "#   continuous distributions = $n_udcont\n"; 
print "#   discrete distributions   = $n_uddiscr\n";

# --- Get list of special generators ----------------------------------------

my %urcont;
my %urdiscr;

my $urmax = "Function";
foreach my $ur (sort @ur_files) {
    my $urcall = $ur;
    $urcall =~ s/\.Rd\s*$//;
    $urmax = $urcall if length($urcall) > length($urmax);
    open UR, "$man_dir/$ur" 
	or die "Cannot open file '$man_dir/$ur' for reading: $!";
    while (<UR>) { 
	chomp;
	next unless $_ =~/^\s*\[Special Generator]\s+\-\-\s+(.*)$/;
	my $tmp = $1;
	$tmp =~ /^Sampling Function:\s*(.*)\s*\.\s*\%\%\s*(Continuous|Discrete)\s*$/
	    or die "Format error in '$_'";
	my $distr = $1;
	my $type = $2;
	if ($type eq "Continuous") {
	    $urcont{$urcall} = $distr;
	}
	if ($type eq "Discrete") {
	    $urdiscr{$urcall} = $distr;
	}
	last;
    }
}

my $n_urcont  = scalar keys %urcont;
my $n_urdiscr = scalar keys %urdiscr;

print "# Special generators\n";
print "#   continuous distributions = $n_urcont\n"; 
print "#   discrete distributions   = $n_urdiscr\n";

# --- Print list of distributions into vignette -----------------------------

my $udheader = 
    "\\begin{tabbing}\n" .
    "\t\\hspace*{1em}\n" .
    "\t\\= \\code{$udmax}~~\\=\\ldots~~\\=  Distribution \\kill\n" .
    "\t\\> \\emph{Function} \\> \\> \\emph{Distribution} \\\\[1ex]\n";
my $udbottom = 
    "\\end{tabbing}\n";

my $udcont_list = "";
foreach my $ud (sort keys %udcont) {
    $udcont_list .= "\t\\> \\code{$ud}\t\\> \\ldots \\> $udcont{$ud} \\\\\n";
}

my $uddiscr_list = "";
foreach my $ud (sort keys %uddiscr) {
    $uddiscr_list .= "\t\\> \\code{$ud}\t\\> \\ldots \\> $uddiscr{$ud} \\\\\n";
}

open DISTR, ">$ud_Vignette"
    or die "Cannot open file '$ud_Vignette' for writing: $!";
print DISTR
    "\\paragraph{Continuous Univariate Distributions ($n_udcont)}\n\n" .
    $udheader . $udcont_list . $udbottom . "\n";
print DISTR
    "\\paragraph{Discrete Univariate Distributions ($n_uddiscr)}\n\n".
    $udheader . $uddiscr_list . $udbottom . "\n";
close DISTR;

# --- Print list of special generators into vignette ------------------------

my $urheader = 
    "\\begin{tabbing}\n" .
    "\t\\hspace*{1em}\n" .
    "\t\\= \\code{$urmax}~~\\=\\ldots~~\\=  Distribution \\kill\n" .
    "\t\\> \\emph{Function} \\> \\> \\emph{Distribution} \\\\[1ex]\n";
my $urbottom = 
    "\\end{tabbing}\n";

my $urcont_list = "";
foreach my $ur (sort keys %urcont) {
    $urcont_list .= "\t\\> \\code{$ur}\t\\> \\ldots \\> $urcont{$ur} \\\\\n";
}

my $urdiscr_list = "";
foreach my $ur (sort keys %urdiscr) {
    $urdiscr_list .= "\t\\> \\code{$ur}\t\\> \\ldots \\> $urdiscr{$ur} \\\\\n";
}

open DISTR, ">$ur_Vignette"
    or die "Cannot open file '$ur_Vignette' for writing: $!";
print DISTR
    "\\paragraph{Continuous Univariate Distributions ($n_urcont)}\n\n" .
    $urheader . $urcont_list . $urbottom . "\n";
print DISTR
    "\\paragraph{Discrete Univariate Distributions ($n_urdiscr)}\n\n".
    $urheader . $urdiscr_list . $urbottom . "\n";
close DISTR;

# --- Print list of distributions into distributions Rd file ----------------

my $udheader = 
    "  \\tabular{lcl}{ \n" .
    "    \\emph{Function} \\tab \\tab \\emph{Distribution} \\cr\n";
my $udbottom = 
    "  }\n";

my $udcont_list = "";
foreach my $ud (sort keys %udcont) {
    $udcont_list .= "    \\code{\\link{$ud}} \\tab \\ldots \\tab $udcont{$ud} \\cr\n";
}

my $uddiscr_list = "";
foreach my $ud (sort keys %uddiscr) {
    $uddiscr_list .= "    \\code{\\link{$ud}} \\tab \\ldots \\tab $uddiscr{$ud} \\cr\n";
}

my $list = 
    "  Continuous Univariate Distributions ($n_udcont):\n\n".
    $udheader . $udcont_list . $udbottom . "\n" .
    "  Discrete Distributions ($n_uddiscr):\n\n".
    $udheader . $uddiscr_list . $udbottom . "\n";

open DISTRIBUTIONS, "$distributions_Rd_file"
    or die "Cannot open file '$distributions_Rd_file' for reading: $!";
my $distributions;
while (<DISTRIBUTIONS>) {
    $distributions .= $_;
}
close DISTRIBUTIONS; 

my $begin = "  %% -- begin: list of distributions --\s*\n";
my $end   = "  %% -- end: list of distributions --\s*\n";

$distributions =~ s/($begin)(.*?)($end)/$1$list$3/s
    or die "Cannot find marker for list of distributions";

open DISTRIBUTIONS, ">$distributions_Rd_file"
    or die "Cannot open file '$distributions_Rd_file' for writing: $!";
print DISTRIBUTIONS $distributions;
close DISTRIBUTIONS;

# --- Print list of distributions into special generators Rd file -----------

my $urheader = 
    "  \\tabular{lcl}{ \n" .
    "    \\emph{Function} \\tab \\tab \\emph{Distribution} \\cr\n";
my $urbottom = 
    "  }\n";

my $urcont_list = "";
foreach my $ur (sort keys %urcont) {
    $urcont_list .= "    \\code{\\link{$ur}} \\tab \\ldots \\tab $urcont{$ur} \\cr\n";
}

my $urdiscr_list = "";
foreach my $ur (sort keys %urdiscr) {
    $urdiscr_list .= "    \\code{\\link{$ur}} \\tab \\ldots \\tab $urdiscr{$ur} \\cr\n";
}

my $list = 
    "  Continuous Univariate Distributions ($n_urcont):\n\n".
    $urheader . $urcont_list . $urbottom . "\n" .
    "  Discrete Distributions ($n_urdiscr):\n\n".
    $urheader . $urdiscr_list . $urbottom . "\n";

open GENERATORS, "$generators_Rd_file"
    or die "Cannot open file '$generators_Rd_file' for reading: $!";
my $generators;
while (<GENERATORS>) {
    $generators .= $_;
}
close GENERATORS; 

my $begin = "  %% -- begin: list of distributions --\s*\n";
my $end   = "  %% -- end: list of distributions --\s*\n";

$generators =~ s/($begin)(.*?)($end)/$1$list$3/s
    or die "Cannot find marker for list of distributions";

open GENERATORS, ">$generators_Rd_file"
    or die "Cannot open file '$generators_Rd_file' for writing: $!";
print GENERATORS $generators;
close GENERATORS;

# --- end -------------------------------------------------------------------

exit 0;

# ---------------------------------------------------------------------------

