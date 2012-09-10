#! /usr/bin/perl
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Run this script in the top-level Runuran directory.
# The script updates version and date in
#   DESCRIPTION
#   man/Runuran-package.Rd
#   inst/doc/inputs/version.tex
# ---------------------------------------------------------------------------

use strict;
use Getopt::Std;

# --- Constants -------------------------------------------------------------

my $vignette_file = "./inst/doc/inputs/version.tex";
my $package_Rd_file = "./man/Runuran-package.Rd";

# --- Usage -----------------------------------------------------------------

sub usage {
    my $progname = $0;
    $progname =~ s#^.*/##g;
        
    print STDERR <<EOM;

usage: $progname -u [-i]

  -u ... update 'version.tex' and 'Runuran-package.Rd'
  -i ... increment version number and update date in 
         file 'DESCRIPTION' first

EOM

    exit -1;
}

# --- Read command line options ---------------------------------------------

my %opts;
getopts('ui', \%opts) or usage();
my $update = $opts{'u'};
my $increment = $opts{'i'};

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

# get data
$description =~ m/^.*\nVersion:\s*(.*?)\s*\n/s 
    or die "Cannot find Version";
my $version = $1;

$description =~ m/^.*\nDate:\s*(.*?)\s*\n/s 
    or die "Cannot find Date";
my $date = $1;
my $longdate;

# --- Increment version number -----------------------------------------------

my @mnames = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );

if ($increment) {
    (my $major, my $minor, my $rev) = split (/\./, $version, 3);
    $version = "$major.$minor.".($rev+1);

    my $year = 1900 + (localtime(time))[5];
    my $month = sprintf("%02d", 1 + (localtime(time))[4]);
    my $day = sprintf("%02d", (localtime(time))[3]);

    $date = "$year-$month-$day";
    $longdate = "$mnames[$month-1] $day, $year";
}
else {
    (my $year, my $month, my $day) = split (/\-/, $date, 3);
    $longdate = "$mnames[$month-1] $day, $year";
}

# --- Print data -------------------------------------------------------------

if ($increment) {
    print "increment version number and date\n";
}
else {
    print "use version number and date as-is\n";
}
print "Version = '$version'\n";
print "Date = '$date'  ($longdate)\n";

# --- Update Vignette -------------------------------------------------------

open VIGNETTE, ">$vignette_file"
    or die "Cannot open $vignette_file for writing: $!";
print VIGNETTE "\\date{Version $version -- $longdate}\n";
close VIGNETTE;

# --- Update DESCRIPTION ----------------------------------------------------

if ($increment) {
    $description =~ s/\nVersion:\s*.*?\n/\nVersion: $version\n/
	or die "DESCRIPTION: Cannot find field 'Version:'";
    $description =~ s/\nDate:\s*.*?\n/\nDate: $date\n/
	or die "DESCRIPTION: Cannot find field 'Date:'";
    open DESC, ">DESCRIPTION"
	or die "Cannot open file 'DESCRIPTION' for writing: $!";
    print DESC $description;
    close DESC;
}

# --- Update Runuran-package -------------------------------------------------

open PACKAGE, "$package_Rd_file"
    or die "Cannot open file 'man/Runuran-package.Rd' for reading: $!";
my $package;
while (<PACKAGE>) {
    $package .= $_;
}
close PACKAGE; 

$version = sprintf("%-14s", $version); 
$package =~ s/\n(\s*Version:\s*\\tab\s+).*?\\cr\n/\n$1$version\\cr\n/
	or die "man/Runuran-package.Rd: Cannot find field 'Version:'";

$date = sprintf("%-14s", $date); 
$package =~ s/\n(\s*Date:\s*\\tab\s+).*?\\cr\n/\n$1$date\\cr\n/
	or die "man/Runuran-package.Rd: Cannot find field 'Date:'";

open PACKAGE, ">$package_Rd_file"
    or die "Cannot open file 'man/Runuran-package.Rd' for writing: $!";
print PACKAGE $package;
close PACKAGE;

# --- end -------------------------------------------------------------------

exit 0;

# ---------------------------------------------------------------------------
