#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Temp qw(tempdir);

$| = 1;
$ENV{GZIP} = '-4';

my %opts = (mode  => 'rta',
	    lanes => 2,
	    tiles => 2,
	    read  => [],
	    filter => 0,
	    prb => 0,
	    clean => 1);

GetOptions(\%opts,
	   "mode=s", "lanes=i", "tiles=i", "read=i@", "filter!", "clean!",
	   "prb!")
    || usage();

my $nreads = @{$opts{read}};
if ($nreads < 1 || $nreads > 3) { usage(); }

if ($opts{mode} !~ /^(1.3|ipar|rta)$/) { usage(); }

my $tmp = tempdir(CLEANUP => $opts{clean});

my @fake_cmd = ('./fake_data',
		'-d', $tmp,
		'-m', $opts{mode},
		'-l', $opts{lanes},
		'-t', $opts{tiles},
		map { ('-r', $_) } @{$opts{read}});
if ($opts{prb}) {
    push(@fake_cmd, '-p');
}

my $bustard_dir;
print "Running : @fake_cmd\n";
open(my $fake, '-|', @fake_cmd)
    || die "Couldn't open pipe to @fake_cmd: $!\n";
while (<$fake>) {
    print $_;

    if (/^(?:Made|Reused) basecalls directory : (.*)\n/) {
	$bustard_dir = $1;
    }
}
close($fake) || die "Error running @fake_cmd\n";

unless ($bustard_dir) {
    die "Couldn't find bustard directory from fake_data output.\n";
}

my @qseqs = glob("$bustard_dir/s_*_qseq.txt");

my $srf = "$tmp/test.srf";

my @i2s_extra;
my @check_extra;

if ($nreads > 2) {
    push(@i2s_extra,
	 '-use_bases',"Y$opts{read}->[0],I$opts{read}->[1],Y$opts{read}->[2]");
}

if ($opts{filter}) {
    push(@i2s_extra, '-filter-bad-reads');
    push(@check_extra, '-f');
}

if ($opts{prb}) {
    push(@i2s_extra, '-logodds');
}

my @i2s_cmd = ('../illumina2srf/illumina2srf',
	       '-clobber', '-o', $srf,
	       '-nse', '-sig',
	       '-bustard-dir', $bustard_dir,
	       '-b',
	       @i2s_extra,
	       @qseqs);

print "Running : @i2s_cmd\n";
system(@i2s_cmd) && die "Error running @i2s_cmd\n";

my @check_cmd = ('./check_srf', '-b', $bustard_dir, @check_extra, $srf);

print "Running : @check_cmd\n";
system(@check_cmd) && die "Error running @check_cmd";

exit;

sub usage {
    print "Usage: $0 -read <len> [-read <len> ...]\n";
    print "Options:\n";
    print " -mode 1.3|ipar|rta       Pipeline type to test\n";
    print " -lanes <n>               Number of lanes to make\n";
    print " -tiles <n>               Number of tiles to make\n";
    print " -read <n>                Add read of length <n> bases\n";
    print " -[no]filter              Do [Don't] test filter bad reads\n";
    print " -[no]prb                 Do [Don't] test .prb files\n";
    print " -[no]clean               ";
    print "There must be between 1 and 3 reads.\n";
    exit(1);
}
