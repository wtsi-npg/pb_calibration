#!/usr/local/bin/perl
#
# check a PB calibartion table
#
# generate some stats

use strict;
use warnings;
use English qw(-no_match_vars);
use Carp;
use Getopt::Long;
use FileHandle;

sub usage;
sub initialise;
sub process;

my $opts = initialise;

process;

exit;

# ----------------------------------------------------------------------
sub usage {

  ## no critic

  print STDERR "\n";
  print STDERR "check_caltable.pl [-verbose] caltable\n";
  print STDERR "\n";
  print STDERR "    options:\n";
  print STDERR "\n";
  print STDERR "    --help         print this message and quit\n";
  print STDERR "\n";

  ## use critic

}

# ----------------------------------------------------------------------
sub initialise {

  my %opts = ( rename => 0, m_min => 25, q_min => 20, j_min => 15 );
  my $rc = GetOptions(\%opts, 'rename!', 'm_min=i', 'q_min=i', 'j_min=i', 'verbose', 'help');
  if ( ! $rc) {
    print {*STDERR} "\nerror in command line parameters\n" or croak 'print failed';
    usage;
    exit 1;
  }

  if (exists $opts{help}) {
    usage;
    exit;
  }

  if (@ARGV != 1) {
    usage;
    exit;
  }

  return \%opts;

}

# ----------------------------------------------------------------------
# calculate the following stats
# m  maximum quality
# q  overall quality
# j  spread

sub process {

    my $file = $ARGV[0];

    my $fh = FileHandle->new();
    $fh->open("<$file") or croak "$file:\n$ERRNO";

    my @cols;
    while (my $line = <$fh>) {
        chomp($line);
        my @F = split(/\t/, $line);
        my ($tf, $ef, $qf);
        if ($#F == 8) {
            ($tf, $ef, $qf) = (4, 5, 8);
        }else{
            die "Invalid calibration table";
        }
        push(@cols,[$F[$tf],$F[$ef],$F[$qf]]);
    }
    close $fh;

    my ($t, $e, $m, $q, $j) = (0, 0, 0.0, 0.0, 0.0);
    map {$t+=$_->[0]; $e+=$_->[1]; $m=$_->[2] if $_->[2] > $m} @cols;
    $q = -10.0*log($e/$t) / log(10.0);
    map {$j +=($_->[2]-$q)**2 * $_->[0] / $t} @cols;
    printf("$file\tm=$m\tq=$q\tj=$j\n") if $opts->{"verbose"};

    my $flag = 0;
    if( $m < $opts->{"m_min"} ) {
        printf("$file\tmaximum quality %d < %d\n", $m, $opts->{"m_min"});
        $flag = 1;
    }
    if( $q < $opts->{"q_min"} ) {
        printf("$file\toverall quality %d < %d\n", $q, $opts->{"q_min"});
        $flag = 1;
    }
    if( $j < $opts->{"j_min"} ) {
        printf("$file\tquality spread %d < %d\n", $j, $opts->{"j_min"});
        # currently we just monitor this
        ####$flag=1;
    }
    if( $flag && $opts->{"rename"}) {
        printf("Renaming $file to $file.ignore\n");
        rename $file,$file.".ignore" or croak 'failed to rename file $file';
    }

}
