#!/usr/local/bin/perl -w


use strict;
use English qw( -no_match_vars ) ;
use Carp;
use Getopt::Long;
use FileHandle;

main();
0;

sub main { 
   
    my $opts = initialise(); 
   
    my $input       = $ARGV[0];
    my $output_good = $ARGV[1];
    my $output_bad  = $ARGV[2];

    print qq{filtering fastq file $input\n};
    
    filter_NNN_fastq($input, $output_good, $output_bad, $opts);
    
}
  
sub filter_NNN_fastq {
    
    my ($input, $output_good, $output_bad, $opts)=@_;

    my $input_fh = FileHandle->new();
    $input_fh->open("<$input") or croak "Error opening input file $input:\n$ERRNO";
    my $good_fh = FileHandle->new();
    $good_fh->open(">$output_good") or croak "Error opening output file $output_good:\n$ERRNO";
    my $bad_fh = FileHandle->new();
    $bad_fh->open(">$output_bad") or croak "Error opening output file $output_good:\n$ERRNO";

    my ($count_all, $count_good, $count_bad) = (0, 0, 0);

    my ($state, $count_n, $output) = (0, 0, "");
    while (my $line = <$input_fh>) {

	$state++;

        if( $state == 1){
            croak "Invalid sequence header $line" unless $line =~ m/^\@\S+$/;
            $output = $line;
        }elsif( $state == 2 ){
            ####croak "Invalid sequence $line" unless $line =~ m/^[ACGTN]+$/;
            $output .= $line;
            my @sequence = split(//, $line);
            map{ $count_n++ if $_ eq 'N' } @sequence;
        }elsif( $state == 3 ){
            croak "Invalid quality header $line" unless $line =~ m/^\+$/;
            $output .= $line;
        }elsif( $state == 4 ){
            ####croak "Invalid quality values $line\n" unless $line =~ m/^[\x21-\x85]+$/;
            $output .= $line;

            $count_all++;
            if( $count_n < $opts->{threshold} ){
                $count_good++;
                $good_fh->print($output);
            }else{
                $count_bad++;
                $bad_fh->print($output);
            }
            ($state, $count_n, $output) = (0, 0, "");
        }
    }

    printf("Fraction N - filtered %2.3f \n", 100.0*$count_good/$count_all);
}

sub initialise{

    my $usage = qq[

    Usage:
            filter_NNN_reads.pl [-help] [-threshold <threshold>] input output1 output2

    ];

    my $options = { threshold => 3 };

    my $rc = GetOptions($options,
			'threshold:i',
                        'help'
   		       );

    die("$usage\n") if( !$rc || $options->{help} );
    die("$usage\n") unless scalar(@ARGV) == 3;

    return $options;
}
