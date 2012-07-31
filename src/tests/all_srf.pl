#!/usr/bin/perl
use strict;
use warnings;

$| = 1;

my $srcdir = $ENV{srcdir} || '.';

my @tests = ([qw(srf_test.pl -mode 1.3 -read 37 -read 37 )],
	     [qw(srf_test.pl -mode ipar -read 37 -read 37 )],
	     [qw(srf_test.pl -mode rta -read 37 -read 37 )],
	     [qw(srf_test.pl -mode 1.3 -filter -read 37 -read 37 )],
	     [qw(srf_test.pl -mode ipar -filter -read 37 -read 37 )],
	     [qw(srf_test.pl -mode rta -filter -read 37 -read 37 )],
	     [qw(srf_test.pl -mode 1.3 -read 120 -read 6 -read 120)],
	     [qw(srf_test.pl -mode ipar -read 120 -read 6 -read 120)],
	     [qw(srf_test.pl -mode rta -read 120 -read 6 -read 120)],
	     [qw(srf_test.pl -mode 1.3 -filter -read 12 -read 2 -read 12)],
	     [qw(srf_test.pl -mode ipar -filter -read 12 -read 2 -read 12)],
	     [qw(srf_test.pl -mode rta -filter -read 12 -read 2 -read 12)],
	     [qw(srf_test.pl -mode 1.3 -filter -prb -read 5)],
	     [qw(srf_test.pl -mode ipar -filter -prb -read 5)],
	     [qw(srf_test.pl -mode rta
		 -read 510 -read 6 -read 10 -lanes 1 -tiles 1)]);

my $passed = 0;
my $failed = 0;

open(my $oldout, '>&', \*STDOUT) || die "Couldn't dup STDOUT: $!\n";
open(my $olderr, '>&', \*STDERR) || die "Couldn't dup STDERR: $!\n";

for (my $i = 1; $i <= @tests; $i++) {
    my $test = $tests[$i - 1];
    printf $oldout "%2d %-70.70s ", $i, "@$test";
    open(STDOUT, '>', "test$i.out") || die "Couldn't open test$i.out : $!\n";
    open(STDERR, '>&', \*STDOUT) || die "Couldn't dup STDOUT: $!\n";
    my $res = system {"$srcdir/$test->[0]"} @$test;
    close(STDERR) || die "Error writing to STDERR: $!\n";
    open(STDERR, '>&', $olderr) || die "Error restoring STDERR: $!\n";
    close(STDOUT) || die "Error writing to STDOUT: $!\n";
    open(STDOUT, '>&', $oldout) || die "Error restoring STDOUT: $!\n";

    if ($res) {
	$failed++;
    } else {
	$passed++;
    }

    printf $oldout "%s\n", $res ? "FAIL" : "PASS";
}

my $all = $passed + $failed;
my $msg = $failed ? "$failed of $all tests failed" : "All $all tests passed";
my $dashes = '=' x length($msg);
print $oldout "$dashes\n$msg\n$dashes\n";

exit($failed ? 1 : 0);
