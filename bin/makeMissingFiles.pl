#!/software/bin/perl

# given a cif, dif, bcl, scl or stats file check the other files exist
# if not create the missing files. Will also check a filter file exists

use strict;
use warnings;
use English qw(-no_match_vars);
use Carp;
use Getopt::Long;
use File::Spec;
use FileHandle;
use POSIX qw(ceil);

sub process;
sub usage;
sub initialise;

my $opts = initialise;

process;

exit;

# ----------------------------------------------------------------------
sub usage {

  ## no critic

  print STDERR "\n";
  print STDERR "makeMissingFiles.pl [-verbose] [-check] [-filter <filter>] [-cif <template>] [-bcl] [-nocif] [-nodif] [-noscl] <cif|dif|bcl|scl|stats>\n";
  print STDERR "\n";
  print STDERR "        -check           just report missing files\n";
  print STDERR "\n";
  print STDERR "        -olb             only check for the files required for OLB i.e. cif files\n";
  print STDERR "        -hiseqx          only check for the files created on a hiseqx i.e. bcl.gz but not cif, dif or stats\n";
  print STDERR "\n";
  print STDERR "        -cif <template>  if cif file is missing make one based on this template\n";
  print STDERR "\n";
  print STDERR "        -filter <filter> if filter file does not exist make this one\n";
  print STDERR "        -bcl             if bcl file is missing make one\n";
  print STDERR "\n";
  print STDERR "        -nocif           do not check for missing cif files\n";
  print STDERR "        -nodif           do not check for missing dif files\n";
  print STDERR "        -noscl           do not check for missing scl files\n";
  print STDERR "\n";

}

# ----------------------------------------------------------------------
sub process {

  my $file = $ARGV[0];
  $file = File::Spec->rel2abs($file);

  my ($intensities, $lane, $cycle, $tile);
  if( $file =~ m/\.cif$/ ){
    die "Invalid cif file $file\n" unless ($intensities, $lane, $cycle, $tile) = ($file =~ m/^(.+)\/L00(\d)\/(C\d+\.1)\/s_\2_(\d+)\.cif$/);
  } elsif( $file =~ m/\.dif$/ ){
    die "Invalid dif file $file\n" unless ($intensities, $lane, $cycle, $tile) = ($file =~ m/^(.+)\/L00(\d)\/(C\d+\.1)\/s_\2_(\d+)\.dif$/);
  } elsif( $file =~ m/\.bcl(.gz)?$/ ){
    die "Invalid bcl file $file\n" unless ($intensities, $lane, $cycle, $tile) = ($file =~ m/^(.+)\/BaseCalls\/L00(\d)\/(C\d+\.1)\/s_\2_(\d+)\.bcl(.gz)?$/);
  } elsif( $file =~ m/\.scl$/ ){
    die "Invalid scl file $file\n" unless ($intensities, $lane, $cycle, $tile) = ($file =~ m/^(.+)\/BaseCalls\/L00(\d)\/(C\d+\.1)\/s_\2_(\d+)\.scl$/);
  } elsif( $file =~ m/\.stats$/ ){
    die "Invalid stats file $file\n" unless ($intensities, $lane, $cycle, $tile) = ($file =~ m/^(.+)\/BaseCalls\/L00(\d)\/(C\d+\.1)\/s_\2_(\d+)\.stats$/);
  }else{
    print {*STDERR} "\nUnexpected file $file\n" or croak 'print failed';
    usage;
    exit 1;
  }

  my $ptile = sprintf("%04d", $tile);

  my ($cif, $dif, $bcl, $scl, $stats, $filter);
  $cif    = "$intensities/L00${lane}/${cycle}/s_${lane}_${tile}.cif";
  $dif    = "$intensities/L00${lane}/${cycle}/s_${lane}_${tile}.dif";
  $bcl    = "$intensities/BaseCalls/L00${lane}/${cycle}/s_${lane}_${tile}.bcl";
  $scl    = "$intensities/BaseCalls/L00${lane}/${cycle}/s_${lane}_${tile}.scl";
  $stats  = "$intensities/BaseCalls/L00${lane}/${cycle}/s_${lane}_${tile}.stats";
  $filter = "$intensities/BaseCalls/L00${lane}/s_${lane}_${ptile}.filter";

  my $fh = FileHandle->new();
  my $data;

  $/ = undef;

  # a cif file must exist unless the -nocif or the -cif or the hiseqx option was specified

  unless (exists($opts->{'nocif'}) || exists($opts->{'cif'}) || exists($opts->{'hiseqx'})) {
    if (exists($opts->{'check'})) {
      warn "Missing cif file $cif\n" unless -e $cif;
    }else{
      die "Missing cif file $cif\n" unless -e $cif;
    }
  }

  if (exists($opts->{'olb'})) {
    exit if (exists($opts->{'check'}));
  }else{
    # hiseqx bcl files are gzipped
    $bcl .= ".gz" if exists($opts->{'hiseqx'});

    # a bcl file must exist unless the -bcl option was specified
    if (exists($opts->{'check'}) || exists($opts->{'bcl'})) {
      warn "Missing bcl file $bcl\n" unless -e $bcl;
    }else{
      die "Missing bcl file $bcl\n" unless -e $bcl;
    }

    # a filter file must exist unless the -filter option was specified
    if (exists($opts->{'filter'})) {
      warn "Missing filter file $filter\n" unless -e $filter;
    }else{
      die "Missing filter file $filter\n" unless -e $filter;
    }

    # hiseqx we don't expect any other files
    unless (exists($opts->{'hiseqx'})) {
      print "Missing dif file $dif\n" unless ( $opts->{'nodif'} || -e $dif );
      print "Missing scl file $scl\n" unless ( $opts->{'noscl'} || -e $scl );
      print "Missing stats file $stats\n" unless -e $stats;
    }

    exit if (exists($opts->{'check'}));
  }

  exit if (exists($opts->{'hiseqx'}));

  unless ( $opts->{'nocif'} || -e $cif ){
    print "Writing zero intensity cif file $cif\n";

    my $intensities = $opts->{'cif'};
    print "Reading intensity file header\n";
    $fh->open("<$intensities") or croak "$intensities:\n$ERRNO";
    $data = <$fh>;
    $fh->close();
    my (@header) = unpack("C13", $data);
    my $dataType = $header[4];
    my $num_entries = $header[9] | $header[10] << 8 | $header[11] << 16 | $header[12] << 24;

    if (exists($opts->{'verbose'})) {
      printf("dataType    = $dataType\n");
      printf("num_entries = $num_entries\n");
    }

    die "Invalid cycle $cycle\n" unless $cycle =~ m/^C(\d+)\.1$/;
    my $firstCycle = $1;

    $header[5] = ($firstCycle      ) & 0xFFFF;
    $header[6] = ($firstCycle >> 8 ) & 0xFFFF;

    $data = pack("C13", @header);
    my @data = (0) x (4 * $num_entries);
    my %format = ("1" => "c*", "2" => "s*", "4" => "l*");
    $data .= pack($format{$dataType}, @data);

    # pad intensity file to a multiple of 4096
    my $s = length($data);
    $s = 4096 * (int($s / 4096) + 1) if $s % 4096;
    my $p = $s - length($data);
    $data .= pack($format{1}, (0) x $p) if $p > 0;

    $fh->open(">$cif") or croak "$cif:\n$ERRNO";
    $fh->write($data);
    $fh->close();
  }

  exit if (exists($opts->{'olb'}));

  # if a filter file exists read it otherwise make one
  if ( -e $filter ){
    print "Reading filter file $filter\n";
  }else{
    makeFilterFile($filter);
  }

  $fh->open("<$filter") or croak "$filter:\n$ERRNO";
  $data = <$fh>;
  $fh->close();

  my ($check, $version) = unpack("ll", $data);
  if( $check == 0 ){
    if( $version == 3 ){
      $data = substr($data,8);
    }else{
      croak "Unexpected version: $version\n";
    }
  }
  my ($nclusters_filtered, @filter) = unpack("lC*", $data);

  my $filtered = 0;
  map { $filtered++ if $_ & 1 } @filter;

  if (exists($opts->{'verbose'})) {
    printf("nclusters_filtered = $nclusters_filtered\n");
    printf("filtered = $filtered\n");
  }

  my $nclusters;
  my @bases = ();
  my @quals = ();
  if ( -e $bcl ){
    unless ( $opts->{'nobcl'} ){
      print "Reading bcl file $bcl\n";
      $fh->open("<$bcl") or croak "$bcl:\n$ERRNO";
      $data = <$fh>;
      $fh->close();

      my @data = ();
      ($nclusters, @data) = unpack("lC*", $data);

      unless ($nclusters == $nclusters_filtered) {
        die "Inconsistent bcl and filter files nclusters = $nclusters expected $nclusters_filtered\n";
      }

      if (exists($opts->{'verbose'})) {
        printf("nclusters = $nclusters\n");
      }
      map{ push(@bases, $_ & 3); push(@quals, $_ >> 2) } @data;
    }
  }else{
    if( $opts->{'rebasecall'} && -e $dif ) {
      print "Re-basecalling intensity file $dif\n";
      print "All qualities set to 2 except for bases called as N which must have a quality of 0\n";

      $fh->open("<$dif") or croak "$dif:\n$ERRNO";
      $data = <$fh>;
      $fh->close();

      my (@header) = unpack("C13", $data);
      my $magic = substr($data,0,4);
      my $dataType = $header[4];
      my $firstCycle = $header[5] | $header[6] << 8;
      my $num_cycles = $header[7] | $header[8] << 8;
      my $num_entries = $header[9] | $header[10] << 8 | $header[11] << 16 | $header[12] << 24;
      $data = substr($data,13);
      my %format = ("1" => "c*", "2" => "s*", "4" => "l*");
      my (@dif_data) = unpack($format{$dataType}, $data);

      unless ($num_entries == $nclusters_filtered) {
        die "Inconsistent dif and filter files nclusters = $num_entries expected $nclusters_filtered\n";
      }

      $nclusters = $nclusters_filtered;
      $data = pack("l", $nclusters);
      my @bcl_data = ();
      for (my $cluster=0; $cluster<$nclusters; $cluster++) {
        my ($max, $base) = (-99999, 0);
        for (my $channel=0; $channel<4; $channel++) {
          if( $dif_data[$num_entries * $channel + $cluster] > $max ) {
             $max = $dif_data[$num_entries * $channel + $cluster];
             $base = $channel;
          }
        }
        # where we have a basecall set the quality to 2 (B)
        my $qual = ($max ? 2 : 0);
        push(@bases, $base);
        push(@quals, $qual);
        push(@bcl_data, $base | ($qual << 2));
      }
      $data .= pack("C*", @bcl_data);
    }else{
      print "Writing zero bcl file $bcl\n";

      $nclusters = $nclusters_filtered;
      $data = pack("l", $nclusters);
      @bases = (0) x $nclusters;
      @quals = (0) x $nclusters;
      my @data = (0) x $nclusters;
      $data .= pack("C*", @data);
    }

    $fh->open(">$bcl") or croak "$bcl:\n$ERRNO";
    $fh->write($data);
    $fh->close();
  }

  unless ( $opts->{'nodif'} || -e $dif ){
    print "Writing zero intensity dif file $dif\n";

    # check the bcl file contains ALL N's <=> qualities ALL 0
    my $count = 0;
    map {$count++ if $_ > 0} @quals;
    die "Not all qualities are zero\n" if $count;

    print "Reading cif file header\n";
    $fh->open("<$cif") or croak "$cif:\n$ERRNO";
    $data = <$fh>;
    $fh->close();
    my (@header) = unpack("C13", $data);
    my $magic = substr($data,0,4);
    my $dataType = $header[4];
    my $firstCycle = $header[5] | $header[6] << 8;
    my $num_cycles = $header[7] | $header[8] << 8;
    my $num_entries = $header[9] | $header[10] << 8 | $header[11] << 16 | $header[12] << 24;

    unless ($num_entries == $nclusters_filtered) {
      die "Inconsistent cif and filter files nclusters = $num_entries expected $nclusters_filtered\n";
    }

    if (exists($opts->{'verbose'})) {
      printf("dataType    = $dataType\n");
      printf("num_entries = $num_entries\n");
    }

    $data = substr($data,0,13);
    my @data = (0) x (4 * $num_entries);
    my %format = ("1" => "c*", "2" => "s*", "4" => "l*");
    $data .= pack($format{$dataType}, @data);

    # pad intensity file to a multiple of 4096
    my $s = length($data);
    $s = 4096 * (int($s / 4096) + 1) if $s % 4096;
    my $p = $s - length($data);
    $data .= pack($format{1}, (0) x $p) if $p > 0;

    $fh->open(">$dif") or croak "$dif:\n$ERRNO";
    $fh->write($data);
    $fh->close();
  }

  unless ( $opts->{'noscl'} || -e $scl ){
    print "Writing zero scl file $scl\n";

    # check the bcl file contains ALL N's <=> qualities ALL 0
    my $count = 0;
    map {$count++ if $_ > 0} @quals;
    die "Not all qualities are zero\n" if $count;

    my $nbytes = POSIX::ceil($nclusters / 4);

    $data = pack("l", $nclusters);
    my @data = (0) x ($nbytes);
    $data .= pack("C*", @data);

    $fh->open(">$scl") or croak "$scl:\n$ERRNO";
    $fh->write($data);
    $fh->close();
  }

  unless ( -e $stats ){
    print "Writing stats file $stats\n";

    my $all = 0;
    my @all = (0,0,0,0);
    my @call = (0,0,0,0);
    my $num_all = 0;
    my @num_call = (0,0,0,0,0);

    if ( -e $dif ) {
      print "Reading intensity file $dif\n";
      $fh->open("<$dif") or croak "$dif:\n$ERRNO";
      $data = <$fh>;
      $fh->close();

      my (@header) = unpack("C13", $data);
      my $magic = substr($data,0,4);
      my $dataType = $header[4];
      my $firstCycle = $header[5] | $header[6] << 8;
      my $num_cycles = $header[7] | $header[8] << 8;
      my $num_entries = $header[9] | $header[10] << 8 | $header[11] << 16 | $header[12] << 24;
      $data = substr($data,13);
      my %format = ("1" => "c*", "2" => "s*", "4" => "l*");
      my (@data) = unpack($format{$dataType}, $data);

      if (exists($opts->{'verbose'})) {
        printf("dataType    = $dataType\n");
        printf("firstCycle  = $firstCycle\n");
        printf("num_cycles  = $num_cycles\n");
        printf("num_entries = $num_entries\n");
      }

      unless ($num_entries == $nclusters_filtered) {
        die "Inconsistent dif and filter files nclusters = $num_entries expected $nclusters_filtered\n";
      }

      print "Calculating stats\n";
      for (my $cluster=0; $cluster<$nclusters; $cluster++) {
        next unless $filter[$cluster] & 1;
        if( $quals[$cluster] ){
          # I've moved this inside the if( $quals[$cluster] ){}, so
          # $all and $all[] are only modifed when $num_all is
          # incremented. This shouldn't make any difference to the
          # stats since bases with quality=0 usually have all
          # intensities zero
          for (my $channel=0; $channel<4; $channel++) {
            $all += $data[$num_entries * $channel + $cluster];
            $all[$channel] += $data[$num_entries * $channel + $cluster];
          }
          $num_all++;
          $call[$bases[$cluster]] += $data[$num_entries * $bases[$cluster] + $cluster];
          $num_call[$bases[$cluster]]++;
        }else{
          $num_call[4]++;
        }
      }
      if( $num_all ) {
        $all /= ( 4 * $num_all);
        for (my $channel=0; $channel<4; $channel++) {
          $all[$channel] /= $num_all;
          next unless $num_call[$channel];
          $call[$channel] /= $num_call[$channel];
        }
      }

      $cycle = $firstCycle - 1;
    }else{
      print "Setting all stats to zero\n";

      $num_call[4] = $nclusters_filtered;

      die "Invalid cycle $cycle\n" unless $cycle =~ m/^C(\d+)\.1$/;
      my $firstCycle = $1;

      $cycle = $firstCycle - 1;
    }

    my @num_all = ($num_all) x 4;

    if (exists($opts->{'verbose'})) {
      printf("cycle   = $cycle\n");
      printf("all     = $all\n");
      printf("All_A   = $all[0]\n");
      printf("All_C   = $all[1]\n");
      printf("All_G   = $all[2]\n");
      printf("All_T   = $all[3]\n");
      printf("Call_A  = $call[0]\n");
      printf("Call_C  = $call[1]\n");
      printf("Call_G  = $call[2]\n");
      printf("Call_T  = $call[3]\n");
      printf("num_A   = $num_call[0]\n");
      printf("num_C   = $num_call[1]\n");
      printf("num_G   = $num_call[2]\n");
      printf("num_T   = $num_call[3]\n");
      printf("num_X   = $num_call[4]\n");
      printf("num_all = $num_all\n");
    }

    $data = pack("ldd4d4l5l4", $cycle,$all,@all,@call,@num_call,@num_all);

    $fh->open(">$stats") or croak "$bcl:\n$ERRNO";
    $fh->print($data);
    $fh->close();
  }

}

# ----------------------------------------------------------------------
sub initialise {

  my %opts;
  my $rc = GetOptions(\%opts, 'verbose', 'check', 'olb', 'hiseqx', 'rebasecall', 'cif=s', 'filter=s', 'bcl', 'nocif', 'nodif', 'noscl', 'help');
  if ( ! $rc) {
    print {*STDERR} "\nerror in command line parameters\n" or croak 'print failed';
    usage;
    exit 1;
  }

  if (exists($opts{'help'})) {
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
sub makeFilterFile {
  my ($filter) = @_;

  my $threshold = 0.6;

  print "Making filter file $filter\n";
  die "Filter file already exists\n" if -e $filter;

  my ($intensities, $lane, $tile);
  die "Invalid filter file $filter\n" unless ($intensities, $lane, $tile) = ($filter =~ m/^(.*\/Data\/Intensities).*\/s_(\d)_(\d+)\.filter$/);
  $tile =~ s/^0+//;

  my $nclusters_filtered = 0;
  my @filter = ();

  for (my $cycle=1; $cycle<=25; $cycle++) {
    my $dif = "${intensities}/L00${lane}/C${cycle}.1/s_${lane}_${tile}.dif";
    print "Reading dif file $dif\n";

    my $fh = FileHandle->new();
    $fh->open("<$dif") or croak "$dif:\n$ERRNO";
    my $data = <$fh>;
    $fh->close();

    my (@header) = unpack("C13", $data);
    my $dataType = $header[4];
    my $num_entries = $header[9] | $header[10] << 8 | $header[11] << 16 | $header[12] << 24;

    if( $cycle == 1 ){
      $nclusters_filtered = $num_entries;
      @filter = (0) x $num_entries;
    }else{
      unless ($num_entries == $nclusters_filtered) {
        die "Inconsistent dif file $num_entries expected $nclusters_filtered\n";
      }
    }

    $data = substr($data,13);
    my %format = ("1" => "c*", "2" => "s*", "4" => "l*");
    my (@data) = unpack($format{$dataType}, $data);

    for (my $cluster=0; $cluster<$num_entries; $cluster++) {
      my @intensities = ();
      for (my $channel=0; $channel<4; $channel++) {
        push(@intensities, $data[$num_entries * $channel + $cluster]);
      }
      @intensities = sort {$b<=>$a} @intensities;
      my $total = $intensities[0] + $intensities[1];
      my $chastity = $total ? ($intensities[0] / $total) : 0.0;
      $filter[$cluster]++ if $chastity < $threshold;
    }
  }

  @filter = map { $_ > 1 ? 0 : 1 } @filter;

  my $data = "";
  if( $filter =~ m/L00\d/ ){
    my ($check, $version) = (0, 3);
    $data .= pack("ll", $check, $version);
  }
  $data .= pack("l", $nclusters_filtered);
  $data .= pack("C*", @filter);

  my $fh = FileHandle->new();
  $fh->open(">$filter") or croak "$filter:\n$ERRNO";
  $fh->write($data);
  $fh->close();

  exit;
}

