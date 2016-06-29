#!/usr/bin/perl
# mgf editor
# author: Dima Ischenko
# version: 0.2
# date start: 19.09.2014
# current date: 16.02.2016

use strict;
use warnings;
use List::Util qw/shuffle/;

sub phelp {
  print "\n mgf editor v 0.2.\n";
  print " Author: Dima Ischenko.\n";

  print "\n";

    print "\t in <id file> <mgf>         :  subset spectra from mgf that in list.\n";
    print "\t ni <id file> <mgf>         :  subset spectra from mgf that not in list.\n";
    print "\t id <mgf> [start]           :  convert TITLE to numbers [start = 1].\n";
    print "\t seq <id_seq file> <mgf>    :  add SEQ field to mgf from seqfile (int tab format).\n";
    print "\t seqc <di_seq_charge> <mgf> :  add SEQ and CHARGE field to mgf from file in tab format.\n";
    print "\t mseq <mascot.csv> <mgf>    :  add SEQ fiield to mgf from mascot result csv file\n";
    print "\t                             * \"pep_scan_title\" and \"pep_seq\" columns required.\n";
    print "\t split <mgf> <N>            :  split mgf into parts by N spectra in each part.\n";
    print "\t rnd <mgf> <N>              :  select random N spectra from mgf.\n";
    print "\t rnds <mgf> <k>             :  select k random spectra for each peptide.\n";
    print "\t filt <mgf>                 :  filter mgf. remove empty lines and comments.\n";
    print "\t inseq <seq file> <mgf>     :  subset spectra from mgf that has seq in list.\n";
    print "\t niseq <seq file> <mgf>     :  subset spectra from mgf that hasn't seq in list.\n";

    print "\n";
    exit;
}

if (@ARGV == 0) {
  phelp;
}

# in ni command
if ($ARGV[0] eq "in" || $ARGV[0] eq "ni") {
  if (@ARGV != 3) {
    phelp;
  }

  my $inpp = $ARGV[0];

  # read titles in hash
  my %titles = ();
  open my $f_in, '<', $ARGV[1] or die $!;

  while(<$f_in>) {
      $_ =~ s/\n|\r//g;
      $titles{$_} = 1;
  }

  close $f_in;

  open my $mgf, '<', $ARGV[2] or die $!;
  # check mgf spectra titles for existing in list
  my $ex = 0;
  while(<$mgf>) {
    $_ =~ s/\n|\r//g;
    if ($_ =~ /BEGIN IONS/) {
      $_ = <$mgf>;
      $_ =~ s/\n|\r//g;
      my @el = split(/=/, $_);
      
      if (($inpp eq "in" && exists $titles{$el[1]}) || ($inpp eq "ni" && ! exists $titles{$el[1]})) {
        $ex = 1;
        print "BEGIN IONS\n";
        print $_."\n";
      } else {
        $ex = 0;
      }
    } else {
    if ($ex == 1) {
      print $_."\n";
    }
    } 
  }

close $mgf;
} elsif ($ARGV[0] eq "id") {
  phelp if (@ARGV < 2);
  my $sid = 1;
  if (@ARGV == 3) {
    $sid = $ARGV[2];
  } elsif (@ARGV > 3) {
    phelp;
  }

  open my $mgf, '<', $ARGV[1] or die $!;
  while(<$mgf>) {
    if ($_ =~ /TITLE/) {
      $_ =~ s/=.*/=$sid/;
      $sid++;
    }
    print $_;
  }
  close $mgf;

} elsif ($ARGV[0] eq "seq") {
  if (@ARGV != 3) {
    phelp;
  }

  # read seq
  my %ids = ();
  open my $ps, '<', $ARGV[1] or die $!;
  while(<$ps>) {
    $_ =~ s/\n|\r//g;
    my @el = split(/\t/, $_);
    $ids{$el[0]} = $el[1];
  }
  close $ps;

  open my $mgf, '<', $ARGV[2] or die $!;
  while(<$mgf>) {
    if ($_ =~ /TITLE/) {
      $_ =~ s/\n|\r//g;
      my @el = split(/=/, $_);
      print $_."\n";
      if (exists $ids{$el[1]}) {
        print "SEQ=".$ids{$el[1]}."\n";
      }
    } elsif ($_ !~ /SEQ=/) {
      print $_;
    }
  }
close $mgf;
} elsif ($ARGV[0] eq "seqc") {
  if (@ARGV != 3) {
    phelp;
  }

  # read seq and charge
  my %ids = ();
  my %idc = ();
  open my $ps, '<', $ARGV[1] or die $!;
  while(<$ps>) {
    $_ =~ s/\n|\r//g;
    my @el = split(/\t/, $_);
    $ids{$el[0]} = $el[1];
    $idc{$el[0]} = $el[2];
  }
  close $ps;

  open my $mgf, '<', $ARGV[2] or die $!;
  while(<$mgf>) {
    if ($_ =~ /TITLE/) {
      $_ =~ s/\n|\r//g;
      my @el = split(/=/, $_);
      print $_."\n";
      if (exists $ids{$el[1]}) {
        print "SEQ=".$ids{$el[1]}."\n";
        print "CHARGE=".$idc{$el[1]}."+\n";
      }
    } elsif ($_ !~ /SEQ=/) {
      print $_;
    }
  }
close $mgf;
} elsif ($ARGV[0] eq "mseq") {
  if (@ARGV != 3) {
    phelp;
  }
  # line by line flag
  my $lbl = 0;
  # read seq
  my %ids = ();
  # hash for column names
  my %clmns = ();
  # line by line
  open my $ps, '<', $ARGV[1] or die $!;
  while(<$ps>) {
    $_ =~ s/\n|\r//g;
    # read column names
    if (!$lbl && $_ =~ /pep_scan_title/) { 
      $lbl = 1;
      my @el = split(/,/, $_);
      for (my $i = 0; $i < @el; $i++) {
        $clmns{$el[$i]} = $i;
      }
      next;
    }
    next if (!$lbl);
    my @el = split(/,/, $_);
    $el[$clmns{"pep_scan_title"}] =~ s/\"//g;
    $el[$clmns{"pep_scan_title"}] =~ s/~/\"/g;
    $ids{$el[$clmns{"pep_scan_title"}]} = $el[$clmns{"pep_seq"}];
  }
  close $ps;

  open my $mgf, '<', $ARGV[2] or die $!;
  # check mgf spectra titles for existing in list
  my $ex = 0;
  while(<$mgf>) {
    $_ =~ s/\n|\r//g;
    if ($_ =~ /BEGIN IONS/) {
      $_ = <$mgf>;
      $_ =~ s/\n|\r//g;
      my @el = split(/=/, $_);
      
      if (exists $ids{$el[1]}) {
        $ex = 1;
        print "BEGIN IONS\n";
        print $_."\n";
        print "SEQ=".$ids{$el[1]}."\n";
      } else {
        $ex = 0;
      }
    } else {
    if ($ex == 1) {
      print $_."\n";
    }
    } 
  }
close $mgf;
} elsif ($ARGV[0] eq "split") {
  if (@ARGV != 3) {
    phelp;
  }

  my $mgfn = $ARGV[1];
  my $pn = 0;
  my $N = $ARGV[2];

  open my $mgf, '<', $mgfn or die $!;
  my $fout;
  my $inc = 0;

  my $ions = $N;
  while (<$mgf>) {
    if ($ions >= $N) {
          $ions = 0;
          close $fout if $inc;
          $pn++;
          open $fout, '>', $mgfn."_part_".$pn.".mgf" or die $!;
        $inc = 1;
      }

      $_ =~ s/\n|\r//g;
      if ($_ =~ 'END') {
          $ions++;
      }
      print $fout $_."\n";
  }

  close $fout;
  close $mgf;
} elsif ($ARGV[0] eq "rnd") {
  if (@ARGV != 3) {
    phelp;
  }

  my $N = $ARGV[2];
  open my $nl, '<', $ARGV[1] or die $!;
  my $ns = 0;
  while (<$nl>) { $ns++ if $_ =~ /TITLE/};
  close $nl;
  if ($ns < $N) {
    print "N more than total number of spectra in file\n";
    exit;
  }

  # random spectra
  my %rs = map {$_ => 1} (shuffle 1..$ns)[0..($N - 1)];
  open my $mgf, '<', $ARGV[1] or die $!;

  my $cs = 0;
  while (<$mgf>) {
    if ($_ =~ /BEGIN ION/) {
      $cs++;
      if (exists $rs{$cs}) {
        print $_;
        while (<$mgf>) {
          print $_;
          last if ($_ =~ /END ION/);
        }
      }
    }
  }
  close $mgf;

}
elsif ($ARGV[0] eq "rnds") {
  if (@ARGV != 3) {
    phelp;
  }

  my $k = $ARGV[2];

  my %seqsp = ();
  my $id = 0;

  open my $nl, '<', $ARGV[1] or die $!;
  while (<$nl>) { 
    chomp $_;
    if ($_ =~ /TITLE/) {
      my @el = split(/=/, $_);
      $id = $el[1];
    }
    if ($_ =~ /SEQ/) {
      my @el = split(/=/, $_);
      my $sq = $el[1];
      if (!(exists $seqsp{$sq})) {
        $seqsp{$sq} = ();
        push(@{$seqsp{$sq}}, $id);
      } else {
        push(@{$seqsp{$sq}}, $id);
      }
    }
  }
  close $nl;

  # subset k spectra for each seq
  my %titles = ();
  foreach my $skey (keys %seqsp) {
    my @sps = @{$seqsp{$skey}};

    if (@sps <= $k) {
      for (my $i = 0; $i < @sps; $i++) {
        $titles{$sps[$i]} = 1;
      }
    } else {
      my @sind = (shuffle 0..(@sps - 1))[0..($k - 1)];
      for (my $i = 0; $i < @sind; $i++) {
        $titles{$sps[$sind[$i]]} = 1;
      }
    }
  }

  for my $sk (keys %titles) {
    #print $sk."\n";
  }

  open my $mgf, '<', $ARGV[1] or die $!;
  # check mgf spectra titles for existing in list
  my $ex = 0;
  while(<$mgf>) {
    $_ =~ s/\n|\r//g;
    if ($_ =~ /BEGIN IONS/) {
      $_ = <$mgf>;
      $_ =~ s/\n|\r//g;
      my @el = split(/=/, $_);
      
      if (exists $titles{$el[1]}) {
        $ex = 1;
        print "BEGIN IONS\n";
        print $_."\n";
      } else {
        $ex = 0;
      }
    } else {
    if ($ex == 1) {
      print $_."\n";
    }
    } 
  }
  close $mgf;
}
elsif ($ARGV[0] eq "filt") {
  if (@ARGV != 2) {
    phelp;
  }
  open my $mgf, '<', $ARGV[1] or die $!;
  while(my $line = <$mgf>) {
    $line =~ s/\n|\r//g;
    next if ($line =~ /^#/ || $line eq "");
    if ($line =~ /(PEPMASS=.*)[ |\t].*/) {
      $line = $1;
    }   
    
    $line =~ s/\t/ /g;
    print $line."\n";
  }
  close $mgf;
} elsif ($ARGV[0] eq "inseq" || $ARGV[0] eq "niseq") {
  if (@ARGV != 3) {
    phelp;
  }
  my $inpp = $ARGV[0];

  # read seq in hash
  my %seq = ();
  open my $f_in, '<', $ARGV[1] or die $!;

  while(<$f_in>) {
      $_ =~ s/\n|\r//g;
      $seq{$_} = 1;
  }

  close $f_in;

  open my $mgf, '<', $ARGV[2] or die $!;
  # check mgf spectra seq for existing in list
  my $ex = 0;
  while(<$mgf>) {
    $_ =~ s/\n|\r//g;
    if ($_ =~ /BEGIN IONS/) {
      my $title = <$mgf>;
      $_ = <$mgf>;
      $_ =~ s/\n|\r//g;
      my @el = split(/=/, $_);
      
      if (($inpp eq "inseq" && exists $seq{$el[1]}) || ($inpp eq "niseq" && ! exists $seq{$el[1]})) {
        $ex = 1;
        print "BEGIN IONS\n";
        print $title;
        print $_."\n";
      } else {
        $ex = 0;
      }
    } else {
    if ($ex == 1) {
      print $_."\n";
    }
    } 
  }

close $mgf;
  }
else {
  phelp;
}
