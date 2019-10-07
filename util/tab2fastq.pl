#!/usr/bin/env perl

use strict;
use warnings;


if (@ARGV < 1) {
        print STDERR "<fastq tab file>\n";
        print STDERR "Will output a FASTQ file\n";
        exit;
}


open IN, $ARGV[0];
my $numLines = 0;
while (my $line = <IN>) {
        $numLines++;
        chomp $line;
        my @data = split("\t", $line);
        print "@",$data[0],"\n",$data[1],"\n+\n",'I' x length($data[1]),"\n";
}
close IN;
print STDERR "\tFound $numLines sequences\n";
