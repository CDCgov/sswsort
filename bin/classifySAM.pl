#!/usr/bin/env perl

# Filename:         classifySAM
# Description:      Classifies gene segments / genomes by Smith-Waterman scores.
#
# Date dedicated:   2022-07-19
# Author:           Samuel S. Shepard, Centers for Disease Control and Prevention
#
# Citation:         Unpublished
#
# =============================================================================
#
#                            PUBLIC DOMAIN NOTICE
#
#  This source code file or script constitutes a work of the United States
#  Government and is not subject to domestic copyright protection under 17 USC §
#  105. This file is in the public domain within the United States, and
#  copyright and related rights in the work worldwide are waived through the CC0
#  1.0 Universal public domain dedication:
#  https://creativecommons.org/publicdomain/zero/1.0/
#
#  The material embodied in this software is provided to you "as-is" and without
#  warranty of any kind, express, implied or otherwise, including without
#  limitation, any warranty of fitness for a particular purpose. In no event
#  shall the Centers for Disease Control and Prevention (CDC) or the United
#  States (U.S.) government be liable to you or anyone else for any direct,
#  special, incidental, indirect or consequential damages of any kind, or any
#  damages whatsoever, including without limitation, loss of profit, loss of
#  use, savings or revenue, or the claims of third parties, whether or not CDC
#  or the U.S. government has been advised of the possibility of such loss,
#  however caused and on any theory of liability, arising out of or in
#  connection with the possession, use or performance of this software.
#
#  Please provide appropriate attribution in any work or product based on this
#  material.

use strict;
use warnings;
use English qw(-no_match_vars);
use Carp qw(croak);

my $hitThreshold       = 800;
my $normScoreThreshold = 0;
my $scoreThreshold     = 0;
my $lengthMinimum      = 0;
my $factor             = 2;

my ( $useMatches,, $interleavedPairs, $altInput ) = ( 0, 0, 0 );
my $lengthMedians = q{};

use Getopt::Long;
GetOptions(
            'use-matches|M'          => \$useMatches,
            'interleaved-pairs|P'    => \$interleavedPairs,
            'alternate-input|A'      => \$altInput,
            'hit-threshold|T=i'      => \$hitThreshold,
            'length-medians|L=s'     => \$lengthMedians,
            'score-minimum|S=i'      => \$scoreThreshold,
            'norm-score-minimum|N=f' => \$normScoreThreshold,
            'length-minimum|G=i'     => \$lengthMinimum
);

if ( scalar @ARGV != 2 ) {
    die("Usage:\n\tperl $PROGRAM_NAME <sam> <fasta>\n\n");
}

# Get Length Medians
my %geneThresholds = ();
if ( $lengthMedians ne q{} ) {
    my $IN;
    open( $IN, '<', $lengthMedians ) or die("Cannot open $lengthMedians for reading.\n");
    local $RS = "\n";
    while ( my $line = <$IN> ) {
        chomp($line);
        my ( $gene, $length ) = split( "\t", $line );
        if ( $length ne q{} ) {
            $geneThresholds{$gene} = $length * $factor;
        }
    }
    close $IN or croak("Cannot open $lengthMedians: $OS_ERROR\n");
}

# FUNCTIONS #
sub countMatch($) {
    my ($cig) = @_;
    my $count = 0;
    while ( $cig =~ /(\d+)([MIDNSHP])/gsmx ) {
        if ( $2 eq ' M ' ) {
            $count += $1;
        }
    }
    return $count;
}

#############

# Open fasta file
local $RS = '>';
my %not_found = ();
my %lengths   = ();
open( my $IN, ' < ', $ARGV[1] ) or croak("$PROGRAM_NAME ERROR: Cannot open $ARGV[1].\n");
while ( my $fasta_record = <$IN> ) {
    chomp($fasta_record);
    my @lines = split( /\r\n|\n|\r/smx, $fasta_record );
    my $id    = shift(@lines);
    if ( defined $id && $id ne q{} ) {
        my $sequence = uc( join( q{}, @lines ) );
        $sequence =~ tr/N//d;
        $lengths{$id}   = length($sequence);
        $not_found{$id} = 1;
    }
}
close $IN or croak("Cannot close $ARGV[1]: $OS_ERROR");

# Open SAM file
local $RS = "\n";

my $previous             = q{};
my $header               = q{};
my $pair                 = 0;
my %scoreByQuery         = ();
my %annotByQuery         = ();
my %strandByQuery        = ();
my %annotsAboveThreshold = ();

open( my $SAM, ' < ', $ARGV[0] ) or die("Cannot open $ARGV[0] for reading.\n");
while ( my $line = <$SAM> ) {
    if ( substr( $line, 0, 1 ) eq '@' ) {
        if ( $line ne $previous ) {
            $header .= $line;
        }
        $previous = $line;
        next;
    }

    chomp($line);

    my ( $qname, $flag, $rn, $pos, $mapq, $cigar, $mrnm, $mpos, $isize, $seq, $qual, $AS );
    if ($altInput) {
        ( $qname, $flag, $rn, $pos, $cigar, $AS ) = split( "\t", $line );
    } else {
        ( $qname, $flag, $rn, $pos, $mapq, $cigar, $mrnm, $mpos, $isize, $seq, $qual, $AS ) = split( "\t", $line );
    }

    # Find strand
    my $strand = ( $flag & 16 ) ? '-' : '+';

    # Reference to annotation
    $rn =~ s/#.+$//smx;
    if ( $rn =~ m/\s*(.+?)\{([^}]+)\}\s*$/smx ) {
        $rn = $2;
    }

    if ($interleavedPairs) {
        $qname = $qname . '_' . ( $pair % 2 );
        $pair++;
    }

    if ( $cigar eq '*' ) {
        next;
    }

    my $score = 0;
    if ($useMatches) {
        $score = countMatch($cigar);
    } elsif ( $AS =~ /AS:\w:(\d+)/smx ) {
        $score = $1;
    } else {
        print STDERR "Warning, using simple matches as back-up: $qname\n";
        $score = countMatch($cigar);
    }

    if ( $score >= $hitThreshold ) {
        $annotsAboveThreshold{$qname}{$rn} = $pos;
    }

    if ( !defined( $scoreByQuery{$qname} ) || $scoreByQuery{$qname} < $score ) {
        $scoreByQuery{$qname}  = $score;
        $annotByQuery{$qname}  = $rn;
        $not_found{$qname}     = 0;
        $strandByQuery{$qname} = $strand;
    }
}
close $SAM or croak("Cannot close file $ARGV[0]: $OS_ERROR");

foreach my $key ( keys(%not_found) ) {
    if ( $not_found{$key} ) {
        print STDOUT $key, "\tUNRECOGNIZABLE\t0\t$lengths{$key}\t*\n";
    } else {

        # Re-work required to sort by query position
        my @annots = sort( keys( %{ $annotsAboveThreshold{$key} } ) );

        if ( scalar @annots > 1 && $hitThreshold > 0 ) {
            print STDOUT $key, "\t*", join( '+', @annots ), "\t-1\t", $lengths{$key}, "\t*\n";
        } else {
            my ( $score, $annot, $strand ) = ( $scoreByQuery{$key}, $annotByQuery{$key}, $strandByQuery{$key} );
            my $length    = $lengths{$key};
            my $normScore = $length > 0 ? $score / $length : 0;
            my $aster     = q{};

            if ( ( $length < $lengthMinimum ) || ( $score < $scoreThreshold && $normScore < $normScoreThreshold ) ) {
                print STDOUT $key, "\tUNRECOGNIZABLE\t$score\t$length\t$strand\n";
            } else {
                my $gene = ( split( '_', $annot ) )[1];
                if ( defined $gene && defined $geneThresholds{$gene} && $length >= $geneThresholds{$gene} ) {
                    $aster = '*';
                }
                print STDOUT $key, "\t", $aster, $annot, "\t", $score, "\t", $length, "\t", $strand, "\n";
            }
        }
    }
}
