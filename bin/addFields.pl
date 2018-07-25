#!/usr/bin/env perl

# Filename:         addFields
# Description:      Adds constant fields to table.
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

my ( $program, $extra, $version, $module ) = ( q{}, q{}, q{}, q{} );
my ( $file, $lineNumber );
use Getopt::Long;
GetOptions(
            'program|P=s'   => \$program,
            'extra|E=s'     => \$extra,
            'version|V=s'   => \$version,
            'module|M=s'    => \$module,
            'file|F'        => \$file,
            'line-number|L' => \$lineNumber
);

my $d = "\t";
my ( $theFile, $line ) = ( q{}, q{} );
my $prefix = $program ne q{} ? $program . $d : q{};
$prefix .= $version ne q{} ? $version . $d : q{};
$prefix .= $module ne q{}  ? $module . $d  : q{};
my $suffix = $extra ne q{} ? $d . $extra : q{};

if ( !defined $file ) {
    if ( defined $lineNumber ) {
        my $i = 0;
        while ( $line = <> ) {
            chomp($line);
            $i++;
            print STDOUT $i, $d, $prefix, $line, $suffix, "\n";
        }
    } else {
        while ( $line = <> ) {
            chomp($line);
            print STDOUT $prefix, $line, $suffix, "\n";
        }
    }
} else {
    foreach $theFile (@ARGV) {
        my $IN;
        open( $IN, '<', $theFile ) or die("Cannot open $theFile for reading: $OS_ERROR\n");
        if ( defined $lineNumber ) {
            my $i = 0;
            while ( $line = <$IN> ) {
                chomp($line);
                $i++;
                print STDOUT $theFile, $d, $i, $d, $prefix, $line, $suffix, "\n";
            }
        } else {
            while ( $line = <$IN> ) {
                chomp($line);
                print STDOUT $theFile, $d, $prefix, $line, $suffix, "\n";
            }
        }
        close($IN) or die("Could not close $theFile: $OS_ERROR\n");
    }
}
