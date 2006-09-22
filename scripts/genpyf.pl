#!/usr/bin/perl -w
# File: genpyf.pl
#
# This file is part of the SpanLib library.
# Copyright (C) 2006  Stephane Raynaud
# Contact: stephane dot raynaud at gmail dot com
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

use strict;

# Inputs
my ($wrapper, $pyf, $precision, $simple, $double) = @ARGV;
if($double eq "") {
	print "Usage: genpyf.pl <f90_wrapper_file_name> <pyf_file_name> {simple|double} <uint:simple> <uit:double>\n";
	exit;
}

# Inits
open(WRAPPER, $wrapper);
open(PYF, "> $pyf");
print PYF "! -*- Mode: f90 -*-\n\n";
my $inside = 0;
my $line;
my $iPrecision = ($precision eq "double")?$double:$simple;

# Loop
while(<WRAPPER>){
	if(/^[\s\t]*(subroutine.*)$/i) {
		if($1 =~ /(.*)&[ \t]*$/) {
			print PYF "$1";
			$line = <WRAPPER>;
			$line =~ s/&(.*)/$1/;
			print PYF $line;
		} else {
			print PYF "$1\n";
		}
		while(<WRAPPER>){
			if(/external/i) {
				$inside = 1;
			} elsif($inside==1){
				if(/(real|integer)/i) {
					$_ =~ s/\(wp\)/($iPrecision)/;
					if($_ =~ /(.*)&[ \t]*$/) {
						print PYF $1;
						$line = <WRAPPER>;
						$line =~ s/&(.*)/$1/;
						print PYF $line;
					} else {
						print PYF $_;
					}
				} elsif(/^[\s\t]*$/) {
					$inside=0;
					last;
				}
			}
		}
	} elsif(/^[\s\t]*(end subroutine.*)$/i) {
		print PYF "$1\n\n";
	}
}
close(PYF);
close(WRAPPER);



