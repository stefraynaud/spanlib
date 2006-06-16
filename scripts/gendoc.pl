#!/usr/bin/perl -w
# File: gendoc.pl
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
use File::Basename;

# Inputs
my ($xmldir, $f90_library, $f90_example) = @ARGV;

# Basic declarations
my (@partNames, $f2xPersonalised, $partName, %parts, $line, $subroutineName, $arguments, $inside);
my ($vars, $type, $intent, $optional, $allocatable, $var, $name, @argDesc, $desc, $opt, %argRef, %arg);
my ($tmp, $tmp2, $i);

# Look
my %cols = (
	'name' => {
		'class' => 'b',
		'width' => '15%'
	},
	'type' => {
		'class' => 's',
		'width' => '15%'
	},
	'long_name' => {
		'class' => '',
		'width' => '705'
	}
);
my %f2xCssColors = (
	'personalised'	=> 'b',
	'comments'		=> 'e i',
	'digits' 		=> 'r',
	'strings'		=> 'o',
	'attributes'	=> 'gr',
	'functions'		=> 'cy',
	'routines'		=> 'vi'
);
my %f2xTags = (
	'functions'		=> 'dot_product|trim|eoshift|modulo|ssyev|not|allocated|allocate|spread|present|sum|matmul|sqrt|transpose|deallocate|pack|unpack|present',
	'routines'		=> 'print|module|contains|interface|logical|subroutine|if|then|else|endif|enddo|do|end|close|call|program|integer|real|open|write|character|use|implicit',
	'attributes'	=> 'procedure|parameter|optional|intent|allocatable|none'
);

# xml header
my $xmlHeader = "<?xml version=\"1.0\"?>\n<!DOCTYPE article PUBLIC \"-//KDE//DTD DocBook XML V4.2-Based Variant V1.1//EN\"\n\"/usr/share/sgml/docbook/xml-dtd-4.2-1.0-17.2/docbookx.dtd\">\n\n";
	
# Files
my %xmlFiles = (
	'subroutines'	=> "$xmldir/doc_subroutines.xml",
	'library'		=> "$xmldir/f90_library.xml",
	'example'		=> "$xmldir/f90_example.xml"
);



# Useful subroutines
sub gen_init_args {
	my $type = shift ;
	print XML_SUBROUTINES "\t<simplesect>\n\t\t<title>$type arguments</title>\n\t\t<informaltable pgwide=\"1\" label=\"\">\n\t\t\t<tgroup cols=\"3\">\n";
	
	my $i = 1;
	for $tmp ('name','type','long_name') {
		print XML_SUBROUTINES "\t\t\t\t<colspec colname=\"col$i\" colnum=\"$i\" width=\"$cols{$tmp}{'width'}\"/>\n";
		$i++;
	}
	print XML_SUBROUTINES "\t\t\t\t<tbody>\n";
}
sub f90toxml {
	my $code = shift;
	my $comment;
	$code =~ s/\&/\&amp;/g;
	$code =~ s/</\&lt;/g;
	$code =~ s/>/\&gt;/g;
	$comment = "";
	$comment = "<phrase role=\"$f2xCssColors{'comments'}\">$2<\/phrase>" if $code =~ s/^([^!]*)(!.*)$/$1/;
	$code =~ s/(['"][^'"]*['"])/<phrase role="$f2xCssColors{'strings'}">$1<\/phrase>/g;
	$code =~ s/\b(\d+)\b/<phrase role="$f2xCssColors{'digits'}">$1<\/phrase>/g;
	$code =~ s/\b($f2xPersonalised)\b/<phrase role="$f2xCssColors{'personalised'}">$1<\/phrase>/g;
	foreach $type ('attributes', 'functions', 'routines') {
		$code =~ s/\b($f2xTags{$type})\b/<phrase role="$f2xCssColors{$type}">$1<\/phrase>/g;
	}
	$code =~ s/(\n)/$comment$1/;
	return $code;
}

# First loop to get subroutine names
$f2xPersonalised="";
open(F90_LIBRARY,$f90_library);
while(<F90_LIBRARY>){
	if(/^[\s\t]*subroutine\s+([\w]+)\b/i) {
		$f2xPersonalised .= "|" if $f2xPersonalised ne "";
		$f2xPersonalised .= "$1";
	}
}
close(F90_LIBRARY);

# Generate the subroutine help file by parsing the f90 sources
open(XML_LIBRARY,"> $xmlFiles{'library'}");
open(XML_SUBROUTINES,"> $xmlFiles{'subroutines'}");
open(F90_LIBRARY,$f90_library);
print XML_LIBRARY "$xmlHeader<programlisting>";
print XML_SUBROUTINES "\t<sect1 id=\"doc_f90subs\">\n\t\t<title>F90 subroutines</title>\n";
while(<F90_LIBRARY>){

	# F90 to XML
	print XML_LIBRARY f90toxml($_);
	
	# A subroutine starts here
	if(/^[\s\t]*subroutine\s+([\w]+)\s*\(([^\)]+)\)/i){
		undef %parts;
		undef @argDesc;
		undef %argRef;
		@partNames = ();
		$partName = "";
		$subroutineName = "$1";
		$arguments = "$2";

		# Initialize arguments description whatever they are
		$arguments =~ s/[\t ]*//g;
		$i = 0;
		for $tmp (split(",",$arguments)) {
			$argRef{$tmp} = $i;
			$argDesc[$i]{'name'} = "$tmp";
			$argDesc[$i]{'type'} = "";
			$argDesc[$i]{'long_name'} = "";
			$argDesc[$i]{'optional'} = 0;
			$i++;
		}
	
		# Parse header
		$inside=0;
		HEADERLOOP: while(<F90_LIBRARY>) {
			if(/^[ \t]*!/) {
				$inside = 1;
				if(/! ([ \w]+):$/) {
					$partName = "$1";
					push(@partNames, "$1") if $partName ne "Title";
					$parts{$partName} = "";	
				} elsif(/!\t(.+)$/ && $partName ne "") {
					$line = $1;
					if ($line =~ /^\- ([^ :]+):(.+)/){
						$tmp = $argRef{"$1"};
						$argDesc[$tmp]{'long_name'} = "$2";
					} elsif($parts{$partName} eq "") {
						$parts{$partName} = "$line"
					} else {
						$parts{$partName} = "$parts{$partName}\n$line";
					}
				}
			} elsif($inside==1) {
				$inside = 0;
				last HEADERLOOP;
			}
		}

		# Parse declarations of external subroutine arguments
		$inside = 0;
		DECLLOOP: while(<F90_LIBRARY>) {
			if(/^[\s\t]*! External/) {
				$inside = 1;
			} elsif($inside==1){
				if(/(real|integer),[ \t]*intent\((in|out|inout)\)[^:]*::(.+)$/i) {
					$desc = "[${2}put,$1]:";
					$vars = "$3";
					$opt = /optional/;
					$vars =~ s/ +//g;
					$vars =~ s/(\([:,]+\))?//g;
					for $var (split(',',$vars)) {
						$tmp = $argRef{$var};
						$argDesc[$tmp]{'type'} = $desc;
						$argDesc[$tmp]{'optional'} = $opt;
						$arguments =~ s/([(, ]+)$var([), ]*)/$1$var=$var$2/ if $opt == 1
						
					}
				} elsif(/^[ \t]*$/){
					$inside=0;
					last DECLLOOP;
				}
			}
		}

		# Now, generate the xml block
		if ($partName ne "") {
			print XML_SUBROUTINES "<sect2 id=\"$subroutineName\">\n\t<title>$parts{'Title'}: <literal role=\"$f2xCssColors{'personalised'}\">$subroutineName</literal></title>\n"; #\t\t<para role=\"sh3\">\n\n";
			$arguments =~ s/(^|[,=\s]+)(\w+)\b(?!=)/$1<phrase role=\"$f2xCssColors{'personalised'}\">$2<\/phrase>/g;
			$arguments =~ s/,/, /g;
			# Usage part
			print XML_SUBROUTINES "\t<simplesect>\n\t\t<title>Usage</title>\n\t\t<programlisting>call <phrase role=\"$f2xCssColors{'personalised'}\">$subroutineName</phrase>($arguments)</programlisting>\n\t</simplesect>\n";
			# Other parts
			for $partName (@partNames){
				# Dependencies part
				if($partName eq "Dependencies"){
					if($parts{$partName} !~ /(LAPACK|BLAS)/){
						$parts{$partName} =~ s/([\w]+)/<link linkend=\"$1\">$1<\/link>/g;
						$parts{$partName} = "\t\t\t<phrase role=\"$f2xCssColors{'personalised'}\"><literal>$parts{$partName}</literal></phrase>";
					} else {
						$parts{$partName} = "<literal>$parts{$partName}</literal>";
					}
				}
				# Title of the part
				print XML_SUBROUTINES "\t<simplesect>\n\t\t<title>$partName</title>\n\t\t<para>\n$parts{$partName}\n\t\t</para>\n" if $parts{$partName} ne "";
				# Description part
				if($partName eq "Description"){
					print XML_SUBROUTINES "\t</simplesect>\n";
					$tmp2 = "Necessary";
					gen_init_args($tmp2);
					my $n = @argDesc;
					for ($i=0;$i<$n;$i++) {
						if($argDesc[$i]{'optional'}==1 && $tmp2 eq "Necessary") {
							print XML_SUBROUTINES "\t\t\t\t</tbody>\n\t\t\t</tgroup>\n\t\t</informaltable>\n\t</simplesect>\n";
							$tmp2 = "Optional";
							gen_init_args($tmp2);
						}
						print XML_SUBROUTINES "\t\t\t\t\t<row>\n";
						for $tmp ('name','type','long_name') {
							print XML_SUBROUTINES "\t\t\t\t\t\t<entry";
							print XML_SUBROUTINES " role=\"$cols{$tmp}{'class'}\"" if $cols{$tmp}{'class'} ne "";
#							print XML_SUBROUTINES " width=\"$cols{$tmp}{'width'}%\">";
							print XML_SUBROUTINES ">";
							print XML_SUBROUTINES "<varname>" if ($tmp eq 'name');
							print XML_SUBROUTINES $argDesc[$i]{$tmp};
							print XML_SUBROUTINES "</varname>" if ($tmp eq 'name');
							print XML_SUBROUTINES "</entry>\n";
						}
						print XML_SUBROUTINES "\t\t\t\t\t</row>\n";
					}
					print XML_SUBROUTINES "\t\t\t\t</tbody>\n\t\t\t</tgroup>\n\t\t</informaltable>\n\t</simplesect>\n";
				}
				print XML_SUBROUTINES "\t</simplesect>\n" if $partName eq "Dependencies";
			}
			print XML_SUBROUTINES "</sect2>\n";
#			print XML_SUBROUTINES "\t\t</para>\n\t</simplesect>\n";
		}
		
	}
	
}
print XML_SUBROUTINES "\t</sect1>\n";
print XML_LIBRARY "</programlisting>\n";
close(F90_LIBRARY);
close(XML_LIBRARY);
close(XML_SUBROUTINES);

# Example source code
open(XML_EXAMPLE,"> $xmlFiles{'example'}");
open(F90_EXAMPLE,$f90_example);
print XML_EXAMPLE $xmlHeader;
print XML_EXAMPLE "<programlisting>";
while(<F90_EXAMPLE>){print XML_EXAMPLE f90toxml($_);}
print XML_EXAMPLE "</programlisting>\n";
close(XML_EXAMPLE);
close(F90_EXAMPLE);

##################################
#my $toc = "\t\t<ul>\n";
#my $docdir = dirname($doc_main);

#	my ($comment, $indent, $newline, %f2xCssColors, @subNames);

#	$code =~ s/ /\&nbsp;/g;
#	$code =~ s/\t/\&nbsp;\&nbsp;\&nbsp;/g;
#	$newline = ($code =~ s/\n//);

#	$code = "$code<br \/>\n" if $newline;

#			$toc = "$toc\t\t\t<li><code role=\"$f2xCssColors{'personalised'}\"><a href=\"#$subroutineName\">$subroutineName</a></code></li>\n";

#$toc = "$toc\t\t</ul>\n";
