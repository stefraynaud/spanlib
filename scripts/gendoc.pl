#!/usr/bin/env perl
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

####################################################################
# Initialisations
####################################################################

# Inputs
#my ($xmldir, $f90_library, $f90_example, $python_module) = @ARGV;
my ($f90_library, $f90_example, $python_module, $python_example1, $python_example2) = @ARGV;

# Basic declarations
my (@partNames, $partName, %parts, $line, $subroutineName, $arguments, $inside, $html);
my ($currentArg, %outputLongNames, $inputArguments, $outputArguments, $isInputArg, $class, $prefix, $indent, $redo, $subroutineTitle, @specs, $insideDescription, $subroutineDescription);
my ($vars, $type, $intent, $optional, $allocatable, $var, $name, @argDesc, $desc, $opt, %argRef, %arg);
my ($tmp, $tmp2, $i, @aTmp);

# Sections
my $upSect = 2;
my $loSect = $upSect + 1;

# Look
my %f2xTags = (
	'functions'		=> 'dot_product|trim|eoshift|modulo|ssyev|not|allocated|allocate|spread|present|sum|matmul|sqrt|transpose|deallocate|pack|unpack|present',
	'routines'		=> 'print|module|contains|interface|logical|subroutine|if|then|else|endif|enddo|do|end|close|call|program|integer|real|open|write|character|use|implicit',
	'attributes'	=> 'procedure|parameter|optional|intent|allocatable|none',
	'personalised'	=> ''
);
my %p2xTags = (
	'functions'		=> 'int|list|id|open|raw_input',
	'commands'		=> 'import|for|in|if|else|is|def|return|print|try|except|elsif',
	'special'		=> 'none|self|True',
	'personalised'	=> ''
);
my $p2xInsideComment = 0;

# Files
my %xmlFiles = (
	'f90_subroutines'		=> "doc_f90_sub_inc.xml",
	'python_subroutines'	=> "doc_pyt_sub_inc.xml",
	'f90_library'			=> "src_f90_lib_inc.xml",
	'f90_example'			=> "src_f90_exa_inc.xml",
	'python_module'		=> "src_pyt_mod_inc.xml",
	'python_example1'		=> "src_pyt_ex1_inc.xml",
	'python_example2'		=> "src_pyt_ex2_inc.xml"
);




####################################################################
# Useful subroutines
####################################################################

####################################################################
# xml header
sub gen_xml_header {
	my $xml_type = shift;
	return "<?xml version=\"1.0\"?>\n<!DOCTYPE $xml_type PUBLIC \"-//OASIS//DTD DocBook XML V4.2//EN\"\n\"http://www.oasis-open.org/docbook/xml/4.2/docbookx.dtd\">\n\n";
}


####################################################################
# Generic header for subroutine inputs (and outputs)
sub gen_init_args {
	my $title = shift ;
	my $oldSec = shift;
	my $tmphtml = '';
	$tmphtml .= "\t\t</itemizedlist>\n\t</simplesect>\n"
		if $oldSec ne "";
	$tmphtml .= "\t<simplesect>\n\t\t<title>$title</title>\n\t\t<itemizedlist>\n";
	return $tmphtml;
}

####################################################################
# Entry for arguments
sub gen_entry_args {
	my ($i) = @_ ;
	# Default values
	if($argDesc[$i]{'long_name'} =~ s/\[default: *([^\]]*)\]//) {
		my $default = $1;
		$argDesc[$i]{'spec'} =~ s/\]/default:$default]/
			if($argDesc[$i]{'spec'} !~ s/default:([^\],]*)/default:$default/);
	}
	$argDesc[$i]{'spec'} !~ s/default:([^\],]*)\b/default:<phrase role="HLpersonalised">$1<\/phrase>/;
	# Build
	my $tmphtml.= "\t\t\t<listitem>\n";
	$tmphtml.= "\t\t\t\t<varname role=\"ARGname\">$argDesc[$i]{'name'}    </varname>\n";
	$tmphtml.= "\t\t\t\t<option><phrase role=\"ARGspec\">$argDesc[$i]{'spec'}    </phrase></option> \n";
	$tmphtml.= "\t\t\t\t :: <emphasis>$argDesc[$i]{'long_name'}</emphasis>\n";
	$tmphtml.= "\t\t\t</listitem>\n";
	return $tmphtml;
}

####################################################################
# Syntax highlighting

# * Fortran 90
sub f90toxml {
	my $code = shift;
	my $comment;
	$code =~ s/\&/\&amp;/g;
	$code =~ s/</\&lt;/g;
	$code =~ s/>/\&gt;/g;
	$comment = "";
	$comment = "<phrase role=\"HLcomments\">$2</phrase>" if $code =~ s/^([^!]*)(!.*)$/$1/;
	$code =~ s/(['"][^'"]*['"])/<phrase role="HLstrings">$1<\/phrase>/ig;
	$code =~ s/\b(\d+)\b/<phrase role="HLdigits">$1<\/phrase>/ig;
	foreach $type ('personalised','attributes', 'functions', 'routines') {
		$code =~ s/\b($f2xTags{$type})\b/<phrase role="HL$type">$1<\/phrase>/ig;
	}
	$code =~ s/(\n)/$comment$1/;
	return $code;
}

# * Python
sub pytoxml {
	my @code;
	$code[0] = shift;
	my $comment = "";
	$code[1] = "";
	$code[0] =~ s/\&/\&amp;/g;
	$code[0] =~ s/</\&lt;/g;
	$code[0] =~ s/>/\&gt;/g;
	# One line long comment
	if($p2xInsideComment == 0 && $code[0] =~ s/^(.*)(""".*""")(.*)$/$1/) {
		$comment = $2;
		$code[1] = $3;
	# Begining of long comment
	} elsif($p2xInsideComment == 0 && $code[0] =~ s/^(.*)(\"\"\".*)/$1/) {
		$comment = $2;
		$p2xInsideComment = 1;
	# End of long comment
	} elsif($p2xInsideComment == 1 && $code[0] =~ s/^(.*\"\"\")(.*)$//) {
		$comment = $1;
		$code[1] = $2;
		$p2xInsideComment = 0;
	# Inside long comment
	} elsif($p2xInsideComment == 1) {
		$comment = $code[0];
		$code[0] = "";
		$code[1] = "";
	# Short comment
	} elsif($code[0] =~ s/^([^#]*)(#.*)$/$1/) {
		$comment = $2;
	}
	# Format the comment
	$comment = "<phrase role=\"HLcomments\">$comment</phrase>" if $comment ne "";
	# Format codes
	foreach my $i (0, 1) {
		$code[$i] =~ s/(['"][^'"]*['"])/<phrase role="HLstrings">$1<\/phrase>/ig;
		$code[$i] =~ s/\b(\d+)\b/<phrase role="HLdigits">$1<\/phrase>/ig;
		foreach $type ('personalised','special', 'functions', 'commands') {
			$code[$i] =~ s/\b($p2xTags{$type})\b/<phrase role="HL$type">$1<\/phrase>/ig;
		}
	}
	# Pack
	my $mycode = "$code[0]$comment$code[1]";
	$mycode =~ s/\n//;
	return "$mycode\n";
}




####################################################################
# Fortran library
####################################################################

####################################################################
# First loop to get subroutine names and convert to xml
open(F90_LIBRARY,$f90_library);
open(XML_F90_LIBRARY,"> $xmlFiles{'f90_library'}");
print XML_F90_LIBRARY gen_xml_header('programlisting')."<programlisting>";
while(<F90_LIBRARY>){
	# Subroutines
	if(/^[\s\t]*subroutine\s+([\w]+)\b/i) {
		$f2xTags{personalised} .= "|" if $f2xTags{personalised} ne "";
		$f2xTags{personalised} .= "$1";
	}
	# F90 to XML
	print XML_F90_LIBRARY f90toxml($_);
}
close(F90_LIBRARY);
print XML_F90_LIBRARY "</programlisting>\n";
close(XML_F90_LIBRARY);

####################################################################
# Generate the subroutine help file by parsing the f90 sources
open(XML_F90_SUBROUTINES,"> $xmlFiles{'f90_subroutines'}");
open(F90_LIBRARY,$f90_library);
$html = "\t<sect$upSect id=\"doc_f90_sub\">\n\t\t<title>F90 subroutines</title>\n";
while(<F90_LIBRARY>) {

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
			$argDesc[$i]{'intent'} = "";
			$argDesc[$i]{'long_name'} = "";
			$argDesc[$i]{'optional'} = 0;
			$i++;
		}

		# Parse header
		$inside=0;
		while(<F90_LIBRARY>) {
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
				last;
			}
		}

		# Parse declarations of external subroutine arguments
		$inside = 0;
		while(<F90_LIBRARY>) {
			if(/^[ \t]*! External/) {
				$inside = 1;
			} elsif($inside==1){
				if(/(real|integer),[ \t]*intent\((in|out|inout)\)[^:]*::(.+)$/i) {
					# Attributes
					$desc = "[intent:${2}put, type:$1]";
					$vars = "$3";
					$opt = /optional/;
					# Splitted lines
					while($vars =~ s/&[\t ]*$//) {
						$tmp = <F90_LIBRARY>;
						last if ! $tmp;
						$tmp =~ s/^[\t ]*\&//;
						$vars .= $tmp;
					}
					# Clean
					$vars =~ s/\n//g;
					$vars =~ s/ +//g;
					while($vars =~ s/\([^()]+\)//){next;}
					# Loop on variables
					for $var (split(',',$vars)) {
						$tmp = $argRef{$var};
						$argDesc[$tmp]{'spec'} = $desc;
						$argDesc[$tmp]{'optional'} = $opt;
						$arguments =~ s/([(, ]+)$var([), ]*)/$1$var=$var$2/ if $opt == 1

					}
				} elsif(/^[ \t]*$/){
					$inside=0;
					last;
				}
			}
		}

		# Now, generate the xml block
		if ($partName ne "") {
			# Title
			$html.= "<sect$loSect id=\"$subroutineName\">\n";
			$html.= "\t<title>$parts{'Title'}: <literal>$subroutineName</literal></title>\n";

			# Usage
			$arguments =~ s/(^|[,=\s]+)(\w+)\b(?!=)/$1<phrase role="HLpersonalised">$2<\/phrase>/g;
			$arguments =~ s/,/, /g;
			$html.= "\t<simplesect>\n";
			$html.= "\t\t<title>Usage</title>\n";
			$html.= "\t\t<programlisting>call <phrase role=\"HLpersonalised\">$subroutineName</phrase>($arguments)</programlisting>\n";
			$html.= "\t</simplesect>\n";

			# Other parts
			for $partName (@partNames){
				# Dependencies part
				if($partName eq "Dependencies"){
					if($parts{$partName} !~ /(LAPACK|BLAS)/){
						$parts{$partName} =~ s/([\w]+)/<link linkend=\"$1\">$1<\/link>/g;
						$parts{$partName} = "\t\t\t<phrase role=\"HLpersonalised\"><literal>$parts{$partName}</literal></phrase>";
					} else {
						$parts{$partName} = "<literal>$parts{$partName}</literal>";
					}
				}
				# Title of the part
				if($parts{$partName} ne "") {
					$html.= "\t<simplesect>\n";
					$html.= "\t\t<title>$partName</title>\n";
					$html.= "\t\t<para>\n$parts{$partName}\n";
					$html.= "\t\t</para>\n";
				}
				# Close description part
				if($partName eq "Description") {
					$html.= "\t</simplesect>\n";
					# Add input arguments
					$tmp = "";
					for my $i (0 .. $#{@argDesc}) {
						# Necessary or optional?
						if ($argDesc[$i]{'optional'} == 1) {
							$tmp2 = "Optional";
						} else {
							$tmp2 = "Necessary";
						}
						$html .= gen_init_args("$tmp2 arguments",$tmp) if $tmp ne $tmp2;
						$tmp = $tmp2;
						# One argument
						$html .= gen_entry_args($i);
					}
					$html.= "\t\t</itemizedlist>\n";
					$html.= "\t</simplesect>\n";
				}
#				$html.= "\t</simplesect>\n";
				$html.= "\t</simplesect>\n" if $partName eq "Dependencies";
			}
			$html.= "</sect$loSect>\n";
		}

	}

}
$html.= "\t</sect$upSect>\n";
print XML_F90_SUBROUTINES $html;
close(F90_LIBRARY);
close(XML_F90_SUBROUTINES);




####################################################################
# Python
####################################################################

####################################################################
# First loop to get subroutine names
open(PYTHON_MODULE,$python_module);
open(XML_PYTHON_MODULE,"> $xmlFiles{'python_module'}");
print XML_PYTHON_MODULE gen_xml_header('programlisting')."<programlisting>";
while(<PYTHON_MODULE>){
	# Subroutines
	if(/^[\s\t]*(class|def)\s+(\w+)\b/i && $2 ne "clean" && $2 ne "__init__") {
		$p2xTags{personalised} .= "|" if $p2xTags{personalised} ne "";
		$p2xTags{personalised} .= "$2";
	}
	# Python to xml
	print XML_PYTHON_MODULE pytoxml($_);
}
close(PYTHON_MODULE);
print XML_PYTHON_MODULE "</programlisting>\n";
close(XML_PYTHON_MODULE);

####################################################################
# Generate the subroutine help file by parsing the python sources
open(XML_PYTHON_SUBROUTINES,"> $xmlFiles{'python_subroutines'}");
$html = "\t<sect$upSect id=\"doc_pyt_sub\">\n\t\t<title>Python functions</title>\n";
$class = "";
$prefix = 'spanlib';
$redo = 0;
open(PYTHON_MODULE,$python_module);
READ_PYTHON_MODULE:
while(<PYTHON_MODULE>){

	$redo=0;

	# A blass starts here
	if(/^[\s\t]*class\s+([\w]+)\s*\(([^\)]+)\)/i){
		$class = "$1";

	# A function or class starts here
	} elsif(/^([\s\t]*)def\s+([\w]+)\s*\(([^\)]+)\)/i) {

		undef @argDesc;
		undef %argRef;
		undef %outputLongNames;

		next if $2 eq "clean";
		$subroutineName = "$2";
		$subroutineTitle = "";
		$subroutineDescription = "";
		$insideDescription = 0;
		$inputArguments = "$3";
		$outputArguments = "";
		$indent = $1;

		# Initialize argument descriptions whatever they are
		$inputArguments =~ s/[\t ]*//g;
		$inputArguments =~ s/\bself\b,?//;
		$i = 0;
		for $tmp (split(",",$inputArguments)) {
			# Optional argument (with default value)?
			$argDesc[$i]{'spec'} = '[intent:input';
			if($tmp =~ s/(.*)=(.*)/$1/) {
				$argDesc[$i]{'spec'} .= ", default:$2]";
			} else {
				$argDesc[$i]{'spec'} .= ']';
			}
			$argRef{$tmp} = $i;
			$argDesc[$i]{'name'} = "$tmp";
			$argDesc[$i]{'long_name'} = "";
			$i++;
		}

		# Parse header
		$inside=0;
		$currentArg="";
		while(<PYTHON_MODULE>) {
			if($inside == 0 && /^[\s\t]+\"\"\"(.*)\"\"\"[\s\t]*$/) {
				# Short header with a title only
				$subroutineTitle = $1;
				last;
			} elsif($inside == 0 && /^[\s\t]+\"\"\"(.*)/) {
				# Begining of a long header
				$subroutineTitle = $1;
				$inside = 1;
			} elsif($inside==1 && /\"\"\"[\s\t]*$/) {
				# End of a long header
				$inside = 0;
				$insideDescription = 0;
				last;
			} elsif($inside==1) {
				# Description of variable (long_name)
				if(/^[\t\s]+([\w\*]+)[\t\s]*::([^:]*)$/) {
					# Input or output?
					$isInputArg = 0;
					for my $key (keys(%argRef)) {
						if($key eq $1) {
							$isInputArg = 1;
							last;
						}
					}
					if($isInputArg == 1)  {
						# Input arguments
						$tmp = $argRef{"$1"};
						$argDesc[$tmp]{'long_name'} = "$2";
						$currentArg=$tmp;
					} else {
						# Outputs
						$outputLongNames{$1} = "$2";
						$currentArg="$1";
					}
				} elsif(/Description:::/){
					$insideDescription = 1;
				} elsif(/:::[\t\s]*/) {
					# An end of argument description
					$currentArg="";
					$insideDescription = 0;
				} elsif("$currentArg" ne "") {
					# Append to the argument description
					if($currentArg =~ /^\d+$/) {
						# Inputs
						$argDesc[$currentArg]{'long_name'} .= $_;
					} else {
						# Outputs
						$outputLongNames{$currentArg} .= $_;
					}
				} elsif($insideDescription == 1) {
					$subroutineDescription .= $_;
				}
			}
		}

		# Search for the return statement
		while(<PYTHON_MODULE>) {
			if((/^$indent\w/ && ! /^#/) || /^[\s\t]*(def|class)/) {
				# We are at the end of __init__
				if($subroutineName eq "__init__") {
					$argRef{"&lt;${class}_object&gt;"} = $i;
					$argDesc[$i]{'name'} = "&lt;${class}_object&gt;";
					$argDesc[$i]{'spec'} = "[intent:output]";
					$argDesc[$i]{'long_name'} = "Object created for further analysis";
					$i++;
				}
				$redo = 1 if /^[\s\t]*(def|class)/;
				last;
			} elsif(/[\s\t]+return (.*)/) {
				# Normal return
				$outputArguments = "$1";
				# Remove spaces
				$outputArguments =~ s/[\t ]*//g;
				# Loop on names
				@aTmp = ();
				for $tmp (split(",",$outputArguments)) {
					# Cleaning
					#$tmp =~ s/(\w+)\.\w+\(\'.+\'\)/$1/ or
					$tmp =~ s/\w+\.\w+\((.+)\)/$1/;
					push(@aTmp, $tmp);
					$argRef{$tmp} = $i;
					$argDesc[$i]{'name'} = "$tmp";
					$argDesc[$i]{'spec'} = "[intent:output]";
					$argDesc[$i]{'long_name'} = "";
					for my $key (keys(%outputLongNames)) {
						if($key eq $tmp) {
							$argDesc[$i]{'long_name'} = $outputLongNames{$tmp};
							last;
						}
					}
					$i++;
				}
				$outputArguments = join(",",@aTmp);
				last;
			}
		}

		# Generate the xml bloc

		# * Title
		$outputArguments =~ s/\b(\w+)\b/<phrase role="HLpersonalised">$1<\/phrase>/g;
		$outputArguments =~ s/,/, /g;
		if ($subroutineName eq "__init__") {
			$subroutineName = $class;
			$outputArguments = "<phrase role=\"HLpersonalised\">&lt;${class}_object&gt;</phrase>"
		} elsif($class ne "") {
			$prefix = "&lt;${class}_object&gt;";
		}
		$html .= "<sect$loSect id=\"${prefix}.$subroutineName\">\n";
		$html .= "\t<title>$subroutineTitle: <literal role=\"HLpersonalised\">${prefix}.$subroutineName</literal></title>\n";

		# * Usage
		$inputArguments =~ s/\b(\w+)([,=\s]+|$)/<phrase role="HLpersonalised">$1<\/phrase>$2/g;
		$inputArguments =~ s/,/, /g;
		$html .=  "\t<simplesect>\n";
		$html .= "\t\t<title>Usage</title>\n";
		$html .= "\t\t<programlisting>$outputArguments = <phrase role=\"HLpersonalised\">$prefix.$subroutineName</phrase>($inputArguments)</programlisting>\n";
		$html .= "\t</simplesect>\n";

		# * Description
		if($subroutineDescription ne "") {
			$html .=  "\t<simplesect>\n";
			$html .= "\t\t<title>Description</title>\n";
			$html .= "\t\t<para>$subroutineDescription</para>\n";
			$html .= "\t</simplesect>\n";
		}

		# Necessary and optional input arguments
		# and outputs
		$tmp = "";
		foreach my $i (0 .. $#{@argDesc}) {
			# Define the type
			if($argDesc[$i]{'spec'} =~ /input/) {
				if($argDesc[$i]{'spec'} =~ /default/) {
					$tmp2 = "Optional arguments"
				} else {
					$tmp2 = "Necessary arguments" ;
				}
			} else {
				$tmp2 = "Outputs" ;
			}
			$html .= gen_init_args($tmp2,$tmp) if $tmp ne $tmp2;
			$tmp = $tmp2 ;
			# One argument
			$html .= gen_entry_args($i);
		}
		$html .= "\t\t</itemizedlist>\n";
		$html .= "\t</simplesect>\n";
		$html .= "</sect$loSect>\n";

		# We finsihed on a def or class statement
		# so we need to parse it again.
		redo if $redo == 1;
	}
}
close(PYTHON_MODULE);
$html .= "\t</sect$upSect>\n";
print XML_PYTHON_SUBROUTINES $html;
close(XML_PYTHON_SUBROUTINES);



####################################################################
# Example source codes
####################################################################

####################################################################
# Fortran
open(XML_EXAMPLE,"> $xmlFiles{'f90_example'}");
open(F90_EXAMPLE,$f90_example);
print XML_EXAMPLE gen_xml_header('programlisting');
print XML_EXAMPLE "<programlisting>";
while(<F90_EXAMPLE>){print XML_EXAMPLE f90toxml($_);}
print XML_EXAMPLE "</programlisting>\n";
close(XML_EXAMPLE);
close(F90_EXAMPLE);

####################################################################
# Python
# 1)
open(XML_EXAMPLE,"> $xmlFiles{'python_example1'}");
open(PYTHON_EXAMPLE,$python_example1);
print XML_EXAMPLE gen_xml_header('programlisting');
print XML_EXAMPLE "<programlisting>";
while(<PYTHON_EXAMPLE>){print XML_EXAMPLE pytoxml($_);}
print XML_EXAMPLE "</programlisting>\n";
close(XML_EXAMPLE);
close(PYTHON_EXAMPLE);
# 2)
open(XML_EXAMPLE,"> $xmlFiles{'python_example2'}");
open(PYTHON_EXAMPLE,$python_example2);
print XML_EXAMPLE gen_xml_header('programlisting');
print XML_EXAMPLE "<programlisting>";
while(<PYTHON_EXAMPLE>){print XML_EXAMPLE pytoxml($_);}
print XML_EXAMPLE "</programlisting>\n";
close(XML_EXAMPLE);
close(PYTHON_EXAMPLE);

