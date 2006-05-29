#!/usr/bin/perl -w
# Spanlib - Spectral Analysis Library
# Copyright (C) 2006  Stephane Raynaud
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
# 
# Contact: stephane dot raynaud at gmail dot com
#
# This script generates the documentation of pcamssa.f90

use strict;
use File::Basename;

# Inputs
my ($doc_top, $doc_mid, $doc_bot, $src_main, $src_example, $doc_main, $doc_code) = @ARGV;

# Declarations
my (@partNames, $partName, %parts, $line, $subroutineName, $arguments, $inside);
my ($vars, $type, $intent, $optional, $allocatable, $var, $name, @argDesc, $desc, $opt, %argRef, %arg);
my ($tmp, $tmp2, $i);
my %cols = (
	'name' => {
		'class' => 'b',
		'width' => '15'
	},
	'type' => {
		'class' => 's',
		'width' => '15'
	},
	'long_name' => {
		'class' => '',
		'width' => '70'
	}
);
my $toc = "\t\t<ul>\n";
my $docdir = dirname($doc_main);
#my $doc_toc = "$docdir/fortran_toc.html";
my $doc_sub = "$docdir/fortran_sub.html";
my $doc_exa = "$docdir/fortran_exa.html";
#print "$doc_sub\n";
#exit;

# Useful subroutines
sub gen_init_args {
	my $type = shift ;
	print SUBFILE "<h4>$type arguments</h4>\n\n<table width=\"100%\">\n";
}
sub f90tohtml {
	my $code = shift;
	my ($comment, $indent, $newline);
	$code =~ s/\&/\&amp;/g;
	$code =~ s/</\&lt;/g;
	$code =~ s/>/\&gt;/g;
	$code =~ s/ /\&nbsp;/g;
	$code =~ s/\t/\&nbsp;\&nbsp;\&nbsp;/g;
	$newline = ($code =~ s/\n//);
	$comment = "";
	$comment = "<span class=\"e i\">$2<\/span>" if 	$code =~ s/^([^!]*)(!.*)$/$1/;
	$code =~ s/['"]([^'"]*)['"]/"<span class="o">$1<\/span>"/g;
	$code =~ s/(\b|&nbsp;)(\d+)(\b|&nbsp;)/$1<span class="r">$2<\/span>$3/g;
	$code =~ s/(\b|&nbsp;)(pcamssa|pca|pcarec|mssa|mssarec|phasecomp|diasym)(\b|&nbsp;)/$1<span class="b">$2<\/span>$3/g;
	$code =~ s/(\b|&nbsp;)(procedure|parameter|optional|intent|allocatable|none)(\b|&nbsp;)/$1<span class="gr">$2<\/span>$3/g;
	$code =~ s/(\b|&nbsp;)(dot_product|trim|eoshift|modulo|ssyev|not|allocated|allocate|spread|present|sum|matmul|sqrt|transpose|deallocate|pack|unpack|present)(\b|&nbsp;)/$1<span class="cy">$2<\/span>$3/g;
	$code =~ s/(\b|&nbsp;)(print|module|contains|interface|logical|subroutine|if|then|else|endif|enddo|do|end|close|call|program|integer|real|open|write|character|use|implicit)(\b|&nbsp;)/$1<span class="vi">$2<\/span>$3/g;
	$code = "$code$comment";
	$code = "$code<br \/>\n" if $newline;
	return $code;
}

# Generate the subroutine help file by parsing the f90 sources
open(SUBFILE,"> $doc_sub");
open(F90FILE,$src_main);
while(<F90FILE>){
	# A subroutine starts here
	if(!/python/ && /^[ \t]*subroutine\s+([\w]+)\s*\(([^\)]+)\)/i){
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
		HEADERLOOP: while(<F90FILE>) {
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
		DECLLOOP: while(<F90FILE>) {
			if(/^[ \t]*! External/) {
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

		# Now, generate the xhtml block
		if ($partName ne "") {
			$toc = "$toc\t\t\t<li><code class=\"b\"><a href=\"#$subroutineName\">$subroutineName</a></code></li>\n";
			print SUBFILE "<h3>$parts{'Title'}: <code><a name=\"$subroutineName\" class=\"b\">$subroutineName</a></code></h3>\n<div class=\"sh3\">\n\n";
			$arguments =~ s/(^|[,=\s]+)(\w+)\b(?!=)/$1<span class=\"b\">$2<\/span>/g;
			$arguments =~ s/,/, /g;
			print SUBFILE "<h4>Usage:</h4>\n<code><p>call <span class=\"b\">$subroutineName</span>($arguments)</p></code>\n\n";
			for $partName (@partNames){
				if($partName eq "Dependencies"){
					$parts{$partName} =~ s/([\w]+)/<a href=\"#$1\">$1<\/a>/g;
					$parts{$partName} = "<span class='b'><code>$parts{$partName}</code></span>";
				}
				print SUBFILE "<h4>$partName</h4>\n<p>\n$parts{$partName}\n</p>\n\n" if $parts{$partName} ne "";
				if($partName eq "Description"){
					$tmp2 = "Necessary";
					gen_init_args($tmp2);
					my $n = @argDesc;
					for ($i=0;$i<$n;$i++) {
						if($argDesc[$i]{'optional'}==1 && $tmp2 eq "Necessary") {
							print SUBFILE "</table>\n\n";
							$tmp2 = "Optional";
							gen_init_args($tmp2);
						}
						print SUBFILE "\t<tr>\n";
						for $tmp ('name','type','long_name') {
							print SUBFILE "\t\t<td";
							print SUBFILE " class=\"$cols{$tmp}{'class'}\"" if $cols{$tmp}{'class'} ne "";
							print SUBFILE " width=\"$cols{$tmp}{'width'}%\">";
							print SUBFILE "<code>" if ($tmp eq 'name');
							print SUBFILE $argDesc[$i]{$tmp};
							print SUBFILE "</code>" if ($tmp eq 'name');
							print SUBFILE "</td>\n";
						}
						print SUBFILE "\t</tr>\n";
					}
					print SUBFILE "</table>\n\n";
				}
			}
			print SUBFILE "</div>\n<br />\n\n";
		}
		
	}
	
}
$toc = "$toc\t\t</ul>\n";
close(F90FILE);
close(SUBFILE);

# Now build the main file
open(MAINFILE,"> $doc_main");

# 1) Top
open(FILE,"$doc_top");
while(<FILE>){print MAINFILE $_;}
close(FILE);

# 2) TOC
print MAINFILE $toc;

# 3) mid
open(FILE,$doc_mid);
while(<FILE>){
	$_ =~ s/<DOC_CODE>/$doc_code/;
	print MAINFILE $_;
}
close(FILE);

# 4) subroutines
open(FILE,$doc_sub);
while(<FILE>){print MAINFILE $_;}
close(FILE);

# 5) Example with syntax colorization and indentation (an attempt)
print MAINFILE "<h2><a name=\"exam\"></a>4. Example <a href=\"#top\" class=\"top\">[Top]</a></h2>\n\n<code><p>\n\n";
open(FILE,$src_example);
while(<FILE>){print MAINFILE f90tohtml($_);}
print MAINFILE "</p></code>\n\n";

# 6) bottom
open(FILE,$doc_bot);
while(<FILE>){print MAINFILE $_;}
close(FILE);
print MAINFILE "\n\n<hr>\n<p class=\"i\">Document generated by Perl (<a href=\"gendoc.pl\"><code>gendoc.pl</code></a>)</p>\n\n</body>\n";
close(MAINFILE);

# Generate the html version of our f90 code
open(FILE,"> $doc_code");
print FILE "<head>\n\t<title>pcamssa.f90</title>\n\t<link rel=\"stylesheet\" type=\"text/css\" href=\"fortran.css\">\n</head>\n\n<body>\n<code>\n\n";
open(CODE,$src_main);
while(<CODE>){print FILE f90tohtml($_);}
close(CODE);
print FILE "</code>\n\n</body>";
close(FILE);
