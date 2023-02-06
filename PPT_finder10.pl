#!/usr/bin/perl

# Script for detection of polypyrimidine tract (PPT) in intronic sequence.
# The script takes as standard in (stdin) a table of intron sequences. The 50 nt closest to the
# intron 3' splice site is enough. The table should be tab delimited and have two columns
# 1) name, 2) sequence
# The script will detect the PPT conforming to the rules set out in:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2991248/
## maximizing for length:
## 1) Both 3' and 5' ends must be pyrimidines;
## 2) No more than two contiguous purines are allowed;
## 3) Every purine segment (length L<3) must be surrounded by at least 4L pyrimidines (this forces
## a minimum pyrimidine content greater than 2/3) distributed in a way that both upstream and
## downstream pyrimidine segments are of length greater or equal to L;
## 4) T(GT)n stretches are allowed;
## 5) Minimum length of 9nt, uridine content greater or equal to 5.
## 6) If multiple qualifying PPT are found only report the one closest to 3' splice site (3'ss)

# Made Feb 2021 by Morten T Venoe.


use strict;
use warnings FATAL => qw ( all );

my $name = "A";
my $seq = "A";
my @first_matches;
my $match;
my $final_match;
my $match1;
my $match2;
my $loc = 1000;
my $ppt;
my $ppt_loc = "1000";
my $ppt_loc_hit;
my $A_ppt;
my $T_ppt;
my $C_ppt;
my $G_ppt;
my $ppt_score;


# Regex

while (<>){
	chomp;
	($name, $seq) = split /\t/, $_;
	print "$name\t";
	# Regex finding all matches that start and end with T or C, allowing A or G if max two adjesant
	while ( $seq =~ /([TC]+([AG]{1,2}[TC]+)*)/gi ) {
	
	# Check if there is one or more matches
	if ( !defined($1) ) {
	} else {
		$match = $1;
		# Count number of T in match
		my $T_count = () = $match =~ /T/gi;
		# Require length min 9 or number of T more than 4 in match
		if ( length($match) > 8 || $T_count > 4) {

			# T(GT)n sequences are allowed even though this would violate other rules. To allow these sequenced the G
			# in T(GT)n G is edited to y so it is preserved in later filters. Three T(GT)n are allowed
			if ( $match =~ /(.*?T)((?:GT)+)$/i ) {
                                my $gt = $2;
                                my $one = $1;
                                $gt =~ s/G/y/ig;
                                $match = "$one$gt";
                        }
			if ( $match =~ /(.*?T)((?:GT)+)(.+?)$/i ) {
				my $gt = $2;
				my $one = $1;
				my $three = $3;
				$gt =~ s/G/y/ig;
				$match = "$one$gt$three";
			}
                        if ( $match =~ /(.*?T)((?:GT)+)(.+?)$/i ) {
                                my $gt = $2;
                                my $one = $1;
                                my $three = $3;
                                $gt =~ s/G/y/ig;
                                $match = "$one$gt$three";
                        }
                        if ( $match =~ /(.*?T)((?:GT)+)(.+?)$/i ) {
                                my $gt = $2;
                                my $one = $1;
                                my $three = $3;
                                $gt =~ s/G/y/ig;
                                $match = "$one$gt$three";
                        }

			# Remove matches where single purines are flanked by < 4 pyrimidines
			# The idea here is to find A or G that violate this role and substitute them with x for later removal
			# These matches are not simply removed now because other subsetions of these matches might be real PPT
			if ( $match =~ /(.+)([yAG][TC]{1})([AG]{1})([TC]{1,2}[yAG])(.+)/i ) {
				$match = "$1$2x$4$5";
			}
			if ( $match =~ /(.+)([yAGx][TC]{1,2})([AG]{1})([TC]{1}[yAGx])(.+)/i ) {
                                $match = "$1$2x$4$5";
                        }
			if ( $match =~ /^([TC]{1})([AG]{1})([TC]{1,2}[yAGx])(.+)/i ) {
                                $match = "$1x$3$4";
                        }
			if ( $match =~ /^([TC]{1,2})([AG]{1})([TC]{1}[yAGx])(.+)/i ) {
                                $match = "$1x$3$4";
                        }
			if ( $match =~ /(.+)([yAGx][TC]{1})([AG]{1})([TC]{1,2})$/i ) {
                                $match = "$1$2x$4";
                        }
			if ( $match =~ /(.+)([yAGx][TC]{1,2})([AG]{1})([TC]{1})$/i ) {
                                $match = "$1$2x$4";
                        }
			# Remove matches where double purines are flanked by < 8 pyrimidines with at least 2 pyrimidines at either end
			# As above A or G that violate this role are substituted with x for later removal
			if ( $match =~  /(.+)([yAGx][TC]{1})([AGx]{2})([TC]{1,50}[yAGx])(.+)/i ) {
                                $match = "$1$2x$4$5";
                        }
			if ( $match =~  /(.+)([yAGx][TC]{2})([AGx]{2})([TC]{1,5}[yAGx])(.+)/i ) {
                                $match = "$1$2x$4$5";
                        }
			if ( $match =~  /(.+)([yAGx][TC]{3})([AGx]{2})([TC]{1,4}[yAGx])(.+)/i ) {
                                $match = "$1$2x$4$5";
                        }
			if ( $match =~  /(.+)([yAGx][TC]{4})([AGx]{2})([TC]{1,3}[yAGx])(.+)/i ) {
                                $match = "$1$2x$4$5";
                        }
			if ( $match =~  /(.+)([yAGx][TC]{5})([AGx]{2})([TC]{1,2}[yAGx])(.+)/i ) {
                                $match = "$1$2x$4$5";
                        }
			if ( $match =~  /(.+)([yAGx][TC]{6,50})([AGx]{2})([TC]{1}[yAGx])(.+)/i ) {
                                $match = "$1$2x$4$5";
                        }


                        if ( $match =~  /^([yAGx][TC]{1})([AGx]{2})([TC]{1,50}[yAGx])(.+)/i ) {
                                $match = "$1x$3$4";
                        }
                        if ( $match =~  /^([yAGx][TC]{2})([AGx]{2})([TC]{1,5}[yAGx])(.+)/i ) {
                                $match = "$1x$3$4";
                        }
                        if ( $match =~  /^([yAGx][TC]{3})([AGx]{2})([TC]{1,4}[yAGx])(.+)/i ) {
                                $match = "$1x$3$4";
                        }
                        if ( $match =~  /^([yAGx][TC]{4})([AGx]{2})([TC]{1,3}[yAGx])(.+)/i ) {
                                $match = "$1x$3$4";
                        }
                        if ( $match =~  /^([yAGx][TC]{5})([AGx]{2})([TC]{1,2}[yAGx])(.+)/i ) {
                                $match = "$1x$3$4";
                        }
                        if ( $match =~  /^([yAGx][TC]{6,50})([AGx]{2})([TC]{1}[yAGx])(.+)/i ) {
                                $match = "$1x$3$4";
                        }


                        if ( $match =~  /(.+)([yAGx][TC]{1})([AGx]{2})([TC]{1,50}[yAGx])$/i ) {
                                $match = "$1$2x$4";
                        }
                        if ( $match =~  /(.+)([yAGx][TC]{2})([AGx]{2})([TC]{1,5}[yAGx])$/i ) {
                                $match = "$1$2x$4";
                        }
                        if ( $match =~  /(.+)([yAGx][TC]{3})([AGx]{2})([TC]{1,4}[yAGx])$/i ) {
                                $match = "$1$2x$4";
                        }
                        if ( $match =~  /(.+)([yAGx][TC]{4})([AGx]{2})([TC]{1,3}[yAGx])$/i ) {
                                $match = "$1$2x$4";
                        }
                        if ( $match =~  /(.+)([yAGx][TC]{5})([AGx]{2})([TC]{1,2}[yAGx])$/i ) {
                                $match = "$1$2x$4";
                        }
                        if ( $match =~  /(.+)([yAGx][TC]{6,50})([AGx]{2})([TC]{1}[yAGx])$/i ) {
                                $match = "$1$2x$4";
                        }


			# Do the initial Regex again to find  all matches after violating A and G have been substituted to x
			while ( $match =~ /([TC]+([AGy]{1,2}[TC]+)*)/gi ) {
			$final_match = $1;

               		# Count number of T in match
                	my $T_count = () = $final_match =~ /T/gi;
			if ( length($final_match) > 8 || $T_count > 4) {

				$loc = length($seq) - index($seq, $final_match) - length($final_match) - 1;

				#Check each candidate PPT for this gene. Only keep the qualifying PPT which is closest to the intron 3'SS
				# Also, substitute y for G to complete T(GT)n sequences
				if ( $loc < $ppt_loc ) {
					$ppt = $final_match;
					$ppt =~ s/y/G/ig;
					$ppt_loc_hit = $loc;

				}
				$ppt_loc = $loc;
			}
			}
		}
	}


	}
	if ( defined($ppt) ) {
		# Count number of A,T,C,G in PPT
		$A_ppt = () = $ppt =~ /A/gi;
		$T_ppt = () = $ppt =~ /T/gi;
		$C_ppt = () = $ppt =~ /C/gi;
		$G_ppt = () = $ppt =~ /G/gi;

		# Calculate PPT score
		$ppt_score = ($A_ppt*(-2) + $C_ppt*2 + $G_ppt*(-2) + $T_ppt*3);
	
		# Print the PPT score, the PPT distance to intron 3'SS and PPT sequence
		print "$ppt_score\t$ppt_loc_hit\t$ppt\n";

		# Undef the ppt variable to make sure this ppt is not erronously assigned to next gene if next gene has no qualifying PPT
		undef($ppt);
	} else {
		print "0\t0\tNA\n";
	}
}

