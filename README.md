# PPT_finder10

Script for detection of polypyrimidine tract (PPT) in intronic sequence

The script takes as standard in (stdin) a table of intron sequences. The 50 nt closest to the
intron 3' splice site is enough. The table should be tab delimited and have two columns
1) name, 2) sequence

The script will detect the PPT conforming to the rules set out in:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2991248/
maximizing for length:
1) Both 3' and 5' ends must be pyrimidines;
2) No more than two contiguous purines are allowed;
3) Every purine segment (length L<3) must be surrounded by at least 4L pyrimidines (this forces
a minimum pyrimidine content greater than 2/3) distributed in a way that both upstream and
downstream pyrimidine segments are of length greater or equal to L;
4) T(GT)n stretches are allowed;
5) Minimum length of 9nt, uridine content greater or equal to 5.
6) If multiple qualifying PPT are found only report the one closest to 3' splice site (3'ss)
