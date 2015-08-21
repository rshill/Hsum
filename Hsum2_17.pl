#!/usr/bin/perl
use Parallel::ForkManager;
use Getopt::Long;

#This program is designed to identify blocks of homozygosity for genome wide affymetrics and illumina  SNP genotyping.

my $version = 17;
my %columns = ("markers", 29, "cM", 32, "bp", 34);   #column number for number of markers, cM size, and bp size in final output.

my $help;
my $percenterrorrate = 5;
my $allowedhet = 1000;
my $maxcon = 2;
my $xatend = 10;
my $yatend = 3;
my $units1 = "cM";
my $units2 = "bp";
my $units3 = "markers";
my $cmcutoff = 0.5;
my $bpcutoff = 1;
my $markercutoff = 25;
my $dir = ".";
my $threads = 1;

GetOptions ( 
	'help' => \$help,
	't:i' => \$threads,
	'p:f'  => \$percenterrorrate,
	'h:i'  => \$allowedhet,
	'c:i'  => \$maxcon,
	'x:i'  => \$xatend,
	'y:i'  => \$yatend,
	'1:s'  => \$units1,
	'2:s'  => \$units2,
	'3:s'  => \$units3,
	'cmcut:f' => \$cmcutoff,
	'bpcut:i' => \$bpcutoff,
	'markcut:i' => \$markercutoff,
	'd:s'  => \$dir,
);

if ($help == 1) {
	print "This script will read a genotype data file and output a list of runs of homozygosity.\n";
	print "It will read through all of the text files in the indicated directory, and generate one output file per input.\n";
	print "The input file must be a tab delimited text file with columns: marker name, dbSNP rs number, chromosome, physical position, genetic position, genotype score, and genotype call.\n";
	print "Genotype call must be AA, AB, BB, or NO.\n";
	print "Command line options:\n";
	print "--help, print this informaiton and quits.\n";
	print "-t, number of threads to use for processing.\n";
	print "-p, % of non-homozygous SNPs allowed in homozygous block, default 5.\n";
	print "-h, max number of heterozygous SNPs allowed in a homozygous block, default 1000 (not used)\n";
	print "-c, max number of consecutive heterozygous SNPs allowed, default 2.\n"; 
	print "-x, number of SNPs at the end of a block to monitor for excess heterozygous snps, default 10\n";
	print "-y, number of heterozygous SNPs allowed in the monitored SNPs at the end of a block, default 3\n";
	print "-1, first criteria to use for sorting output, can be 'markers', 'bp', or 'cM', default cm.\n";	
	print "-2, second criteria to use for sorting output, must be one of the options above, and cannot be the same given for -1, default bp.\n";	
	print "--cmcut, cutoff used for smallest cM block reported, default 0.5.\n";
	print "--bpcut, cutoff used for smallest bp block reported, default 1.\n";
	print "--markcut, cutoff used for smallest number of markers in a block, default 25\n";
	print "-d, locaton of the genotype files, default \"\\.\"\n";
	print "\n\n";
	die;
}

#check if sorting units are correct and define third sort
if ($units1 !~ /(markers)|(cM)|(bp)/) {
	die "error: -1 must be markers, cM, or bp\n";
} elsif ($units2 !~ /(markers)|(cM)|(bp)/) {
	die "error: -2 must be markers, cM, or bp\n";
} elsif ($units1 eq $units2) {
	die "error: -1 and -2 must be different\n";
} else {
	my @ops = qw(markers cM bp);
	my @ops2 = grep { $units1 ne $_ } @ops;
	my @ops3 = grep { $units2 ne $_ } @ops2;
	$units3 = $ops3[0];
}

my $mark1 = $columns{$units1};
my $mark2 = $columns{$units2};
my $mark3 = $columns{$units3};
my %unitsort = ("markers", $markercutoff, "cM", $cmcutoff, "bp", $bpcutoff );
my $cutoff = $unitsort{$units1};
my $cutoff2 = $unitsort{$units2};
my $cutoff3 = $unitsort{$units3};
print "$cutoff, $cutoff2, $cutoff3\n";

$time = time;
my $method = "Hsum2_CO-${cutoff}${units1}_ME-${percenterrorrate}_MH-${allowedhet}_MC-${maxcon}_BE-${yatend}-${xatend}";
#get list of files
opendir DIR, "$dir" || die "cannot open $dir, $!\n";
my @textfiles = grep /.+GT\.txt$/, readdir DIR;
closedir DIR;

my $pm = new Parallel::ForkManager($threads);

#cycle through and read files one at a time
foreach my $file (@textfiles) {
	$pm -> start and next;
	open (DATAFILE, "$file") || die "cannot open $file, $!";
	print "working on file $file\n";
	$prevchrom = 1;
	$mem = $error = $operror = $perror = $fperror = $concount = $activeblock = $snpcount = $hetsnpcount = $hetsnptally = $blockendhet = 0;
	$foundfirst = $firsthomo = $firsthet = $activehet = $lasthomo = $lasthet = $ahetsnpcount = $number = 0;
	$icmdist = $ocmdist = $imbdist = $ombdist = 0;
	$snplist = $hetsnplist = "";
	$blockcount = -1;
	@buffer = @blocks = @line = @memory = @memory2 = @memory3 = @tmp1 = @prevsnp = @blockend = @sort = %order = ();
	while (($#buffer >= 0) ? ($_ = join "\t", @{$buffer[0]}) : ($_ = <DATAFILE>)) {   #if there is info in @buffer then look at them first, otherwise goto next line
		($#buffer >= 0) && (shift @buffer);
		next unless /^(SNP)|(rs)/;                                                   #ignore non SNP lines
		chop;
		@line = ();                                        
		@line = split /\t/;
		&chromcheck;                                                          #check to see if new snp is on same chromosome as previous snp  
		@prevsnp = ();
		@prevsnp = @line;                                                     #sets @prevsnp equal to line incase it is the last snp on chrom (used in chromcheck)
		my @gthettest = ();
		@gthettest = split(//, $line[6]);                                     #split genotype call to test for homozygosity
		if (($gthettest[0] eq $gthettest[1])||($line[6] =~ /^[Nn][Oo]/)) {    #test to see if two genotypes are the same, so it will id real gt calls as well as AB calls.  Old: do this if the snp is not heterozygous, ie AA, BB, NoCall ... or anything else for that matter
			if ($line[6] =~ /^[Nn][Oo]/) {                                 #if NoCall or NO then do not start a new block, but continue the block if open
				($activeblock == 1) && (&continueblock);
			} else {
				($activeblock == 1) ? (&continueblock) : (&startblock);       #if there is a block open then add to that block, else start a new block
			}
		} else {                                                              #do this if nep is heterozygous
			@memory2 = @line;                                             #@memory2 always equals last heterozygous snp
			($activeblock == 1) && (&checkerror);                         #if there is a block open, check the heterozygous status
		}
	}
	if ($activeblock == 1) {                                                      #if there is an open block after the last snp is read, the close out that block
		($mem == 0) && (@memory3 = @prevsnp) && ($mem = 1);                   #before closing the block, set memory3 to previous snp and mem to 1 unless mem is already 1
		&closeblock;
	}
	
	$a= "";                                                                       #creates unique name for output file based on the sample name
	$b = 0;
	$sample= $file;
	$sample =~ s/(.*)_GT.txt/\1/;
	$newfile = "${sample}_Hsum2${a}.txt";
	$dbfile = "${sample}_dbfile${a}.txt";
	while (-T $newfile) {
		$b++;
		$a = "-$b";
		$newfile = "${sample}_Hsum2${a}.txt";
		$dbfile = "${sample}_dbfile${a}.txt";
		$version = $a;
	}

	open (OUT, ">$newfile") || die "cannot open $newfile";
	open (DB, ">$dbfile") || die "cannot open $dbfile";
	print OUT "Blocks of homozygosity larger than $cutoff $units1 with an allowable error rate of ${percenterrorrate}%, a maximum of $allowedhet heterozyous markers per block, a maximum of $maxcon consecutive heterozygous markers, and no more than $yatend heterozygous markers in the last $xatend markers\n\n";
	print OUT "Block#\t#Markers\tcM flanking\tMb flanking\tcM homozygous\tMb homozygous\tChromosome\t%Heterozygous\t#Heterozygous\tFirst flanking affy\tfirst flanking rs\tfirst flanking cM\tfirst flanking Mb\tfirst affy\tfirst rs\tfirst cM\tfirst Mb\tlast affy\tlast rs\tlast cM\tlast Mb\tlast flanking affy\tlast flanking rs\tlast flanking cM\tlast flanking Mb\tList of genotypes\n";
	foreach ($c = 0; $c < @blocks; $c++) {
		$order{$c} = $blocks[$c][$mark1];
		$order2{$c} = $blocks[$c][$mark2];
		$order3{$c} = $blocks[$c][$mark3];
	}
	@sort = reverse sort {($order{$a} <=> $order{$b}) || ($order2{$a} <=> $order2{$b}) || ($order3{$a} <=> $order3{$b})} keys(%order);  
	foreach ($d = 0; $d < @sort; $d++) {
		$number = $d +1;
		($blocks[$sort[$d]][$mark1] < $cutoff) && last;
		($blocks[$sort[$d]][$mark2] < $cutoff2) && next;    #cutoff list for second criteria
		($blocks[$sort[$d]][$mark3] < $cutoff3) && next;    #cutoff list for third criteria
		$outline = join "\t",$number,@{$blocks[$sort[$d]]}[29,32,34,31,33,2,30,35,0,1,4,3,7,8,11,10,14,15,18,17,21,22,25,24,28];
		print OUT "$outline\n";
		print DB "$sample\t$method\t$time\t$version\t$outline\n";
	}
	close OUT;
	close DB;
	$pm -> finish; #exit from child proces
}

$pm -> wait_all_children;  #wait here untill all of the childreh are processed
print "all children have returned home\n";

sub startblock {
	$blockcount++;                                                                #increase block count by 1, it is preset to -1 so first block is 0 
	$activeblock = 1;                                                             #set avctiveblock toggle to 1, open block
	@memory = @line;                                                              #$memory will be last homozygous marker in block, it gets reset with each homogygous snp
	if (@memory2 == ()) {
		@{$blocks[$blockcount]} = (@line, @line);                             #if first snp in file then memory2 will be empty
	}else{
		@{$blocks[$blockcount]} = (@memory2, @line);                          #start new block in @blocks, will get all info from the last heterozygous snp (@memory2) -- or current snp if first on chrom -- and current snp (@line)
	}
	$snplist = "";                                                                #reset snplis
	$snplist = "$line[0],$line[6]";                                               #enter snp name and genotype into $snplist
	$snpcount = 1;                                                                #start counting the number of snps in the block at 1
	@blockend = ();
	push @blockend, [@line];                                                        #put line into @blockend to monitor the last $xatend snps for heterozygous snps
	$blockendhet = 0;                                                             #set heterozyous count for the last $xatend snps
}


sub continueblock {
	if ($mem == 1) {                                                              #test if there is an active heterzygous list that needs to be added to the current block before adding the current line
		$snplist .= $hetsnplist;                                              #adds the list of heterozygus snps to the overall list of snps
		$snpcount += $concount;                                               #adds the heterozygous snps to the overall count of snps
		$hetsnpcount = $hetsnptally;                                          #adds the heterozygous snps to the total count of heterozygus snps
		$operror = $perror;                                                   #set old percent error ($operror) equal to percent error ($perror)
		$hetsnplist = "";                                                     #reset $hetsnplist because it will concnatinate
		$concount = 0;                                                        #reset $concount
		$mem = 0;                                                             #reset mem toggel to 0, no longer in het loop
	}
	@memory = @line;                                                              #@memory must be last homozygous snp
	$snplist .= ",$line[0],$line[6]";                                             #add new snp name and genotype to list
	$snpcount++;                                                                  #increment snp count
	if ($#blockend == ($xatend -1)) {                                              #if @blockend is full, then remove the first line and add the new line
		@tmp2 = @{$blockend[0]};
		shift @blockend;
		push @blockend, [@line];
		my @gthettest = ();
		@gthettest = split(//, $tmp2[6]);
		(($gthettest[0] eq $gthettest[1])||($tmp2[6] =~ /^[Nn][Oo]/)) || ($blockendhet--);
		@tmp2 = ();
	} else {                                                                      #if @blockend is not full, then add new line
		push @blockend, [@line];
	}
}


sub checkerror {
	($mem == 0) && (@memory3 = @line) && ($operror = $perror);                    #if first het snp then set @memory3 equal to first heterozygous snp and preserves percent error
	$error++;                                                                     #increment total number of heterozygous snps
	$concount++;                                                                  #increment number of consecutive heterozygous snp
	$hetsnptally++;                                                               #total number of heterozygous snps in block
	$perror = ($error / ($snpcount + $concount)) * 100;                           #calculates percent error
	if ($#blockend >= ($xatend -1)) {                                              #if @blockend is full, then remove the first line and add the new line
		@tmp2 = @{$blockend[0]};
		shift @blockend;
		push @blockend, [@line];
		my @gthettest = ();
		@gthettest = split(//, $tmp2[6]);
		(($gthettest[0] eq $gthettest[1])||($tmp2[6] =~ /^[Nn][Oo]/)) && ($blockendhet++);                              #if the line removed is homozygous then add 1 to the $blockendhet because new line is heterozygous
		@tmp2 = ();
	} else {                                                                      #if @blockend is not full, then add new line
		push @blockend, [@line];
		$blockendhet++;                                                       #add 1 to $blockendhet because new line is heterozygous
	}	
	if ($blockendhet >= $yatend) {
		$foundfirst = $activatehet = $ahetsnpcount = 0;
		$firsthomo = $firsthet = $#blockend;
		for my $d (0..$#blockend) {                                              #look through lines in array of the last $xatend SNPs and find the fist heterozygous snp and last heterozgous and homozygous snps
			my @gthettest = ();
			@gthettest = split(//, $blockend[$d][6]);
			if (($gthettest[0] eq $gthettest[1])||($blockend[$d][6] =~ /^[Nn][Oo]/)) {
				$activatehet = 1;                                     #found first homozygous snp, now start looking for next heterozygous snp
				$lasthomo = $d;
			}else{
				$lasthet = $d;
				($foundfirst == 1) && ($ahetsnpcount++);
			}
			if (($foundfirst == 0) && (($gthettest[0] ne $gthettest[1]) && ($blockend[$d][6] !~ /^[Nn][Oo]/)) && ($activatehet == 1)){       #v16:fixed bug here, was set to look for homozygous snp, not het snp  #don't start looking for first heterozygus marker unless a homozygous markers has already been found    
				$firsthet = $d;
				$firsthomo = ($d - 1);
				$foundfirst = 1;
				$ahetsnpcount++;0
			}
		}
		$snplist .= $hetsnplist;                                              #add the last heterozygous snp data on before adjusting it
		$snpcount += $concount;      
		$hetsnpcount = $hetsnptally; 
		$hetsnplist = "";            
		$concount = 0;               
		$mem = 0;                    
		@memory = @{$blockend[$firsthomo]};                                      #set @memory equal to the first homozygous snp in endblock
		@line = @{$blockend[$firsthet]};                                         #set @line equal to the first heterozygous snp in endblock
		$snplist =~ s/(.*$blockend[$firsthomo][0],$blockend[$firsthomo][6]).*/\1/;   #adjust $snplist to eliminate all snps after the last homozygous snp
		$snpcount -= ($#blockend - $firsthomo);                               #adjust snpcount to remove extra snps
		$hetsnpcount -= $ahetsnpcount;                                        #adjust hetnspcount to remove extra heterozygous snps
		@buffer = ();
		if ($firsthet < $#blockend) {                                         #if there are lines left in the blockend, then they need to be put into an array (@buffer) and sent back through incase there is a block start in there
			for $e ($firsthet..$#blockend) {
				push @buffer, [@{$blockend[$e]}, "\n"];
			}
		}
		&closeblock;
	} elsif ($perror > $percenterrorrate) {                                       #moved this and maxcon and allowedhet to after yatend.     #if the new percent error exceedes user limit then close block
		&closeblock;
	} elsif ($concount > $maxcon) {                                               #if the new consecuive heterozygus snp count exceeds user limit then close block
		&closeblock;
	} elsif ($hetsnptally > $allowedhet) {                                        #if the new total number of snps exceeds user limit then close block
		&closeblock;
	} else {                                                                      #if the new heterozygus snp doesnot execced limit then  
		$mem = 1;                                                             #set $mem to 1 (open het  block),
		$hetsnplist .= ",$line[0],$line[6]";                                  #add snp name and genotype to hetsnplist
	}
}


sub closeblock {
	($mem == 1) && (@line = @memory3);                                             #if comming from checkerror then need to set @line equal to first heterozygous snp
	$fperror = ($hetsnpcount / $snpcount) * 100;
	push @{$blocks[$blockcount]}, @memory, @line, $snplist, $snpcount, $fperror;   #push onto block last homozygous snp (@memory), last heterozygous snp (@line), snp list, snp count, final percent error
	$icmdist = $blocks[$blockcount][18] - $blocks[$blockcount][11];                #calculates size of inner block in cM
	$ocmdist = $blocks[$blockcount][25] - $blocks[$blockcount][4];                 #calculates size of outer block in cM
	$imbdist = $blocks[$blockcount][17] - $blocks[$blockcount][10];                #caluclates size of inner block in Mb
	$ombdist = $blocks[$blockcount][24] - $blocks[$blockcount][3];                 #calculates size of outer block in Mb
	push @{$blocks[$blockcount]}, $icmdist, $ocmdist, $imbdist, $ombdist, $hetsnpcount;          #pushes inner and outter block distance in cM and Mb onto blocks
	$activeblock = $snpcount = $mem = $error = $operror = $perror = 0;             #clear variables that cycle with each block
	$concount = $icmdist = $ocmdist = $imbdist = $ombdist = 0;
	$blockendhet = $hetsnpcount = $hetsnptally = 0;
	$snplist = $hetsnplist = "";                                                   #clear snplist
	@memory = @memory3 = ();                                                       #clear memory arrays, but not @memory2 because next snp could be start of next block
}


sub chromcheck {
	if ($prevchrom == $line[2]) {                                                   #if new snp is on same chromosome as previous snp, do nothing
		return;
	} elsif ($activeblock == 1) {                                                   #if on a new chromosome and there is an open block then
		@memory2 = @line;                                                       #set memory2 to current line, this way if this snp is homozygous and start of new block, the last heterozygous snp before block will actually be first snp on chromosome
		@tmp1 = @line;                                                          #&closeblock will change the value of line to prevsnp, need to reset line to current line before returning
		($mem == 0) && (@memory3 = @prevsnp) && ($mem = 1);                     #if heterozygous memory is not in use, then last snp on previous chromosome was homozygous, thus need to set @memory3 equal to that snp and $mem equal to 1 so it will close properly (else it will try to use current line on wrong chromosome 
		&closeblock;                                                            #close the open block then return to the top
		@line = @tmp1;
		$prevchrom = $line[2];
		@tmp1 = ();
		return;                                                                 #after closing block for block on previous chromosome, the proceed to top and analyze current line
	} else {                                                                        #if new chromosome but no open block:
		$prevchrom = $line[2];                                                  #reset $prevchrom
		@memory2 = @line;                                                       #set memory2 to current line, this way if this snp is homozygous and start of new block, the last heterozygous snp before block will actually be first snp on chromosome
		return;
	}
}

