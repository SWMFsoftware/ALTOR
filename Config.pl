#!/usr/bin/perl -i
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
use strict;

our $Component = "PC";
our $Code = "ALTOR";
our $MakefileDefOrig = 'src/Makefile.def';
our @Arguments = @ARGV;

my $config     = "share/Scripts/Config.pl";
if(-f $config){
    require $config;
}else{
    require "../../$config";
}

# Variables inherited from share/Scripts/Config.pl
our %Remaining; # Unprocessed arguments
our $ERROR;
our $WARNING;
our $Help;
our $Verbose;
our $Show;
our $ShowGridSize;
our $NewGridSize;
our $NewGhostCell;

&print_help if $Help;

# Equation and user module variables
my $Src         = 'src';


# Grid size variables
my $NameSizeFile1 = "srcBATL/BATL_size.f90";
my $NameSizeFile = "$Src/PIC_ModSize.f90";
my $GridSize;
my ($nI, $nJ, $nK, $nPType, $MaxBlock);

# Read previous grid size, equation and user module
&get_settings;

foreach (@Arguments){
    if(/^-s$/)                {$Show=1;                        next};

    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}

&set_grid_size if $NewGridSize and $NewGridSize ne $GridSize;
# Show grid size in a compact form if requested
print "Config.pl -g=$nI,$nJ,$nK,$nPType,$MaxBlock",
    ,"\n" if $ShowGridSize and not $Show;


my $Settings = &current_settings; print $Settings if $Show;


exit 0;

#############################################################################

sub get_settings{

    # Read size of the grid from $NameSizeFile
    open(FILE, $NameSizeFile1) or die "$ERROR could not open $NameSizeFile\n";
    while(<FILE>){
	next if /^\s*!/; # skip commented out lines
        $nI=$1           if /\bnI\s*=\s*(\d+)/i;
	$nJ=$1           if /\bnJ\s*=\s*(\d+)/i;
	$nK=$1           if /\bnK\s*=\s*(\d+)/i;
        $MaxBlock=$1     if /\bMaxBlock\s*=\s*(\d+)/i;
    }
    close FILE;
    open(FILE, $NameSizeFile) or die "$ERROR could not open $NameSizeFile\n";
    while(<FILE>){
	next if /^\s*!/; # skip commented out lines
	$nPType=$1        if /\bnPType\s*=\s*(\d+)/i;
    }
    close FILE;

    die "$ERROR could not read MaxBlock from $NameSizeFile\n" 
	unless length($MaxBlock);


    $GridSize = "$nI,$nJ,$nK,$nPType,$MaxBlock";

}

#############################################################################

sub set_grid_size{

    $GridSize = $NewGridSize;

    if($GridSize=~/^[1-9]\d*,[1-9]\d*,[1-9]\d*,[0-9]\d*,[1-9]\d*
       $/x){
	($nI,$nJ,$nK,$nPType,$MaxBlock)= split(',', $GridSize);
    }elsif($GridSize){
	die "$ERROR -g=$GridSize should be ".
	    #"4". 
	    " positive integers separated with commas\n";
    }

 
    print "Writing new grid size $GridSize into ".
	"$NameSizeFile ...\n";

    @ARGV = ($NameSizeFile);
    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(nPType\s*=[^0-9]*)(\d+)/$1$nPType/i;
	print;
    }
    @ARGV = ($NameSizeFile1);
    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(MaxBlock\s*=[^0-9]*)(\d+)/$1$MaxBlock/i;
	s/\b(nI\s*=[^0-9]*)(\d+)/$1$nI/i;
	s/\b(nJ\s*=[^0-9]*)(\d+)/$1$nJ/i;
	s/\b(nK\s*=[^0-9]*)(\d+)/$1$nK/i;
	print;
    }
}




sub current_settings{

    $Settings = 
	"Number of cells in a block        : nI=$nI, nJ=$nJ, nK=$nK\n";
    $Settings .= 
	"Number of different types of particles        : nPType=$nPType\n";
    $Settings .= 
	"Max. number of blocks/PE          : MaxBlock=$MaxBlock\n";

}

#############################################################################

sub print_help{

    print "
Additional options for ALTOR/Config.pl:

-g=NI,NJ,NK,MAXBLK   
                Set grid size. NI, NJ and NK are the number of cells 
                in the I, J and K directions, respectively. 
                MAXBLK is the maximum number of blocks per processor.


Examples for ALTOR/Config.pl:


Set block size to 8x8x8, number of blocks to 400",
":

    Config.pl -g=8,8,8,400",
"

Show settings for BATSRUS:

    Config.pl -s
\n";
    exit 0;


}

