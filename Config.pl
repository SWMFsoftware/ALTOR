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
my $NameBatlFile = "srcBATL/BATL_size.f90";
my $NameSizeFile = "$Src/PC_ModSize.f90";
my $GhostCell;
my $GridSize;
my $Hybrid;
my $NewHybrid;
my $ReadHybrid;
my ($nI, $nJ, $nK, $MaxBlock, $nPType, $nElectronMax);


# Read previous grid size, equation and user module
&get_settings;

foreach (@Arguments){
    if(/^-s$/)                {$Show=1;                        next};
    if(/^-ng=(.*)$/i)         {$NewGhostCell=$1; next};
    if(/^-hybrid=(.*)$/i)         {$ReadHybrid=$1; next};
    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}
$NewHybrid=0;
$NewHybrid=1 if($ReadHybrid eq "Y");

&set_grid_size if ($NewGridSize and $NewGridSize ne $GridSize)
    or            ($NewGhostCell and $NewGhostCell ne $GhostCell)
    or            ($ReadHybrid and $NewHybrid ne $Hybrid);
# Show grid size in a compact form if requested
print "Config.pl -g=$nI,$nJ,$nK,$MaxBlock,$nPType,$nElectronMax",
    ,"\n" if $ShowGridSize and not $Show;


my $Settings = &current_settings; print $Settings if $Show;


exit 0;

#############################################################################

sub get_settings{

    # Read size of the grid from $NameSizeFile
    open(FILE, $NameBatlFile) or die "$ERROR could not open $NameBatlFile\n";
    while(<FILE>){
	next if /^\s*!/; # skip commented out lines
        $nI=$1           if /\bnI\s*=\s*(\d+)/i;
	$nJ=$1           if /\bnJ\s*=\s*(\d+)/i;
	$nK=$1           if /\bnK\s*=\s*(\d+)/i;
	$nPType=$1       if /\bnKindParticle\s*=\s*(\d+)/i;
	$GhostCell=$1    if /\bnG\s*=\s*(\d)/;	
    }
    close FILE;
    open(FILE, $NameSizeFile) or die "$ERROR could not open $NameSizeFile\n";
    while(<FILE>){
	next if /^\s*!/; # skip commented out lines
        $MaxBlock=$1     if /\bMaxBlock\s*=\s*(\d+)/i;
	$nElectronMax=$1        if /\bnElectronMax\s*=\s*(\d+)/i;
	$Hybrid=$1        if /\bnHybrid\s*=\s*(\d+)/i;
    }
    close FILE;

    die "$ERROR could not read MaxBlock from $NameSizeFile\n" 
	unless length($MaxBlock);


    $GridSize = "$nI,$nJ,$nK,$MaxBlock,$nPType,$nElectronMax";

}

#############################################################################

sub set_grid_size{

    $GridSize = $NewGridSize;

    if($GridSize=~/^[1-9]\d*,[1-9]\d*,[1-9]\d*,[0-9]\d*,[0-9]\d*,[0-9]\d*
       $/x){
	($nI,$nJ,$nK,$MaxBlock,$nPType,$nElectronMax)= split(',', $GridSize);
    }elsif($GridSize){
	die "$ERROR -g=$GridSize should be ".
	    #"4". 
	    " positive integers separated with commas\n";
    }
    $GhostCell = $NewGhostCell if $NewGhostCell;
    $Hybrid    = $NewHybrid    if $ReadHybrid;

 
    print "Writing new grid size $GridSize into ".
	"$NameSizeFile ...\n";

    @ARGV = ($NameSizeFile);
    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(nElectronMax\s*=[^0-9]*)(\d+)/$1$nElectronMax/i;
	s/\b(MaxBlock\s*=[^0-9]*)(\d+)/$1$MaxBlock/i;
	s/\b(nHybrid\s*=[^0-9]*)(\d+)/$1$Hybrid/i;
	print;
    }
    @ARGV = ($NameBatlFile);
    while(<>){
	if(/^\s*!/){print; next} # Skip commented out lines
	s/\b(nKindParticle\s*=[^0-9]*)(\d+)/$1$nPType/i;
	s/\b(nI\s*=[^0-9]*)(\d+)/$1$nI/i;
	s/\b(nJ\s*=[^0-9]*)(\d+)/$1$nJ/i;
	s/\b(nK\s*=[^0-9]*)(\d+)/$1$nK/i;
	s/\b(nG\s*=[^0-9]*)\d/$1$GhostCell/i;
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
    $Settings .= 
	"Max. number of electrons          : nElectronMax=$nElectronMax\n";
    $Settings .= 
	"Use hybrid scheme                 : nHybrid=$Hybrid\n";

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

