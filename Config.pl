#!/usr/bin/perl -i
use strict;

our $Component = "PC";
our $Code = "ALTOP";
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

&print_help if $Help;

# Equation and user module variables
my $Src         = 'src';


# Grid size variables
my $NameSizeFile = "$Src/PIC_ModSize.f90";
my $GridSize;
my ($nX, $nY, $nZ, $nPType, $MaxBlock);

# Read previous grid size, equation and user module
&get_settings;

foreach (@Arguments){
    if(/^-s$/)                {$Show=1;                        next};

    warn "WARNING: Unknown flag $_\n" if $Remaining{$_};
}

&set_grid_size if $NewGridSize and $NewGridSize ne $GridSize;
# Show grid size in a compact form if requested
print "Config.pl -g=$nX,$nY,$nZ,$nPType,$MaxBlock",
    ,"\n" if $ShowGridSize and not $Show;


my $Settings = &current_settings; print $Settings if $Show;


exit 0;

#############################################################################

sub get_settings{

    # Read size of the grid from $NameSizeFile
    open(FILE, $NameSizeFile) or die "$ERROR could not open $NameSizeFile\n";
    while(<FILE>){
	next if /^\s*!/; # skip commented out lines
        $nX=$1           if /\bnX\s*=\s*(\d+)/i;
	$nY=$1           if /\bnY\s*=\s*(\d+)/i;
	$nZ=$1           if /\bnZ\s*=\s*(\d+)/i;
        $MaxBlock=$1     if /\bMaxBlock\s*=\s*(\d+)/i;
	$nPType=$1        if /\bnPType\s*=\s*(\d+)/i;
    }
    close FILE;

    die "$ERROR could not read MaxBlock from $NameSizeFile\n" 
	unless length($MaxBlock);


    $GridSize = "$nX,$nY,$nZ,$nPType,$MaxBlock";

}

#############################################################################

sub set_grid_size{

    $GridSize = $NewGridSize;

    if($GridSize=~/^[1-9]\d*,[1-9]\d*,[1-9]\d*,[0-9]\d*,[1-9]\d*
       $/x){
	($nX,$nY,$nZ,$nPType,$MaxBlock)= split(',', $GridSize);
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
	s/\b(MaxBlock\s*=[^0-9]*)(\d+)/$1$MaxBlock/i;
	s/\b(nX\s*=[^0-9]*)(\d+)/$1$nX/i;
	s/\b(nY\s*=[^0-9]*)(\d+)/$1$nY/i;
	s/\b(nZ\s*=[^0-9]*)(\d+)/$1$nZ/i;
	s/\b(nPType\s*=[^0-9]*)(\d+)/$1$nPType/i;
	print;
    }
}




sub current_settings{

    $Settings = 
	"Number of cells in a block        : nX=$nX, nY=$nY, nZ=$nZ\n";
    $Settings .= 
	"Number of different types of particles        : nPType=$nPType\n";
    $Settings .= 
	"Max. number of blocks/PE          : MaxBlock=$MaxBlock\n";

}

#############################################################################

sub print_help{

    print "
Additional options for BATSRUS/Config.pl:

-g=NI,NJ,NK,MAXBLK   
                Set grid size. NI, NJ and NK are the number of cells 
                in the I, J and K directions, respectively. 
                MAXBLK is the maximum number of blocks per processor.


Examples for BATSRUS/Config.pl:


Set block size to 8x8x8, number of blocks to 400",
":

    Config.pl -g=8,8,8,400",
"

Show settings for BATSRUS:

    Config.pl -s
\n";
    exit 0;


}

