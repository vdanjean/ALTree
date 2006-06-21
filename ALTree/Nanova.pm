package ALTree::Nanova;

use strict;

BEGIN {
    use Exporter   ();
    use vars       qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
    
    # set the version for version checking
    #$VERSION     = 1.00;
    # if using RCS/CVS, this may be preferred
    #$VERSION = do { my @r = (q$Revision: 153 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r }; # must be all one line, for MakeMaker
    
    @ISA         = qw(Exporter);
    @EXPORT      = qw(); #(&func1 &func2 &func4);
    %EXPORT_TAGS = ( );     # eg: TAG => [ qw!name1 name2! ],
    
    # your exported package globals go here,
    # as well as any optionally exported functions
    @EXPORT_OK   = qw();
}
use vars      @EXPORT_OK;

use ALTree::Utils qw(erreur);
use Data::Dumper;


# This function transforms the tree structure into the matrix used by the library NAnova 

sub Tree2mat
{
    my $present_node = shift; 
    my $vect = shift; # structure transitoire: chemin de la racine à une feuille
    my $mat = shift; # The matrix which is filled by the function
    
    if ($present_node->NbChildren()==0)  {
	push (@{$mat}, $vect);
    } else {
	my $number=0;
	for my $child ($present_node->GetChildrenList()) {  
	    push(@{$vect}, $number);
	    $number ++;
	    Tree2mat($present_node, $vect, $mat);
	}
    }
}

sub WriteMat
{
    my $mat = shift;
    
    for (my $i=0; $i<=$#$mat; $i++)  {
	foreach my $elem (@{$mat->[$i]}) {
	    print $elem, "\t";
	}
	print "\n";
    }
}
 
1;
