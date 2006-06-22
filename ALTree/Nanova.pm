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
    my @vect=(); # dernier chenin parcouru
    my @mat;
    my $height=$present_node->{"height"};

    print STDERR "heigh=", $present_node->{"height"}, "\n";
    for (my $i=0; $i<$height; $i++) {
	push(@vect, -1);
    }

    my $tree2mat;
    $tree2mat = sub {
	my $present_node = shift;

	if ($present_node->NbChildren()==0)  {
	    for (my $i=$present_node->{"level"}; $i<$height; $i++) {
		$vect[$i]++;
	    }
	    my @tab=@vect;
	    push (@mat, \@tab);
	} else {
	    for my $child ($present_node->GetChildrenList()) {  
		$vect[$present_node->{"level"}]++;
		$tree2mat->($child);
	    }
	}
    };
    $tree2mat->($present_node);
    return \@mat;
}

# Fille the various tabular necessary for NAnaova
sub FillTableaux
{

    my $present_node = shift;
    my $values = shift;
    my $groups = shift;
    my $nb_term = shift;

    if ($present_node->NbChildren()==0)  {
	$presentNode->GetQuantiList();
	push @{$values}, @{$present_node->GetQuantiList()};
	push @{$groups}, $present_node->NbQuanti(); 
	$nb_term++;
    } else {
	for my $child ($present_node->GetChildrenList()) {  
	    FillTableaux($child, $value, $groups, $nb_term);
	}
    }


}

sub WriteMat
{

    my $mat = shift;

    print STDERR Dumper($mat);

    return;
    for (my $i=0; $i<=$#$mat; $i++)  {
	foreach my $elem (@{$mat->[$i]}) {
	    print $elem, "\t";
	}
	print "\n";
    }
}
 
1;
