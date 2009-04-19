package ALTree::Utils;
use Pod::Usage;

use strict;

BEGIN {
    use Exporter   ();
    use vars       qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
    
    # set the version for version checking
    #$VERSION     = 1.00;
    # if using RCS/CVS, this may be preferred
    #$VERSION = do { my @r = (q$Revision$ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r }; # must be all one line, for MakeMaker
    
    @ISA         = qw(Exporter);
    @EXPORT      = qw(); #(&func1 &func2 &func4);
    %EXPORT_TAGS = ( );     # eg: TAG => [ qw!name1 name2! ],
    
    # your exported package globals go here,
    # as well as any optionally exported functions
    @EXPORT_OK   = qw(&erreur &internal_error);
}
use vars      @EXPORT_OK;


sub erreur
{
    my $msg=shift;
    my $use_pod=shift;

    if (not defined($use_pod) or $use_pod) {
	pod2usage("Error: ".$msg);
    } else {
	print STDERR "Error: ".$msg;
	exit 1;
    }
}

sub internal_error
{
    my $msg=shift;

    die("Internal error: $msg\n".
	"Please, report this bug (with all is needed to reproduce it) to:\n".
	"Claire.Bardel\@univ-lyon1.fr\n");
}

1;

package DataType;
use constant SNP => 0;
use constant DNA => 1;

package PhylProg;
use constant PHYLIP => 0;
use constant PAUP => 1;
use constant PAML => 2;

package ANC;
use constant Rooted => 0;
use constant Unrooted => 1;
use constant OutGroup => 2;

no strict;

@Name=("rooted using an ancestral sequence",
       "unrooted", "rooted using an outgroup");
