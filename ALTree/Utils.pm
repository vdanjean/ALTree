package ALTree::Utils;
use Pod::Usage;

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
	"bardel\@vjf.inserm.fr\n");
}

1;

package DataType;
use constant SNP => 0;
use constant DNA => 1;

package PhylProg;
use constant PHYLIP => 0;
use constant PAUP => 1;
use constant PAML => 2;
