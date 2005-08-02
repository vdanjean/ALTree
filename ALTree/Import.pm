package ALTree::Import;

use strict;
use ALTree::Sens;
use ALTree::Tree;
use ALTree::Foret;
use ALTree::Node;
use ALTree::SitePerTree;
use ALTree::SitePerForet;
use ALTree::SiteSensPerForet;

BEGIN {
    use Exporter   ();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    
    # set the version for version checking
    #$VERSION     = 1.00;
    # if using RCS/CVS, this may be preferred
    $VERSION = do { my @r = (q$Revision$ =~ /\d+/g); 
		    sprintf "%d."."%02d" x $#r, @r }; # must be all one line, for MakeMaker
    
    @ISA         = qw(Exporter);
    @EXPORT      = qw();
    #%EXPORT_TAGS = ( );     # eg: TAG => [ qw!name1 name2! ],
    
    # your exported package globals go here,
    # as well as any optionally exported functions
    @EXPORT_OK   = qw();
}
#our @EXPORT_OK;

INIT { 
    
}
 
# exported package globals go here
#our $Var1;
#our %Hashit;

# non-exported package globals go here
#our @more;
#our $stuff;

# initialize package globals, first exported ones
#$Var1   = '';
#%Hashit = ();

# then the others (which are still accessible as $Some::Module::stuff)
#$stuff  = '';
#@more   = ();

# all file-scoped lexicals must be created before
# the functions below that use them.

# file-private lexicals go here
#my $priv_var    = '';
#my %secret_hash = ();

# here's a file-private function as a closure,
# callable as &$priv_func;  it cannot be prototyped.
#my $priv_func = sub {
#    # stuff goes here.
#};

# make all your functions, whether exported or not;
# remember to put something interesting in the {} stubs
#sub func1      {}    # no prototype
#sub func2()    {}    # proto'd void
#sub func3($$)  {}    # proto'd to 2 scalars

# this one isn't exported, but could be called!
#sub func4(\%)  {}    # proto'd to 1 hash ref

END { }       # module clean-up code here (global destructor)

##################################################################
# Fonctions de seuil du chi2 (pré-calcul puis stockage dans un tableau)

1;
