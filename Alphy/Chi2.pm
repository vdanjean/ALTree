package chi2;

use strict;
use Alphy::CUtils;

BEGIN {
    use Exporter   ();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    
    # set the version for version checking
    #$VERSION     = 1.00;
    # if using RCS/CVS, this may be preferred
    $VERSION = do { my @r = (q$Revision$ =~ /\d+/g); 
		    sprintf "%d."."%02d" x $#r, @r }; # must be all one line, for MakeMaker
    
    @ISA         = qw(Exporter);
    @EXPORT      = qw(NON_SIGNIFICATIF SIGNIFICATIF &chi2_significatif &definition_p_chi2 &reech_chi2);
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

use constant NON_SIGNIFICATIF => 0;
use constant SIGNIFICATIF => 1;
use constant PERM => 1000;
##################################################################
# Fonctions de seuil du chi2 (pré-calcul puis stockage dans un tableau)

 my ($chi2_p)="chi2_p doit être initialisé !";
#y ($chi2_p)=0.05;
my ($test_prop_p)="test_prop_p doit être initialisé !";
my ($chi2_seuil)=[];
# test de significativité. Doit retourner vrai ou faux.
sub chi2_significatif {
    my ($ddl) = shift;
    my ($chi2) = shift;

    if ($ddl < 1) {
	warn "chi[$ddl] asked...\n";
    }

    if (not defined($chi2_seuil->[$ddl])) {
	#my $c=`critchi2 $chi2_p $ddl`+0; # Verif que les 2 appels sont équivalents
	$chi2_seuil->[$ddl]=CUtils::critchi($chi2_p, $ddl);
	#if ($c != $$chi2_seuil[$ddl]) {
	  # print STDERR "Critchi2 : $c != $$chi2_seuil[$ddl]\n";
	#}
	#print "chi2_seuil[$ddl]= $$chi2_seuil[$ddl]\n";
	#warn "Seuil chi2 non défini pour ddl $ddl";
	#return 0;
    }
    return ($chi2 > $chi2_seuil->[$ddl]);
}

sub definition_p_chi2
{
    my($p)=shift;
    my($pprop)=shift;
    if (defined $p) {
	$chi2_p=$p;
    }
    if (defined $pprop) {
	$test_prop_p=$pprop;
    }
}

sub chi2_fisher_significatif
{
    my($pvalue)=shift;

    return ($pvalue < $chi2_p);
}

##################################################################
# Rééchantillonnage

# Quelques variables globales pour aller plus vite (même si c'est à éviter
# en général)
my(@fils_c);
my(@fils_m);
my($compteur);
my($sum_malade);
my($sum_controle);
my($sum_total);
my($nb_fils);


my(@th_c, @th_m);
my($clades);
sub calcule_chi2
{
    my($i, $chi2, $temp);
    $chi2=0;
    for ($i=0; $i < $nb_fils; $i++){
	$temp=($fils_c[$i]-$th_c[$i]);
	$chi2+=($temp)*($temp)/$th_c[$i];

	$temp=($$clades[$i]-$fils_c[$i]-$th_m[$i]);
	$chi2+=($temp)*($temp)/$th_m[$i];
    }
    #print "Chi2 : $chi2\n";
    return $chi2;
}

sub reech_chi2
{
    $sum_malade=shift;
    $sum_controle=shift;
    $sum_total=$sum_malade+$sum_controle;
    $nb_fils=shift;
    my($chi2_reel)=shift;
    $clades=shift;
    my($p_val)=0;
    my($i, $k);
    #my($alea);

# nb_fils correspond en fait a tous les groupes sur lesquelles il faut faire
# le chi2.
# Cet ensemble de groupe a au total: $sum_malade et $sum_controle individus
# (respectivement malades et controles)
# clades est une référence sur un tableau contenant les effectifs globaux de
# chaque clade

    
    #print "Reechantillonage chi2 : ddl : ", ($nb_fils-1), " M: $sum_malade C: $sum_controle\n   ";
    #print "Chi2 réel : $chi2_reel\n";
    #print "Clades: ";
    #for $i (@{$clades}) { print "$i "; } print "\n";

    # Pré calcul des effectifs théoriques
    for ($i=0; $i < $nb_fils; $i++){
	$th_c[$i]=($sum_controle*$$clades[$i])/($sum_total);
	$th_m[$i]=($sum_malade*$$clades[$i])/($sum_total);
    }

    my($clade, $alea, $malades, $controles);
    my($nb_tests)=PERM;
    for ($k=1;$k<=$nb_tests; $k++){
	$malades=$sum_malade;
	$controles=$sum_controle;
	for($clade=0; $clade<$nb_fils; $clade++) {
	    $fils_m[$clade]=0;
	    $fils_c[$clade]=0;
	    for($i=0; $i<$$clades[$clade]; $i++) {
		$alea=rand($malades+$controles);
		if ($alea < $malades) {
		    $malades--;
		    $fils_m[$clade]++;
		} else {
		    $controles--;
		    $fils_c[$clade]++;
		}
	    }
	}

	if (calcule_chi2 >= $chi2_reel){
	    $p_val++;
	}
    }
    #DEBUG  print "CHI2=$chi2_reel\n";
    #print"nb de chi2 non calculable (effectif nul)= $compteur\n";
    #DEBUG print "p_val1 = ", $p_val/$nb_tests,"\n";
    # print "chi2_p771=$chi2_p\n";
    return ($p_val/$nb_tests);
}

sub reech_significatif 
{
    my($p_val)=shift;
    my($nb_tests)=PERM;
  #DEBUG  print "Chi2P= $chi2_p\n";
  #DEBUG  print "p=  ", $p_val , "\n";
  #DEBUG  print "test=", $p_val/$nb_tests<=$chi2_p, "\n";
    return ($p_val<=$chi2_p);
}

1;
