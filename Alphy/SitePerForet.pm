package Alphy::SitePerForet;

################################################################
################################################################
####################### SitePerForet  ##########################
################################################################
################################################################

use base qw(Alphy::Base Alphy::Site);
use sort '_quicksort';

sub New { # [classe] site_nb
    my $class=shift;
    my $self={};
    my $site_nb=shift;
    bless($self, $class);
    $self->InitSite($site_nb);
    $self->_init("nb_mut" => 0, @_);
    return $self;
}

sub NewSens {
    my $self=shift;
    my $sens=shift;
    return Alphy::SiteSensPerForet->New($sens, $self);
}

sub mysort {
    my $func=shift;
    my @tab=@_;

    my @tab2=();

    if (scalar(@tab) == 0) {
	return @tab2;
    }
 #   print "Site ", $tab[0]->GetSiteNb(), "\n";
#    foreach my $sens (@tab) {
#	print $sens->GetVi()," ", $sens->GetSensLabel(), "\n";
#    }
    my $size;
    while (($size=scalar(@tab)) > 1) {
	my $min=0;
	for(my $i=1; $i<$size; $i++) {
	    if ($func->($tab[$i], $tab[$min]) < 0) {
		$min=$i;
	    }
	}
	push @tab2, $tab[$min];
#	print "min=$min adding ",$tab[$min]->GetSensLabel(),"\n";
#	print "Old tab\n";
#	foreach my $sens (@tab) {
#	    print $sens->GetVi()," ", $sens->GetSensLabel(), "\n";
#	}
	#my @tab3=splice(@tab, $min, 1);
	my @tab3;
	for (my $i=0; $i<=$#tab; $i++) {
	    if ($i != $min) {
		push @tab3, $tab[$i];
	    }
	}
#	print "NewTab\n";
#	foreach my $sens (@tab3) {
#	    print $sens->GetVi()," ", $sens->GetSensLabel(), "\n";
#	}
	@tab=@tab3;
    }
    push @tab2, $tab[0];
#    print "Site ", $tab[0]->GetSiteNb(), " result\n";
#    foreach my $sens (@tab2) {
#	print $sens->GetVi()," ", $sens->GetSensLabel(), "\n";
#    }
    return @tab2;

}

sub functrie {
    my $a=shift;
    my $b=shift;
    #print "a=$a, b=$b\n";
    #if (not defined($a)) {
	#return 0;
	#print "a=$a b=$b\n";
	#die "beauk\n";
    #}
   # my $value=($b->GetVi() <=> $a->GetVi());
   # print "Comparing ", $b->GetVi(), " and ", $a->GetVi(), " = $value\n";
   # return $value;#
    return ($b->GetVi() <=> $a->GetVi());
}

sub _EnsureViMax {
    my($self)=shift;
    if (not exists ($self->{"V_i_max"})) {
	my @toto=$self->GetSensList();
	#print "site: ", $self->GetSiteNb(), "\n";
	#for my $ref (@toto) {
	    #print "Ref: $ref\n";
	#}
	my @tab_trie=mysort \&functrie, @toto;
	#for my $ref (@toto) {
	    #print "Ref: $ref\n";
	#}

	#my @tab3;
	#for (my $i=0; $i<=$#toto; $i++) {
	#    print "value $i = $toto[$i]\n";
	#    push @tab3, $toto[$i];
	#}

	#my @tab_trie=sort {
	#    $b->GetVi() <=> $a->GetVi();
	#} @tab3;
	
	my($Vimax)=$self->{"V_i_max"}=$tab_trie[0]->GetVi();
	for (my $i=1; $i<=$#tab_trie; $i++) {
	    if ($tab_trie[$i]->GetVi() !=  $Vimax) {
		splice(@tab_trie, $i); # Elimine tout à partir de indice $i
		last;
	    }
	}
	#print "Vi max : ", $Vimax, "\n";
	$self->{"V_i_max_tab"}=\@tab_trie;
    }
}
sub GetViMax {
    my $self=shift;
    $self->_EnsureViMax();
    return $self->{"V_i_max"};
}   
sub GetViMaxSens {
    my($self)=shift;
    my($index)=shift;
    $self->_EnsureViMax();
    return $self->{"V_i_max_tab"}->[$index];
}
#NbViMaxSens() => nombre sitesensperforet
sub NbViMaxSens {
    my($self)=shift; 
    $self->_EnsureViMax();
    return (scalar @{$self->{"V_i_max_tab"}});
}
#GetViMaxSensList() => array de tous les sitessensperforet (avec Vi==ViMax)
sub GetViMaxSensList {
    my($self)=shift; 
    $self->_EnsureViMax();
    return @{$self->{"V_i_max_tab"}}; # Comment ça marche ce truc?
}
 
1;
