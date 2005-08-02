package ALTree::Foret;

################################################################
################################################################
####################### Foret           ########################
################################################################
################################################################

use base qw(ALTree::Base ALTree::SiteCollection);

sub New { # [classe]
    my $class=shift;
    my $self={};
    bless($self, $class);
    $self->InitSiteCollection();
    $self->_init(@_);
    $self->Debug("creating Foret\n");
    return $self;
}

sub AddTree {
    my $self=shift;
    push @{$self->{"trees"}}, @_;
}
sub GetTree {
    my $self=shift;
    my $index=shift;
    return $self->{"trees"}->[$index];
}    
sub GetTreesList {
    my $self=shift;
    return @{$self->{"trees"}};
}

sub ProvideSite {
    my $self=shift;
    my $site_nb=shift;

    if (not $self->HasSiteIndex($site_nb)) {
	$self->AddSite(ALTree::SitePerForet->New($site_nb));
    }
    return $self->GetSite($site_nb);
}

sub CalculVi {
    my $self=shift;

    foreach my $tree ($self->GetTreesList()) {
	foreach my $site ($tree->GetSitesList()) {
	    foreach my $sens ($site->GetSensList()) {
		$self->ProvideSite($site->GetSiteNb())
		    ->ProvideSens($sens->GetSensStruct())
		    ->PlusVi($sens->GetVit());
	    }
	}
    }
}

sub _EnsureViMax {
    my($self)=shift;
    if (not exists ($self->{"V_i_max"})) {
	my @tab_trie=sort {
	    $b->GetViMax() <=> $a->GetViMax()} 
	   $self->GetSitesList();
	$self->{"V_i_max"}=$tab_trie[0]->GetViMax();
	$self->{"V_i_max_tab"}=\@tab_trie;
    }
}
sub GetViMax {
    my $self=shift;
    $self->_EnsureViMax();
    return $self->{"V_i_max"};
}   
sub GetViMaxSite {
    my($self)=shift;
    my($index)=shift;
    $self->_EnsureViMax();
    return $self->{"V_i_max_tab"}->[$index];
}
sub GetViMaxSiteList {
    my($self)=shift; 
    $self->_EnsureViMax();
    return @{$self->{"V_i_max_tab"}}; 
}
sub NbViMaxSite {
    my($self)=shift; 
    $self->_EnsureViMax();
    return (scalar @{$self->{"V_i_max_tab"}});
}

sub _EnsureViMaxSens {
    my($self)=shift;
    if (not exists ($self->{"V_i_max_sens_tab"})) {
	my @tab_trie=sort {
	    $b->GetVi() <=> $a->GetVi()
	    } map { $_->GetSensList(); } $self->GetSitesList();
	#if ($tab_trie[0]->GetVi() != $self->GetViMax()) {
	#    die "Arghh\n";
	#}
	$self->{"V_i_max_sens_tab"}=\@tab_trie;
    }
}
sub GetViMaxSens {
    my($self)=shift;
    my($index)=shift;
    $self->_EnsureViMaxSens();
    return $self->{"V_i_max_sens_tab"}->[$index];
}
sub GetViMaxSensList {
    my($self)=shift; 
    $self->_EnsureViMaxSens();
    return @{$self->{"V_i_max_sens_tab"}}; 
}
sub NbViMaxSens {
    my($self)=shift; 
    $self->_EnsureViMaxSens();
    return (scalar @{$self->{"V_i_max_sens_tab"}});
}

1;
