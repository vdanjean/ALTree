package Alphy::SiteSensPerTree;

################################################################
################################################################
####################### SiteSensPerTree ########################
################################################################
################################################################

use base qw(Alphy::Base Alphy::SiteSens);

# Structure SiteSens
#   "site_struct" -> Site
#   "sens_label" -> String
#   "node_list" -> Array of (Node)
#    "m_it" -> Interger # nb mutation of this chage in the tree
#    "R_it" -> Interger # nb co-mutations of this change with character S
#    "V_it" -> Integer  # (R_it-E_it)/sqrt(E_it)

sub New { # [classe] sens_label site_struct
    my $class=shift;
    my $self={};
    my $sens=shift;
    my $site=shift;
    bless($self, $class);
    $self->InitSiteSens($sens, $site);
    $self->_init("node_list" => [], "R_it" => 0, @_);
    return $self;
}

# Appelé par Node->AddApo
sub _AddNode {
    my $self=shift;
    my $node=shift;
    if (exists($self->{"m_it"})) {
	die "_AddNode called after GetMit";	
    }
    push @{$self->{"node_list"}}, $node;
}
# Appelé par Node->DeleteAllApo
sub _DeleteNode {
    my $self=shift;
    my $node=shift;
    my @new_node_list=grep { $_ != $node } @{$self->{"node_list"}};
    if (scalar(@new_node_list)+1 != $self->NbNodes()) {
	die "Error while removing a node from a SiteSensPerTree";
    }
    $self->{"node_list"}=\@new_node_list;
}
sub NbNodes {
    my $self=shift;

    return scalar (@{$self->{"node_list"}});
}
sub GetNode {
    my $self=shift;
    my $index=shift;
    return $self->{"node_list"}->[$index];
}
sub GetNodesList {
    my $self=shift;
    return @{$self->{"node_list"}};
}

sub SetStep {
    my $self=shift;
    my $step=shift;
    $self->{"Step"}=$step;
}

sub GetStep {
    my $self=shift;
    my $step=shift;
    return $self->{"Step"};
}

sub GetMit {
    my($self)=shift;
    if (not exists($self->{"m_it"})) {
	$self->{"m_it"}=$self->NbNodes();
    }
    return ($self->{"m_it"});
}

sub IncRit {
    my($self)=shift;
    if (exists($self->{"V_it"})) {
	die "IncRit called after GetVit";	
    }    
    $self->{"R_it"}++;
}
sub GetRit {
    my($self)=shift;
    return ($self->{"R_it"});
}

sub GetVit {
    my $self=shift;

    my $Rit=$self->GetRit();
    my $Eit=$self->GetEit();

    if (not exists($self->{"V_it"})) {
	if ($Eit == 0) {
	    $self->{"V_it"}=0;
	} else {
	    $self->{"V_it"}=($Rit-$Eit)/sqrt($Eit);
	}
    }
    return ($self->{"V_it"});
}

sub SetEit {
    my($self)=shift;
    my($value)=shift;
    if (exists($self->{"V_it"})) {
	die "SetEit called after GetVit";	
    }    
    $self->{"E_it"}=$value;
}

sub GetEit {
    my($self)=shift;
    return ($self->{"E_it"});
}
1;
