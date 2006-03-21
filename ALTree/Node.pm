package ALTree::Node;

################################################################
################################################################
####################### Node          ##########################
################################################################
################################################################

use base 'ALTree::Base';

# Structure Node
#   "id" -> String
#   "children" -> Array of (Node)
#   "father" -> Node
#   "apo" -> Hash of ('num_apo' => SiteSens)
#   "br_len" -> Integer
#   "case" -> Integer
#   "control" -> Integer
#   "level" -> Integer 
#   "height" -> Integer
#   "sequence" ->  string # used only with phylip
#   "oldfather" -> Node  # après fusion des branches nulles (anciennes relations)
#   "oldchildren" -> Array of (Node)  # après fusion des branches nulles (anciennes relations de parenté)
#   "label" -> string # nom noeuds apres fusion 3+(4) par exemple 
#   "res_test" -> Array of (TestResults)

#Creation d'une structure
sub New { # [classe] id  
    my $class=shift;
    my $self={};
    my $id=shift;
    bless($self, $class);
    $self->_init("id" => $id, "apo" => {}, 
		 "children" => [], @_);
    $self->Debug("creating Node $id\n");
    return $self;
}

sub GetId {
    my $self=shift;
    return $self->{"id"};
}

sub SetFather {
    my $self=shift;
    my $father=shift;
    $self->{"father"}=$father;
}
sub GetFather {
    my $self=shift;
    return $self->{"father"};
}
sub HasFather {
    my $self=shift;
    return defined($self->{"father"});
}
sub Father {
    my $self=shift;
    my $newfather=shift;
    if (defined($newfather)) {
	$self->SetFather($newfather);
    }
    return $self->GetFather();
}
sub RecordOrigFather {
    my $self=shift;
    my $father=shift;
    if ($self->HasFather()) {
    	$self->{"orig_father"}=$self->GetFather();
    }
}
sub GetOrigFather {
    my $self=shift;
    if (exists($self->{"orig_father"})) {
    	return $self->{"orig_father"};
    } else {
    	return $self->GetFather();
    }
}

sub SetCase {
    my $self=shift;
    my $value=shift;
    $self->{"case"}=$value;
}
sub GetCase {
    my $self=shift;
    return $self->{"case"};
}
sub EraseCase {
    my $self=shift;
    delete($self->{"case"});
}


sub SetControl {
    my $self=shift;
    my $value=shift;
    $self->{"control"}=$value;
}
sub GetControl {
    my $self=shift;
    return $self->{"control"};
}
sub EraseControl {
    my $self=shift;
    delete($self->{"control"});
}

sub EraseQuanti {
    my $self=shift;
    delete($self->{"quanti"});
}

sub SetBrLen {
    my $self=shift;
    my $br_len=shift;
    $self->{"br_len"}=$br_len;
}
sub GetBrLen {
    my $self=shift;
    return $self->{"br_len"};
}
sub HasBrLen {
    my $self=shift;
    return exists($self->{"br_len"});
}

sub AddOldChild {
    my $self=shift;
    push @{$self->{"oldchildren"}}, @_;
}
sub GetOldChildrenList {
    my $self=shift;
    return @{$self->{"oldchildren"}};
}

sub AddChild {
    my $self=shift;
    push @{$self->{"children"}}, @_;
}
sub GetChildrenList {
    my $self=shift;
    return @{$self->{"children"}};
}
sub DeleteChild {
    my $self=shift;
    my $todelete=shift;
    my @newchidren=grep {$_ != $todelete} $self->GetChildrenList();
    $self->{"children"}=\@newchidren;
    return;
#    my $children=$self->{"children"};
#    my($i);
#    for ($i=0; $i<=$#$children; $i++) { 
#	if ($children->[$i]==$todelete) {
#	    splice(@{$children}, $i, 1);
#	    return;
#	}
#    }
}
sub NbChildren {
    my $self=shift;
    return scalar(@{$self->{"children"}});
}
sub GetChild {
    my $self=shift;
    my $index=shift;

    return $self->{"children"}->[$index];
}
sub ForgetChildren {
    my $self=shift;
    $self->{"children"}=[];
}

sub SetSequence {
    my($self)=shift;
    my($sequence)=shift;
    $self->{"sequence"}=$sequence;
}
sub GetSequence {
    my($self)=shift;
    return $self->{"sequence"};
}

sub AddApo {
    my $self=shift;
    my $site_sens=shift;
    my $apo=$site_sens->GetSiteNb();
    if (exists($self->{"apo"}->{$apo})) {
	die "Adding aleardy present apo $apo in ", $self->Name(),"\n";
    }
    $self->Debug("Adding apo $apo in ", $self->Name(),"\n");
    $self->{"apo"}->{$apo}=$site_sens;
    $site_sens->_AddNode($self);
}
sub NbApo {
    my $self=shift;
    return scalar (keys(%{$self->{"apo"}}));
}
sub GetApoIndexList {
    my $self=shift;
    return keys(%{$self->{"apo"}});
}
sub GetApoList {
    my $self=shift;
    return values(%{$self->{"apo"}});
}
sub GetApo {
    my $self=shift;
    my $apo=shift;
    return $self->{"apo"}->{$apo};
}
sub DeleteAllApo {
    my $self=shift;
    foreach my $apo ($self->GetApoList()) {
	$apo->_DeleteNode($self);
    }
    $self->{"apo"}={};
}

##Modifie Claire
sub NbApoStep {
    my $self=shift;
    my($compteur)=0;
    my($apo);
    foreach $apo ($self->GetApoList()) {
	$compteur+=$apo->GetStep();
    }
    return $compteur;
}

sub Name {
    my $self=shift;
    if (defined ($self->{"label"})) {
	return $self->{"label"};
    } else {
	return $self->GetId();
    }
}

#Faire meme chose pour case, control et br_len (HasBrLen), level, high, sequence, label, oldfather
#Faire oldchildren (GetOldChildren, SetOldChildren, ) 


1;
