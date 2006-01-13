package ALTree::Tree;

################################################################
################################################################
####################### Tree            ########################
################################################################
################################################################

use base qw(ALTree::Base ALTree::SiteCollection);

# "nodes" -> Hash of ('id' -> Node)
# "sites" -> Hash of ('site_nb' -> SitePerTree)

sub New { # [classe]
    my $class=shift;
    my $self={};
    bless($self, $class);
    $self->InitSiteCollection();
    $self->_init("nodes" => {}, @_);
    $self->Debug("creating Tree\n");
    return $self;
}

sub AddNode {
    my $self=shift;
    my $node=shift;

    my $id=$node->GetId();
    $self->{"nodes"}->{$id}=$node;
}
sub GetNode {
    my $self=shift;
    my $id=shift;
    if (defined ($self->{"nodes"}->{$id})) {
	return $self->{"nodes"}->{$id};
    } else {
	return undef;
    }
}    
sub HasNodeIndex {
    my $self=shift;
    my $id=shift;
    return exists($self->{"nodes"}->{$id});
}
sub GetNodesIndexList {
    my $self=shift;
    return keys(%{$self->{"nodes"}});
}
sub GetNodesList {
    my $self=shift;
    return values(%{$self->{"nodes"}});
}

#sub _SetOutgroup {
#    my $self=shift;
#    my($outgroup)=shift;
#    $self->{"outgroup"}=$outgroup;
#}
sub GetOutgroup {
     my $self=shift;
     if (not exists($self->{"outgroup"})) {
       FIND: {
	   foreach my $clef ($self->GetNodesIndexList()) {
	       if ($clef eq "H000") {
		   $self->{"outgroup"}=$self->GetNode($clef);
		   last FIND;
	       }
	   }
	   die "No outgroup corresponding to id=H000 identified in the tree\n";
       }
     }
     return ($self->{"outgroup"});
}

sub SetNbBrNonNulle {
    my $self=shift;
    my($value)=shift;
    $self->{"nb_br_non_nulle"}=$value;
} 
sub GetNbBrNonNulle {
    my $self=shift;
    return ($self->{"nb_br_non_nulle"});
}

sub _SetRoot {
    my $self=shift;

    my(@roots);
    #print STDERR "nodes: \n";
    foreach my $node ($self->GetNodesList()) {
	if (not $node->HasFather()) {
	    push @roots, $node;
	    #print STDERR "root: ", $node->Name(), "\n";
	}	    
	#else { print STDERR "node: ", $node->Name(), " father ", $node->GetFather()->Name(), "\n"; }
    }
    #Verification that there is only one root
    if (scalar(@roots) !=1) {
	die "error in the tree: ", scalar(@roots), " roots for the tree!\n";
    }
    $self->{"root"}=$roots[0];
} 
sub GetRoot {
    my $self=shift;
    if (not exists($self->{"root"})) {
	$self->_SetRoot();
    }
    return ($self->{"root"});
}

sub ChangeRoot {
    my $self=shift;
    my $newroot=shift;

    if ($newroot == $self->GetRoot()) {
	return;
    }
    my $oldfather=$newroot->GetFather();
    if (not defined($oldfather)) {
	die ("Non root node of the tree has no father\n"
	     ."Do you use a node of the correct tree ?");
    }
    $self->ChangeRoot($oldfather);
    
    # Some integrity tests...
    if ($oldfather->NbApo() != 0) {
	die "Root has apomorphies !";
    }
    if ($oldfather->GetBrLen() != 0) {
	die "Root has non null BrLen !";
    }
    
    foreach my $apo ($newroot->GetApoList()) {
	$oldfather->AddApo($apo->GetSensRev());
    }
    $newroot->DeleteAllApo();
    #print "Setting BRLen to ",$newroot->GetBrLen()," for ", $oldfather->Name()," from ", $newroot->Name(),"\n";
    $oldfather->SetBrLen($newroot->GetBrLen());
    $newroot->SetBrLen(0);

    $oldfather->SetFather($newroot);
    $newroot->SetFather(undef);
    $newroot->AddChild($oldfather);
    $oldfather->DeleteChild($newroot);

    $self->{"root"}=$newroot;
}

1;
