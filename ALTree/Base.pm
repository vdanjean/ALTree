package ALTree::Base;

###########################################
########  MAIN DATA STRUCTURES    #########
###########################################

sub _init {
    my $self=shift;
    if (@_) {
        my %extra = @_;
        @$self{keys %extra} = values %extra;
    }
}

sub Debug {
    my $self=shift;

    #print STDERR @_;
}

1;
