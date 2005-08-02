package ALTree::Sens;

################################################################
################################################################
####################### Sens          ##########################
################################################################
################################################################

use base 'ALTree::Base';

#Creation d'une structure
sub New { # [classe] sens  
    my $class=shift;
    my $self={};
    my $sens=shift;
    bless($self, $class);
    
    if ($sens !~ /^\s*([0-9a-zA-Z?]+)\s*[-=_]*[>]\s*([0-9a-zA-Z?]+)\s*$/) {
	die "Invalid sens : $sens\n";
    }

    $self->_init("start" => $1, "end" => $2, "switch" => 0, @_);
    return $self;
}

#Creation d'une structure
sub NewRev { # [classe] sens  
    my $class=shift;
    my $self={};
    my $sens=shift;
    bless($self, $class);
    
    $self->_init("start" => $sens->{"end"}, 
		 "end" => $sens->{"start"}, "switch" => 0, @_);
    return $self;
}

sub _makeLabel {
    my $self=shift;
    my $start=shift;
    my $end=shift;

    if ($self->{'switch'}) {
	return $end."-->".$start;
    } else {
	return $start."-->".$end;
    }
}

sub GetLabel {
    my $self=shift;
    return $self->_makeLabel($self->{"start"}, $self->{"end"});
}

sub GetLabelRev {
    my $self=shift;
    return $self->_makeLabel($self->{"end"}, $self->{"start"});
}

sub Switch {
    my $self=shift;

    $self->{'switch'} = 1-$self->{'switch'};
}

1;
