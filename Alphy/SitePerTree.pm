package Alphy::SitePerTree;

################################################################
################################################################
####################### SitePerTree   ##########################
################################################################
################################################################

use base qw(Alphy::Base Alphy::Site);
use Alphy::SiteSensPerTree;

# Structure SitePerTree
#   "site_nb" -> Integer
#   "CI" -> Float (=1/nb_mut)
#   "sens_struct" -> Hash of ('sens_label' -> SiteSens)
#   "nb_mut" -> Integer

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
    return Alphy::SiteSensPerTree->New($sens, $self);
}

sub SetCI {
    my $self=shift;
    my $CI=shift;
    $self->{"CI"}=$CI;
}
sub GetCI {
    my $self=shift;
    if (exists($self->{"CI"})) {
	return $self->{"CI"};
    } else {
	my $nb_mut=$self->GetNbMut();
	if ($nb_mut != 0) {
	    $self->SetCI(1/$nb_mut);
	    return 1/$nb_mut;
	}
	die "Unable to provide CI\n";
    }
}

sub IncNbMut {
    my $self=shift;
    if (exists($self->{"CI"})) {
	die "IncNbMut called after GetCI or SetCI";
    }
    $self->{"nb_mut"}++;
}

1;
