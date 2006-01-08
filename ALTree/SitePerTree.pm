package ALTree::SitePerTree;

################################################################
################################################################
####################### SitePerTree   ##########################
################################################################
################################################################

use base qw(ALTree::Base ALTree::Site);
use ALTree::SiteSensPerTree;

# Structure SitePerTree
#   "site_nb" -> Integer
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
    return ALTree::SiteSensPerTree->New($sens, $self);
}

1;
