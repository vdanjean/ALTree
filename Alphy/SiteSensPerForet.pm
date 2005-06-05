package Alphy::SiteSensPerForet;

################################################################
################################################################
####################### SiteSensPerForet #######################
################################################################
################################################################

use base qw(Alphy::Base Alphy::SiteSens);

# Structure SiteSensPerForet
#   "site_struct" -> Site
#   "sens_label" -> String
#   "V_i" -> Float

sub New { # [classe] sens_label site_struct
    my $class=shift;
    my $self={};
    my $sens=shift;
    my $site=shift;
    bless($self, $class);
    $self->InitSiteSens($sens, $site);
    $self->_init("V_i" => 0, @_);
    return $self;
}

sub PlusVi {
    my $self=shift;
    my $value=shift;
    $self->{"V_i"} += $value;
}
sub GetVi {
    my $self=shift;
    return $self->{"V_i"};
}    

