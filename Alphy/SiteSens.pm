package Alphy::SiteSens;

################################################################
################################################################
####################### SiteSens        ########################
################################################################
################################################################

use base 'Alphy::Base';

# Structure SiteSens
#   "site_struct" -> Site
#   "sens" -> Sens
#   "rev" -> SiteSens

sub InitSiteSens { # [Obj] Sens_label SiteSens
    my $self=shift;
    my $sens=shift;
    my $site=shift;

    $self->_init("sens" => $sens, "site_struct" => $site, @_);
    $self->Debug("creating SiteSens $sens\n");
}

sub GetSensLabel {
    my $self=shift;
    return $self->{"sens"}->GetLabel();
}

sub GetSensStruct {
    my $self=shift;
    return $self->{"sens"};
}

sub GetSite {
    my $self=shift;
    return $self->{"site_struct"};
}

sub GetSiteNb {
    my $self=shift;
    return $self->GetSite()->GetSiteNb();
}

sub LinkRev {
    my $siteSens=shift;
    my $siteSensRev=shift;

    $siteSens->{'rev'}=$siteSensRev;
    $siteSensRev->{'rev'}=$siteSens;
}

sub GetSensRev {
    my $self=shift;

    if (not defined($self->{'rev'})) {
	die "LinkRev not called\n";
    }
    return $self->{'rev'};
}    

1;
