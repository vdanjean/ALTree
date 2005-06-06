package Alphy::SiteCollection;

################################################################
################################################################
####################### SiteCollection  ########################
################################################################
################################################################

use base 'Alphy::Base';

sub InitSiteCollection {
    my $self=shift;
    $self->_init("sites" => {});
}

sub AddSite {
    my $self=shift;
    my $site=shift; 
    my $site_nb;
    $site_nb=$site->GetSiteNb();
    $self->{"sites"}->{$site_nb}=$site;
}
sub GetSite {
    my $self=shift;
    my $site_nb=shift;
    #my $site=$self->{"sites"}->{$site_nb};
    #if (not defined($site)) {
    #die "The site number $site_nb does not exist";
    #}
    return $self->{"sites"}->{$site_nb};
}
sub HasSiteIndex {
    my $self=shift;
    my $site_nb=shift;
    return exists($self->{"sites"}->{$site_nb});
}
sub GetSitesList {
    my $self=shift;
    return values(%{$self->{"sites"}});
}
1;
