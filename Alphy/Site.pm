package Alphy::Site;

################################################################
################################################################
####################### Site          ##########################
################################################################
################################################################

use base 'Alphy::Base';
use Alphy::Sens;

# Structure Site
#   "site_nb" -> Integer
#   "sens_struct" -> Hash of ('sens_label' -> SiteSens)

sub InitSite { # [obj] $site
    my $self=shift;
    my $site_nb=shift;
    $self->_init("site_nb" => $site_nb, "sens_struct" => {});
    $self->Debug("creating Site $site_nb\n");
}

sub GetSiteNb {
    my $self=shift;
    return $self->{"site_nb"};
}

sub HasSensIndex {
    my $self=shift;
    my $sens=shift;

    return exists($self->{"sens_struct"}->{$sens->GetLabel()});
}
sub AddSens {
    my $self=shift;
    my $sens=shift;

    my($ref_site_sens)=$self->NewSens($sens);
    $self->{"sens_struct"}->{$sens->GetLabel()}=$ref_site_sens;
    my $sensRev=Alphy::Sens->NewRev($sens);
    my($ref_site_sens_rev)=$self->NewSens($sensRev);
    $self->{"sens_struct"}->{$sensRev->GetLabel()}=$ref_site_sens_rev;
    Alphy::SiteSens::LinkRev($ref_site_sens, $ref_site_sens_rev);
}
sub GetSens {
    my $self=shift;
    my $sens=shift;

    return $self->{"sens_struct"}->{$sens->GetLabel()};
}
sub ProvideSens {
    my $self=shift;
    my $sens=shift;
    if (not $self->HasSensIndex($sens)) {
	$self->AddSens($sens);
    } # creation du hash ref_site_sens et d'une ref dessus
    return $self->GetSens($sens);
}
sub GetSensIndexList {
    my $self=shift;
    return keys(%{$self->{"sens_struct"}});
}
sub GetSensList {
    my $self=shift;
    return values(%{$self->{"sens_struct"}});
}

sub NewSens {
    my $self=shift;
    my $sens=shift;
    
    die "This method needs to be overwriten\n";
}

1;
