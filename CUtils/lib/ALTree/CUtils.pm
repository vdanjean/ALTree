package ALTree::CUtils;

use 5.008;
use strict;
use warnings;
use Carp;

require Exporter;
use AutoLoader;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use ALTree::CUtils ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	bilateral
	critchi
	left
	pochisq
	right
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

our $VERSION = '1.1';

sub AUTOLOAD {
    # This AUTOLOAD is used to 'autoload' constants from the constant()
    # XS function.

    my $constname;
    our $AUTOLOAD;
    ($constname = $AUTOLOAD) =~ s/.*:://;
    croak "&ALTree::CUtils::constant not defined" if $constname eq 'constant';
    my ($error, $val) = constant($constname);
    if ($error) { croak $error; }
    {
	no strict 'refs';
	# Fixed between 5.005_53 and 5.005_61
#XXX	if ($] >= 5.00561) {
#XXX	    *$AUTOLOAD = sub () { $val };
#XXX	}
#XXX	else {
	    *$AUTOLOAD = sub { $val };
#XXX	}
    }
    goto &$AUTOLOAD;
}

require XSLoader;
XSLoader::load('ALTree::CUtils', $VERSION);

# Preloaded methods go here.

# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

ALTree::CUtils - Perl extension for blah blah blah

=head1 SYNOPSIS

  use ALTree::CUtils;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for ALTree::CUtils, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.

=head2 Exportable functions

  double bilateral(double a, double b, double c, double d)
  double critchi(double p, int df)
  double left(double a, double b, double c, double d)
  double pochisq(double x, int df)
  double right(double a, double b, double c, double d)
  hash_ref double_permutation(int nb_sample, int nb_chi2, array_ref data)
    data : array of nb_sample*nb_chi2 doubles
           chi2_0_of_sample_0, chi2_1_of_sample_0, ..., chi2_0_of_sample_1, ...
    return keys:
      pmin => double
      chi2 => array of doubles

=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Claire Bardel, E<lt>Claire.Bardel@univ-lyon1.frE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2005 by Claire Bardel

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.4 or,
at your option, any later version of Perl 5 you may have available.


=cut
