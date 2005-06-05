use 5.008004;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'alphy',
    PMLIBDIRS         => [ 'Alphy' ],
    DIR               => [ 'CUtils' ],
    VERSION_FROM      => 'alphy', # finds $VERSION
    PREREQ_PM         => {}, # e.g., Module::Name => 1.1
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'alphy', # retrieve abstract from module
       AUTHOR         => 'Claire Bardel <bardel@vjf.inserm.fr>') : ()),
    EXE_FILES         => [qw (alphy alphy-add-S alphy-convert)],
);
