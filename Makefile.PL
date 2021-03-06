use 5.008004;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'altree',
    PMLIBDIRS         => [ 'ALTree' ],
    DIR               => [ 'CUtils' ],
    VERSION_FROM      => 'progs/altree',
    depend	      => { Makefile => '$(VERSION_FROM)' },
    LICENSE           => 'perl',
    PREREQ_PM         => {}, # e.g., Module::Name => 1.1
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'progs/altree', # retrieve abstract from module
       AUTHOR         => 'Claire Bardel <Claire.Bardel@univ-lyon1.fr>') : ()),
    EXE_FILES         => [qw (progs/altree progs/altree-add-S progs/altree-convert)],
);
