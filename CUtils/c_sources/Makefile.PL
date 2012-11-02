
use ExtUtils::MakeMaker;
$Verbose = 1;
WriteMakefile(
    NAME   => 'ALTree::CUtils',
    SKIP   => [qw(all static static_lib dynamic dynamic_lib)],
    clean  => {'FILES' => 'libaltree-cutils$(LIB_EXT)'},
    INC    => '-Wall -I. -Ic_sources',
);

sub MY::top_targets {
        '
all :: static

pure_all :: static

static ::       libaltree-cutils$(LIB_EXT)

libaltree-cutils$(LIB_EXT): $(O_FILES)
	$(AR) cr libaltree-cutils$(LIB_EXT) $(O_FILES)
	$(RANLIB) libaltree-cutils$(LIB_EXT)

';
}

