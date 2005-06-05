
use ExtUtils::MakeMaker;
$Verbose = 1;
WriteMakefile(
    NAME   => 'Alphy::CUtils',
    SKIP   => [qw(all static static_lib dynamic dynamic_lib)],
    clean  => {'FILES' => 'libalphy-cutils$(LIB_EXT)'},
);

sub MY::top_targets {
        '
all :: static

pure_all :: static

static ::       libalphy-cutils$(LIB_EXT)

libalphy-cutils$(LIB_EXT): $(O_FILES)
	$(AR) cr libalphy-cutils$(LIB_EXT) $(O_FILES)
	$(RANLIB) libalphy-cutils$(LIB_EXT)

';
}

