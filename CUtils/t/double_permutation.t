# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl CUtils.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More tests => 3;
use Test::Deep;
BEGIN { use_ok('ALTree::CUtils') };

#########################

# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

$res=ALTree::CUtils::DoublePermutation(2,3, [0, 0.1]);
is($res, undef, "undef when bad args");

$res=ALTree::CUtils::DoublePermutation(10, 2, 
[12.1, 432.2, 123.3, 15.4, 3453.5, 34253.6, 12.7, 23.8, 23.9, 10, 
11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
);

$p=0.00001;
cmp_deeply($res, { "pmin" => num(0.2, $p),
		  "chi2" => [ num(0.8, $p), num(0.1, $p)],
		  "distrib_pmin" => 
		      [num(0.1, $p), num(0.1, $p), num(0  , $p), num(0.2, $p),
		       num(0.2, $p), num(0.8, $p), num(0.6, $p), num(0.5, $p),
		       num(0.4, $p), num(0.3, $p)]
		}, "double_permutation");

