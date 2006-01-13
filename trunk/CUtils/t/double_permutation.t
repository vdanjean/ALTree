# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl CUtils.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More tests => 3;
BEGIN { use_ok('ALTree::CUtils') };

#########################

# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

$res=ALTree::CUtils::double_permutation(2,3, [0, 0.1]);
is($res, undef, "undef when bad args");

$res=ALTree::CUtils::double_permutation(10, 2, 
[12.1, 432.2, 123.3, 15.4, 3453.5, 34253.6, 12.7, 23.8, 23.9, 10, 
11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
);

is_deeply($res, { "pmin" => 0.2,
		  "chi2" => [ 0.8, 0.1] }, "double_permutation");

