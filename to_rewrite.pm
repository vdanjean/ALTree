
package main;
use strict;
use diagnostics;
use warnings;

use ALTree::Chi2 ();
use ALTree::Import;
use ALTree::Utils qw(erreur);
use ALTree::Input qw(PrepareTree);
#use Newchi2treeUtils;
use TamuAnova;

sub parcours_nosplit_chi2split
{
    my($tabnodes_a_traiter)=shift;
    my($prolonge)=shift;
    my($splitmode)=shift;
    my($node_ecriture)=shift;
    my($sign_util)=shift; # vaut 1 si on a besoin de la significativité, 0 sinon
    my($node, $child, @tab_noeuds_suivants);
    my($val)=0;
    my($test, $p_val);
    my($test_results);
 
    $test_results->{"ddl"}=scalar(@{$tabnodes_a_traiter})-1; # Nb branches -1
    # $test n'est valable que si $sign_util est à YES
    ($p_val, $test)=CalculChi2($tabnodes_a_traiter, $test_results->{"ddl"},
			       $test_results, $sign_util );
    $test_results->{"node_teste"}=$node_ecriture;
    push (@{$node_ecriture->{"res_test"}}, $test_results);
    $test_results->{"level"}=scalar(@{$node_ecriture->{"res_test"}})-1;
    
    if ($sign_util== SignUtil::YES && $test==1 && $splitmode == SplitMode::CHI2SPLIT) { # sign et que on on est en chi2split
	foreach $node (@{$tabnodes_a_traiter}) {
	    if (NbFils($node) != 0) {
		my @children=$node->GetChildrenList();
		parcours_nosplit_chi2split(\@children, 
					   $prolonge, $splitmode, $node);
	    }
	}
    } elsif ($sign_util== SignUtil::NO || $test==0 || $splitmode == SplitMode::NOSPLIT) { # ou alors on est en nosplit
	foreach $node (@{$tabnodes_a_traiter}) {
	    if (NbFils($node) != 0) {
		$val=1;
		foreach $child ($node->GetChildrenList()) {
		    push (@tab_noeuds_suivants, $child);
		}
	    } else {
		if ($prolonge == 1) {
		    push (@tab_noeuds_suivants, $node);
		}
	    }
	}
	if ($val==1) {
	    parcours_nosplit_chi2split(\@tab_noeuds_suivants, 
				       $prolonge, $splitmode, $node_ecriture, $sign_util);
	} else {
	    return;
	}
    }
}

sub ParcoursQuanti
{
    my($tabnodes_a_traiter)=shift;
    my($prolonge)=shift;
    my($splitmode)=shift;
    my($node_ecriture)=shift;
    my($sign_util)=shift; # vaut 1 si on a besoin de la significativité, 0 sinon
    my($node, $child, @tab_noeuds_suivants);
    my($val)=0;
    my($test, $res_anova);
    my($test_results);
    
#    $test_results->{"ddl"}=scalar(@{$tabnodes_a_traiter})-1; # Nb branches -1
    my @valeurs;
    my @facteurs;
    my $i=0;
#DEBUG    print STDERR "TTTT ", scalar (@{$tabnodes_a_traiter}), "\n";
    foreach $node (@{$tabnodes_a_traiter}) {
	$i++;
	foreach my $case (@{$node->{"quanti"}}) {
	    push (@valeurs, $case->[0]);
	    push (@facteurs, $i);
	}
    }
    my  $nb_factors=$i;
    $test_results->{"nb_facteurs"}=$nb_factors;
# DEBUG    print STDERR "node ";
#    for (my $i=0; $i<=$#valeurs; $i++) {
#	print STDERR " $valeurs[$i]";
#	print STDERR " ($facteurs[$i])";
#    }
#    print STDERR "\n";
	
	
#    if ($sign_util==SignUtil::YES) {
    ($res_anova, $test)=CalculAnovaOneWay($tabnodes_a_traiter, \@valeurs, \@facteurs, $test_results, $sign_util, $nb_factors );
#    } elsif ($sign_util==SignUtil::NO) { 
#	($res_anova)=CalculAnovaOneWay($tabnodes_a_traiter, \@valeurs, \@facteurs, $test_results, $sign_util, $nb_factors);
#   }
    $test_results->{"node_teste"}=$node_ecriture;
    push (@{$node_ecriture->{"res_test"}}, $test_results);
    $test_results->{"level"}=scalar(@{$node_ecriture->{"res_test"}})-1;
    
    if ($sign_util== SignUtil::YES && $test==ALTree::Chi2::SIGNIFICATIF 
	&& $splitmode == SplitMode::CHI2SPLIT) { 
	# sign et que on on est en chi2split
	foreach $node (@{$tabnodes_a_traiter}) {
	    if (NbFils($node) != 0) {
		my @children=$node->GetChildrenList();
		ParcoursQuanti(\@children, 
			       $prolonge, $splitmode, $node);
	    }
	}
    } elsif ($sign_util== SignUtil::NO || $test==ALTree::Chi2::NON_SIGNIFICATIF 
	     || $splitmode == SplitMode::NOSPLIT) {
	# ou alors on est en nosplit
	foreach $node (@{$tabnodes_a_traiter}) {
	    if (NbFils($node) != 0) {
		$val=1;
		foreach $child ($node->GetChildrenList()) {
		    push (@tab_noeuds_suivants, $child);
		}
	    } else {
		if ($prolonge == 1) {
		    push (@tab_noeuds_suivants, $node);
		}
	    }
	}
	if ($val==1) {
	    ParcoursQuanti(\@tab_noeuds_suivants, 
			   $prolonge, $splitmode, $node_ecriture, $sign_util);
	} else {
	    return;
	}
    } else {
	die("Arghhh");
    }
}

sub CalculChi2
{
    my($tabnodes_a_traiter)=shift;
    my($ddl)=shift; 
    my($test_results)=shift;
    my($sign_util)=shift;
    my($chi2, $chi2invalid, $error, $sum_control, $sum_case);
    my($significatif)=ALTree::Chi2::NON_SIGNIFICATIF;
    my($p_value);
    
    ($chi2, $chi2invalid, $error, $sum_control, $sum_case)= 
	ClassicalChi2($tabnodes_a_traiter);
    if ($error != 0) {
	# TODO: A vérifier : est-ce OK de mettre $significatif à 0
	# la valeur est utilisée au retour de cette fonction
#	$significatif=0;
	if ($error == 1) {
	    $test_results->{"texte"}=
		"No cases,  ($sum_control controls)";
	    if ($sign_util==SignUtil::YES) {
		$significatif=ALTree::Chi2::NON_SIGNIFICATIF;
	    }
	} elsif ($error == 2) {
	    $test_results->{"texte"}="No controls: only $sum_case cases";
	    if ($sum_case>=Seuil::ONLY_CASE) {
		 if ($sign_util==SignUtil::YES) {
		     $significatif=1;
		     $test_results->{"sign"}=ALTree::Chi2::SIGNIFICATIF;
		 }
	    } else {
		if ($sign_util==SignUtil::YES) {
		    $significatif=0;
		    $test_results->{"sign"}=ALTree::Chi2::NON_SIGNIFICATIF;
		}	
	    }
	    #$test_results->{"sign"}=ALTree::Chi2::NON_SIGNIFICATIF;
	} elsif ($error == 4) {
	    $test_results->{"texte"}="Only one clade";
	    if ($sign_util==SignUtil::YES) {		
		$significatif=0;
		$test_results->{"sign"}=ALTree::Chi2::NON_SIGNIFICATIF;
	    }
	    # Manque plein de trucs par rapport à la fonction dans chi2tree...
	} else {
	    die "invalid error $error\n";
	}
    } else {
	if ($chi2invalid !=0) {
	    $test_results->{"warning"}="Small sample size correction used";   
	    # J'ai pas compté dans combien de branches...
	    if ($ddl == 1) {
		$p_value=ALTree::CUtils::bilateral($tabnodes_a_traiter->[0]->{"case"},
					   $tabnodes_a_traiter->[0]->{"control"},
					   $tabnodes_a_traiter->[1]->{"case"},
					   $tabnodes_a_traiter->[1]->{"control"});
		if ($sign_util==SignUtil::YES) {
		    $significatif=ALTree::Chi2::chi2_fisher_significatif($p_value);
		}
	    } else {
		my(@clades, $node);
		foreach $node (@{$tabnodes_a_traiter}) {
		    push @clades, ($node->{"case"} + $node->{"control"});
		    # remplit un tableau contenant les effectifs
		    # totaux des différents clades à utiliser dans le
		    # réechantillonnage
		}
		($p_value)= ALTree::Chi2::reech_chi2($sum_case, $sum_control,
					     $ddl+1, $chi2, \@clades);
		$test_results->{"warning"}.=" ($p_value)";
		if ($sign_util==SignUtil::YES) {
		    $significatif= ALTree::Chi2::reech_significatif ($p_value);
		    if ($significatif != 
			ALTree::Chi2::chi2_significatif($ddl, $chi2)) {
			$test_results->{"warning"}.=" Result has changed !";
		    }
		}
	    }
	} else {
	    if ($sign_util==SignUtil::YES) {
		$significatif=ALTree::Chi2::chi2_significatif($ddl, $chi2);
	    }
	    #my $p=`pochisq $chi2 $ddl`+0; # Verif que les 2 appellent 
	                                   #bien la même chose!
	    $p_value=ALTree::CUtils::pochisq($chi2,$ddl);
	    #if ($p != $p_value) {
	    #print STDERR "pochisq: $p != $p_value !\n";
	    #}
	}
	if ($sign_util==SignUtil::YES) {
	    if ($significatif) {
		$test_results->{"sign"}=ALTree::Chi2::SIGNIFICATIF;
	    #$test_results->{"texte"}.="significatif";
	    } else {
		$test_results->{"sign"}=ALTree::Chi2::NON_SIGNIFICATIF;
		# $test_results->{"texte"}.="non significatif";
	    }
	}
	$test_results->{"chi2"}=$chi2;
	$test_results->{"p_val"}=$p_value;
	#$test_results->{"texte"}.=" [p_value_chi2=$p_value]";
    }
#    if ($sign_util==SignUtil::YES) {
	return ($p_value, $significatif);
#    } else {
#	return ($p_value);
#    }
}

sub CalculAnovaOneWay
{
    my $tabnodes_a_traiter=shift;
    my $valeurs=shift;
    my $facteurs=shift;
    my $test_results =shift;
    my $sign_util=shift;
    my $nb_factors=shift;

    my $significatif=ALTree::Chi2::NON_SIGNIFICATIF;
    my $res_anova;
    if (scalar (@{$tabnodes_a_traiter}) < 2) {
	$test_results->{"texte"}=
	    "Only one category";
	if ($sign_util==SignUtil::YES) {
	    $significatif=ALTree::Chi2::NON_SIGNIFICATIF;
	}
    } else {
	if (scalar (@{$valeurs}) != scalar (@{$facteurs})) { 
	    erreur("Error in the anova data: the number of values ", scalar @${valeurs}, " and the number of factors ", scalar @{$facteurs}, " should be the same\n");
	} else {
	    $res_anova=TamuAnova::anova($valeurs, $facteurs,  $nb_factors);
	    #DEBUG   print STDERR $nb_factors, " ", scalar(@{$valeurs}), "\n";
	    $test_results->{"F"}=$res_anova->{"F"};
	    $test_results->{"p_val"}= $res_anova->{"p"};
	    #  $test_results->{"warning"}.=" ($p_value)";
	    if ($sign_util==SignUtil::YES) {
		$significatif = ALTree::Chi2::chi2_fisher_significatif($res_anova->{"p"});
	    }
	    
	}
    }
    if ($sign_util==SignUtil::YES) {
	if ($significatif) {
	    $test_results->{"sign"}=ALTree::Chi2::SIGNIFICATIF;
	    #$test_results->{"texte"}.="significatif";
	} else {
	    $test_results->{"sign"}=ALTree::Chi2::NON_SIGNIFICATIF;
	    # $test_results->{"texte"}.="non significatif";
	}
    }
    return ($res_anova, $significatif);

#DEBUG if (defined ($res_anova)) {
#DEBUG    print STDERR "RESANOVA ", $res_anova->{"F"}, "  ",  $res_anova->{"p"}, "\n";
}

sub InfosAffichees
{
    my($node)=shift;
    my($mode)=shift;
    my($chaine)=Name($node);
    my($lbl_test)=0;
    my $test;
    
    $chaine.=" (LEVEL: ".$node->{"level"}.")";
    if ($mode==1 || $mode == 2) { # Affiche ou pas les case/control
	$chaine.=" case/control:".$node->{"case"}."/".$node->{"control"};
    }
    if ($mode==3 || $mode == 4) {
	$chaine.= sprintf " mean:%.2f",$node->{"moyenne"};
    }
    if (1) { # affiche les apomorphies
	$chaine.="\n";
	foreach my $apo ($node->GetApoList()) {
	    $chaine.= ("  Site: ".$apo->GetSiteNb." Sens: ".$apo->GetSensLabel()."\n");
	}
    }
    $chaine.="\n";
    if ($mode==1 || $mode == 2) { # affiche ou pas les ddl
	if (defined $node->{"res_test"}) {
	    for $test (@{$node->{"res_test"}}) {
		$chaine.= sprintf "[%d] ddl=%d", 
		$test->{"level"}, $test->{"ddl"};
		if ($test->{"ddl"} > 0) {
		    $chaine.= sprintf " chi2=%.2f p_value_chi2=%.3g",
		    $test->{"chi2"}, $test->{"p_val"};
		    # TODO : ça arrive quand on a que des malades ou témoins
		    # dans les clades...
		    if (not defined($test->{"chi2"})) {
			print "chi2 for ", Name($node),
			"(", $test->{"ddl"}, ")", "\n";
		    }
		    if (not defined($test->{"p_val"})) {
			print "p_val for ", Name($node), 
			"(", $test->{"ddl"}, ")", "\n";
		    }
		    
		    if ($mode ==2) {
			if (defined($test->{"sign"})) {
			    if ($test->{"sign"} == ALTree::Chi2::NON_SIGNIFICATIF) {
				$chaine .= " (non significatif)";
			    } elsif ($test->{"sign"} == ALTree::Chi2::SIGNIFICATIF) {
				$chaine .= " (significatif)";
			    } else {
				ALTree::Utils::internal_error("unknown value ".
							      $test->{"sign"});
			      }
			}		    
			if (defined($test->{"texte"})) {
			    $chaine .= "\n".$test->{"texte"};
			}
			if (defined($test->{"warning"})) {
			    $chaine .= "\n".$test->{"warning"};
			}
		    }
		}
		$chaine.="\n";
	    }
	}
    } elsif ($mode == 3 || $mode ==4) {
	if (defined $node->{"res_test"}) {
	    for $test (@{$node->{"res_test"}}) {
		$chaine.= sprintf "[%d] nb_fact=%d", 
		$test->{"level"}, $test->{"nb_facteurs"};
		if ($test->{"nb_facteurs"} > 1) {
		    $chaine.= sprintf " F=%.2f p_value=%.3g",
		    $test->{"F"}, $test->{"p_val"};
		    # TODO : ça arrive quand on a que des malades ou témoins
		    # dans les clades...
		    if (not defined($test->{"F"})) {
			print "F for ", Name($node),
			"(", $test->{"nb_facteurs"}, ")", "\n";
		    }
		    if (not defined($test->{"p_val"})) {
			print "p_val for ", Name($node), 
			"(", $test->{"nb_facteurs"}, ")", "\n";
		    }
		    
		    if ($mode == 4) {
			if (defined($test->{"sign"})) {
			    if ($test->{"sign"} == ALTree::Chi2::NON_SIGNIFICATIF) {
				$chaine .= " (non significatif)";
			    } elsif ($test->{"sign"} == ALTree::Chi2::SIGNIFICATIF) {
				$chaine .= " (significatif)";
			    } else {
				ALTree::Utils::internal_error("unknown value ".
							      $test->{"sign"});
			      }
			}		    
			if (defined($test->{"texte"})) {
			    $chaine .= "\n".$test->{"texte"};
			}
			if (defined($test->{"warning"})) {
			    $chaine .= "\n".$test->{"warning"};
			}
		    }
		}
		$chaine.="\n";
	    }
	}
    }
    return($chaine); 
    
}

1;
