package ALTree::Input;

use strict;
use ALTree::Utils qw(erreur);
use Data::Dumper;
###########################################
########  Read file fonctions     #########
###########################################


#Read the file correspond.txt and put haplotype ID, nb case and nb
#controle in a hash: $correspondance{$HaploID}->{"case"} and
#$correspondance{$HaploID}->{"controle"}

sub ReadCorrespond 
{
    my($name_correspond) =shift;
    my($ligne, @tableau);
    my(%correspondance);
    open (CORRESP, '<', $name_correspond) || die "Unable to open file $name_correspond: $!\n";
    while ($ligne=<CORRESP>) {
	chomp($ligne);
	if ($ligne =~ /^\s*$/) {
	    next;
	}
	# On peut mettre '#' ou ';' pour introduire une ligne de commentaire
	my $di="#";
	if ($ligne =~ /^\s*[$di;]/) {
	    next;
	}
	@tableau=split(/\s+/, $ligne);
	if ($#tableau != 2) {
	    ALTree::Utils::erreur("I do not find 3 columns in $name_correspond at".
		   " line $.:\n$ligne\nPlease, check the syntax.\n",0);
	}

	my $id=shift @tableau;
	if ($id !~ /^[a-zA-Z0-9]+$/) {
	    ALTree::Utils::erreur("The haplotype name '$id' contains unauthorized".
		   " characters\nin $name_correspond at".
		   " line $.:\n$ligne\nPlease, check the syntax.\n",0);
	}
	
	for(my $i=0; $i<2; $i++) {
	    $_=shift @tableau;
	    my $value=0;
	    my $type;
	    if (/^([a-zA-Z]+)_?([0-9]+)$/) {
		$type=$1;
		$value += $2;
	    } else {
		ALTree::Utils::erreur("I do not understand '$_' in $name_correspond at".
		       " line $.:\n$ligne\nPlease, check the syntax.\n", 0);
	    }
	    $_=$type;
	    if (/^(m(alade(s)?)?)|(case(s)?)$/) {
		$type="case";
	    } elsif (/^c(ontrols?)?$/) {
		$type="control";
	    } else {
		ALTree::Utils::erreur("I do not understand '$_' in $name_correspond at".
		       " line $.:\n$ligne\nPlease, check the syntax.\n", 0);
	    }
	    if (defined($correspondance{$id}->{$type})) {
		ALTree::Utils::erreur("For the second time, I read the number of".
		       " ${type}s of '$id'\nin file '$name_correspond' at".
		       " line $.:\n$ligne\nPlease, check the syntax.\n", 0);
	    }
	    $correspondance{$id}->{$type}=$value;
	}
    }
    #my($clefs);
    #DEBUG
    #foreach $clefs (keys %correspondance) {
    #print "$clefs case: ", $correspondance{$clefs}->{"case"}, "\n";
    #	print "$clefs, control: ",$correspondance{$clefs}->{"control"}, "\n";
    #}
    return(\%correspondance);
}

# return :
# Hash of
#   "filename" => String
#   "trees" =>
#     Array (size number of trees) of
#       Hash of
#         "index" => Int
#         "line" => Int
#         "file" => ref of first hash
#         "tab_longbranche" =>
#            Array of (son, father, branch length)
#
#         "tab_seq" =>
#           Hash of
#             son => Array of (father, yes|no|maybe, .....1..0..1..)
#
#         "tab_infoapo" =>
#            Array of (father, son, apo
# number, CI, nb steps, direction of change)
#         "nb_br_non_nulle" =>
#            Int
#         ["outgroup" => String(leaf name)]
#         ["ancestor" => String(sequence)]
#         "tree" => Tree

sub ReadInputFile1
{
    my($input_file)=shift;
    my($phylo_program)=shift;
    my($datatype)=shift;
    my($ancetre)=shift;
    my($identifiant);
    my(%file)=("filename" => $input_file);

    if ($phylo_program == PhylProg::PAUP) {
	$identifiant="Tree number";
    } elsif  ($phylo_program == PhylProg::PHYLIP) { 
	$identifiant="requires a total of";
    } else {
	$identifiant = "TREE # ";
    }
    my $tab_arbres;
    if ($phylo_program == PhylProg::PAUP) {
	$tab_arbres=ReadPAUP($input_file);
    } elsif  ($phylo_program == PhylProg::PHYLIP) { 
    	if ($datatype != DataType::SNP) {
	    erreur("Datatype DNA not yet implemented\n", 0);
	}
	$tab_arbres=ReadPHYLIP($input_file);
    } else {
	$tab_arbres=ReadPAML($input_file);
    }
    $file{"trees"}=$tab_arbres;

    my $index=0;
    for my $arbre (@{$tab_arbres}) {
	$arbre->{"index"}=$index++;
	$arbre->{"file"}=\%file;
	$arbre->{"tree"}=readTreeOld($phylo_program, $arbre, 
				     $datatype, $ancetre);
    }
    $file{"nb_trees"}=$index;
    print STDERR "read done\n";
    return (\%file);
}

sub ReadPAUP
{
    my $file=shift;
    my($ligne);
    my @trees;
    
    open (INPUT, '<', $file) || die "Unable to open file $file: $!\n";
  TREE: 
    {
	my(@tab_longbranche, @tab_infoapo);
	my($nb_br_non_nulle)=0;
	my(%tree);
	
      FIND_TREE:
	{
	    while ($ligne=<INPUT>) {
		last FIND_TREE if ($ligne =~ /^\s+Node\s+to node/);
	    }
	    # Fin du fichier
	    last TREE;
	}
	$tree{"line"}=$.;
	    
      READ_LONGBRANCHE:
	{
	    my $pos = index($ligne, "to node");
	    while ($ligne=<INPUT>) {
		if ($ligne =~ /-------------------------/) {
		    next;
		}		
		if ($ligne =~ /^\s*Sum/) {
		    last READ_LONGBRANCHE;
		}
		chomp($ligne);
		my $son = substr($ligne,0,$pos-1);
		#   print STDERR "son1:$son\n";
		my $other_infos = substr($ligne,$pos);
		$other_infos =~ s/^\s*//;
		#   print STDERR "otherinfos:$other_infos\n";
		if ($son =~ /^\s*([0-9]+)\s*$/) {
		    $son = $1;
		} elsif ($son =~ /\s*([^\s](.*[^\s])?)\s+\([0-9]*\)(\*)?\s*$/) {
		    $son = $1;
		    if (defined($3)) {
			if (defined($tree{"outgroup"})) {
			    erreur("I found a second outgroup '$son'\n".
				   "in file '$file' at line $.\n".
				   "(the first outgroup was".
				   " '$tree{outgroup}')\n", 0);
			}
			$tree{"outgroup"} = $son;
		    }
		} else {
		    erreur("Sorry, I am unable to understand:\n$ligne\n".
			   "in file '$file' at line $.\n".
			   "(bad branch '$son' while reading branch".
			   " lengths)\n", 0);
		}
		my @tableau = split(/\s+/,  $other_infos);
		if ($#tableau != 3) {
		    erreur("Sorry, I am unable to understand:\n$ligne\n".
			   "in file '$file' at line $.\n".
			   "(bad number of columns while reading branch".
			   " lengths)\n", 0);
		}
		# We add (son, father, branch length)
		push @tab_longbranche, [$tableau[0], $son, $tableau[1]];
		if ($tableau[2] != 0) {
		    $nb_br_non_nulle++;
		}
	    }
	    next TREE;
	}
      FIND_APO:
	{
	    while ($ligne=<INPUT>) {
		if ($ligne =~ /^\s+Branch\s+Character\s+Steps\s+CI\s+Change/) {
		    last FIND_APO;
		}
	    }
	    next TREE;
	}
      READ_APO:
	{
	    my $position = index($ligne, "Character");
	    my $son;
	    my $father;

	    while ($ligne=<INPUT>) {
		chomp($ligne);
		if ($ligne =~ /^\s*$/) {
		    last READ_APO;
		}
		if ($ligne =~ /^\s*------------------/) {
		    next;
		}
		my $name_haplo = substr($ligne,0,$position-1);
		if ($name_haplo =~ /^\s*$/) {
		    if (not defined($son)) {
			erreur("Sorry, I am unable to understand:\n$ligne\n".
			       "in file '$file' at line $.\n".
			       "(no branch while reading".
			       " apomorphies)\n", 0);
		    }
		} elsif ($name_haplo =~ /^\s*(root|node_[0-9]+) --> (.*[^\s])\s*$/){
		    $father=$1;
		    $son=$2;
		    $father=~ s/^node_//;
		    $son=~ s/^node_//;
		} else {
			erreur("Sorry, I am unable to understand:\n$ligne\n".
			       "in file '$file' at line $.\n".
			       "(bad branch while reading".
			       " apomorphies)\n", 0);
		}


		my $infos = substr($ligne,$position);
		#print STDERR "infos:$infos\n";
		my @tab_infos=split(/\s+/, $infos);
		if ($#tab_infos != 5) {
		    erreur("Sorry, I am unable to understand:\n$ligne\n".
			   "in file '$file' at line $.\n".
			   "(bad number of columns while reading".
			   " apomorphies)\n", 0);
		}		

		my $sens_chgt = join(' ', $tab_infos[3], "-->", $tab_infos[5]);

		# (father, son, apo number, CI, nb steps, direction of change)
		push @tab_infoapo, [$father, $son, $tab_infos[0],
				$tab_infos[2], $tab_infos[1],
				ALTree::Sens->New($sens_chgt)];

	    }
	}
	$tree{"tab_longbranche"} = \@tab_longbranche;
	$tree{"tab_infoapo"} = \@tab_infoapo;
	$tree{"nb_br_non_nulle"} = $nb_br_non_nulle;
	push @trees, \%tree;
	redo TREE;
    } continue {
	erreur("Unexpected end of file while reading a tree\n".
	       "in file '$file' at line $.\n", 0);
    }
    close(INPUT);
    #print STDERR Dumper(@trees);
    return \@trees;
}

sub ReadPHYLIP
{
    my $file=shift;
    my($ligne);
    my @trees;
    
    open (INPUT, '<', $file) || die "Unable to open file $file: $!\n";
  TREE: 
    {
	my(@tab_longbranche);
	my(%tab_seq);
	my(%tree);
	
      FIND_TREE:
	{
	  FIND_ANCESTOR:
	    {
		while ($ligne=<INPUT>) {
		    last FIND_TREE if ($ligne =~ /^From\s+To\s+Any\s+Steps/);
		    last FIND_ANCESTOR if ($ligne =~ /^best guesses of ancestral states:/);
		    }
		# Fin du fichier
		last TREE;
	    }
	    my $ancestor="";
	  READ_ANCESTOR:
	    {
		while ($ligne=<INPUT>) {
		    chomp($ligne);
		    if ($ligne =~ /.*0\!\s+([01?\s]+)$/) {
			my @tab=split(/\s+/, $1);
			$ancestor.=join('', @tab);
		    } elsif ($ligne =~ /^[\s]*$/) {
			last READ_ANCESTOR;
		    } elsif ($ligne =~ /^[0-9\s]*$/) {
		    } elsif ($ligne =~ /^\s*\*-+\s*$/) {
		    } else {
			erreur("Sorry, I am unable to understand:\n$ligne\n".
			       "in file '$file' at line $.\n".
			       "(while reading ancestor)\n", 0);			
		    }
		}
		next TREE;
	    }
	    $tree{"ancestor"}=$ancestor;
	    #print STDERR "ancestor=$ancestor\n";
	    while ($ligne=<INPUT>) {
		last FIND_TREE if ($ligne =~ /^From\s+To\s+Any\s+Steps/);
	    }
	    next TREE;
	}
	$tree{"line"}=$.;
	    
      READ_LONGBRANCHE:
	{
	    while ($ligne=<INPUT>) {
		chomp($ligne);

		if ($ligne =~ /^\s*([0-9a-zA-Z_]+)\s+([0-9a-zA-Z_]+)\s+(yes|no|maybe)\s+([01.? ]+)\s*$/) {
		    #print STDERR "trouvé! $ligne\n";
		    my($father)=$1;
		    my($son)=$2;
		    my($step2)=$3;
		    my $seq=join('', split(/\s+/, $4));
		    my $step=$seq;
		    $step =~ s/\.//g;
		    $step=length($step);
		    push @tab_longbranche, [$father, $son, $step];
		    if (defined($tab_seq{$son})) {
			erreur("It is the second time I read an information".
			       "for the branch leading to '$son':\n$ligne\n".
			       "in file '$file' at line $.\n".
			       "Please, check your data.\n", 0);
		    }
		    $tab_seq{$son}=[$father, $step2, $seq];
		} else {
		    if (scalar(@tab_longbranche)==0) {
			next;
		    } elsif ($ligne =~ /^\s*$/) {
			last READ_LONGBRANCHE;
		    } else {
			erreur("Sorry, I am unable to understand:\n$ligne\n".
			       "in file '$file' at line $.\n".
			       "(while reading branches)\n", 0);			
		    }
		}
	    }
	    if (scalar(@tab_longbranche)==0) {
		next TREE;
	    }
	}
	$tree{"tab_longbranche"} = \@tab_longbranche;
	$tree{"tab_seq"} = \%tab_seq;
	push @trees, \%tree;
	redo TREE;
    } continue {
	erreur("Unexpected end of file while reading a tree\n".
	       "in file '$file' at line $.\n", 0);
    }
    close(INPUT);
    #print STDERR Dumper(@trees);
    return \@trees;
}

sub ReadPAML 
{
    my $file=shift;
    my($ligne);
    my @trees;
    
    open (INPUT, '<', $file) || die "Unable to open file $file: $!\n";
  TREE: 
    {
	my(@tab_longbranche, @tab_infoapo);
	my($nb_br_non_nulle)=0;
	my(%tree);
	
      FIND_TREE:
	{
	    while ($ligne=<INPUT>) {
		last FIND_TREE if ($ligne =~ /^\s*Summary of changes along branches/);
	    }
	    # Fin du fichier
	    last TREE;
	}
	$tree{"line"}=$.;	
      READ_LONGBRANCHE:
	{
	  READ_BRANCH:
	    while ($ligne=<INPUT>) {
		chomp($ligne);
		if ($ligne =~ /^\s*List of extant and reconstructed sequences/) {
		    last READ_LONGBRANCHE;
		}
		if ($ligne =~ /^\s*$/) {
		    next;
		}
		if ($ligne =~ /^Check root for directions of change./) {
		    next;
		}
		if (($ligne =~ /\s*Branch [0-9]+:\s+([0-9]+)\.\.([0-9]+)\s*$/) || 
		    ($ligne =~ /\s*Branch [0-9]+:\s+([0-9]+)\.\.[0-9]+\s+[(](.*)[)]\s*$/)) {
		    my $father=$1;
		    my $son=$2;
		    my $br_len=0;
		  READ_APO:
		    {
			while ($ligne=<INPUT>) {
			    if ($ligne =~ /^\s*$/) {
				next;
			    }
			    if ($ligne =~ /\s*([0-9]+)\s+([0-9A-Za-z?_-])\s+[0-9.]+\s*->\s*([0-9A-Za-z?_-])\s*([0-9\.]*)?\s*/) {
				# A VERIFIER: NORMALEMENT NE DOIT SE PRODUIRE QU'EN BOUT DE BRANCHES!
				if ($3 eq "?") {
				    next;
				}
				#father
				#son
				# apomorphie number
				# CI
				# nb steps =1 for SNPs 
				# direction of the change 
				$br_len++;
				push @tab_infoapo, [$father, $son, $1, undef,
						    1, ALTree::Sens
						    ->New($2."->".$3)];				
			    } else {
				last READ_APO;
			    }
			}
			next TREE;
		    }

		    # We add (son, father, branch length)
		    push @tab_longbranche, [$father, $son, $br_len];
		    if ($br_len != 0) {
			$nb_br_non_nulle++;
		    }
		    # On boucle : on a déjà trouvé une autre ligne
		    redo READ_BRANCH;
		} else {
		    erreur("Sorry, I am unable to understand:\n$ligne\n".
			   "in file '$file' at line $.\n", 0);
		}
	    }
	    next TREE;
	}
	$tree{"tab_longbranche"} = \@tab_longbranche;
	$tree{"tab_infoapo"} = \@tab_infoapo;
	$tree{"nb_br_non_nulle"} = $nb_br_non_nulle;
	push @trees, \%tree;
	redo TREE;
    } continue {
	erreur("Unexpected end of file while reading a tree\n".
	       "in file '$file' at line $.\n", 0);
    }
    close(INPUT);
    #print STDERR Dumper(@trees);
    return \@trees;
}

sub readTree {
    my $phylo_program=shift;
    my $input_file=shift;
    my $datatype=shift;
    my $ancetre=shift;
    my $number_arbre=shift;

    return $input_file->{"trees"}->[$number_arbre]->{"tree"};
}

#use Data::Dumper;
sub readTreeOld {
    my $phylo_program=shift;
    my $tab_arbre=shift;
    my $datatype=shift;
    my $ancetre=shift;
    my $tree;

    my $tab_longbranche=$tab_arbre->{"tab_longbranche"};
    $tree=TreeBuilding($tab_longbranche);
    
    #print STDERR Dumper($tab_arbre);
    #print STDERR Dumper($tree);
    if ($phylo_program == PhylProg::PHYLIP) {
	my $racine=$tree->GetRoot();
	my $ancetre_seq;
	if (not defined($tab_arbre->{"ancestor"})) {
	    if ($ancetre) {
		$ancetre_seq=$ancetre;
	    } else {
		erreur ("You have forgotten the option --anc-seq!\n", 0);
	    }
	} else {
	    $ancetre_seq=$tab_arbre->{"ancestor"};
	    if ($ancetre && $ancetre ne $ancetre_seq) {
		print STDERR "ancestor defined twice, the sequence entered with the -anc-seq option will be ignored\n";
	    }
	}
	BuildInfoApoPhylip($tab_arbre, $ancetre_seq, $racine);
	#FillTreeApo1Phylip($tree, $tab_longbranche, $ancetre_seq);
	#FillTreeApo2Phylip($racine, $racine, $ancetre_seq, $tree);
    }elsif ($phylo_program == PhylProg::PAUP ||
	$phylo_program == PhylProg::PAML) {
	
    } else {
	ALTree::Utils::internal_error("Phylogeny program not supported");
    }
    
    my $tab_infoapo=$tab_arbre->{"tab_infoapo"};
    FillTreeApoInfo($tree, $tab_infoapo);

    if (not defined($tab_arbre->{"has_ambiguity"})) {
	# CheckApoBrlen($tree); Pb avec phylip..
	my $nb_br_non_nulle=$tab_arbre->{"nb_br_non_nulle"};
	$tree->SetNbBrNonNulle($nb_br_non_nulle);

	CheckApoBrlen($tree);
    }

    return $tree;
}

###########################################
#########  BUILDING OF THE TREE  ##########
###########################################


# Build the tree and add the branch length info
sub TreeBuilding{
    my($tab_longbranche)=shift;
    my($i);
    my($tree)=ALTree::Tree->New();

    #print "TreeBuilding\n";
    for ($i=0;$i<=$#$tab_longbranche;$i++) {
	my($pere_id, $fils_id, $long_br); # variables intermédiaires pour 
	# lisibilite du prog
	    
	$pere_id=$tab_longbranche->[$i]->[0];
	$fils_id=$tab_longbranche->[$i]->[1];
	$long_br=$tab_longbranche->[$i]->[2];

	my $pere;
	my $fils;
	if (not $tree->HasNodeIndex($fils_id)) {
	    $fils=ALTree::Node->New($fils_id); #creation de la structure Node
	    $tree->AddNode($fils);     # puis, on la met dans tree
	} else {
	    $fils=$tree->GetNode($fils_id);
	}
	if (not $tree->HasNodeIndex($pere_id)) {
	    $pere=ALTree::Node->New($pere_id); #creation de la structure Node
	    $tree->AddNode($pere);     # puis, on la met dans tree
	} else {
	    $pere=$tree->GetNode($pere_id);
	}
	$pere->AddChild($fils);
#	print "arbre{pere}->{children}->[0]->{id} ", $arbre->{$pere}->{"children"}->[0]->{"id"}, "\n"; 
	if ($fils->HasFather()) {
	    die ($fils->Name()." already have a father: ".$fils->GetFather()->Name().
		 " and I would like ".$pere->Name()." as father\n");
	} else {
	    $fils->SetFather($pere);
	}
	if ($fils->HasBrLen()) {
	    die($fils->Name()." already have a branch length: ".$fils->GetBrLen()."\n");
	} else {
	    $fils->SetBrLen($long_br);
	}
    }
    return($tree);
}

sub FillTreeApoInfo
{
    my($tree)=shift;
    my($tab_infoapo)=shift;
    my($i);
 
    #print "FillTreeApoInfo\n";
    for ($i=0;$i<=$#$tab_infoapo;$i++) {
	my($pere_id, $fils_id, $apo_num, $apo_CI, $apo_steps, $apo_sens)
	    = @{$tab_infoapo->[$i]};
	#print "branche pere($pere) -> fils($fils)\n";
	if (not $tree->HasNodeIndex($pere_id) or
	    not $tree->HasNodeIndex($fils_id)) {
	    my $id;
	    if (not $tree->HasNodeIndex($pere_id)) {
		$id=$pere_id;
	    } else {
		$id=$fils_id;
	    }
	    erreur("For the apomorphy $apo_num ('$pere_id' -> '$fils_id'),".
		   " the node name\n'$id' has not been found in the tree\n".
		   "Please, check your data.", 0);
	}
	my $pere=$tree->GetNode($pere_id);
	#print $pere->Name();
	my $fils=$tree->GetNode($fils_id);
	#print STDERR "fils =  $fils_id, pere=$pere_id\n";
	if ($fils->GetFather() != $pere) {
	    erreur("I read the apomorphy $apo_num ('$pere_id' -> '$fils_id').\n".
		   "However, the node '$fils_id' has '".
		   $fils->GetFather()->Name()."' as father in the tree.\n".
		   "Please, check your data.", 0);
	}
	
	my $site;
	if (not $tree->HasSiteIndex($apo_num)) { 
	    $site=ALTree::SitePerTree->New($apo_num);
	    $tree->AddSite($site);
	    if (defined $apo_CI) { # For PAML trees, CI is not defined
		$site->SetCI($apo_CI);
	    }
	} else {
	    $site=$tree->GetSite($apo_num);
	}

	my($ref_site_sens)=$site->ProvideSens($apo_sens);
	$ref_site_sens->SetStep($apo_steps); # For PAML, $aposteps=1

	$fils->AddApo($ref_site_sens); # lie arbre et hash_site_sens
    }
}

sub BuildInfoApoPhylip
{
    my($tab_arbre)=shift;
    my($ancetre_seq)=shift;
    my($racine)=shift;
    
    my @tab_infoapo;
    my $nb_br_non_nulle=0;
    my $has_ambiguity=0;
    my($tab_seq)=$tab_arbre->{"tab_seq"};

    my $set_seq;

    my $root=$racine->GetId();

    $set_seq=sub {
	my $son=shift;

	print STDERR "managing $son\n";
	my $branch=$tab_seq->{$son};
	if (scalar(@{$branch})==4) {
	    return;
	}

	my $father=$branch->[0];
	my $fatherseq;
	
	if (not defined($tab_seq->{$father})) {
	    if ($father ne $root) {
		erreur("The node '$son' seems to be a son of '".
		       $father."'\nbut I cannot get info about this branch ".
		       "in the tree number ".$tab_arbre->{"index"}.
		       "\nthat begins at line ".$tab_arbre->{"line"}.
		       " in file '".$tab_arbre->{"file"}->{"filename"}.
		       "'\nPlease, check your data.\n", 0);
	    }
	    $fatherseq=$ancetre_seq;
	} else {
	    if (scalar(@{$tab_seq->{$father}}) == 3) {
		$set_seq->($father);
	    }
	    $fatherseq=$tab_seq->{$father}->[3];
	}
	my $childseq=$branch->[2];
	my $longueur=length($childseq);

	if (length($fatherseq) != $longueur) {
	    erreur("Reconstructed sequences have not the same length".
		   " at node '".$father."'\n($fatherseq) and at node '".
		   $son."' ($childseq)\n", 0);
	}
	my($br_len)=0;
	my($newchildseq)="";
	for (my $i=0; $i<$longueur;$i++) {
	    my($fathersite)= substr($fatherseq,$i,1);
	    my($childsite)= substr($childseq, $i,1);
	    if ($childsite eq ".") {
		$childsite=$fathersite;
		
	    } elsif ($fathersite ne $childsite) {
		if ($fathersite eq '?' || $childsite eq '?') {
		    $has_ambiguity++;
		} 
		$br_len++;
		my($apo_num)=$i+1;
		my($apo_sens)=ALTree::Sens->New($fathersite."->".$childsite);
		push @tab_infoapo, [$father, $son, $apo_num,
				    undef, # Pas de CI dans les données avec PHYLIP
				    1, # Pas de gestion des micro sat.
				    $apo_sens];
	    }
	    $newchildseq.=$childsite;
	}
	if ($br_len != 0) {
	    $nb_br_non_nulle++;
	}
	push @{$branch}, $newchildseq;
    };

    for my $son (keys %{$tab_seq}) {
	$set_seq->($son);
    }

    $tab_arbre->{"tab_infoapo"}=\@tab_infoapo;
    $tab_arbre->{"nb_br_non_nulle"}=$nb_br_non_nulle;
    if ($has_ambiguity) {
	$tab_arbre->{"has_ambiguity"}=$has_ambiguity;
    }
    #delete($tab_arbre->{"tab_seq"});
}


sub CheckApoBrlen
{
    my($tree)=shift;
    my($node);

    $tree->GetRoot()->SetBrLen(0); # to prevent "uninitialized value"

    foreach $node ($tree->GetNodesList()) {
	if ($node->GetBrLen() != $node->NbApoStep()) { 
            # check if nb_apo correspond to br_len
	    die "Error in the tree: branch length= ", $node->GetBrLen(),
	    " but ", $node->NbApoStep(),
	    " apomorphies are defined for node ", $node->Name(), "\n";
	}
    }	
}

1;
