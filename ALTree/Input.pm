package ALTree::Input;

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
	    if (exists($correspondance{$id}->{$type})) {
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

sub ReadInputFile1
{
    my($input_file)=shift;
    my($phylo_program)=shift;
    my($datatype)=shift;
    my($ancetre)=shift;
    my($ligne);
    my($identifiant);
    
    if ($phylo_program == PhylProg::PAUP) {
	$identifiant="Tree number";
    } elsif  ($phylo_program == PhylProg::PHYLIP) { 
	$identifiant="requires a total of";
    } else {
	$identifiant = "TREE # ";
    }
    open (INPUT, '<', $input_file) || die "Unable to open file $input_file: $!\n";
    my($indice)=0;
    my(@tab_arbres)=(); # contain ref on a tab containing lines for 1 tree 
    while ($ligne=<INPUT>) {
	chomp($ligne);
	if ($ligne=~/$identifiant\s*[0-9.]*/) {
	    $indice++;
	}
	push (@{$tab_arbres[$indice]}, $ligne);
    }
    for my $arbre (@tab_arbres) {
	if ($arbre != $tab_arbres[0]) {
	    $arbre=readTreeOld($phylo_program, $arbre, 
			       $datatype, $ancetre);
	}
    }
    shift @tab_arbres;
    return (\@tab_arbres);
}


sub ReadPAUP
{
    my($tab_arbre)=shift;
    my($ligne, @tableau);
    my($marqueur)=0;
    my(@tab_longbranche, @tab_infoapo);
    my($i)=0;
    my($j)=0;
    my($nb_br_non_nulle)=0;
    my ($position)=0;
    my($pos)=0;
    foreach $ligne (@{$tab_arbre}) {
	
	if ($ligne =~ /^\s+Node\s+to node/) {
	    $marqueur=1;
	}
	if ($ligne =~ /^\s*Sum/) {
	    # print "marqueur=$marqueur\n";
	    $marqueur=0;
	}
	if ($ligne =~ /^Apomorphy lists:/) {
	    $marqueur=2;
	    next;
	} 
	
	if ($marqueur==1) {
	   # $ligne =~ s/^\s+//;
	    #if ($ligne =~ /root/) {
	#	next;
	 #   }
	    if ($ligne =~ /---------------------------/) {
		next;
	    }
	    if ($ligne =~ /Node\s+to node/) {
		$pos = index($ligne, "to node"); 
		#print STDERR "pos = $pos\n";
		next;
	    }	
	    my $son = substr($ligne,0,$pos-1);
	 #   print STDERR "son1:$son\n";
	    my $other_infos = substr($ligne,$pos);
	    $other_infos =~ s/^\s*//;
	 #   print STDERR "otherinfos:$other_infos\n";
	    $son =~ /\s*([^\(\)]*\S)\s\(?[0-9]*\)?\s*$/; #Revoir avec Vince
	    $tab_longbranche[$i]->[1]= $1; #son
	    my @tableau = split(/\s+/,  $other_infos);
	    if ($#tableau != 3) {
		die "error: not 4 columns, $#tableau columns\n";
	    } else { 
		$tab_longbranche[$i]->[0]= $tableau[0]; #father
		$tab_longbranche[$i]->[2]= $tableau[1]; #branch length
		if ($tableau[2] != 0) {
		    $nb_br_non_nulle++;
		}
	    }
	    $i++;
	}

	my($modificateur)=0;
	my($sens_chgt);
	
	if ($marqueur==2) {
	    my $has_sens=0;
	    if ($ligne =~ /^\s*$/ || $ligne =~ /^\s*------------------/) {
		next;
	    } 
	    
	    if ($ligne =~ /^\s+Branch/) {
		$position = index($ligne, "Character");
	#	print STDERR "position = $position\n";
		next;
	    }
	    #print STDERR $ligne, "\n";
	    my $name_haplo = substr($ligne,0,$position-1);
	    my $infos = substr($ligne,$position);
	    #print STDERR "infos:$infos\n";
	    my @tab_infos=split(/\s+/, $infos);
	    $tab_infoapo[$j]->[2]=$tab_infos[0]; # apomorphie number
	    $tab_infoapo[$j]->[3]=$tab_infos[2]; # CI
	    $tab_infoapo[$j]->[4]=$tab_infos[1]; # nb steps 
	    $sens_chgt = join(' ', $tab_infos[3], $tab_infos[4], $tab_infos[5]);
	    if ($name_haplo =~ /^\s*((root|node_).*\S)\s*$/) { 
		$name_haplo= $1;
		my @tab_name = split(/ --> /, $name_haplo);
	#	print STDERR "tab0:$tab_name[0]\n";
		$tab_name[0]=~ s/node_//g;
	#	print STDERR "2tab0:$tab_name[0]\n";
		$tab_name[1]=~ s/node_//g;
		$tab_infoapo[$j]->[0]=$tab_name[0]; # father
		#print STDERR "tab1:$tab_name[1]\n";
	#	$tab_name[1] =~ s/\s/_/g;	
		#print STDERR "2tab1:$tab_name[1]\n";
		$tab_infoapo[$j]->[1]=$tab_name[1]; # son
	    } elsif ($name_haplo =~ /^\s*$/) {
		$tab_infoapo[$j]->[0] = $tab_infoapo[$j-1]->[0]; # father
		$tab_infoapo[$j]->[1] = $tab_infoapo[$j-1]->[1]; #son	
	    } else {
		die "Unknown line $name_haplo\n";
	    }
	    #print STDERR $tab_infoapo[$j]->[0], ":father ", $tab_infoapo[$j]->[1], ":son  ",$tab_infoapo[$j]->[2]=$tab_infos[0], ":aponum ",  $tab_infoapo[$j]->[3]=$tab_infos[2], ":CI ", $tab_infoapo[$j]->[4]=$tab_infos[1], ":nbstep ", $sens_chgt, ":sens\n"; 
	    
	    $has_sens=1;
	    
	    if ($has_sens) {
		$sens_chgt =~ s/=/-/g;
		$tab_infoapo[$j]->[5]=ALTree::Sens->New($sens_chgt); # direction of the change 
		$j++;
	    }
	}
	
    }
    # print "nb_non_nul=$nb_br_non_nulle\n";
    return (\@tab_longbranche, \@tab_infoapo, $nb_br_non_nulle);
}

sub ReadPHYLIP
{
    my($tab_arbre)=shift;
    my($num_arbre)=shift;
    my($data_type)=shift;
    my $ancetre =shift;
    my($ligne, @tableau);
    my($marqueur)=0;
    my(@tab, $tabinfo, $j);
    
    foreach $ligne (@{$tab_arbre}) {
	$ligne =~ s/^\s+//;
    	if ($data_type == DataType::SNP) {
	    #print STDERR "line: $ligne\n";
	    if ($ligne =~ /best guesses of ancestral states:/) {
		$marqueur=1;
		if ($ancetre ne "") {
		    print STDERR "ancestor defined twice, the sequence entered with the -anc-seq option will be ignored\n";
		    $ancetre="";
		}
	    }
	    if ($ligne =~ /^From\s+To\s+Any\s+Steps/) {
		$marqueur=2;
		$j=0;
	    }
	    if ($marqueur==1) {
		if ($ligne =~ /.*0\!\s+[01?\s]+/) {
		    
		    @tab=split(/\s+/, $ligne);
		    shift(@tab);
		    my($i);
		    for ($i=0; $i<=$#tab; $i++) {
			$ancetre.=$tab[$i];
		    }
		}
		#print STDERR "ancetre=$ancetre ($ligne)\n";
	    }
	    if ($marqueur==2) {
	        #print STDERR "Reading line $ligne\n";
		#if ($ligne =~ /root/) {
		#    next;
		if ($ligne =~ /^\s*([0-9a-zA-Z_]+)\s+[0-9a-zA-Z_]+\s+(yes|no|maybe)\s+[01.? ]+\s*$/) {
		    #print STDERR "trouvé! $ligne\n";
		    @tab=split(/\s+/, $ligne);
		    $tabinfo->[$j]->[0]=shift(@tab);
		    $tabinfo->[$j]->[1]=shift(@tab);
		    $tabinfo->[$j]->[2]=shift(@tab);
		    $tabinfo->[$j]->[3]=join('',@tab);
		    $j++;
		}
	    }
	    
	   
	} else { # $data_type == DataType::DNA
	    erreur("Datatype DNA not yet implemented\n", 0);
	    #A faire
	}
    }
    if ($ancetre eq "") {
	erreur ("You have forgotten the option --anc-seq!\n", 0);
    }
    return ($tabinfo, $ancetre);
}

sub ReadPAML 
{
    my($tab_arbre)=shift;
    my($ligne);
    my($i)=-1;
    my($j)=0;
    my($nb_br_non_nulle)=0;
    my($marqueur)=0;
    my(@tab_longbranche, @tab_infoapo);
    my($has_mutation)=0;
    foreach $ligne (@{$tab_arbre}) {
	chomp($ligne);
	if ($ligne =~ /^\s*Summary of changes along branches/) {
	    $marqueur=1;
	    next;
	}
	if ($ligne =~ /^\s*List of extant and reconstructed sequences/) {
	    if ($has_mutation != 0) {
		$nb_br_non_nulle++;
	    }
	    $tab_longbranche[$i]->[2]=$has_mutation; 
	    $marqueur=0;
	    last;
	}

	if ($marqueur==1) {
	    $ligne =~ s/^\s+//;
	    if (($ligne =~ /\s*Branch [0-9]+:\s+([0-9]+)\.\.([0-9]+)\s*$/) || 
		($ligne =~ /\s*Branch [0-9]+:\s+([0-9]+)\.\.[0-9]+\s+[(](.*)[)]\s*$/)) {
		# Début de branche
		#print "Cas 1 ($1 $2)\n";
		$i++;
		if ($has_mutation != 0) {
		    $nb_br_non_nulle++;
		}
		if ($i>0) {
		    $tab_longbranche[$i-1]->[2]=$has_mutation; #branch length for preceeding branch
		}
		$has_mutation=0;
		$tab_longbranche[$i]->[0]= $1; #father
		$tab_longbranche[$i]->[1]= $2; #son
	    #} elsif ($ligne =~ /\s*([0-9]+)\s+([0-9A-Za-z?_-])\s+[0-9.]+\s*->/) {
	    } elsif ($ligne =~ /\s*([0-9]+)\s+([0-9A-Za-z?_-])\s+[0-9.]+\s*->\s*([0-9A-Za-z?_-])\s*([0-9\.]*)?\s*/) {
		
		# A VERIFIER: NORMALEMENT NE DOIT SE PRODUIRE QU'EN BOUT DE BRANCHES!
		if ($3 eq "?") {
		    next;
		}
		$has_mutation++;
		$tab_infoapo[$j]->[0] = $tab_longbranche[$i]->[0]; #father
		$tab_infoapo[$j]->[1] = $tab_longbranche[$i]->[1]; #son
      		$tab_infoapo[$j]->[2]=$1; # apomorphie number
      		#$tab_infoapo[$j]->[3]=; # CI
		$tab_infoapo[$j]->[4]=1; # nb steps =1 for SNPs 
		my($sens_chgt)=$2."->".$3;
		$tab_infoapo[$j]->[5]=ALTree::Sens->New($sens_chgt); # direction of the change 
		$j++;
	    }
	}
    }
    return (\@tab_longbranche, \@tab_infoapo, $nb_br_non_nulle);
}

sub readTree {
    my $phylo_program=shift;
    my $tab_arbres=shift;
    my $datatype=shift;
    my $ancetre=shift;
    my $number_arbre=shift;

    my $tree;
    my($tab_longbranche, $tab_infoapo, $ancetre_seq);

    my($nb_br_non_nulle);
 
    return $tab_arbres->[$number_arbre]->{"tree"};
}

#use Data::Dumper;
sub readTreeOld {
    my $phylo_program=shift;
    my $tab_arbre=shift;
    my $datatype=shift;
    my $ancetre=shift;
    my $tree;
    my($tab_longbranche, $tab_infoapo, $ancetre_seq);

    my($nb_br_non_nulle);
    if ($phylo_program == PhylProg::PAUP) {
	($tab_longbranche, $tab_infoapo, $nb_br_non_nulle)
	    =ReadPAUP($tab_arbre);
    } elsif ($phylo_program == PhylProg::PHYLIP) {
	($tab_longbranche, $ancetre_seq)
	    =ReadPHYLIP($tab_arbre, $datatype, $ancetre);
    } elsif ($phylo_program == PhylProg::PAML) {
	($tab_longbranche, $tab_infoapo, $nb_br_non_nulle)
	    =ReadPAML($tab_arbre);
    } else {
	ALTree::Utils::internal_error("Phylogeny program not supported");
    }
    $tree=TreeBuilding($tab_longbranche);
    
    #print STDERR Dumper($tab_arbre);
    #print STDERR Dumper($tree);
    if ($phylo_program == PhylProg::PAUP ||
	$phylo_program == PhylProg::PAML) {
	$tree->SetNbBrNonNulle($nb_br_non_nulle);
	FillTreeApoInfoPAUP($tree, $tab_infoapo);
	CheckApoBrlen($tree);
    } elsif ($phylo_program == PhylProg::PHYLIP) {
	my $racine=$tree->GetRoot();
	FillTreeApo1Phylip($tree, $tab_longbranche, $ancetre_seq);
	FillTreeApo2Phylip($racine, $racine, $ancetre_seq, $tree);
    } else {
	ALTree::Utils::internal_error("Phylogeny program not supported");
    }
   # CheckApoBrlen($tree); Pb avec phylip..
    
    my %infos=("tree" => $tree);
    return \%infos;
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
	#print "branche pere($pere) -> fils($fils)\n";
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
	    die ($fils->Name()." already have a father: ".$fils->GetFather()->Name()."\n");
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

sub FillTreeApoInfoPAUP
{
    my($tree)=shift;
    my($tab_infoapo)=shift;
    my($i);
 
    #print "FillTreeApoInfo\n";
    for ($i=0;$i<=$#$tab_infoapo;$i++) {
	my($pere_id, $fils_id, $apo_num, $apo_CI, $apo_steps, $apo_sens)
	    = @{$tab_infoapo->[$i]};
	#print "branche pere($pere) -> fils($fils)\n";
	if (not $tree->HasNodeIndex($pere_id) || 
	    not $tree->HasNodeIndex($fils_id)) {
	    die "unknown node: $pere_id or $fils_id\n";
	}
	my $pere=$tree->GetNode($pere_id);
	#print $pere->Name();
	my $fils=$tree->GetNode($fils_id);
	#print STDERR "fils =  $fils_id, pere=$pere_id\n";
	if ($fils->GetFather() != $pere) {
	    die ("Inconsistant data while analysing apomophy informations:\n".
		 "=> node '$fils_id' is not the son of '$pere_id', but of '",
		 $fils->GetFather()->Name()."'\n"		 
		 );
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

sub FillTreeApo1Phylip # Put sequence in structure node 
{
    my($tree)=shift;
    my($tabinfo)=shift;
    my($ancestor)=shift;
    my($fils);
    my($sequence);
  
    for (my $i=0; $i<=$#$tabinfo; $i++) {
	$fils=$tabinfo->[$i]->[1];
	$sequence=$tabinfo->[$i]->[3];
#	print "fils=$fils, seq=$sequence\n";
	$tree->GetNode($fils)->SetSequence($sequence);
    }
    
    $tree->GetRoot()->SetSequence($ancestor);
}

sub FillTreeApo2Phylip
{
    my($present_node)=shift;
    my($racine)=shift;
    my($ancetre_seq)=shift;
    my($tree)=shift;
    my($fatherseq);
    my($childseq);
    my($child);
    my($newchildseq);
    
    #print "noeud: ", $present_node->{"id"}, "\n";
    if ($present_node eq $racine) {
	$fatherseq=$ancetre_seq;
    } else {
	$fatherseq=$present_node->GetFather()->GetSequence();
    }
    $childseq=$present_node->GetSequence();
    my($longueur)=length($childseq);
    if (length($fatherseq) != $longueur) {
	die "Error: reconstructed sequences have not the same length at note ",
	$present_node->GetFather()->Name(), " ($fatherseq) and at node ",
	$present_node->Name(), " ($childseq)\n";
    } 
    my($br_len)=0;
    for (my $i=0; $i<$longueur;$i++) {
	my($fathersite)= substr($fatherseq,$i,1);
	my($childsite)= substr($childseq, $i,1);
	if ($childsite eq ".") {
	    $childsite=$fathersite;
	    
	}
	elsif ($fathersite ne $childsite) {
	    $br_len++;
	    my($apo_num)=$i+1;
	    my($apo_sens)=ALTree::Sens->New($fathersite."->".$childsite);
	    my $site;
	    if (not $tree->HasSiteIndex($apo_num)) { 
       		$site=ALTree::SitePerTree->New($apo_num);
		$tree->AddSite($site);
	    } else {
		$site=$tree->GetSite($apo_num);
	    }
	    $site->IncNbMut();

	    my($ref_site_sens)=$site->ProvideSens($apo_sens);
	    
	    $present_node->AddApo($ref_site_sens); # lie arbre et hash_site_sens
	}
	$newchildseq.=$childsite; #rajoute
	$present_node->SetSequence($newchildseq); # rjouté
    }
    $present_node->SetBrLen($br_len);
    foreach $child ($present_node->GetChildrenList()) { 
	FillTreeApo2Phylip($child, $racine, $ancetre_seq, $tree);
    }
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
