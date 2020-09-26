#!/bin/sh
USAGE="USAGE: $0 [-f] [-hg] [-com] [-p] species1 species2 G1.el G2.el OBOfile.obo gene2go [alignFile(s)]

Notes:
    - default: for each GO term, print p-values compared to random alignment based on Poisson distribution.
    - '-f' means: for each GO term, list its freuencies in G1, G2, and the alignment
    - '-hg' means: for each GO term, print p-values according to the HypeGeometric distribution (WARNING: slow!)
    - '-com' means: for each GO term, print p-values according to our novel exact Combinatiorial method
    - '-p' means: for each pair of proteins in the alignment, list the GO terms they share
    - networks MUST be in edgelist format, with file names that end in '.el'
    - obo file name MUST end in '.obo'
"
die() { echo "$USAGE
FATAL: $@">&2; exit 1
}
FREQ=0
HyperGeo=0
ExactComb=0
LIST_PROTEINS=0

while true; do
    case "$1" in
    -p) LIST_PROTEINS=1; shift;;
    -f) FREQ=1; shift;;
    -hg) HyperGeo=1; shift;;
    -com) ExactComb=1; shift;;
    -*) die "unknown option '$1'";;
    *) break ;;
    esac
done

[ $# -ge 7 ] || die "not enough args"
s1=$1
s2=$2
shift 2

TMPDIR=`mktemp -d`
trap "/bin/rm -rf $TMPDIR" 0 1 2 3 15
( BioGRIDname $s1; BioGRIDname $s2) > $TMPDIR/names.txt

hawk '
    BEGIN{OBO_ARGIND=1e30; LIST_PROTEINS='$LIST_PROTEINS'; HyperGeo='$HyperGeo'; ExactComb='$ExactComb';FREQ='$FREQ'}
    ARGIND==1{ # get tax IDs
	ASSERT(NF>=3,"tax ID from BioGRIDname line is too short: "$0);
	name[FNR]=$1;for(i=3;i<=NF;i++)tax[FNR][$i]=taxIDs[$i]=tax2net[$i]=FNR
	next
    }
    index(FILENAME,".el")+2==length(FILENAME) { # network files
	++deg[ARGIND-1][$1]
	++deg[ARGIND-1][$2]
	nodes[$1]=nodes[$2]=1
	next;
    }

    ############################# START READING OBO FILE ######################
    index(FILENAME,".obo")+3==length(FILENAME) { # OBO file
	gsub("!.*$","")  # delete all comments
	if(/^id:/){id=$2;OBO_ID[$2]=1}
	if(/^is_a:/){OBO_P[id][$2]=1}
	if(/^name:/){gsub("name: ","");gsub(" ","_");OBOname[id]=$0}
	if(/^namespace:/){OBO_NS[id]=$2}
	if(/^alt_id:/||/^consider:/){OBOmap[$2][id]=1}
	if(/^is_obsolete:/){OBO_Ob[id]=1}
	if(/^.Typedef/){nextfile} # this marks the end of the actual heirarchy definition
	next;
    }
    # Recursively find the root starting at GO term g, assigning ancestors along the way (do not include self as ancestor)
    function OBOtraceRoot(g,g1,		p,i,n){
	if(g1!=g) OBOancestors[g][g1]=1;
	if(isarray(OBO_P[g1])) { # This is only false once we reach one of the three roots of the GO hierarchy
	    for(p in OBO_P[g1]){
		OBOtraceRoot(g,p)
	    }
	}
    }
    ENDFILE{if(index(FILENAME,".obo")+3==length(FILENAME)){ # parse the OBO file
	numNets = length(deg);
	OBO_ARGIND = ARGIND
	# Post-process the OBO file:
	for(g in OBO_ID) if(!OBO_Ob[g]) # if not obsolete
	    OBOtraceRoot(g,g)
    }}
    ############################# END OF READING OBO FILE ######################


    ############################# START READING GENE2GO FILE ######################
    ARGIND==OBO_ARGIND+1 && !/^#/ && ($1 in taxIDs)&&($2 in nodes) {
    # Read the gene2go file
	net=tax2net[$1]; p=$2; g=$3;
	# Ensure all of the following are set to 1
	pGO[net][p][g] = GOp[net][g][p] = GOobserved[net][g] = allGO[g] = 1
	next
    }
    ENDFILE{ if(ARGIND==OBO_ARGIND+1) {
	# process GOfreqs of the ancestors of the ones observed in the gene2go file
	for(net in GOobserved) {
	    for(g in GOobserved[net]) { # for every GO term that annotates SOME protein in this network...
		if(g in OBOancestors) { # there can be old GO terms in gene2go that are not in a newer OBO file
		    for(g1 in OBOancestors[g]){ # for each ancestor GO term...
			# For all the proteins that g annotates, record the fact that the ancestor g1 also annotates them.
			for(p in GOp[net][g]){
			    pGO[net][p][g1] = GOp[net][g1][p] = 1
			}
		    }
		}
	    }
	}
	# Now that all proteins have been updated with OBO traces through all GO terms, count GO frequencies.
	for(net in GOobserved)for(g in GOp[net]) {
	    GOfreq[net][g] = length(GOp[net][g])
	}
	if(0)for(net in GOobserved) {
	    for(g in GOobserved[net]) {
		GOfreq[net][g]+=GOobserved[net][g] # "+=" since it may be non-zero b/c ancestor of another protein.
		if(g in OBOancestors) { # there can be old GO terms in gene2go that are not in a newer OBO file
		    for(g1 in OBOancestors[g]){
			GOfreq[net][g1] += GOobserved[net][g] # all ancestors get a frequency bump from g.
			# For all the proteins that g annotates, record the fact that the ancestor g1 also annotates them.
			for(p in GOp[net][g]){
			    ++pGO[net][p][g1]; ++GOp[net][g1][p]
			}
		    }
		}
	    }
	}
    }}
    ############################# END PROCESSING GENE2GO FILE & ANCESTOR GO TERMS ######################


    ###################### ALIGNMENT FILES ##############################
    ARGIND>=OBO_ARGIND+2{
	if(_BEGINFILE) {
	    _BEGINFILE=0;
	    print "alignment_BEGINFILE",ARGIND,FILENAME
	    delete numAligPairs
	}
	ASSERT(numNets == 2, "Sorry, can only work with 2-network alignfiles at the moment");
	ASSERT($1 in deg[1], "line "FNR" col1 of alignment contains unkown node "$1);
	ASSERT($2 in deg[2], "line "FNR" col2 of alignment contains unkown node "$2);
	e=0 # entropy
	p=1 # p-value
	ASSERT(NF == numNets,"expecting "numNets" columns in alignfile "FILENAME", line "FNR": "$0);
	for(i=1;i<=numNets;i++)if(!($i in pGO[i])) next; # assuming fixed align file, not arbitrary clusters
	if(LIST_PROTEINS){
	    ASSERT(NF==2,"oops, LIST_PROTEINS can only currently handle exactly 2 proteins...")
	    printf "%s,%s",$1,$2
	}
	for(g in pGO[1][$1])if(g in pGO[2][$2]){
	    if(LIST_PROTEINS) printf "\t%s",g
	    pg=GOfreq[1][g]*GOfreq[2][g]/length(deg[1])/length(deg[2])
	    p*=pg
	    e+=-pg*log(pg)
	    ++numAligPairs[g]
	}
	if(LIST_PROTEINS) print ""
    }

    {_BEGINFILE=0} # first line of any file turns it off
    BEGINFILE{printf "*** %s ***\n", FILENAME;_BEGINFILE=1}

    # Poisson measure:
    # For any given protein in network G, the probability it is annotated with GO term g is GOfreq[G][g]/n.
    # Thus, for a random alignment, a pair of proteins are *both* annotated with g with probability
    # (GOfreq[1][g]/n1) * (GOfreq[2][g]/n2).
    # Given that there are n1 pairs in the alignment, the expected number of *pairs* that are annotated with g is thus
    # n1 * (the above), and the n1s cancel out, giving GOfreq[1][g]*GOfreq[2][g]/n2 == lambda
    # Thus, for random alignments, the number of pairs annotated with g is Poisson distributed with mean lambda.
    # If we are looking at the concatenation of k different alignment files, then multiply lambda by numAligFiles.
    # Except we may not actually have n1 pairs because some will have *no* annotations, so use numProteins in gene2go
    ENDFILE {if(ARGIND>=OBO_ARGIND+2) {
	print "Statistics for ENDFILE "FILENAME

	totalExpected=0;
	totalPairs=0
	n1=length(deg[1])
	n2=length(deg[2])

	if(ExactComb){
	    Exact_A=logAlignSearchSpace(n1,n2); # this is a constant, compute it only once
	    printf "Computed logAlignSearchSpace(%d,%d)=%g\n",n1,n2,Exact_A
	}

	for(g in allGO) if(g in GOfreq[1] && g in GOfreq[2])
	{
	    printf "%s", g
	    if(FREQ) {
		printf " GOfreq %d %d %d", GOfreq[1][g],GOfreq[2][g],numAligPairs[g]
	    }
	    expected = GOfreq[1][g]*GOfreq[2][g]/n2 #*numAligFiles;
	    totalExpected += expected
	    totalPairs+=numAligPairs[g]
	    pValuePoisson = Poisson1_CDF(expected,numAligPairs[g]);
	    if(pValuePoisson == 0 && numAligPairs[g]>expected) pValuePoisson=1e-320 # recover from underflow
	    printf " PoissonTail( %g %d )= %g", expected, numAligPairs[g], pValuePoisson;

	    if(HyperGeo) {
		# number of pairs that share GO term g, indep of alignment, in prep for hyperGeom
		hyper_k= numAligPairs[g];
		hyper_n= n1 #*numAligFiles
		ASSERT((g in GOfreq[1] && g in GOfreq[2]),g " not in one of the networks");
		ASSERT(GOfreq[1][g]>0,"GOfreq[1]["g"] is "GOfreq[1][g]);
		ASSERT(GOfreq[2][g]>0,"GOfreq[2]["g"] is "GOfreq[2][g]);
		l1=MIN(GOfreq[1][g], GOfreq[2][g]);
		l2=MAX(GOfreq[1][g], GOfreq[2][g]);
		ASSERT(l1>0,"l1");
		ASSERT(l2>0,"l2");
		hyper_K= l1*l2 - l1*(l1-1)/2 #*numAligFiles
		hyper_N= n1*n2 - n1*(n1-1)/2 #*numAligFiles
		printf " HyperGeom( %d %d %d %d )=",g,hyper_k,hyper_n,hyper_K,hyper_N;fflush("")
		#pValueHyperGeom = HyperGeomTail(hyper_k,hyper_n,hyper_K,hyper_N) # This is WAY slow!
		# logHyperGeoTail is much faster:
		pValueHyperGeom = exp(logHyperGeomTail(hyper_k,hyper_n,hyper_K,hyper_N))
		printf " %g", pValueHyperGeom
		#printf "diff %g %g %g\n", pValuePoisson, pValueHyperGeom, ABS(pValuePoisson-pValueHyperGeom)/MIN(pValuePoisson,pValueHyperGeom)
	    }
	    if(ExactComb){
		Exact_MU=logCountGOtermAlignments(n1,n2,GOfreq[1][g],GOfreq[2][g],numAligPairs[g]);
		printf " ExactComb (%g %g) = %g",Exact_MU,Exact_A,exp(Exact_MU-Exact_A);fflush("")
	    }
	    print ""
	}
	printf "GO total expected Poisson mean is %g, got %d, \n", totalExpected, totalPairs;
	#printf "log10(p-value) %g\n", Log10Poisson1_CDF(totalExpected,totalPairs);
    }}' $TMPDIR/names.txt "$@"
