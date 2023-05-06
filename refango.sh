#!/bin/sh
USAGE="USAGE: $0 [-eK] [-wK] [-f] [-hg] [-pK] [-ebm] species1 species2 G1.el G2.el OBOfile.obo gene2go [alignFile(s)]
PURPOSE: evaluate the functional significance of a network alignment using GO terms.
Networks MUST be in edgelist format, with file names that end in '.el'; obo file name MUST end in '.obo'.
Default behavior: for each GO term, print p-values compared to random alignment based on (1) Poisson distribution, and
    (2) our novel exact Combinatiorial method (see our paper 'Exact p-values for network alignments')
speciesN: use our (admittedly odd) abbreviations, eg HS for HSapiens, RN for RNorvegicus, etc.
Note: Errors have a severity from 0 (not severe) to 9 (Fatal). Severity < 9 errors can be converted to warnings, or ignored.
OPTIONS:
-eK means: only halt if severity >= K (Default: K=0, ie., halt on all errors)
-wK means: only warn if severity >= K (Default: K=0, ie., warn on all errors that didn't cause a halt)
-f means: for each GO term, list its frequencies in G1, G2, as well as the alignment
-hg means: for each GO term, print p-values according to the HypeGeometric distribution (WARNING: slow!)
-pK means: for each pair of proteins in the alignment, list their GO terms based on the value of K:
    K=s (shared) means print comma-separated list of GO terms SHARED by the two proteins
    K=a (all) means print the above shared list, a tab, then GO terms of protein 1, a tab, then GO terms of protein 2.
-ebm means: combine the p-values across GO terms using the Empirical Brown's Method."

# Functions
die(){ (echo "$USAGE"; echo "FATAL ERROR: $@")>&2; exit 1; }
warn(){ (echo "WARNING: $@")>&2; }
not(){ if eval "$@"; then return 1; else return 0; fi; }
newlines(){ awk '{for(i=1; i<=NF;i++)print $i}' "$@"; }
parse(){ awk "BEGIN{print $@}" </dev/null; }

# generally useful Variables
BASENAME=`basename "$0" .sh`
[ $BASENAME = "$BASENAME" ] || die "something weird with filename in '$BASENAME'"

# Temporary Filename + Directory (both, you can use either, note they'll have different random stuff in the XXXXXX part)
TMP=`mktemp /tmp/$BASENAME.XXXXXX`
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
 trap "/bin/rm -rf $TMP $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

PATH="`dirname "$0"`:$PATH"
export PATH
FREQ=0
HyperGeo=0
ExactComb=1
LIST_PROTEINS=0
ERRORLEVEL=0
WARNLEVEL=0
EBM=false

while true; do
    case "$1" in
    -[Ee][Bb][Mm]) EBM=true; LIST_PROTEINS=2; FREQ=1;; # turn on additional options needed for EBM to work
    -e[0-9]) ERRORLEVEL=`echo $1 | sed 's/^-e//'`;
	[ "$ERRORLEVEL" -ge 0 -a "$ERRORLEVEL" -le 9 ] || die "ERRORLEVEL must be 0 through 9 inclusive, not '$ERRORLEVEL'";;
    -w[0-9]) WARNLEVEL=`echo $1 | sed 's/^-w//'`;
	[ "$WARNLEVEL" -ge 0 -a "$WARNLEVEL" -le 9 ] || die "WARNLEVEL must be 0 through 9 inclusive, not '$WARNLEVEL'";;
    -p*) case "$1" in
	    -p|-ps) LIST_PROTEINS=1;;
	    -pa) LIST_PROTEINS=2;;
	    *) die "Unknown option '$1'";;
	esac
	;;
    -f) FREQ=1;;
    -hg) HyperGeo=1;;
    -*) die "unknown option '$1'";;
    *) break ;;
    esac
    shift
done

[ $# -ge 7 ] || die "too few args"
s1=$1
s2=$2
shift 2

( BioGRIDname $s1; BioGRIDname $s2) > $TMPDIR/names.txt

hawk '
    BEGIN{OBO_ARGIND=1e30; LIST_PROTEINS='$LIST_PROTEINS'; HyperGeo='$HyperGeo'; ExactComb='$ExactComb';FREQ='$FREQ'}

    # Expected number of aligned orthologs in a random alignment of G1 and G2 with n1,n2 nodes and h common ortholog pairs.
    function ExpectedPairedOrthologs(h,n1,n2, hg,k) {hg=0;for(k=0;k<h;k++)hg+=(h-k)/(n1*n2-k*(n1+n2-k));return hg}

    function logAlignSearchSpace(n1,n2, result){
	if(n1 in _memLogAligSS && n2 in _memLogAligSS[n1]) return _memLogAligSS[n1][n2];
	ASSERT(n1>=0&&n2>=0,"AligSearchSpace: (n1,n2)=("n1","n2") cannot be negative");
	if(n1>n2) result=logAlignSearchSpace(n2,n1); else result = logFact(n2)-logFact(n2-n1);
	return (_memLogAligSS[n1][n2] = result);
    }
    function AlignSearchSpace(n1,n2){return exp(logAlignSearchSpace(n1,n2))}

    function CountGOtermAlignments(n1,n2,l1,l2,k,      ll,M,U,mu,muMin,muMax) {
	ASSERT(n1<=n2, "Sorry, shared probability of GO terms requires n1<=n2");
	ASSERT(k>=0, k" must be greater than zero in CountGoTermAligs");
	ll=MIN(l1,l2); # lower and upper lambdas
	if(k>ll) return 0;
	if(ll==0)return (k==0?AlignSearchSpace(n1,n2):0);
	if(l2==n2)return (k==l1?AlignSearchSpace(n1,n2):0);
	if(l1>n2-l2) # There are more annotated pegs than unannotated holes; at least l1-(n2-l2) anopegs *must* match
	    if (k < l1-(n2-l2)) return 0;
	M=choose(l1,k) * choose(l2,k) * fact(k); # aligning the k matched pairs
	U=0
	muMin=MAX(0,(n1-k)-(n2-l2));
	muMax=MIN(n1-l1,l2-k);
	for(mu=muMin;mu<=muMax;mu++){ # sum over possible values for numAnnotatedPegs aligning to l2-k annotated holes.
	    AS1=AlignSearchSpace(mu,l2-k); # aligning annot. pegs to unannot. holes
	    AS2=AlignSearchSpace(n1-l1-mu,n2-l2-(l1-k)); # remaining unannot pegs aligned to unannot holes
	    choices = choose(n1-l1,mu) * choose(n2-l2,l1-k)
	    U += choices * AS1 * AS2
	}
	return M * fact(l1-k)*U;
    }

    # Below is just the logarithmic version of the above to handle much bigger numbers.
    function logCountGOtermAlignments(n1,n2,l1,l2,k,      ll,M,U,Utmp,mu,muMin,muMax) {
	ASSERT(n1<=n2, "Sorry, shared probability of GO terms requires n1<=n2");
	ASSERT(k>=0, k" must be greater than zero in logCountGoTermAligs");
	ll=MIN(l1,l2); # lower and upper lambdas
	if(k>ll) return log(0);
	if(ll==0)return (k==0?logAlignSearchSpace(n1,n2):log(0));
	if(l2==n2)return (k==l1?logAlignSearchSpace(n1,n2):log(0));
	if(l1>n2-l2) # There are more annotated pegs than unannotated holes; at least l1-(n2-l2) anopegs *must* match
	    if (k < l1-(n2-l2)) return log(0);
	M=logChoose(l1,k) + logChoose(l2,k) + logFact(k);
	muMin=MAX(0,(n1-k)-(n2-l2));
	muMax=MIN(n1-l1,l2-k);
	U=0
	for(mu=muMin;mu<=muMax;mu++){ # sum over possible values for numAnnotatedPegs aligning to l2-k annotated holes.
	    # do NOT try any memoization here; too many parameters means not enough repeats -> actually SLOWER
	    AS1=logAlignSearchSpace(mu,l2-k); # aligning annot. pegs to unannot. holes
	    AS2=logAlignSearchSpace(n1-l1-mu,n2-l2-(l1-k)); # remaining unannot pegs aligned to unannot holes
	    Utmp = AS1+AS2 + logChoose(n1-l1,mu)
	    U=LogSumLogs(U,Utmp);
	}
	U += logChoose(n2-l2,l1-k)
	return  M + logFact(l1-k)+U;
    }

    function CHECK(level,cond,str) {if(level>='$ERRORLEVEL')ASSERT(cond,str);else if(level>='$WARNLEVEL')WARN(cond,str)}
    ARGIND==1{ # get tax IDs
	CHECK(9,NF>=3,"tax ID from BioGRIDname line is too short: "$0);
	name[FNR]=$1; #printf "name[%d]=%s[%s]", FNR, $1,$2
	for(i=3;i<=NF;i++){++tax2net[$i][FNR]; #printf "\t%d",$i
	}
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
	#print "Finished parsing OBO file; numNets is",numNets,"OBO_ARGIND is",OBO_ARGIND
	for(g in OBO_ID) if(!OBO_Ob[g]) # if not obsolete
	    OBOtraceRoot(g,g)
    }}
    ############################# END OF READING OBO FILE ######################


    ############################# START READING GENE2GO FILE ######################
    ARGIND==OBO_ARGIND+1 && !/^#/ && ($1 in tax2net)&&($2 in nodes) {
    # Read the gene2go file
	p=$2; g=$3; # Ensure all of the following are set to 1
	for(net in tax2net[$1]){
	    pGO[net][p][g] = GOp[net][g][p] = GOobserved[net][g] = 1
	    #printf "pGO[%d][%s][%s]=1\n",net,p,g
	}
	next
    }
    ENDFILE{ if(ARGIND==OBO_ARGIND+1) {
	# print "Finished parsing gene2go file"
	# process GOfreqs of the ancestors of the ones observed in the gene2go file
	for(net in GOobserved) {
	    for(g in GOobserved[net]) { # for every GO term that annotates SOME protein in this network...
		if(g in OBOancestors) { # there can be old GO terms in gene2go that are not in a newer OBO file
		    for(g1 in OBOancestors[g]){ # for each ancestor GO term...
			# For all the proteins that g annotates, record the fact that the ancestor g1 also annotates them.
			for(p in GOp[net][g]){
			    ++pGO[net][p][g1]; ++GOp[net][g1][p];
			}
		    }
		}
	    }
	}
	# Now that all proteins have been updated with OBO traces through all GO terms, count GO frequencies.
	for(net in GOobserved)for(g in GOp[net]) {
	    allGO[g] = 1;
	    GOfreq[net][g] = length(GOp[net][g])
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
	CHECK(9,numNets == 2, "Sorry, can only work with 2-network alignfiles at the moment");
	CHECK(5,$1 in deg[1], "line "FNR" col1 of alignment contains unknown node "$1);
	CHECK(5,$2 in deg[2], "line "FNR" col2 of alignment contains unknown node "$2);
	e=0 # entropy
	p=1 # p-value
	CHECK(8,NF == numNets,"expecting "numNets" columns in alignfile "FILENAME", line "FNR": "$0);
	for(i=1;i<=numNets;i++)if(!($i in pGO[i])) next; # assuming fixed align file, not arbitrary clusters
	if(LIST_PROTEINS){
	    CHECK(8,NF==2,"oops, LIST_PROTEINS can only currently handle exactly 2 proteins...")
	    printf "%s,%s\t",$1,$2; GOlist="";
	}
	for(g in pGO[1][$1]){
	    if(g in pGO[2][$2]){
		++numAligPairs[g]
		if(LIST_PROTEINS) GOlist = GOlist sprintf("%s,",g);
		pg=GOfreq[1][g]*GOfreq[2][g]/length(deg[1])/length(deg[2])
		p*=pg
		e+=-pg*log(pg)
	    }
	}
	if(LIST_PROTEINS){
	    if(LIST_PROTEINS == 2) {
		if(!GOlist) GOlist="0"; # there are no shared GO terms
		GOlist = GOlist "\t"; for(g in pGO[1][$1]) GOlist=GOlist sprintf("%s,",g);
		GOlist = GOlist "\t"; for(g in pGO[2][$2]) GOlist=GOlist sprintf("%s,",g);
	    }
	    gsub(",\t","\t",GOlist); gsub(",$","",GOlist); # get rid of trailing commas
	    print GOlist; # could be empty if LIST_PROTEINS < 2, and print with a newline regardless
	}
    }

    {_BEGINFILE=0} # first line of any file turns it off
    BEGINFILE{printf "*** ARGIND %d FILENAME %s ***\n", ARGIND, FILENAME;_BEGINFILE=1}

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

	# Print the header for all GO terms
	printf "GO_term"
	if(FREQ) {
	    printf " GOfreq shared G1 G2", numAligPairs[g], GOfreq[1][g], GOfreq[2][g]
	}
	printf " PoissTail( l k )= p-value\n";
	for(g in allGO)
	{
	    if(!(g in GOfreq[1])) GOfreq[1][g]=0;
	    if(!(g in GOfreq[2])) GOfreq[2][g]=0;
	    printf "%s", g
	    if(FREQ) {
		printf " GOfreq %d %d %d", numAligPairs[g], GOfreq[1][g],GOfreq[2][g]
	    }
	    expected = GOfreq[1][g]*GOfreq[2][g]/n2 #*numAligFiles;
	    totalExpected += expected
	    totalPairs+=numAligPairs[g]
	    pValuePoisson = Poisson1_CDF(expected,numAligPairs[g]);
	    if(pValuePoisson == 0 && numAligPairs[g]>expected) pValuePoisson=1e-320 # recover from underflow
	    printf " ( %g %d )= %g", expected, numAligPairs[g], pValuePoisson;

	    if(HyperGeo) {
		# number of pairs that share GO term g, indep of alignment, in prep for hyperGeom
		hyper_k= numAligPairs[g];
		hyper_n= n1 #*numAligFiles
		CHECK(5,(g in GOfreq[1] && g in GOfreq[2]),g " not in one of the networks");
		CHECK(5,GOfreq[1][g]>0,"GOfreq[1]["g"] is "GOfreq[1][g]);
		CHECK(5,GOfreq[2][g]>0,"GOfreq[2]["g"] is "GOfreq[2][g]);
		l1=MIN(GOfreq[1][g], GOfreq[2][g]);
		l2=MAX(GOfreq[1][g], GOfreq[2][g]);
		CHECK(5,l1>0,"l1");
		CHECK(5,l2>0,"l2");
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
		printf " ExactComb (%d %d %d) = %g",GOfreq[1][g],GOfreq[2][g],numAligPairs[g],logPrint(Exact_MU-Exact_A);fflush("")
	    }
	    print ""
	}
	printf "GO total expected Poisson mean is %g, got %d, \n", totalExpected, totalPairs;
	#printf "log10(p-value) %g\n", Log10Poisson1_CDF(totalExpected,totalPairs);
    }}' $TMPDIR/names.txt "$@" | tee $TMPDIR/GOeval.out

if $EBM; then
    echo "Please wait while we run the Empirical Brown's Method to combine the above p-values... may take some time... " >&2
    echo "(Note: there is a C-compiled version of ebm in my github repo for libwayne, in libwayne/tests/ebm.c)" >&2
    gawk '/^GO:/{pV[$1]=$NF}/\tGO:/{n=split($2,a,",");for(i=1;i<=n;i++)GOpair[a[i]][$1]=pair[$1]=1}END{printf "GOterm\tpValue"; for(p in pair)printf "\t%s",p; print ""; for(g in GOpair){printf "%s\t%g",g,pV[g]; for(p in pair)printf "\t%d",1*GOpair[g][p]; print ""}}' $TMPDIR/GOeval.out | tee $TMPDIR/ebm.in | ebm.sh -v || trap "" 0 1 2 3 15
fi
