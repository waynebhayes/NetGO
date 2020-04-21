#!/bin/sh
# Sept. 2019
# NetGO.awk (c) Wayne B. Hayes (whayes@uci.edu). 
# Yes, we call this NetGO.awk even though technically it's a Bourne Shell script. Most of the "meat" of the implementation
# is in the awk at the bottom of this Bourne shell script.

exeDIR=`dirname "$0"`
PATH="$exeDIR:$PATH" # needed for "hawk" (Hayes awk)

USAGE="USAGE: $0 [-verbose] [-L|-R] OBOfile.obo gene2goFile alignFile[s]

    -L: 'Lenient'. The default behavior is what we call 'Dracanion' in the paper, which
    insists that a GO term must annotate every protein in a cluster for it to count.
    The Lenient option gives a GO term a weight per-cluster that is scaled by the
    number of proteins it annotates (so long as it's more than 1).

    -R: Don't compute NetGO at all; instead, assume the alignFiles are *Resnik* outputs, and
    compte the NetGO-weighted Resnik score.  The first word is ignored (usually 'max', meaning
    the Maximum variant of Resnik), the 2nd and 3rd columns are the proteins aligned, and the
    4th column is the Resnik score.

    OBOfile.obo: obo file from the Gene Ontology website; should be same version as gene2go.

    gene2goFile: a standard-format gene2go file downloaded from the GO consortium's
    website. For now we only use columns 2 (protein name) and 3 (GO term). We ignore
    the species, which is a simplification since about 20% of proteins appear in two
    species, but only 2% appear in more than 2. Should be same version as OBO file above.

    alignFile: each line consists of a cluster of any number of proteins; proteins can
    appear in more than one cluster. Note it must use the same protein naming convention
    as your gene2go file.

"
die() { echo "$USAGE" >&2; echo "$@" >&2; exit 1
}

# References:
# Main paper this script is based on is the 2019 paper by Hayes, "NetGO: a data-driven approach..."
# Older paper with more justification: Hayes and Mamano 2017, doi: 10.1093/bioinformatics/btx716 "SANA NetGO: a combinatorial approach..."

# IMPLEMENTATION DECISIONS
# Though the math for NetGO is straightforward, there are are many implementation decisions.
# Few of these are carved in stone.
#
# First, for now WE IGNORE SPECIES IN THE gene2go FILE.
# Most proteins (80%) occur in only one species, so only 20% occur in more than one species, and only about 2% occur in more
# than two species. It is a good bet (though not certain) that the *exact* same protein will perform a highly similar set of
# functions in whatever species it appears in. So GO terms valid for one species are likely to be valid for the same protein
# in another species.
# Furthermore, the species-taxonomic ID (which is the first columns of the gene2go file) is not unique. For example, at last
# count S.cerevisiae has 321 (!!!) different strains, and thus over 300 taxonomic IDs. Even mouse had over a dozen different
# taxonomic IDs. Dealing with this complexity would be a mess.
# For the above reasons, for now we ignore the first column (taxonomic ID) of the gene2go file, and assume than the pair
# (gene,GO) is enough to determine GO counts, GO frequencies, and protein annotations.

# Second, protein names are always in flux, and some proteins have vastly more annotations than others. Even if you use the
# same naming convention in your alignment file as the appears in the gene2go file (which is necessary for this script to work),
# there will likely be some proteins in your alignment that do not appear in the gene2go file, and some proteins in the gene2go
# file for your species that are not in your alignment file. Finally, the gene2go file is HUGE, and we don't want to read the
# whole thing in every time. Thus, this script needs to know YOUR list of proteins (gotten from your alignment/cluster file)
# before it reads the gene2go file; it will only record annotation information for the proteins you need.

# The differing sets of proteins (your alignment file vs. those in the gene2go file) also brings up the question of "how many
# proteins", and "how many GO terms", are there, when we need to know their count? The GO term frequency (kappa_g in the paper)
# is based on how many times we see the GO term g *among the proteins in your alignment file*---and recall that we also ignore
# species ID, so this count can fluctuate a bit depening on your input data.

# CLUSTER SIZE
# How many proteins are in your cluster? This depends on how we count. If the exact same protein name appears twice in the
# same cluster, does it count for 1, or 2? That is, do we remove duplicates? This may depend on what you want to do.
# Certainly seeing the same protein $p$ a second time does not add to the set of GO annotations in this cluster after we've
# seen $p$ before, so should it change the *count* of these GO terms? 
# Things to consider in this question: though you may be tempted to remove duplicate proteins, this means you don't get any
# points for correctly putting both copies of the same protein into the same cluster *if they came from different species*;
# this is a significant drawback because correctly aligning a protein to itself between species is a pretty beneficial thing
# to do, and getting no points for it sucks. Furthermore, we claim in the paper that aligning a network to itself should,
# of course, result in a NetGO score of (very close to) 1; if you remove duplicates then aligning a network to itself results
# in every cluster having only one protein, and a NetGO score of exactly zero.
# For these reasons, we have decided to allow duplicates. The number of times a protein $u$ appears in a cluster is kept
# in the variable M[u] (M=membership set).

# CODING CONVENTIONS
# AWK is an error prone language. I like it primarily because it's supremely portable (e., any Linux distribution), it's
# lightning fast (usually faster than Java or Python), and it can be used on the Unix command line to do simple (or complex)
# tasks. Variables have dynamic type, and there are only 3 types: double precision, string, and array. All arrays are
# associative (ie., dictionaries), and an element of an array can be another array, and different elements of the same array
# can have different type.
#
# To implement "sets" in AWK, I use A[x]=1 to mean that x is in the set A; x can be a number, string, etc. Sometimes
# I interpret it as a multi-set, where the value is an integer rather than Boolean. such as the membership array
# M[] which contains string elements that are names of proteins that occur in the cluster currently being processed.
# (The multi-set is how we keep track of multiple occurances of the same protein in the same cluster.)
# Note also that in AWK, if a variable (even an array element) doesn't exist and you try to access it, it will be created
# and then returned with value "" (empty string) or zero, depending on context. THIS IS DANGEROUS, because it means you should
# not check for the existence of a member in an array by asking for its value, because each one you ask for, if it does not
# exist, creates a new (useless) member of the array. This can quickly cause the memory usage of awk to baloon. Instead,
# use the syntax (x in A) to check if A[x] exists, befory trying to access it. The (x in A) syntax does not create A[x] if
# it doesn't already exist.

# GLOBAL VARIABLES:
# By default all variables in AWK are global. The only exceptions are function parameters (see below).
# The primary data variables (global) in this program are set in the ARGIND program blocks at the bottom of the AWK script:
#
# CA is the "cluster alignment". Contents: CA[line_number][column_number]=proteinName
# That is, CA[L] is the cluster of proteins that occured on line L of the input cluster (align) file. The number of protein
# names, including duplicates, that occured on that line is length(CA[line_number]).
#
# pC: protein-to-cluster multi-set; tells us which cluster(s) a protein appeared in, and its multiplicity in that cluster
# Contains: pC[p][L]=count means protein p occurs in cluster (line number) L, count times.
# length(pC) is thus the number of unique protein names that occured in your alignment file;
# length(pC[p]) is the number of clusters # it appears in.
#
# pGO: the protein-to-GO multiset map
# pGO[p][g]=count is the number of times p is annotated by g (including OBO-inferred annotations) across all species,
#    but only for proteins p in the alignment
# FOR NOW WE IGNORE SPECIES ID in the gene2go file. This means that there's a possibility that p has been annotated with g in
# another species, but not the one you're aligning. This will be (a) rare, and (b) probably not an entirely invalid assumption.
#
# GOp: GO-to-protein multiset (inverse relationship of the above, but a multiset instead of just a set)
# GOp[g][p]=count is the number of times GO term g (including OBO-inferred annotations) annotates protein p across all species,
#    but only for proteins p in the input alignment.
# This usually means that the GO term occurs with protein p in the gene2go file many times, each with a different evidence code.
# This is useful information to have (eg., more evidence codes makes the annotation more certain to be correct).
#
# EMPTY SETS
# In several places, we assume that the length of an array defines how many relationships it has. The problem is that,
# since arrays are actually associative, if the relationship does not exist at all, then the element will be not be zero,
# but will instead not exist. This breaks some parts of the code. To fix it, we sometimes create an empty array, to specify
# "this relationship does not exist". Unfortunately there is no way in AWK to create an empty array. The only way to do it
# is to bring the array into existence by accessing an (as-yet non-existant) member; this creates the array, and the one
# new member.  Then you just delete that one member, leaving an empty array.
# The syntax to do this (seen below in several places) is:
# A[0]=1; delete A[0]; # create dummy element A[0], then delete it, leaving A as an array with zero elements.

# awk allows no local variables other than parameters, but it also allows you to declare more parameters than are passed.
# Thus the convention is to declare local variables as extra parameters (with extra whitespace or a newline after the
# true parameter list).

RESNIK=0
DRACONIAN=1
VERBOSE=0
case "$1" in
-verbose|-V) VERBOSE=1;shift;;
esac
case "$1" in
-R*) RESNIK=1; shift;;
-L*) DRACONIAN=0;shift;;
-*) die "unknown option '$1'";;
esac
[ $# -ge 3 ] || die "expecting at least 3 arguments: OBOfile.obo, gene2goFile, and at least one clusterAlignFile"

OBO=$1; shift
(fgrep -q 'id: GO:' $OBO && fgrep -q 'is_a: GO:' $OBO) || die "first argument must be a go.obo file"

GENE2GO=$1; shift
[ `cat $GENE2GO | head -10 | grep -c '	GO:'` -ge 9 ] || die "2nd argument must be the gene2go file"

for i
do
    hawk '
    # Return the "knowledge" (ie., specificity) of a single GO term g.
    function K_g(g){if(g in GOp) return 1/length(GOp[g]); else return 0}

    # Return the sum of specifities of a set of GO terms.
    function K_gset(T,	g,sum){sum=0;for(g in T)sum+=K_g(g); return sum}

    # Return the "knowledge" level of a protein p, which is just K_gset(T), where T is its sets of GO annotations.
    function K_p(p,	sum){sum=0;if(p in pGO)return K_gset(pGO[p]); else return 0}

    # Knowledge in a simple pairwise alignment in which A[u]=v
    # This function was mostly just used for testing simple cases early on; it is superceded by K_AC below.
    function K_A2(A,
	sum,u,v,Tuv){sum=0;for(u in A)if(u in pGO){v=A[u];if(v in pGO){SetIntersect(Tuv,pGO[u],pGO[v]); sum+=K_gset(Tuv)}}; return sum}

    # Knowledge in a general alignment of protein clusters
    # Note technically it would be nice to have a K_C function (knowledge inside a cluster) but AWK makes this hard, so
    # instead we just compute K_C inline, inside this function.
    function K_AC(C,
	sum,numClusterFields,cl,i,u,g,T,M,K_C){
	sum=0;
	for(cl=1;cl<=length(C);cl++){
	    K_C=0;
	    if(VERBOSE) printf "Cluster %d numProteins %d\n",cl,length(C[cl])
	    numClusterFields=length(C[cl]);
	    delete M; M[0]=1;delete M[0]; # initialize M to empty list; M=set of protein members;
	    delete T; T[0]=1;delete T[0]; # initialize T to empty list; T=set of GO terms across the cluster;
	    if(numClusterFields>1){ # have to match more than one protein to be interesting (handle duplicates below)
		# Now go through the cluster, removing elements from T that do not occur in subsequent columns (DRACONIAN)
		# and incremently updating the count T[g] if the GO term annotates more than one protein.
		for(i=0;i<numClusterFields;i++){
		    u=C[cl][i]
		    ASSERT(u!="-"&&u!="_"&&u!="NA","INTERNAL ERROR: invalid protein got into K_AC");
		    if(i==0) { # initialize T to the GO terms of first protein in the cluster:
			u=C[cl][0];
			++M[u] # multi-set: keep track if protein u occurs more than once.
			for(g in pGO[u])++T[g];
		    } else {
			++M[u];
			if(DRACONIAN) {
			    for(g in T)if(!(g in pGO[u]))delete T[g] # cumulative set intersection of common GO terms
			} else {
			    for(g in T)++T[g]; # cumulative union of common GO terms
			}
		    }
		    if(VERBOSE){
			printf "\t%s (%d GOs { ", u,length(T)
			for(g in T) printf "%s(%d) ",g,GOfreq[g]
			printf "} K(%s)=%g)\n",u,K_gset(T)
		    }
		}
	    }
	    if(length(T)>0 && # if this cluster has any annotations...
		(length(M)>1 || # and it has more than one protein...
		    (length(M)==1 && M[u]>1))) # or the same protein multiple times...
	    {
		if(DRACONIAN) K_C += K_gset(T)
		else
		    for(g in T)if(T[g]>1) # g must annotate more than 1 protein to get counted
			K_C += T[g]*K_g(g)/numClusterFields
	    }
	    if(VERBOSE){
		printf "\t ClusterCommonGOs {"
		for(g in T)printf " %s",g
		printf " } K_C=%g\n",K_C
	    }
	    sum+=K_C
	}
	return sum
    }

    # Knowledge in a alignment of a bunch of clusters
    function sim_A2(A){return K_A2(A)/length(GOfreq)}

    #Old 2-column alignment files, obsolete but simple to test
    #ARGIND==1{C[$1]=$2; C_[$2]=$1}
    #ARGIND==2{if($2 in C || $2 in C_){++GOfreq[$3];++pGO[$2][$3];++GOp[$3][$2]}}

    BEGIN{DRACONIAN='$DRACONIAN';VERBOSE='$VERBOSE';RESNIK='$RESNIK'}
    #BEGINFILE{printf "Reading file %d %s\n",ARGIND,FILENAME}
    #ENDFILE{printf "Finished reading file %d %s\n",ARGIND,FILENAME}

    #Clusters version
    # CA[][]=cluster alignment; pC[p] = clusters this protein is in.
    ARGIND==1{
	startCol=1;endCol=NF
	if(RESNIK){
	    ASSERT(NF==5,"Expecting Resnik files to have exactly 5 columns, including leading NAF");
	    NAF[FNR]=$1
	    startCol=3;endCol=4
	    R[FNR]=$NF
	}
	n=0;
	for(i=startCol;i<=endCol;i++)if($i!="_"&&$i!="-"&&$i!="NA"){CA[FNR][n++]=$i;++pC[$i][FNR]}
    }

    # Read the OBO file
    ARGIND==2{gsub("!.*$","")  # delete all comments
	if(/^id:/){id=$2;OBO_ID[$2]=1}
	if(/^is_a:/){OBO_P[id][$2]=1}
	if(/^name:/){gsub("name: ","");gsub(" ","_");OBOname[id]=$0}
	if(/^namespace:/){OBO_NS[id]=$2}
	if(/^alt_id:/||/^consider:/){OBOmap[$2][id]=1}
	if(/^is_obsolete:/){OBO_Ob[id]=1}
	if(/^.Typedef/){nextfile} # this marks the end of the actual heirarchy definition
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
    ENDFILE{ if(ARGIND==2){
	    # Post-process the OBO file:
	    for(g in OBO_ID) if(!OBO_Ob[g]) # if not obsolete
		OBOtraceRoot(g,g)
    }}

    # Read the gene2go file
    ARGIND==3 && !/^#/ {
	# In the following, we only link the protein p to the GO term g if p is in the alignment.
	# Note that for algorithms that produce tiny alignments on large networks (ie., ignoring many proteins)
	# this could make frequent, vague GO terms look more specific than they really are.
	p=$2; if(p in pC) {
	    g=$3; ++GOfreqOBS[g]; # observed GO frequency in the gene2go file among proteins in the alignment
	    ++pGO[p][g]; ++GOp[g][p]
	}
    }
    ENDFILE{ if(ARGIND==3) {
	# loop over the GO terms in the gene2go file that are associated with proteins in the alignment,
	# increasing their ancestor GOfreqs as well
	for(g in GOfreqOBS) { # loop over GO terms from gene2go that annotate proteins in our alignment
	    #ASSERT(g in GOp, "Oops, GO term "g" is not in GOp, which should not happen");
	    GOfreq[g]+=GOfreqOBS[g] # Use "+=" since it may be non-zero due to being an ancestor of another protein.
	    if(g in OBOancestors) { # there can be old GO terms in gene2go that are not in a newer OBO file
		for(g1 in OBOancestors[g]){
		    GOfreq[g1] += GOfreqOBS[g] # all the ancestor GO terms get a frequency bump from the annotation of g.
		    # For all the proteins that g annotates, record the fact that the ancestor g1 also annotates them.
		    for(p in GOp[g]){
			#ASSERT(p in pC,"Oops, protein "p" is in GOp but not pC")
			++pGO[p][g1]; ++GOp[g1][p]
		    }
		}
	    }
	    #else printf "WARNING: %s not found in OBO file\n", $3 > "/dev/fd/2"
	}
	delete GOfreqOBS # not needed anymore
    }}

    END{
	if(length(pGO) < 0.01 * length(pC)) {
	    print "Warning: Fewer than 1% of the proteins in your alignment have GO annotations. This typically happens\n" \
		"when your alignment files use a different gene/protein naming convention than the gene2go file.\n" \
		"(eg., your alignment uses Uniprot, Ensembl, or Official Names but gene2go uses BioGRID IDs, etc.)" >"/dev/fd/2"
	}
	for(p in pC) {
	    if(!(p in pGO)){pGO[p][0]=1;delete pGO[p][0]} #proteins with no GO terms need pGO[p] explicit empty list.
	    for(g in pGO[p]) sumGOp+=length(pC[p])*K_g(g) # total number of (equivalent) GO terms used in this alignment
	}
	if(RESNIK) {
	    for(line in R){
		Kmax=0;
		Kmin=1e30;
		for(col in CA[line]){Kmin=MIN(Kmin,K_p(CA[line][col]));Kmax=MAX(Kmax,K_p(CA[line][col]))}
		Kweight += Kmin
		Resnik+=R[line]*Kmin
		if(VERBOSE) {
		    printf "%d\t%s\t%s\t%s\t(%g,%g)\n",NAF[line],CA[line][0],CA[line][1],R[line],Kmin,Kmax
		}
	    }
	    printf "%s: numPairs %d numP %d sumGO %d GOcorpus %d Rweight(A) %g WeightedResnik %g\n", ARGV[1], length(CA), length(pGO), sumGOp, length(GOfreq), Kweight, Resnik/Kweight
	}
	else {
	    know=K_AC(CA);
	    printf "%s: numClus %d numP %d sumGO %d GOcorpus %d K(A) %g score %g\n", ARGV[1], length(CA), length(pGO), sumGOp, length(GOfreq), know, know/sumGOp
	}
    }' "$i" $OBO $GENE2GO
done
