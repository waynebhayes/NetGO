#!/bin/sh
USAGE="USAGE: $0 {-combine | -eval | -create numSamples {plus the 5 params required by NetGO2int.sh}
    use -create to create the data files
    use -combine to combine a bunch of them into one
    use -eval to evaluate one or a combined set after-the-fact"
die() { (echo "$FATAL ERROR: $@"; echo "$USAGE"; exec ./NetGO2int.sh) >&2; exit 1
}

# Note: PHI = empirically observed count for the tuplet [n1 n2 g l1 l2 k]; phi (scalar) is used when reading input
case "$1" in
-*create) [ "$#" -eq 7 ] || die "Need exactly 7 parameters to create"
    N="$2"
    shift 2
    ./NetGO2int.sh "$@" | ./ra.eval $N | # output is: g l1 l2 k, for each nonzero k in each of N alignments
	exec hawk '/^[NG]/{print;next}
	NF==4{g=$1;l1=$2;l2=$3;k=$4;++PHI[g][l1][l2][k];next}
	{ # we only get here if NF!=4 so at this point assertion failure is guaranteed
	    ASSERT(NF==4,"oops, NF!=4 during create, line "FNR": "$0);
	}
	END{for(g in PHI)for(l1 in PHI[g])for(l2 in PHI[g][l1])for(k in PHI[g][l1][l2])
		if(1*PHI[g][l1][l2][k]>0) # if the frequency of (g,l1,l2,k) > 0, print it
		    print g,l1,l2,k,PHI[g][l1][l2][k];
	}'
    ;;
-*combine) shift; # take a bunch of the above files of random sample alignments and merge them into one larger sample file.
    exec hawk '/^N/{N+=$2;next}
	/^G/{gsub("^G",""); # remove the letter "G" to extract the integer (1 or 2) after it
	    G=$1;
	    if(G in n)
		ASSERT($2==n[G],"Wrong nodeCount in "FILENAME":"$0", previously G"G" had "n[G]);
	    else
		n[G]=$2;
	    next # skip to the next line of input
	}
	NF==5{ # we expect exactly 5 columns, otherwise error
	    g=$1;l1=$2;l2=$3;k=$4;phi=$5; # phi = frequency we saw g,l1,l2,k (n1,n2 are implicit)
	    PHI[g][l1][l2][k]+=phi;
	    next # skip to the next line to avoid wasting CPU on the ASSERT below
	}
	{ # we only get here if NF!=5, so once we get here we are guaranteed an assertion failure
	    ASSERT(NF==5,"oops, NF!=5 during combine at "FILENAME" line "FNR":"$0);
	}
	END{printf "N %d\nG1 %d\nG2 %d\n",N,n[1],n[2];
	    for(g in PHI)for(l1 in PHI[g])for(l2 in PHI[g][l1])for(k in PHI[g][l1][l2])
		if(1*PHI[g][l1][l2][k]>0) print g,l1,l2,k,PHI[g][l1][l2][k];
	}' "$@"
    ;;
-*eval) shift; # take ONE sample file on the standard input, and evaluate theory vs. experiment
    [ $# = 0 ] || die "not expecting any args after '-eval'"
    exec hawk '/^N/{N=$2;next}/^G/{n[++G]=$2;next}
	NF==5 {
	    g=$1; l1=$2; l2=$3; k=$4; phi=$5;
	    predict=exp(logCountGOtermAlignments(n[1],n[2],l1,l2,k)-logAlignSearchSpace(n[1],n[2]));
	    observed=phi/N;
	    printf "%28s\t%.9f %.9f %10s %.9f\n",
		$0, # print the whole input line for debugging purposes
		predict, observed, "", # print the two values, plus more whitespace
		predict/observed # finally, print the ratio of the two
	    next
	}
	{ASSERT(NF==5,"oops, expecting exactly 5 columns on --eval")}' "$@"
    ;;
*) die "first argument must be '-create' or '-eval'";;
esac
