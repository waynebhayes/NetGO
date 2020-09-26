#!/bin/sh
USAGE="USAGE: $0 {-combine | -eval | -create numSamples {plus the 5 params required by NetGO2int.sh}
    use -create to create the data files
    use -combine to combine a bunch of them into one
    use -eval to evaluate one or a combined set after-the-fact"
die() { (echo "$FATAL ERROR: $@"; echo "$USAGE"; exec ./NetGO2int.sh) >&2; exit 1
}

case "$1" in
-*create) [ "$#" -eq 7 ] || die "Need exactly 7 parameters to create"
    N="$2"
    shift 2
    ./NetGO2int.sh "$@" | ./ra.eval $N | # output is: g l1 l2 k, for each nonzero k in each of N alignments
	exec hawk '/^[NG]/{print;next}NF==4{g=$1;l1=$2;l2=$3;k=$4;++shared[g][l1][l2][k];next}
	{ASSERT(NF==4,"oops, NF!=4 during create");}
	    END{for(g in shared)for(l1 in shared[g])for(l2 in shared[g][l1])for(k in shared[g][l1][l2])
		    if(1*shared[g][l2][l2][k]>0) print g,l1,l2,k,shared[g][l2][l2][k];
	    }'
    ;;
-*combine) shift;
    exec hawk '/^N/{N+=$2;next}
	/^G/{gsub("^G",""); G=$1;
	if(G in n) ASSERT($2==n[G],"Wrong nodeCount in "FILENAME":"$0", previously G"G" had "n[G]);else n[G]=$2; next}
	NF==5{g=$1;l1=$2;l2=$3;k=$4;count=$5; shared[g][l1][l2][k]+=count; next}
	{ASSERT(NF==5,"oops, NF!=5 during combine at "FILENAME" line "FNR":"$0);}
	END{printf "N %d\nG1 %d\nG2 %d\n",N,n[1],n[2];
	    for(g in shared)for(l1 in shared[g])for(l2 in shared[g][l1])for(k in shared[g][l1][l2])
		if(1*shared[g][l1][l2][k]>0) print g,l1,l2,k,shared[g][l2][l2][k];
	}' "$@"
    ;;
-*eval) shift;
    [ $# = 0 ] || die "not expecting any args after '-eval'"
    exec hawk '/^N/{N=$2;next}/^G/{n[++G]=$2;next}
	NF==5 {
	    g=$1; l1=$2; l2=$3; k=$4; count=$5;
	    predict=exp(logCountGOtermAlignments(n[1],n[2],l1,l2,k)-logAlignSearchSpace(n[1],n[2]));
	    actual=count/N;
	    printf "%18s %10s %.9f %.9f %10s %.9f\n",$0,"",predict,actual,"",predict/actual
	    next
	}
	{ASSERT(NF==5,"oops, expecting exactly 5 columns on --eval")}' "$@"
    ;;
*) die "first argument must be '-create' or '-eval'";;
esac
