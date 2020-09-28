#!/bin/sh
USAGE="$0 species1Name species1Name G1.el G2.el gene2go
Purpose: Integer-ize the node and GO term info (no edges) for usage in C program to generate zillions of random alignments
Output: to the standand output, a file of format below.  NOTE: on output, nodes and GO terms are numbered 1..n, not 0..n-1.
maxGO   # total number of GO terms across both networks
G1 n1   # the string G1 followed by number of nodes in G1
1 numGO1 <list of GO terms>
2 numGO2 <list of GO terms>
...
n1 numGOn2 <list of GO termS>
G2 n2
1 numGO1 <list of GO terms>
...
n2 numGOn2 ..."
die() { (echo "$USAGE"; echo "FATAL ERROR: $@")>&2; exit 1
}
[ $# -eq 5 ] || die "missing arguments"
gene2GO="$5"; [ -f "$5" ] || die "can't find file '$5'"

tax1=`if echo "$1" | grep '^[0-9]*$' >/dev/null; then echo $1; else BioGRIDname $1 | cut -f3- | newlines | while read t; do if fgrep -x "$t" "$gene2GO.taxIDs"; then break; fi; done; fi`
[ -n "$tax1" -a "$tax1" -gt 0 ] || die "Can't find taxonomic ID for species '$1'"
tax2=`if echo "$2" | grep '^[0-9]*$' >/dev/null; then echo $2; else BioGRIDname $2 | cut -f3- | newlines | while read t; do if fgrep -x "$t" "$gene2GO.taxIDs"; then break; fi; done; fi`
[ -n "$tax2" -a "$tax2" -gt 0 ] || die "Can't find taxonomic ID for species '$2'"

G1="$3"; [ -f "$G1" ] || die "can't find file '$G1'"
G2="$4"; [ -f "$G2" ] || die "can't find file '$G2'"

exec hawk 'BEGIN{
    tax2G['$tax1']=1;tax2G['$tax2']=2;
    for(i=1;i<=2;i++){
	pName2int[i][-1]=0;
	delete pName2int[i][-1];
    }
}

# Edgelist files--but we only extract the nodes.
ARGIND==1||ARGIND==2{
    for(i=1;i<=2;i++)if(!($i in pName2int[ARGIND])){
	pName2int[ARGIND][$i]=++maxV[ARGIND];
	pName[ARGIND][maxV[ARGIND]]=$i;
    }
}

ARGIND==3 && ($1 in tax2G) {
    G=tax2G[$1];
    if($2 in pName2int[G]) {
	nodeInt=pName2int[G][$2];
	if(!($3 in GO2int)){GO2int[$3]=++maxGO; GO[maxGO]=$3}
	++GpGO[G][nodeInt][GO2int[$3]]
    }
}
END{
    Gmap[1]=1; Gmap[2]=2;
    if(maxV[1]>maxV[2]){Gmap[1]=2;Gmap[2]=1}
    printf "maxGO %d\n", maxGO
    for(G=1;G<=2;G++){
	printf "G%d %d\n",G, maxV[Gmap[G]]
	for(i=1;i<=maxV[Gmap[G]];i++){
	    ASSERT(i==pName2int[Gmap[G]][pName[Gmap[G]][i]]);
	    printf "%d",i;
	    if(i in GpGO[Gmap[G]]){
		printf " %d", length(GpGO[Gmap[G]][i]);
		for(j in GpGO[Gmap[G]][i])printf " %d", j;
	    } else
		printf " 0";
	    print ""
	}
    }
}' "$G1" "$G2" "$gene2GO"

