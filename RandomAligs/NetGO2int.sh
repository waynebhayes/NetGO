#!/bin/sh
USAGE="$0 species1Name species1Name G1.el G2.el gene2go
Purpose: Integer-ize the node and GO term info (no edges) for usage in C program to generate zillions of random alignments"
die() { (echo "$USAGE"; echo "FATAL ERROR: $@")>&2; exit 1
}
[ $# -eq 5 ] || die "Expecting exactly 5 arguments"
tax1=`BioGRIDname $1 | cut -f3`
tax2=`BioGRIDname $2 | cut -f3`
G1="$3"; [ -f "$G1" ] || die "can't find file '$G1'"
G2="$4"; [ -f "$G2" ] || die "can't find file '$G2'"
gene2GO="$5"; [ -f "$5" ] || die "can't find file '$5'"

grep -l "^$tax1	" "$gene2GO" >/dev/null || die "$gene2GO file doesn't contain taxonomic ID '$tax1'"
grep -l "^$tax2	" "$gene2GO" >/dev/null || die "$gene2GO file doesn't contain taxonomic ID '$tax2'"

hawk 'BEGIN{
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
	if(!($3 in GO2int)){GO2int[$3]=++maxGO; GO[maxGO]=$3}
	++GpGO[G][pName2int[G][$2]][GO2int[$3]]
    }
}
END{
    Gmap[1]=1; Gmap[2]=2;
    if(maxV[1]>maxV[2]){Gmap[1]=2;Gmap[2]=1}
    for(G=1;G<=2;G++){
	printf "G%d %d\n",G,maxV[Gmap[G]];
	for(i=1;i<=maxV[Gmap[G]];i++){
	    ASSERT(i==pName2int[Gmap[G]][pName[Gmap[G]][i]]);
	    printf "%d",i;
	    if(i in GpGO[Gmap[G]])for(j in GpGO[Gmap[G]][i])printf " %d",j;
	    print ""
	}
    }
}' "$G1" "$G2" "$gene2GO"

