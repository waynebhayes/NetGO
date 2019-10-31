#!/bin/sh
KOUNT=0
PRINT=0
case $1 in
-[pP]*) PRINT=1; shift;;
-[kK]*) KOUNT=1; shift;;
esac

awk 'BEGIN{PRINTPATHS='$PRINT';KOUNT='$KOUNT'}
    # Print all paths to the root, starting at g
    function OBORootPaths(g,d,	p,i,n){
	if(d==0) delete _OBOcounts
	else ++_OBOcounts[g]
	if(PRINTPATHS)printf "%s",g
	if(isarray(P[g])) {
	    n=0;
	    for(p in P[g]){
		if(PRINTPATHS){
		    if(n>0)for(i=0;i<d;i++)printf "\t"
		    printf "\t"
		}
		OBORootPaths(p,d+1)
		n++
	    }
	} else if(PRINTPATHS) {
	    print ""
	}
    }

    ARGIND==1{gsub("!.*$","")  # delete all comments
	gsub("GO:","") # when you remove the GO: part it fits in 7 columns so tabs always move just 1 tabstop.
	if(/^id:/){id=$2;ID[$2]=1}
	if(/^is_a:/){P[id][$2]=1}
	if(/^name:/){gsub("name: ","");gsub(" ","_");name[id]=$0}
	if(/^namespace:/){NS[id]=$2}
	if(/^alt_id:/||/^consider:/){map[$2][id]=1}
	if(/^is_obsolete:/){Ob[id]=1}
	if(/^.Typedef/){exit} # this marks the end of the actual heirarchy definition
    }
    END{
	for(g in ID) if(!Ob[g]) { # if not obsolete
	    OBORootPaths(g,0)
	    if(KOUNT){
		for(g1 in _OBOcounts)printf "%d:%s ",_OBOcounts[g1],g1
		print ""
	    }
	    else {
		for(g1 in _OBOcounts)OBOancestors[g][g1]=1
	    }
	}
    }' "$@"
