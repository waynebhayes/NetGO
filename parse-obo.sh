#!/bin/sh
awk '{gsub("!.*$","")  # delete all comments
	gsub("GO:","") # when you remove the GO: part it fits in 7 columns so tabs always move just 1 tabstop.
    }
    /^id:/{id=$2;ID[$2]=1}
    /^is_a:/{P[id][$2]=1}
    /^name:/{gsub("name: ","");gsub(" ","_");name[id]=$0}
    /^namespace:/{NS[id]=$2}
    /^alt_id:/||/^consider:/{map[$2][id]=1}
    /^is_obsolete:/{Ob[id]=1}
    /^.Typedef/{exit} # this marks the end of the actual heirarchy definition

    # Print the path to the root, starting at g
    function RootPaths(g,d,	p,i,n){
	printf "%s",g
	if(isarray(P[g])) {
	    n=0;
	    for(p in P[g]){
		if(n>0)for(i=0;i<d;i++)printf "\t"
		printf "\t"
		RootPaths(p,d+1)
		n++
	    }
	} else {
	    print ""
	}
    }

    END{
	for(g in ID) if(!Ob[g]) { # if not obsolete
	    RootPaths(g,0)
	}
    }' "$@"
