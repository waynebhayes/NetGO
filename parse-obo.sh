#!/bin/sh
awk '{gsub("!.*$","")} # delete all comments
    /^id:/{id=$2;ID[$2]=1}
    /^is_a:/{P[id][$2]=1}
    /^name:/{gsub("name: ","");gsub(" ","_");name[id]=$0}
    /^namespace:/{NS[id]=$2}
    /^alt_id:/||/^consider:/{map[$2][id]=1}
    /^is_obsolete:/{Ob[id]=1}
    /^.Typedef/{exit} # this marks the end of the actual heirarchy definition

    # Print the path STARTING at the root, down to (but not including) g
    function PrintPaths(g, p,d){
	if(isarray(P[g]))for(p in P[g]){
	    d=PrintPaths(p)
	    printf "%s(%d)\t",p,d
	    return d+1
	} else return 0
    }

    END{
	for(g in ID) if(!Ob[g]) { # if not obsolete
	    d=PrintPaths(g)
	    printf "%s(%d)\n",g,d
	}
    }' "$@"
