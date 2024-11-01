ls "$1" | # list old directory containing names
    grep _ | sed 's/_/ /g' | # get names with spaces
    while read name; do
	echo "$name" >&2;
	fgrep "$name" names.dmp | # egrep below remove entries that are not the actual species or a variant
	    egrep -v " x |virus|expression|vector|fusant|killer|particle| DNA |library|hybrid|synthetic|(ex|of) $name| and |graft" |
	    cut -f1 | # get the taxonomic ID
	    sort -un | # sort them numerically--the lowest one is the "main" one
	    tee "$name"
    done
