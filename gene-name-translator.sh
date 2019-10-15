#!/bin/sh
USAGE="USAGE: $0 [speciesRegexp] [InputFile] (either can be empty for stdin)"
die() { echo "$USAGE">&2; echo "$@" >&2; exit 1
}

DATABASE_DIR=/home/wayne/extra2/preserve/BioGRID/3.5.176/IDENTIFIERS.tab
[ -d $DATABASE_DIR ] || die "set DATABASE_DIR in $0"

InputFile=""
SPECIES_FILE=`ls $DATABASE_DIR/* | grep -i "$1"`
[ `echo "$SPECIES_FILE" | wc -l` -eq 1 ] || die "$SPECIES_FILE:
expecting only one species file, but found the above files. Make your regexp more strict."
[ -f "$SPECIES_FILE" ] || die "cannot find species file '$SPECIES_FILE'"
shift
if [ $# -eq 0 ]; then
    InputFile="-"
fi # force awk to read stdin if we are given no arguments

TAB="	"
# The IDENTIFIER file looks like this (species column missing if we use species-specific files)
#id	value	type	species
#1	1	BIOGRID	Arabidopsis thaliana
#1	ArthMr001	SYSTEMATIC NAME	Arabidopsis thaliana
#1	ArthMr001	ORDERED LOCUS	Arabidopsis thaliana
#1	rrn26	OFFICIAL SYMBOL	Arabidopsis thaliana
#1	814566	ENTREZ_GENE	Arabidopsis thaliana
#1	ETG814566	ENTREZ_GENE_ETG	Arabidopsis thaliana
#2	2	BIOGRID	Arabidopsis thaliana
#2	ArthMp007	SYSTEMATIC NAME	Arabidopsis thaliana
#2	ArthMp007	ORDERED LOCUS	Arabidopsis thaliana
awk -F"$TAB" '
    BEGIN{if(FILENAME=="-")print "\nReady to take TAB-SEPARATED 4-column queries (qtype,qname, rtype,rspecies)"}
    ARGIND==1&&FNR>1{ # skip header line
	$0=tolower($0) # only deal in lower-case
	bg=$1; name=$2; type=$3
	TYPE[name][bg]=type
	NAMES[bg][type]=name
    }
    ARGIND>1{
	if(NF!=3){printf "\nexpecting 3 TAB-SEPARATED fields: queryType, queryName, outType. All but name can be a regexp, \".*\" for anything\n" > "/dev/fd/2";next}
	$0=tolower($0)
	qType=$1
	qName=$2
	if(!(qName in TYPE)){printf "sorry, no such symbol found \"%s\"\n",qName > "/dev/fd/2"; next}
	rType=$3
	for(b in TYPE[qName])if(match(TYPE[qName][b],qType)){ # for all BioGRID IDs that have a match to this name,type pair
	    #printf "Symbol \"%s\" of type \"%s\" is BioGRID id \"%s\"\n",qName,TYPE[qName][b],b
	    for(t in NAMES[b])if(match(t,rType))
	    {
		printf "%s\t%s\t",TYPE[qName][b],qName
		printf "%s\t%s\t%s\n", s,t,NAMES[b][t]
	    }
	}
    }' $SPECIES_FILE $InputFile "$@"
