#!/bin/sh

ID_FILE=/home/wayne/extra2/preserve/BioGRID/3.5.176/IDENTIFIERS.tab/IDENTIFIERS-3.5.176.tsv
#ID_FILE=/home/wayne/extra2/preserve/BioGRID/3.5.176/IDENTIFIERS.tab/IDENTIFIERS-short.tsv
echo "You will need 40GB of RAM and 5-10 minutes as we load this huge file:"
(cd `dirname $ID_FILE`; ls -l `basename $ID_FILE`)
LINES=`wc -l < $ID_FILE`

INFILE=
if [ $# -eq 0 ]; then INFILE="-"; fi # force awk to read stdin if we are given no arguments

TAB="	"
# The file looks like this:
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
    BEGINFILE{if(FILENAME=="-")print "\nReady to take TAB-SEPARATED 4-column queries (qtype,qname, rtype,rspecies)";else printf "%s:\n0%", FILENAME >"/dev/fd/2"}
    ARGIND==1&&FNR>1{ # skip header line
	$0=tolower($0) # only deal in lower-case
	bg=$1; name=$2; type=$3; species=$4
	TYPE[name][bg]=type
	NAMES[bg][type][species]=name
	percent=int(100*FNR/'$LINES');
	if(percent!=oldPercent){printf "\r%d%%",percent >"/dev/fd/2";oldPercent=percent}
    }
    ARGIND>1{
	if(NF!=4){printf "\nexpecting 4 TAB-SEPARATED fields: queryType, queryName, outType, species. All but name can be a regexp, \".*\" for anything\n" > "/dev/fd/2";next}
	$0=tolower($0)
	qType=$1
	qName=$2
	if(!(qName in TYPE)){printf "sorry, no such symbol found \"%s\"\n",qName > "/dev/fd/2"; next}
	rType=$3
	species=$4
	for(b in TYPE[qName])if(match(TYPE[qName][b],qType)){ # for all BioGRID IDs that have a match to this name,type pair
	    #printf "Symbol \"%s\" of type \"%s\" is BioGRID id \"%s\"\n",qName,TYPE[qName][b],b
	    for(t in NAMES[b])if(match(t,rType))for(s in NAMES[b][t]) if(match(s,species))
	    {
		printf "%s\t%s\t",TYPE[qName][b],qName
		printf "%s\t%s\t%s\n", s,t,NAMES[b][t][s]
	    }
	}
    }' $ID_FILE $INFILE "$@"
