#!/bin/sh
DATABASE_DIR=/home/wayne/extra1/preserve/BioGRID/raw/4.4.209/IDENTIFIERS.tab
USAGE="PURPOSE: translate between the (bloody annoying) gazillions of naming schemes for genes and proteins,
by leveraging BiGRID's database of these names.

Do an 'ls' on $DATABASE_DIR for the entire list of species; some examples are
	`cd $DATABASE_DIR && ls Hom* Sac* Rat* Mus_* D*_mel* A*_th* S*_pombe | fmt -120`

Commonly used types include: ENTREZ_GENE    OFFICIAL_SYMBOL    SWISS-PROT    SYSTEMATIC_NAME    UNIPROT-ACCESSION.

For example, the SANA/networks/HSapiens.el edge list in the SANA repo uses ENTREZ_GENE names; to convert to SWISS-PROT, use:

    gene-name-translator.sh Homo_sapiens ENTREZ_GENE SWISS-PROT SANA-repo/networks/HSapiens.el

USAGE: $0 [-ID IDENTIFIER_FILE] species queryType outType [InputFile]
    Any of them (except -ID argument) can be regular expressions."

die() { echo "$USAGE">&2; echo "FATAL ERROR:$@" >&2; exit 1
}
[ $# -gt 0 ] || die "Expecting at least 3 arguments"
case "$1" in
-ID) SPECIES_FILE="$2"; shift 2;;
*)
    [ -d $DATABASE_DIR ] || die "set DATABASE_DIR in $0"

    InputFile=""
    SPECIES_FILE=`ls $DATABASE_DIR/* | grep -i "$1"`
    [ `echo "$SPECIES_FILE" | wc -l` -eq 1 ] || die "$SPECIES_FILE:
    expecting only one species file, but found the above files. Make your species regexp more strict."
    ;;
esac
[ -f "$SPECIES_FILE" ] || die "cannot find species file '$SPECIES_FILE'"
qType="$2"
rType="$3"
[ "$qType" = "" -o "$rType" = "" ] && die "qType and rType cannot be empty"
shift 3
if [ $# -eq 0 ]; then
    InputFile="-"
fi # force awk to read stdin if we are given no arguments
TTY=0
if isatty; then TTY=1;fi

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
    BEGIN{qType="'"$qType"'";rType="'"$rType"'"; IGNORECASE=1; if('$TTY'&&ARGV[2]=="-")print "Ready to take name queries"}
    ARGIND==1&&FNR>1{ # skip header line
	bg=$1; name=$2; type=$3
	TYPE[toupper(name)][bg]=type # force UPPERCASE for array indices since IGNORECASE has no effect on them.
	NAMES[bg][toupper(type)]=name
    }
    ARGIND>1{
	delete outCols;
	for(col=1;col<=NF;col++) {
	    qName=$col; QNAME=toupper(qName);
	    if(!(QNAME in TYPE)){
		printf "No symbol \"%s\" of type \"'"$qType"'\" in '"$SPECIES_FILE"'\n",qName > "/dev/fd/2"
		printf "UNKNOWN_SOURCE_TYPE\t%s\tUNKNOWN_DEST_TYPE\tUNKNOWN\n", qName
		next
	    }
	    for(b in TYPE[QNAME])if(match(TYPE[QNAME][b],qType)){ # for all BioGRID IDs that have a match to this name,type pair
		#printf "Symbol \"%s\" of type \"%s\" is BioGRID id \"%s\"\n",qName,TYPE[qName][b],b
		for(t in NAMES[b])if(match(t,rType))
		    if(NF==1) {
			printf "%s\t%s\t",TYPE[QNAME][b],qName
			printf "%s\t%s\n", t,NAMES[b][t]
		    } else {
			if(col in outCols) outCols[col]=outCols[col]"|"
			outCols[col]=outCols[col]NAMES[b][t]
		    }
	    }
	    if(!(col in outCols)) outCols[col]=sprintf("NO_MATCH_%s",qName)
	}
	if(NF>1){
	    printf "%s",outCols[1]
	    for(col=2;col<=NF;col++)
		printf "\t%s",outCols[col]
	    print ""
	}
    }' $SPECIES_FILE $InputFile "$@"
