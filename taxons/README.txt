To get the most recent list of taxonomic IDs, visit
    https://ftp.ncbi.nih.gov/pub/taxonomy
and download "taxdmp.zip". It will contain the file names.dmp. Then take any existing taxon directory (in which the names
are separated by _ rather than spaces), and run the following script ./get-new-taxons.sh.

Then put the resulting files into a new taxons.DATE directory, and to create symbolic links to the names I sometimes use
(eg "Homo sapiens" -> HSapiens), run the following
command line inside the new taxon directory:

ls | sed 's/_/ /' -e 's/ \(.\)/ \1 /' | awk '{F=$1"_"$2$3; new=sprintf("%s%s%s",substr($1,1,1),toupper($2),$3); printf "ln -s \"%s\" %s\n", F,new}'
