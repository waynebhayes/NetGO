#/bin/sh
################## SKELETON: DO NOT TOUCH THESE 2 LINES
BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME [-h] [data file]
PURPOSE: compute holistic p-value of possibly correlated variates using Empirical Brown's Method
    input: m+1 lines, with n+2 columns each.
    The top line is a "header" line; this line is completely ignored, and can be whatever you want
    Each line represents a variable; the first entry is the 'name' (whatever you want, but must be non-empty string); then the p-value; and then the n raw samples of that variable.
    The Empirical Brown's Method (Poole 2016) computes the covariance of all the variables and spits out a
    holistic p-value that is <= product of p-values, accounting approximately for inter-dependencies."

################## SKELETON: DO NOT TOUCH CODE HERE
# check that you really did add a usage message above
USAGE=${USAGE:?"$0 should have a USAGE message before sourcing skel.sh"}
die(){ echo "$USAGE${NL}FATAL ERROR in $BASENAME:" "$@" >&2; exit 1; }
[ "$BASENAME" == skel ] && die "$0 is a skeleton Bourne Shell script; your scripts should source it, not run it"
echo "$BASENAME" | grep "[ $TAB]" && die "Shell script names really REALLY shouldn't contain spaces or tabs"
[ $BASENAME == "$BASENAME" ] || die "something weird with filename in '$BASENAME'"
warn(){ (echo "WARNING: $@")>&2; }
not(){ if eval "$@"; then return 1; else return 0; fi; }
newlines(){ awk '{for(i=1; i<=NF;i++)print $i}' "$@"; }
parse(){ awk "BEGIN{print $*}" </dev/null; }

# Temporary Filename + Directory (both, you can use either, note they'll have different random stuff in the XXXXXX part)
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

#################### END OF SKELETON, ADD YOUR CODE BELOW THIS LINE

while true; do
    case "$1" in
    -h) die "<printing help message only>" ;;
    -*) die "option '$1' not supported";;
    *) break;;
    esac
done

tail -n +2 "$@" > $TMPDIR/input

hawk 'function TransformData(     j) {
	for(j=1;j<=NF;j++) {
	    data[NR][j] = -2*log(StatHistECDF(NR, data[NR][j]));
	}
    }
	NR==1{n=NF-2}
	{ # every line including NR==1
	    ASSERT(NF-2==n,"number of columns must be constant; first column had "n" but this one has "NF-2);
	    name[NR]=$1; pVal[NR]=$2;
	    for(i=1;i<=n;i++) {
		data[NR][i]=sample=$(i+2)
		StatAddSample(NR,sample);
	    }
	    #printf "DATA[%d]",NR>"/dev/stderr";for(i=1;i<=n;i++)printf " %g",data[NR][i]>"/dev/stderr";print "">"/dev/stderr";
	    for(i=1;i<=n;i++){
		data[NR][i] = (data[NR][i] - StatMean(NR))/StatStdDev(NR);
		StatHistAddSample(NR, data[NR][i]);
	    }
	    #printf "NORM[%d]",NR>"/dev/stderr";for(i=1;i<=n;i++)printf " %g",data[NR][i]>"/dev/stderr";print "">"/dev/stderr";
	    StatHistMakeCDF(NR);
	    #printf "ECDF[%d]",NR>"/dev/stderr";for(i=1;i<=n;i++)printf " [%d][%d][%g]=%g",NR,i,data[NR][i],StatHistECDF(NR,data[NR][i])>"/dev/stderr";print "">"/dev/stderr";
	    TransformData(); # no arguments, it always works on the current line
	    #printf "-LOG[%d]",NR>"/dev/stderr";for(i=1;i<=n;i++)printf " %g",data[NR][i]>"/dev/stderr";print "">"/dev/stderr";
	}
    END{
	m=NR
	# Compute covariances across all samples of all input variables.
	for(i=1;i<=m;i++) for(j=i+1;j<=m;j++)for(k=1;k<=n;k++) CovarAddSample(i" "j, data[i][k], data[j][k]);
	# Now perform Empirical Browns Method
	df_fisher = Expected = 2.0*m;
	cov_sum = 0;
	for(i=1;i<=m;i++) { 
	    #printf "row %4d", i > "/dev/stderr"
	    for(j=i+1;j<=m;j++) {
		covar=CovarCompute(i" "j);
		#printf " %.5f", covar > "/dev/stderr";
		cov_sum += covar;
	    }
	    #print "" > "/dev/stderr"
	}
	#printf "cov_sum %g\n",cov_sum > "/dev/stderr";
	Var = 4.0*m+2*cov_sum;
	c = Var/(2.0*Expected);
	df_brown = 2.0*Expected**2/Var;
	if(df_brown > df_fisher) {
	    df_brown = df_fisher
	    c = 1.0
	}
	#printf "c = %g\n", c > "/dev/stderr"
	pProd=1.0;
	x=0; # twice the sum of logs of p-values
	for(i=1;i<=NR;i++) {x += -log(pVal[i]); pProd *= pVal[i]}
	x *= 2;
	#printf("x = %g\n", x) > "/dev/stderr";
	p_brown = exp(logChi2_pair(int(df_brown+0.5), 1.0*x/c))

	ASSERT(p_brown >= pProd,
	    "Oops, something wrong: p_brown should be < product(pVals), but p_brown ="p_brown" product ="pProd);
	printf("%g %g (2nd is product)\n", p_brown, pProd);
    }' $TMPDIR/input
