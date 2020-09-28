#include <stdio.h>
#include "misc.h"
#include "sets.h"
#include "rand48.h"

#define MAX_NODES 22500 // used mostly for testing
#define MAX_GO 22500 // used mostly for testing

// These are mutually exclusive
#define TEST_IO 0
#define TEST_SHUFFLE 0
#define EVALUATE 1
#define MUTUALLY_EXCLUSIVE (TEST_IO+TEST_SHUFFLE+EVALUATE)

static int maxGO, numSamples;

typedef struct _NetGO {
    int n; // numNodes
    int *numGO; // number of GO terms for each node
    int **GOarray; // vector of GO terms for each node.
    SET **GOset; // set of GO terms for each node
    int *GOfreq; // array of size maxGO counting total frequency of all GO terms in this network
} NETGO;

NETGO *ReadNetGO(FILE *fp)
{
    NETGO *N = Malloc(sizeof(NETGO));
    char line[100*BUFSIZ];
    int i, G;
    assert(fgets(line, sizeof(line), fp));
    assert(2==sscanf(line,"G%d %d\n", &G, &(N->n)));
    assert(G==1 || G==2);
    assert(N->n>0);

    N->numGO = (int*)Calloc(sizeof(int),N->n);
    N->GOarray = (int**)Calloc(sizeof(int*),N->n);
    N->GOset = (SET**)Calloc(sizeof(SET*),N->n);
    N->GOfreq = (int*)Calloc(sizeof(int),maxGO);

    for(i=0; i<N->n; i++)
    {
	int ii,j;
	if(1!=fscanf(fp,"%d", &ii))Fatal("scanning ii at i=%d",i);
	if(--ii!=i) Fatal("i is %d, ii is %d",i,ii);
	assert(1==fscanf(fp," %d", &(N->numGO[i])));
	N->GOarray[i] = (int*)Calloc(sizeof(int),N->numGO[i]);
	N->GOset[i] = SetAlloc(maxGO);
	for(j=0;j<N->numGO[i];j++) {
	    int GOterm;
	    assert(1==fscanf(fp," %d", &GOterm));
	    --GOterm; // numbered from 1, make it zero.
	    assert(0 <= GOterm && GOterm < maxGO);
	    N->GOarray[i][j] = GOterm;
	    SetAdd(N->GOset[i], GOterm);
	    ++N->GOfreq[GOterm];
	}
	assert(0==fscanf(fp,"\n"));
    }
    return N;
}


void EvalAlignment(NETGO *netGO[], int A[])
{
    // For each GO term g, numShared[g] is number of protein pairs that shared g in the alignment A[].
    // In the combinatorial analysis in the paper, numShared[g] is simply "k".
    static int g, numShared[MAX_GO];
    for(g=0;g<maxGO;g++) numShared[g]=0;
    int n[2], i[2];
    n[0]=netGO[0]->n;
    n[1]=netGO[1]->n;
    for(i[0]=0;i[0]<n[0];i[0]++) // for each node in in the smaller network
    {
	i[1] = A[i[0]];
	int j, nG[2], whoMin = 0, other; // who has fewer GO terms, and who's the other?
	nG[0] = netGO[0]->numGO[i[0]];
	nG[1] = netGO[1]->numGO[i[1]];
	if(nG[0] > nG[1]) whoMin = 1;
	other=1-whoMin;
	for(j=0;j<nG[whoMin];j++)
	{
	    int GOterm = netGO[whoMin]->GOarray[i[whoMin]][j];
	    if(SetIn(netGO[other]->GOset[i[other]], GOterm)) ++numShared[GOterm];
	}
    }

    //printf("A"); for(g=0;g<n[0];g++)printf(" %d",A[g]); putchar('\n');
    for(g=0;g<maxGO;g++) if(numShared[g]>0) // avoid output when k=0 since output becomes HUGE
	    printf("%d %d %d %d\n", g, netGO[0]->GOfreq[g], netGO[1]->GOfreq[g], numShared[g]);
}

int Shuffle(int A[], int n)
{
    int i;
    for(i=0;i<n-1;i++)
    {
	int who = i+1+(n-i-1)*drand48();
	assert(0 < who && who < n); // >0 because we don't swap with ourselves
	int tmp = A[who]; A[who]=A[i]; A[i]=tmp;
	assert(0 <= tmp && tmp < n);
    }
    return n;
}

int main(int argc, char *argv[])
{
    assert(MUTUALLY_EXCLUSIVE == 1);
    FILE *fp = stdin;
    int seed = GetFancySeed(true);
    srand48(seed);

    assert(argc==2);
    numSamples = atoi(argv[1]);
    assert(numSamples > 0);
    
    assert(1==fscanf(fp,"maxGO %d\n", &maxGO));
    assert(0 <= maxGO && maxGO < MAX_GO);
    NETGO *netGO[2];
    netGO[0] = ReadNetGO(fp);
    netGO[1] = ReadNetGO(fp);

#if TEST_IO // verify by outputting exactly what we read, bit for byte
    int G,j,k;
    printf("maxGO %d\n", maxGO);
    for(G=0;G<=1;G++)
    {
	printf("G%d %d\n", G+1, netGO[G]->n);
	for(j=0;j<netGO[G]->n;j++)
	{
	    printf("%d %d", j+1, netGO[G]->numGO[j]);
	    for(k=0;k<netGO[G]->numGO[j];k++) printf(" %d",netGO[G]->GOarray[j][k]+1);
	    puts("");
	}
    }
#endif

    int i, n1=netGO[1]->n, A[n1];
#if TEST_SHUFFLE
    int j, n0 = netGO[0]->n;
    static int AA[MAX_NODES][MAX_NODES];
    assert(0 <= n0 && n0 <= n1 && n1 <= MAX_NODES);
    for(i=0;i<n1;i++)A[i]=i; // start with the identity mapping
    for(i=0;i<numSamples;i++) {
	Shuffle(A, n1);
	for(j=0;j<n1;j++)++AA[j][A[j]];
    }
    for(i=0;i<n1;i++) for(j=0;j<n1;j++)
	printf("%d %d %d\n", i,j, AA[i][j]);
    fprintf(stderr, "Expect about:  %g counts per bin (%d/%d)\n", numSamples*1.0/(n1+1), numSamples,(n1+1));
#endif

#if EVALUATE
    printf("N %d\n",numSamples);
    for(i=0;i<=1;i++) printf("G%d %d\n", i+1, netGO[i]->n);
    for(i=0;i<n1;i++)A[i]=i; // start with the identity mapping
    for(i=0;i<numSamples;i++) {
	Shuffle(A, n1);
	EvalAlignment(netGO, A);
    }
    //PrintSummary(netGO);
#endif
    return 0;
}
