#include <stdio.h>
#include "misc.h"
#include "sets.h"

static int maxGO;

typedef struct _NetGO {
    int n; // numNodes
    int *numGO; // number of GO terms for each node
    int **GOarray; // vector of GO terms for each node.
    SET **GOset; // set of GO terms for each node
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
	}
	assert(0==fscanf(fp,"\n"));
    }
    return N;
}

int main(int argc, char *argv[])
{
    FILE *fp = stdin;
    assert(1==fscanf(fp,"maxGO %d\n", &maxGO));
    NETGO *netGO[2];
    netGO[0] = ReadNetGO(fp);
    netGO[1] = ReadNetGO(fp);

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
}
