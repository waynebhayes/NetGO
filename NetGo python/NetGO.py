#!/usr/bin/python
import sys
from _collections import defaultdict
class NetGO:
    def __init__(self):
        self.DRACONIAN = True
        self.VERBOSE = False
    def newClusters(self, alignFile):
        self.pC = defaultdict(dict)
        self.CA = defaultdict(list)
        self.encounteredProteins = set()
        with open(alignFile) as alignFileOpen:
            self.get_pC_CA(alignFileOpen)
    def newg2g(self, g2gFile):
        self.GOp = defaultdict(dict)
        self.pGO = defaultdict(dict)
        self.GOfreq = defaultdict(int)
        with open(g2gFile) as g2gFileOpen:
            self.get_pGO_GOp(g2gFileOpen)
    def SetIntersect(self, T1, T2):
        '''
        Takes pGO[protein], so a dict of GO terms, 1 meaning the protein
        has been annotated.
        Returns a dict of GO terms as well.
        '''
        out = {}
        if T1.size()>T2.size():
            for g in T2:
                if g in T1:
                    out[g]=1
        if T1.size()<=T2.size():
            for g in T1:
                if g in T2:
                    out[g]=1
        return out
    
    def K_g(self, g):
        '''
        Takes a GO term as a string, returns K(g) as int
        '''
        if(g in self.GOp):
            return float(1)/float(len(self.GOp[g]))
        else:
            return 0
    def K_gset(self, T):
        '''
        Takes pGO[p], a dict of GO terms with the names as keys, returns K(T)
        '''
        out=0
        for g in T:
            out+=self.K_g(g)
        return out
    def K_p(self, p):
        '''
        Takes a protein name as a string, returns K(p)
        '''
        if p in self.pGO:
            return self.k_gset(self.pGO[p])
        else:
            return 0
    def K_A2(self, A):
        '''
        Takes a pairwise alignment in the format A[u]=v, where u and v are pairs of proteins.
        Returns K(A)
        '''
        out = 0
        for u in A:
            if u in self.pGO:
                v = A[u]
                if v in self.pGO:
                    out+=self.K_gset(self.SetIntersect(self.pGO[u], self.pGO[v]))
        return out
    def K_AC(self, C):
        '''
        C is a dict of clusters in the format {cl:
        '''
        out = 0
        for cl in C:
            if self.VERBOSE:
                print("cluster ", cl, " num of proteins ", len(C[cl]))
            #K_C is the K(C) value for each cluster cl in C
            K_C = 0
            #M is a dict of proteins of the following format: {protein:number of times occurring in cl...}
            #T is a dict of GO terms of the following format: {GO term:number of times occurring in cl... }
            #both are across cluster cl
            M, T = defaultdict(int), defaultdict(int)
            numClusterFields=len(C[cl])
            #only calculate K_C if cluster cl has more than one protein
            if numClusterFields>1:
                '''
                u=C[cl][0]
                #if the first protein, u, is not a placeholder, then add it to M
                if u!="_" or u!="NA":
                    M[u]+=1
                    #then, for each GO term annotating protein u, add it to T
                    for g in self.pGO[u]:
                        T[g]+=1
                '''
                #Now, we iterate over all proteins in cluster cl (skipping the first one, which we already processed). Add them to M.
                for i in range(0, numClusterFields):
                    u=C[cl][i]
                    assert(u!="-" and u!="_" and u!="NA"),"INTERNAL ERROR: invalid protein got into K_AC"
                    #get first protein from cluster, add its GOterms to T 
                    if i==0:
                        u=C[cl][0]
                        M[u]+=1
                        for g in self.pGO[u]:
                            T[g]+=1
                    else:
                        M[u]+=1
                        #If we are in draconian mode:
                        #Check that, for each protein u in M, each GO term g in T annotates u.
                        #If the Go term g does not annotate any one of the proteins, then remove it from T.
                        if self.DRACONIAN:
                            Tremove = set()
                            for g in T:
                                if g not in self.pGO[u]:
                                    Tremove.add(g)
                            for g in Tremove:
                                T.pop(g)
                        #If we are not in draconian mode:
                        #Increment g's entry in T by one.
                        else:
                            for g in T:
                                T[g]+=1
                    if self.VERBOSE:
                        sys.stdout.write("\t%s (%d GOs { " % (u,len(T)))
                        for g in T:
                            sys.stdout.write("%s(%d) "%(g,self.GOfreq[g]))
                        sys.stdout.write("} K(%s)=%g)\n" % (u,self.K_gset(T)))
                    #If cluster cl has any annotations, and either more than one protein, or one protein that occurs more than once...
                if len(T)>0 and (len(M)>1 or (len(M)==1 and M[u]>1)):
                    if self.DRACONIAN:
                        #If we are in draconian mode, just add K_gset(T)
                        K_C+=self.K_gset(T)
                    else:
                        #If we are not in draconian mode, weed out any GO terms that annotate 1 or 0 proteins.
                        for g in T:
                            if T[g]>1:
                                K_C+=T[g]*self.K_g(g)/numClusterFields
            if self.VERBOSE:
                sys.stdout.write("\tClusterCommonGOs {")
                for g in T:
                    sys.stdout.write(g)
                sys.stdout.write("}, K_C=%g" % (K_C))
                sys.stdout.write("\n")
                sys.stdout.flush()
            out+=K_C
        return out
    def sim_A2(self, A):
        return self.K_A2(A/len(self.GOp))
    def get_pC_CA(self, alignFile):
        lineCount = 0
        for cluster in alignFile:
            lineCount+=1
            for protein in cluster.split("\t"):
                if protein!="_" and protein!="NA" and protein!="-":
                    #set pC
                    if lineCount in self.pC[protein.strip()]:
                        self.pC[protein.strip()][lineCount]+=1
                    else:
                        self.pC[protein.strip()][lineCount]=1
                    #set CA
                    self.CA[lineCount].append(protein.strip())
                    self.encounteredProteins.add(protein.strip())
    def get_pGO_GOp(self, g2gFile):
        next(g2gFile)
        for line in g2gFile:
            protein, GOterm = line.split("\t")[1:3]
            self.GOfreq[GOterm]+=1
            if protein in self.encounteredProteins:
                self.pGO[protein][GOterm] = 1
                if protein not in self.GOp[GOterm]:
                    self.GOp[GOterm][protein]=1
                else:
                    self.GOp[GOterm][protein]+=1
if __name__ == "__main__":
    s = NetGO()
    assert(len(sys.argv)>=3),"USAGE=USAGE: $0 [-L] gene2goFile alignFile[s] -L: 'Lenient'. The default behavior is what we call 'Dracanion' in the paper, which insists that a GO term must annotate every protein in a cluster for it to count. The Lenient option gives a GO term a weight per-cluster that is scaled by the number of proteins it annotates (so long as it's more than 1). alignFile: each line consists of a cluster of any number of proteins; proteins can appear in more than one cluster. Note it must use the same protein naming convention as your gene2go file. gene2goFile: a standard-format gene2go file downloaded from the GO consortium's website. For now we only use columns 2 (protein name) and 3 (GO term). We ignore the species, which is a simplification since about 20% of proteins appear in two species, but only 2% appear in more than 2."
    if sys.argv[1]=="-L":
        startVal=3
        geneFileIndex=2
        s.DRACONIAN=False
    elif sys.argv[1]=="-verbose" and sys.argv[1]!="-verbose":
        startVal=3
        geneFileIndex=2
        s.VERBOSE=True
    elif sys.argv[2]=="-L" and sys.argv[1]=="-verbose":
        startVal=4
        geneFileIndex=3
        s.DRACONIAN=False
        s.VERBOSE=True
    else:
        startVal=2
        geneFileIndex=1
    for i in range(startVal, len(sys.argv)):
        s.newClusters(sys.argv[i])
        s.newg2g(sys.argv[geneFileIndex])
        assert(len(s.pGO)>=0.01*len(s.pC)), "Warning: Fewer than 0.1% of the proteins in your alignment have GO annotations. This typically happens\n when your alignment files use a different gene/protein naming convention than the gene2go file. '/dev/fd/2'"
        sumGOp = 0
        count=0
        for p in s.pC:
            if not p in s.pGO:
                s.pGO[p].clear()
                count+=1
            for g in s.pGO[p]:
                sumGOp = sumGOp + len(s.pC[p])*s.K_g(g)
        know = s.K_AC(s.CA)
        print(sys.argv[i],
                ": numClus ", len(s.CA),
                " numP ", len(s.pGO),
                " numGO ", len(s.GOp),
                " GOcorpus ", len(s.GOfreq),
                " k ", know,
                " score ", float(know)/float(sumGOp))






















                    
