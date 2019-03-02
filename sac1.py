from igraph import *
from scipy import spatial
import pandas as pd
import numpy as np
import sys

adjM=[]
adj_copy=[]

def main(a,alpha):
    global adjM
    global adj_copy
    attrL = pd.read_csv('data/fb_caltech_small_attrlist.csv')
    V = len(attrL)
    with open('data/fb_caltech_small_edgelist.txt') as f:
        edges = f.readlines()
    edges = [tuple([int(x) for x in line.strip().split(" ")]) for line in edges]

    g = Graph()
    g.add_vertices(V)
    g.add_edges(edges)

    for col in attrL.keys():
        g.vs[col] = attrL[col]

    def cosineSimilarity(v1, v2, g):
        vec1 = list(g.vs[v1].attributes().values())
        vec2 = list(g.vs[v2].attributes().values())
        #prod=0
        #for i in range(len(vec1)):
        #    prod = prod + vec1[i]*vec2[i]
        prod=1 - spatial.distance.cosine(vec1, vec2)
        return prod

    adjM = np.zeros((V, V))
    for i in range(V):
        for j in range(V):
            adjM[i][j] = cosineSimilarity(i, j, g)
    adj_copy = np.array(adjM)

    def uniform(C):
        uniformC = []
        d = {}
        c = 0
        for assgnmnt in C:
            if assgnmnt in d:
                uniformC.append(d[assgnmnt])
            else:
                uniformC.append(c)
                d[assgnmnt] = c
                c = c + 1
        return uniformC

    def cM(g, C):
        C1 = list(Clustering(C))
        a = 0.0
        for c in C1:
            b = 0.0
            for v1 in c:
                for v2 in C:
                    if (v1 != v2):
                        b = b + adjM[v1][v2]
        b=b/len(c)
        a=a+b
        x=a/(len(set(C))) + g.modularity(C)
        return x

    def contractGraph(a, alpha):

        def phase_1(g, alpha, C):
            V = len(g.vs)
            i=0
            flag = 0

            def Qgain(alpha, C, g, vx, vy):
                Qg1 = g.modularity(C)
                t = C[vx]
                C[vx] = vy
                Qg2 = g.modularity(C)
                C[vx] = t
                z1 = Qg2 - Qg1
                S = 0.0
                indices = [i for i, x in enumerate(C) if x == vy]
                for v in indices:
                    S = S + adjM[vx][v]
                z2 = S / (len(indices) * len(set(C)))

                return (alpha * z1) + ((1 - alpha) * z2)

            Cz = range(V)
            while (flag == 0 and i < 15):
                flag = 1
                for vi in Cz:
                    v_ = -1
                    maxgain = 0.0
                    cc = list(set(C))
                    for index in range(len(cc)):
                        if (C[vi] != cc[index]):
                            dQ = Qgain(alpha, C, g, vi, cc[index])
                            if (dQ > maxgain):
                                maxgain = dQ
                                v_ = cc[index]
                    if (maxgain > 0.0 and v_ != -1):
                        flag = 0
                        C[vi] = v_
                i += 1

            return C

        def phase_2(g, C):
            C1 = uniform(C)
            temp = list(Clustering(C1))
            x=list(set(C))
            l=len(x)
            adjM = np.zeros((l, l))
            for i in range(l):
                for j in range(l):
                    s=0.0
                    for k in temp[i]:
                        for f in temp[j]:
                            s=s+adj_copy[k][f]
                    adjM[i][j]=s
            g.contract_vertices(C1)
            g.simplify(combine_edges=sum)
            return

        V = g.vcount()
        print(V)
        C = phase_1(g, alpha, list(range(V)))
        print('Number of Communities after Phase 1')
        print(len(set(C)))
        C = uniform(C)
        m1 = cM(g, C)
        phase_2(g, C)
        V = g.vcount()
        C2 = phase_1(g, alpha, list(range(V)))
        C2x = uniform(C2)
        c_p_2 = list(Clustering(C2x))
        m2 = cM(g, C)
        C1x = uniform(C)
        c_p_1 = list(Clustering(C1x))

        return m1,m2,c_p_1,c_p_2

    def createFile(clusters):
        file = open("communities_" + str(a[alpha]) + ".txt", 'w+')
        for c in clusters:
            for i in range(len(c)):
                file.write("%s," % c[i])
            file.write('\n')
        file.close()

    m1, m2, c1, c2 = contractGraph(a, alpha)
    if (m1 > m2):
        createFile(c1)
        print('Phase 1 clusters have higher modularity')
        return c1
    else:
        createFile(c2)
        print('Phase 2 clusters have higher modularity')
        return c2

if __name__ == "__main__":
    k = float(sys.argv[1])
    fileMap={0:0, 1:1, 0.5:5}

    main(fileMap, k)

