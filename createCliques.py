#!/usr/bin/env python3

# Dependent imports
from Bio import SeqIO
from joblib import Parallel, delayed
from tqdm import tqdm
import igraph

# Default imports
import argparse
import pathlib
import time

# Add option to include a list of annotations
# Or just process them straight from the fasta?
# Only needs so can label containing cliques by annot ID

def get_arguments():
    # Input
    parser = argparse.ArgumentParser()
    parser.add_argument('-v','--version',action='version', version='1.0')

    opt = parser._action_groups.pop()
    req = parser.add_argument_group('required arguments')


    req.add_argument('-d','--cdhit', help='The cd-hit cluster file (.clstr)', required=True, type=pathlib.Path)
    req.add_argument('-g','--genomes', nargs='+', help='Prokka output of genomes (.ffn)', required=True, type=pathlib.Path)
    #parser.add_argument('--out', help='COG network output file', required=True, type=pathlib.Path)
    opt.add_argument('-j', help='Jaccard inclusion threshold for COG co-occurrence weights, default=0.9', required=False, type=float, default=0.9)
    opt.add_argument('-p','--cpus', help='Number of parallel processes, default=1', required=False, type=int, default=1)
    opt.add_argument('-e','--exclude', help='Genes to exclude (chrom/lonely)', required=False, type=pathlib.Path)
    opt.add_argument('-a','--annot', help='Database of gene annotations and categories', required=False, type=pathlib.Path)

    # Output
    req.add_argument('-c','--cliquesDB', help='Output database of cliques discoveredi (.cliques)', required=True, type=pathlib.Path)
    req.add_argument('-n','--nucDB', help='Output fasta file of the sequences present in the final cliques (.fasta)', required=True, type=pathlib.Path)
    
    parser._action_groups.append(opt)
    args = parser.parse_args()
    return args


def main():
    global pa
    # Get command line arguments
    args = get_arguments()

    # Check input files are valid and out dest.

    # Create presence/absence matrix
    ## Import exclude genes and only add to pa/bg if not in excl
    cpa_out = create_pa(args.cdhit,args.genomes)
    pa = cpa_out[0]
    baseGenes = cpa_out[1]
    seqs = cpa_out[2]

    # (optional) Exclude genes
    if args.exclude:
        excl_out = exclude_genes(pa,baseGenes,args.exclude)
        pa = excl_out[0]
        baseGenes = excl_out[1]


    if args.annot:
        annot_out = getAnnots(args.annot)
        annotCat = annot_out[0]
        topAnnot = annot_out[1]

    # Calculate J network
    net = create_net(baseGenes,args.cpus,args.j)


    
    # Optionally annotate genes and then cliques
    # Annotate by fasta description using regex
    # Annotate by specialist lists (gene name, description)
    # Annotate by hmms



    # Discover cliques
    cl = create_cliques(net, baseGenes, annotCat, topAnnot)

    # Ouptut cliqueDB and seqDB
    output_cliques(seqs,cl,baseGenes,args.cliquesDB,args.nucDB)

    return 0


def getAnnots(annotF):
    annotCat = {}
    topAnnot = 0
    with open (annotF,'r') as annot:
        for line in annot:
            tmp = line.split()
            annotCat[tmp[0]] = [int(tmp[1]),tmp[2]]
            if int(tmp[1])>topAnnot:
                topAnnot = int(tmp[1])
    return (annotCat,topAnnot)






def output_cliques(seqs,cl,baseGenes,cliquesOut,seqOut):

    print("Writing cliques and their sequences...")
    ti = time.time()

    bgs = set(baseGenes)

    # Output sequences of the genes
    with open(seqOut,'w') as seqF:
        for s in seqs:
            if s[0] in bgs:
                seqF.write(">" + s[0] + "\n")
                seqF.write(str(s[1]) + "\n")
    
    # Output the database of cliques
    with open(cliquesOut,'w') as cliqF:
        for c in cl:
            cliqF.write(str(c[0]) + "\t" + "\t".join(c[1]) + "\n")

    print("... completed in (", round(time.time()-ti,2), ") seconds",end="\n\n")

    return 0


def create_cliques(net, baseGenes, annotCat, topAnnot):

    print("Discovering cliques...")
    ti = time.time()

    graph = igraph.Graph()
    graph.add_vertices(sorted(baseGenes))
    graph.add_edges(net)

    cliques = graph.maximal_cliques()
    cl = []

    
    for i,c in enumerate(cliques):

        vert = [graph.vs()[v]['name'] for v in c]
        #lv = len(vert)
        maxAnnot = [0]*(topAnnot+1)

        # Add ID and Desc to ve
        for v in range(len(vert)):
            ve = vert[v]
            if ve in annotCat:
                cat = annotCat[ve][0]
                desc = annotCat[ve][1]
                vert[v] = ve + '(' + str(cat) + ',' + desc + ')'
                maxAnnot[cat] = maxAnnot[cat]+1

        mostAnnot = maxAnnot.index(max(maxAnnot))
        if mostAnnot==0: mostAnnot=1


        cl.append([mostAnnot,vert])#,mostAnnot]

    print("... completed in (", round(time.time()-ti,2), ") seconds",end="\n\n")

    return cl



def exclude_genes(pa,baseGenes,exclF):
    with open(exclF,'r') as excl_raw:
        for l in excl_raw:
            excl = l.rstrip()
            if excl in baseGenes:
                exclI = baseGenes.index(excl)
                baseGenes.pop(exclI)
                pa.pop(exclI)
    return (pa,baseGenes)


def create_pa(cdhit,genomes):

    print("Creating presence/absence matrix...")
    ti = time.time()

    baseGenes = []
    bGid = {}
    
    repGene = {}

    # Get baseGenes from cdhit and mapping
    seqs = []

    tmp = []
    with open(cdhit, 'r') as cdhit_raw:
        for line in cdhit_raw:
            if line[0] == '>':
                tmp.append(["",[]])
            else:
                gene = line.split(' ')[1][1:15]
                if line.find('*')>-1:
                    tmp[-1][0] = gene
                tmp[-1][1].append(gene)

    for i,t in enumerate(tmp):
        rep = t[0]
        other = t[1]
        baseGenes.append(rep)
        bGid[rep] = i
        for o in other:
            if o not in repGene:
                repGene[o] = rep
            else:
                #print("Gene already present in mapping: ",gene)
                # Likely as prokka gave same ID as was same gene
                pass

    # Store presence for each baseGenes, which order fastest? Only one read over genomes

    pa = [set() for _ in baseGenes]

    for i,g in enumerate(tqdm(genomes,ncols=70)):
        with open(g, 'r') as g_raw:
            for record in SeqIO.parse(g_raw, "fasta"):
                if record.id in repGene:
                    pa[bGid[repGene[record.id]]].add(i)
                    seqs.append([record.id,record.seq])
                else:
                    print("Error | gene not present in cd-hit clstr file: ",gene)
                    quit()

    print("... completed in (", round(time.time()-ti,2), ") seconds",end="\n\n")

    return (pa,baseGenes,seqs)


def calculate_j(i,J):
    global pa
    lbj = len(pa)
    psi = pa[i]
    subnet = []
    for j in range(i+1,lbj):
        psj = pa[j]
        iNj = len(psi.intersection(psj))
        iUj = len(psi.union(psj))
        jac = iNj/iUj
        if jac>=float(J): # 2 mutual co-occurrences implied by j>=0.09
            subnet.append([i,j])#,"{0:.4f}".format(jac),str(iNj)))
    return subnet


def create_net(baseGenes,cpus,J):
    global pa
    
    print("Creating co-occurrence network...")
    ti = time.time()

    lbj = len(baseGenes)

    # Iterate through and remove those with <2 observations from pa and baseGenes
    for i in range(len(pa)-1,-1,-1):
        if len(pa[i])<2:
            pa.pop(i)
            baseGenes.pop(i)

    lbj = len(baseGenes)

    net = Parallel(n_jobs=cpus)(delayed(calculate_j)(i,J) for i in tqdm(range(lbj),ncols=70))

    #print_net(net,baseGenes)
    #quit()


    outNet = []
    for subnet in net:
        for row in subnet:
            outNet.append([baseGenes[row[0]],baseGenes[row[1]]])

    print("... completed in (",round(time.time()-ti,2),") seconds",end="\n\n")

    return outNet

# Not necessary later as this script will then use net to make cliques
def print_net(net,baseGenes):

    fOut = "out.net"
    with open(fOut,'w') as fOut_raw:
        fOut_raw.write("Gene1\tGene2\tJ\tg1Ng2\n")
        for subnet in net:
            for row in subnet:
                tmp = list(row)
                tmp[0] = baseGenes[tmp[0]]
                tmp[1] = baseGenes[tmp[1]]
                fOut_raw.write("\t".join(tmp) + "\n")


if __name__ == '__main__':
    main()


