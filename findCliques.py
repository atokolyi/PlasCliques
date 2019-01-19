#!/usr/bin/env python3

# Dependent imports
from Bio import SeqIO
from tqdm import tqdm
from joblib import Parallel, delayed
import multiprocessing

# Default imports
import argparse
import pathlib
import os
import tempfile
import copy
import time


def get_arguments():


    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', action='version', version='1.0')

    opt = parser._action_groups.pop()
    req = parser.add_argument_group('required arguments')


    req.add_argument('-a','--assemblies', help='Bacterial isolate assembly/ies (.fasta)', nargs='+', required=True, type=pathlib.Path)

    req.add_argument('-c','--cliquesDB', help='Folder with clique database (CLDB)', required=True, type=pathlib.Path)

    opt.add_argument('-s','--minCliqueSize', help='Minimum size (integer) of cliques included, default=2', required=False, type=int, default=2)

    opt.add_argument('-m','--match', help="Match file to match isolate backbones against a DB of large plasmids (.match)", required=False, type=pathlib.Path)

    opt.add_argument('-l','--outLevel', help='Integer level of output verbosity, 0: Cliques (default), 1: Previous + categories, 2: Previous + genes in cliques, 3: Previous + gene descriptions', type=int, required=False, default=0)
    opt.add_argument('-b','--bulkOut', help='For matrix output, specify an output file (.grid)', type=pathlib.Path, required=False)
    opt.add_argument('-p','--cpus', help='Number of parallel processes, default=1', type=int, default=1, required=False)

    # Maincliques argument for displaying cliques above a certain size, have default of 1?
    
    parser._action_groups.append(opt)
    args = parser.parse_args()
    return args


def main():
    global cl_M

    # Load the arguments
    args = get_arguments()

    nucDB = str(args.cliquesDB) + '/' + str(args.cliquesDB)[-1] + '.fa'
    cliquesDB = str(args.cliquesDB) + '/' + str(args.cliquesDB)[-1] + '.cliques'
    
    # Load the DBs of cliques and sequences
    db_out = load_dbs(cliquesDB)
    annotCat = db_out[0]
    cliqAnnotID = db_out[1]
    cl_M = db_out[2]
    
    # If matching, load the DB of matches
    if args.match:
        matchDB = load_match(args.match)
    

    # Scan for cliques in each assembly (parallel)
    print("Discovering present cliques...",file=sys.stderr)
    ti = time.time()
    net = Parallel(n_jobs=args.cpus)(delayed(processInput)(i=p,ass=args.assemblies[p],minCliqueSize=args.minCliqueSize,fileOut=args.bulkOut,nucDB=nucDB,match=args.match,cliqAnnotID=cliqAnnotID) for p in tqdm(range(len(args.assemblies)),ncols=80))
    print("... completed in (", round(time.time()-ti,2), ") seconds",file=sys.stderr)

    #if args.bulkOut:
        # Process input   ,p,args.assembly[i],
    #    net = Parallel(n_jobs=args.cpus)(delayed(processInput)(i=p,ass=args.a[p],minCliqueSize=args.minCliqueSize,fileOut=args.bulkOut) for p in tqdm(range(len(args.a)),ncols=80))
    #else:
    #    net = Parallel(n_jobs=args.cpus)(delayed(processInput)(i=p,ass=args.a[p],minCliqueSize=args.minCliqueSize,fileOut=args.bulkOut) for p in tqdm(range(len(args.a))))
    
    # Flush cache and print output success message and total time.

    # Report clique output for each assembly
    output_results(net,args.outLevel,args.bulkOut,args.match,cliqAnnotID,cl_M,annotCat)




def load_match(match):
    matchDB = []
    with open(match,'r') as matchF:
        for line in matchF:
            tmp = line.split()
            if len(tmp)==2:
                name = tmp[0]
                els = list(map(int,tmp[1].split(',')))
                matchDB.append([name,els])
    return matchDB


# Assign IDs to the singletons (streamlined later depending on what we include)
# Also codes for cliques or singletons C1,2,3 S1,2,3
## As currently including excluded and clique genes
## Just making IDs larger than they need to be

def load_dbs(cliques):

    # For each clique, store clique ID > ID+desc of constituent genes
    # To then print as int in heatmap and in colname for mouse over clique
    # So store individual in annotCat + all in cliqueCat


    # Remove annotation capabilities??
    # Only needed for verbose output, maybe doesn't even help
    annotCat = {}

    #annotCat[name] = [cat,desc]

    # Load clique DBs
    cliqAnnotID = [] # Annotation category of cliques, remove??
    cl_M = [] # List of backbones IDs to genes which are in


    with open(cliques,'r') as cliqueDBf:
        for i,line in enumerate(cliqueDBf):

            # if last char of gene is ')', extract info, save to annotCat, add to cliqueCat
            tmp = line.rstrip().split('\t')

            cat = tmp[0]
            cliqAnnotID.append(cat)
            raw_vert = tmp[1:]
            vert = []
            
            for v in raw_vert:
                
                if v[-1]==')':
                    g = v.split('(')
                    ve = g[0]
                    vt = '('.join(g[1:])
                    g = vt[0:-1].split(',')
                    cat = g[0]
                    desc = ",".join(g[1:])
                    annotCat[ve] = [cat,desc]
                    vert.append(ve)
                    
                else:
                    vert.append(v)

            
            cl_M.append(set(vert))

    return (annotCat, cliqAnnotID, cl_M) 




def processInput(i,ass,minCliqueSize,fileOut,nucDB,match,cliqAnnotID):
    
    global cl_M
    # Prepare variables
    cl = copy.deepcopy(cl_M)
    subnet = []

    # Blastn for clique genes
    with tempfile.TemporaryDirectory() as tmp:
        
        cmd = "makeblastdb -in " + str(ass) + "  -dbtype nucl -out " + str(tmp) + "/tmp_db"
        res = os.popen(cmd).read() # Check for success
        cmd = "blastn -db " + str(tmp) + "/tmp_db -query " + str(nucDB) + " -outfmt '6 qseqid sacc pident qlen length sstart send bitscore' | awk '($5/$4)>=0.9 && $3>=90'"
        res = os.popen(cmd).read().rstrip().split('\n') # Check for success

        resX = []
        
    # Remove overlapping genes with lower bit scores
    ## If matches overlap keep those with highest bit score
    ## As we're only considering >90% coverage, some smaller gene forms may appear
    ## This will make sure we just keep the larger/better one
    for rei in res:
        tmp = rei.split('\t')
        if len(tmp)==8:
            tmp[2] = float(tmp[2])
            tmp[3:7] = map(int, tmp[3:7])
            tmp[7] = float(tmp[7])
            tmp.append(1)
            if tmp[5]>tmp[6]: # Make so start is smaller than end
                tmp[-1] = 0
                tmp5 = tmp[5] # If transcribed different direction this should be okay
                tmp[5] = tmp[6]
                tmp[6] = tmp5
            found = False
            for re_j in reversed(range(len(resX))):
                # If either is fully inside the other, keep the one with highest bit
                # rej inside tmp   or    tmp inside rej
                rej = resX[re_j]
                if tmp[-1]==rej[-1] and tmp[1]==rej[1]: # On same strand and contig
                    rej_in_tmp = (tmp[5]<=rej[5] and rej[6]<=tmp[6])
                    tmp_in_rej = (tmp[5]>=rej[5] and rej[6]>=tmp[6])
                    
                    if rej_in_tmp or tmp_in_rej:
                        found = True
                        if tmp[7]>rej[7]: # Keep highest bit score
                            #if rej[0]=='LBDJCPPM_00041': print('\t',rej,'\n\t',tmp)
                            resX.append(tmp)
                            resX.pop(re_j)
            if not found:
                resX.append(tmp)
    

    # Collect present cliques 
    pres = [y[0] for y in resX]
    presS = set(pres)


    bb = []
    for key,val in enumerate(cl):
        # key is bb_ID, val is set of bb elements

        # Iterate through all cliques and see if their els are subset of pres
        # If so count count the min occurance of all el in bb
        # Only want to allow duplicates when matching though, as not so informative otherwise?
        if len(val)>=minCliqueSize:
            if val.issubset(presS):
                occ = min([pres.count(z) for z in val]) 
                if fileOut:
                    bb.append(cliqAnnotID[key])
                else:
                    bb.extend([key+1]*occ) # +1 to make 1 indexing
            else:
                if fileOut:
                    bb.append(0)
    

    # If match, compare against matchDB, else just appnd
    # Now that we have our bbs, compare against match
    if match:
        results = []
        while len(bb)>0:
            most = [-1]*5
            for a,aP in enumerate(matchDB):
                saP = set(aP[1])
                fAbs = sum(el in bb for el in aP[1])
                ratio = fAbs/len(aP[1])
                fIn = 0
                # Exponentially more valued as ratio approaches 1
                if ratio>0: fIn = fAbs**ratio # Amount in common weighted by ratio
                if fIn>most[1]:
                #if ratio>most[3] or (most[3]==1 and ratio==1 and fAbs>most[2]):
                    most = [aP[0],fIn,fAbs,ratio,a]
            if most[1]==0:
                break
            results.append([round(most[1],2),most[2],round(most[3]*100,2),most[0]])
            for el in matchDB[most[4]][1]:
                if el in bb: bb.remove(el)
        results = sorted(results,key=lambda l:l[2], reverse=True)
        subnet = [str(ass),results]
    

    # Not matching, just enumerating cliques
    # Removing those which overlap
    else:
        subnet = [str(ass).split('/')[-1],bb]

    return subnet


def output_results(net,outLevel,fileOut,match,cliqAnnotID,cl_M,annotCat):

    # Terminal colours
    GREEN  = '\33[92m'
    ORANGE = '\33[93m'
    RED  = '\33[91m'
    BOLD = '\033[1m'
    END = '\033[0m'

    if fileOut:
        with open(fileOut, 'w') as raw:
            for subnet in net: # For each isolate in the input
                raw.write(subnet[0].replace('#','_') + "\t")
                els = list(map(str,subnet[1]))
                #els = "\n".join(els)
                els = "\t".join(els)
                raw.write(els)
                raw.write("\n")

    if match:
        for subnet in net:
            print(BOLD + subnet[0].split('/')[-1] + END)
            print(BOLD + "Score\tB+S\t%Found\tPlasmidID" + END)
            for r in [r for r in subnet[1] if r[2]>30]:
                print(r[0],r[1],sep="\t",end="\t")
                COLOUR = RED
                if r[2]>30: COLOUR = ORANGE
                if r[2]>50: COLOUR = GREEN
                print(COLOUR + str(r[2]) + END,end="\t")
                print(r[3])
            print(BOLD + "Predicted plasmids: " + END,len([r for r in subnet[1] if r[2]>50]),end="\n\n")
    else:

        for subnet in net:

            isolate = subnet[0].replace('#','_')
            bbs = list(map(str,subnet[1]))

            if outLevel==0:
                print(isolate,end="\t")
                print(",".join(bbs))
            if outLevel==1:
                print(isolate,end="\t")
                for i in range(len(bbs)):
                    bbs[i] = bbs[i] + "(" + cliqAnnotID[int(bbs[i])] + ")"
                print(",".join(bbs))
            if outLevel==2:
                print(isolate,end="\t")
                for i in range(len(bbs)):
                    bb = int(bbs[i])
                    cat = cliqAnnotID[bb]
                    print("\t",cat,"\t",bb,end="\t",sep="")
                    for ve in cl_M[bb]:
                        if ve in annotCat:
                            print(annotCat[ve][1],end="; ")
                    print()
            if outLevel==3:
                print(isolate)
                for bb in bbs:
                    print("\t",cliqAnnotID[int(bb)],"\t",bb,sep="",end="\n")
                    for ve in cl_M[int(bb)]:
                        if ve in annotCat:
                            print("\t",annotCat[ve][0],ve,annotCat[ve][1],sep="\t")
                        else:
                            print("\t",1,ve,'-',sep="\t")
                    print()


if __name__ == '__main__':
    main()
