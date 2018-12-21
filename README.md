# Creating a database of cliques (category_annot.py)
### Summary
Starting with a collection of genome assemblies, the output of this pipeline is a collection of genes, with groups which are highly co-occurring being clustered into functional groups (a clique). This can then be used as input to `backboneCat.py` to discover their presence in novel isolates. 

An overview of this process is that we extract genes from the assemblies, cluster them by similarity, optionally annotate the clusters, detect the presence of these clusters in the original isolates to create a presence/absence matrix, use this to calculate co-occurrence similarities between genes

### Prerequisites
- python3 (3.6.6)
- Biopython (1.72)
`pip3 install biopython`
- tqdm (4.26.0)
`pip3 install tqdm`
- joblib (0.12.5)
`pip3 install joblib`
- python-igraph (0.7.1.post6)
`pip3 install python-igraph`

#
#
#### 1) Extracting genes with prokka
~~
python runProkka_SLURM.py --contigs completePlasmids/*.fasta --outdir prokkaOut --proteins /vlsci/SG0006/shared/data/klebsiella/references/NTUH-K2044.aa --genus Klebsiella --species pneumoniae --extension .fasta
cat prokkaOut/*/*.ffn > allGenes.ffn
~~

#### 2) Clustering genes with CD-HIT-EST
~~
cd-hit-est -M 16000 -T 8 -G 1 -g 1 -i allGenes.ffn -c 0.9 -s 0.9 -n 8 -o clusteredGenes.ffn
~~

#
#### 3) Create a network of highly co-occurring genes
Make `scripts/clusterAnalyse_set_mp.py` very barebones so easily outputs netj9.tsv

#
#### 4) Annotate the clustered genes (optional)
Arrange your annotations in the format below, with a seperate file for each annotation category as input to the clique finding script.
~~
UniqueGeneID    AnnotationDescription
~~
A list of UniqueGeneIDs which are to be excluded from clique creation (i.e because they are only observed once [plus how to calculate] or map undesireably [chromosomal])

#
#### 5) Create the final database of cliques
~~
~~

# 

#
#
#
# Discovering cliques in isolates (backboneCat.py)
