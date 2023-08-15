# Evolutionary history and Functional inference of the Oskar protein 
- _Leo Blondel, Savandara Besse, Emily L Rivard, Guillem Ylla, Cassandra G Extavour_
- Harvard University (2020)

## Now published in Molecular, Biology and Evolution (2021)
> Leo Blondel, Savandara Besse, Emily L Rivard, Guillem Ylla, Cassandra G Extavour, Evolution of a Cytoplasmic Determinant: Evidence for the Biochemical Basis of Functional Evolution of the Novel Germ Line Regulator Oskar, Molecular Biology and Evolution, Volume 38, Issue 12, December 2021, Pages 5491–5513, https://doi.org/10.1093/molbev/msab284

<br>

## 1. Resources 
### Genome and transcriptome sources 
Available in `./Data/01_Oskar_identification/2019`
- genome_insect_database.csv
- transcriptome_insect_database.csv

### LOTUS, OSK and OSKAR HMM models  
Available in `./Data/Oskar_hmm/V4`
- OSKAR_CONSENSUS.hmm
- LOTUS_CONSENSUS.hmm
- OSK_CONSENSUS.hmm

____________

## 2. Listing and summaries of the provided scripts

### 01_Oskar_identification
#### Genomes and transcriptomes preprocessing
- `1.1_insect_database_builder.py`

Generates the insect and crustacean genome and transcriptome tables containing taxonomic information, bioproject and biosample IDs for all available species in NCBI.

- `1.2_data_downloader.py`

Downloads genomes and transcriptomes from NCBI using the accession numbers from the genome and transcriptome tables generated by 1.1_insect_database_builder.py

- `1.3_run_augustus_training.py`

Generates sbatch scripts to run the Augustus training and creates the necessary input folders from genomic fasta and gff files.

#### Creation of protein sequence databases + Identification of oskar orthologs
- `1.4_run_augustus.py`

Generates sbatch scripts that will launch Augustus in annotation mode over all genomes with an automatic model assignment based on the closest available insect order last common ancestor model.

- `1.5_Oskar_tracker.ipynb`

Collects all putative protein sequences from annotated genomes and transcribed transcriptomes
  - For TSA: Uses Transeq to generate translated sequences from transcriptomes
  - For GCF: Collects protein sequences available in annotated genomes
  - For GCA: Collects generated protein sequences from Augustus de novo annotations <br>
 
Tracks Oskar candidates from all gathered protein sequences using `execute_hmmsearch.py`. Identifies Oskar sequences using LOTUS and OSK HMMs and validates them through cross-validation using the OSKAR HMM. Finally, it filters out the duplicated sequences to generate a final filtered dataset as a fasta files.

_Generated outputs:_
- Available in `./Data/01_Oskar_identification/oskar_tracker_results/`
  - Search metadata result table 
    - search_results.csv
  - Raw sequences
    - long_oskar.fasta
    - oskar_filtered.fasta 
    - oskar_all.fasta
  - Alignments
    - long_oskar.aligned.fasta
    - oskar_all.aligned.fasta
    - oskar_filtered.aligned.fasta
    - oskar_filtered.aligned.LOTUS_domain.fasta
    - oskar_filtered.aligned.OSK_domain.fasta

Generates Table S1. 

<br>

### 02_Oskar_analyses
#### Correlative analysis of assembly quality and absence of oskar
- `2.1_Oskar_discovery_quality.ipynb`

It generates the plots shwown on Figure S1 and statistically compares the distributions of each of the 8 available assembly statistics: Contig and Scaffold N50, Contig and Scaffold L50, contig and Scaffold counts and Number of Contig and Scaffold per genome length.

#### TSA metadata parsing and curation
- `2.3_Oskar_tissues_stages.ipynb`

Collects, identifies and parses the tissues and developmental stages where oskar was found in transcriptomes. The initial collection of metadata is extracted from the description of Biosample projects. Consensus keywords were created for redundant tissues and developmental stages. 
The results of this analysis can be visualized in Figure 3 and Figure S5.  

_Generated outputs:_
- Available in `./Data/02_Oskar_analyses/2.3/`
  - 2.3.1.insect_metadata_biosample_database.csv
  - 2.3.2.oskar_all_tissues_stages.csv
  - 2.3.3.oskar_all_tissues.csv
  - 2.3.4.oskar_all_stages.csv

#### Conservation differences between the Holo and Hemimetabolous Oskar sequences
- `2.4_Oskar_pgc_specification.ipynb`

Classifies and sorts Oskar sequences according to a database of germ cells specification strategies across insects extracted from an upublished literature review. 

_Generated outputs_:
- Available in `./Data/02_Oskar_analyses/2.4/FASTA/`
  - OSKAR_holometabola.fasta
  - LOTUS_holometabola.fasta
  - OSK_holometabola.fasta
  - OSKAR_hemimetabola.fasta
  - LOTUS_hemimetabola.fasta
  - OSK_hemimetabola.fasta

- `2.5_Oskar_dimer_monomer_extraction.ipynb`

Sorts Oskar sequences in different groups based on the LOTUS monomeric or dimeric state. 

_Generated outputs_:
- Available in `../Data/02_Oskar_analyses/2.5/FASTA/`
  - OSKAR_Dimeric_alignment.fasta
  - OSKAR_Monomeric_alignment.fasta

#### Computation of secondary structure conservation

- `2.6_Oskar_lotus_osk_structures.ipynb`

Generates the secondary structure prediction for the two Oskar domains: LOTUS and OSK, using the JPred 4 algorithm. Generates the plots shown in Figure S8.

_Generated outputs_:
- Available in `../Data/02_Oskar_analyses/2.6/`
  - Raw structures (in STRUCTURES folder)
    - LOTUS_raw_structures.faa
    - OSK_raw_structures.faa
    - LOTUS_structures_alignment.fasta
  - Alignments (in STRUCTURES folder)
    - OSK_structures_alignment.fasta
    - STRUCTURES/LOTUS_structures_alignment.trimmed_0.3.fasta
    - STRUCTURES/OSK_structures_alignment.trimmed_0.3.fasta
  - JPRED scores
    - LOTUS_jpred_scores.csv
    - OSK_jpred_scores.csv

#### Phylogenetic inference of Oskar sequences in the Hymenopteran
- `2.7_Oskar_duplication.ipynb`

Runs the phylogenetic reconstruction and figure generation for the hymenopteran oskar duplication seen in Figure 4 and Figure S4. 

#### Dimensionality reduction of oskar alignment sequence space
- `2.8_Oskar_MCA_Analysis.ipynb`

Performs the Multiple Correspondence Analysis (MCA) on the full lenght Oskar, LOTUS and OSK alignments and generates the plots shown in Figure S6.

#### Analyses about the evolution of Oskar 
- `2.9_Oskar_Tree_of_Evolution.ipynb`

Generates Figure 3 where Oskar identification statistics across insect orders were plotted with tissue and developmental stage information. The underlying phylogenetic relationships were extracted from _Misof et al. - Phylogenomics resolves the timing and pattern of insect evolution (Science, 2014); <url>https://doi.org/10.1126/science.1257570</url>_

- `2.10_Long_Oskar_Evolution.ipynb`

Computes the overall alignment occupancy in the Dipteran Oskar alignment (trimmed for a minimum of 10% overall occupancy for any position). Groups the results by Dipteran families and generates Figure S7. 
 
<br>

### 03_Oskar_score_generation

#### Computation of the JSD score
- `score_conservation.py`

Computes the Jensen-Shannon Divergence (JSD) score, a measure of how much information any position in the alignment brings to the overall alignment. We used the initial implementation of JSD score proposed in _Capra & Singh - Predicting functionally important residues from sequence conservation (Bioinformatics, 2007); <url>https://doi.org/10.1093/bioinformatics/btm270</url>_

#### Computation of Oskar conservations scores including the electrostatic and hydrophobic conservation score
- `besse_blondel_conservation_scores.py`

Computes three conservation scores for a given protein alignment and saves the scores in a CSV table
- Valdar score calculates an overall amino acid conservation score by taking account the transition probabilities, stereochemical properties and amino acid frequencies gaps for  each position in the alignment. Implementation based on _Valdar - Scoring residue conservation. (Proteins, 2002); <url>https://doi.org/10.1002/prot.10146</url>_
- Electrostatic conservation score calculates the conservation of electrostatic properties for each position in the alignment
- Hydrophobic conservation score calculates the conservation of hydrophobic properties for  each position in the alignment

#### Computation of the RNA binding affinity score
- `RNABindR_run_predictions.py`

Automatized the prediction of RNA-binding sites by generating requests to the RNABindR web-server using a Fasta alignment as an input.

#### Conservation differences between the Holo and Hemimetabolous Oskar sequences  
- `3.1_LogRatio_Bootstrap.ipynb`

Calculates the conservation bias score from Hemimetabolous and Holometabolous sequences and generates the plots shown in Figure 5.

#### Generation of scores and mapping tables for Oskar visualization  
- `3.2_generate_scores.ipynb`

Computes (if missing) and concatenates all generated conservation scores resulting from `score_conservation.py` and `besse_blondel_conservation_scores.py`. Generates the required mapping table to vizualize the scores on the D. melanogaster LOTUS and OSK structures.

_Generated outputs:_
- Available in `./Data/03_Oskar_scores_generation/CSV/`
  - scores.csv
  - mapping.csv

#### Scripts used for figures 
- `2.2_Oskar_insect_repartition.ipynb`

Generates the plots shown in Figure S1.

- `3.3_Make_Sequence_Logos.ipynb`

Generates all weblogos as shown in Figure 6 and Figure 7 and Figure S8. 

<br>

### 04_Oskar_visualization
#### Visualization of conservation scores
- `Oskar_pymol_visualization.py`

_Required inputs:_
- PDB structures
    - ./04_Oskar_visualization/5A4A.pdb (OSK)
    - ./04_Oskar_visualization/5NT7.pdb (LOTUS-VASA)
    - ./04_Oskar_visualization/2DB3.pdb (VASAR-NA)
- CSV tables
    - ./Data/03_Oskar_scores_generation/CSV/scores.csv
        - all conservation scores saved in a CSV file generated with scripts in `03_Oskar_scores_generation`
    - ./Data/03_Oskar_scores_generation/CSV/mapping.csv
        - mapping table between oskar alignment positions and oskar structure positions

_Guidelines:_
- Run this script in Pymol
- Command line template: SHOW, STRUCTURENAME, SCORE, AREA (option)
    - STRUCTURENAME types:Hydro
        - LOTUSVASA
        - OSK
    - SCORE types : Raw score - LOTUS monomers - LOTUS dimers - Holometabola - Hemimetabola
        - (Conservation) Valdar, Valdarmon, Valdardim, Valdarholo, Valdarhemi
        - (Conservation) JSD, JSDmon, JSDdim, JSDholo, JSDhemi
        - (Electrostaticity) Elec, Elecmon, Elecdim, Elecholo, Elechemi
        - (Hydrophobicity) Hydro, Hydromon, Hydrodim, Hydroholo, Hydrohemi
        - (Valdar Holo / Hemi) Logratio
        - (RNA-binding) RNABindR
    - AREA types :
        - For OSK --> OSK, core, rna-binding
        - For LOTUSVASA --> LOTUS, dimerization-interface, vasa-interface
        - LOTUSRNA function to generate RNA / LOTUSVASA alignment
          
____________

## 3. Software and libraries

| Type     | Name           | Version  | Source                                                          |
|----------|----------------|----------|-----------------------------------------------------------------|
| Software | HMMER          | 3.1.b2   | http://hmmer.org/                                               |
| Software | Pymol          | 1.8.x    | https://pymol.org/                                              |
| Software | rsync          | 3.1.2    | http://rsync.samba.org/                                         |
| Software | Python3        | 3.7      | https://www.python.org/                                         |
| Software | Mrbayes        | 3.2.6    | http://nbisweden.github.io/MrBayes/                             |
| Software | trimal         | 1.2rev59 | http://trimal.cgenomics.org/                                    |
| Software | transeq        | 6.6.0.0  | http://emboss.sourceforge.net/apps/cvs/emboss/apps/transeq.html |
| Software | augustus       | 2.5.5    | http://augustus.gobics.de/                                      |
| Software | JPred4         | 4        | http://www.compbio.dundee.ac.uk/jpred/                          |
| Software | RNABindR       | 2        | ailab1.ist.psu.edu/RNABindR/                                    |
| Software | Inkscape       | 0.92.3   | https://inkscape.org/                                           |
| Library  | jupyter        | 4.4.0    | https://jupyter.org/                                            |
| Library  | ete3           | 3.3.1    | http://etetoolkit.org                                           |
| Library  | pandas         | 0.25.1   | https://pandas.pydata.org/                                      |
| Library  | mca            | 1.0.3    | https://pypi.org/project/mca/                                   |
| Library  | fuzzywuzzy     | 0.17.0   | https://github.com/seatgeek/fuzzywuzzy                          |
| Library  | BeautifulSoup4 | 4.6.3    | https://pypi.org/project/beautifulsoup4/                        |
| Library  | biopython      | 1.74     | https://pypi.org/project/biopython/                             |
| Library  | numpy          | 1.16.2   | https://www.numpy.org/                                          |
| Library  | seaborn        | 0.9.0    | https://seaborn.pydata.org/                                     |
| Library  | matplotlib     | 3.0.0    | https://matplotlib.org/                                         |
| Library  | scipy          | 1.1.0    | https://www.scipy.org/                                          |
| Library  | progressbar    | 3.38.0   | https://github.com/niltonvolpato/python-progressbar/            |

<br>

## 4. Code Maintenance
- Leo Blondel (@Xqua)
- Savandara Besse (@ladyson1806)
