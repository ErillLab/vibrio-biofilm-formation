"""Promoter analysis for biofilm formation regulon in Vibrio cholareae"""
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import motifs
from Bio.Alphabet import IUPAC
import os
import time
import pickle
from Bio.Align.Applications import ClustalwCommandline
import glob
from math import log
from matplotlib import pyplot as plt
import numpy as np
from collections import defaultdict
import collections
import pylab
import itertools
import pandas as pd
import matplotlib

Entrez.email = 'sefa1@umbc.edu'
EVAL_THRESHOLD = 10**-50
DATA_DIR = './data'

def safe_log2(x):
    """Implements log2, but defines log2(0) = 0"""
    return log(x,2) if x > 0 else 0

def reverse_complement(seq):
    """Reverse complement of a sequence"""
    return Seq(seq).reverse_complement().tostring()

def base_freq(seq):
    """Base frequencies of a sequence"""
    epsilon = 1e-10
    bg_dict = dict()
    for b in IUPAC.ambiguous_dna.letters:
        bg_dict[b] = (epsilon + seq.count(b)) / len(seq)
    return bg_dict

def read_locus_tags():
    """Read the list of locus tags."""
    with open('locus_tags.txt') as f:
        locus_tags = f.read().split()
    return locus_tags

def read_promoters():
    """Read list of promoter sequences"""
    with open("vpsRT_promoters.fna") as f:
        promoters = list(SeqIO.parse(f, "fasta"))
    return promoters

def gene_positions(locus_tag):
    """Given a locus tag return a tuple: end position of upstream gene and the
    start position of the gene of interest.
    Start and end positions for genes of interest are hard-coded here.
    """
    d = {'VC0665': (27, 300),
         'VC0916': (39, 300),
         'VC0917': (22, 334),
         'VC0928': (105, 300),
         'VC0930': (37, 300),
         'VC0934': (3, 492),
         'VCA0952': (3, 580),
         'VC1888': (3, 337),
         'VC2647': (3, 361),
         'VC0583': (5, 300),
         'VCA0075': (10, 300),
         'VCA0785': (81, 249)}
    return d[locus_tag]

def read_blast():
    """Parse blast results and return list of blast records."""
    blast_results_file = "blast_results.xml"
    with open(blast_results_file) as results_handle:
        blast_records = list(NCBIXML.parse(results_handle))
    return blast_records

def load_bioproject_pickle():
    cache_file = "bioproject_cache.pkl"
    with open(cache_file) as f:
        cache = pickle.load(f)
    #print "cache pickle file loaded [%d]" % len(cache)
    return cache

def get_bioprojects(accs, cache=None):
    """Return the dictionary of bioproject numbers for all accession numbers"""
    cache_file = "bioproject_cache.pkl"
    if not cache:
        try:
            cache = load_bioproject_pickle()
        except IOError:
            print "can't open the pickle file, creating one."
            cache = {}
    for acc in accs:
        if acc in cache:
            continue            # already exists, skip
        print "bioproject for %s not in cache file, fetching from NCBI." % acc
        time.sleep(0.5)
        handle = Entrez.efetch(db='nuccore', id=acc, retmode='gbwithparts',
                               rettype='text')
        rec = SeqIO.read(handle, 'gb')
        # dbxrefs feature has the following format
        # dbxrefs=['BioProject:PRJNA224116 BioSample:SAMN02469781'])
        cache[acc] = [tok.split(':')[1]
                      for tok in ' '.join(rec.dbxrefs).split(' ')
                      if (tok.startswith('BioProject:') or
                          tok.startswith('Project:'))][0]
    with open(cache_file, 'w') as f:
        pickle.dump(cache, f)
    return dict((acc, cache[acc]) for acc in accs)

def get_bioproject(acc, *args, **kwargs):
    """return the bioproject number for the given accession number"""
    return get_bioprojects([acc], *args, **kwargs)[acc]

def get_species(blast_records):
    """Get the list of species that have blast hits for all queries"""
    unique = lambda xs: list(set(xs))
    # get all genome accessions and get their bioproject numbers
    genome_accessions = []
    for blast_record in blast_records: # go over blast results of each query
        blast_record.alignments = [alignment
                                   for alignment in blast_record.alignments
                                   for hsp in alignment.hsps
                                   if hsp.expect < EVAL_THRESHOLD]
        genome_accessions.extend((aln.accession
                                  for aln in blast_record.alignments))
    bp_dict = get_bioprojects(genome_accessions)
    # find bioproject numbers (equivalently genomes) that have hits for all
    # query sequences
    universal_bps = set(bp_dict.values())
    # filter these set as you iterate over the blast results
    for blast_record in blast_records:
        accessions = set(bp_dict[aln.accession]
                         for aln in blast_record.alignments)
        universal_bps &= accessions
        #print len(universal_bps)
    return list(universal_bps)

def get_org_name(bp):
    """Return the organism name, given the bio-project identifier"""
    try:
        bp = bp.replace('PRJNA', '')
        handle = Entrez.esummary(db='bioproject', id=bp)
        rec = Entrez.read(handle, validate=False)
        return rec['DocumentSummarySet']['DocumentSummary'][0]['Organism_Name']
    except RuntimeError:
        return "Cannot get document summary"

def get_org_names(bioprojects):
    """Given a list of bioproject ids
    return the organism name for all bioprojects."""
    bp_dict = load_bioproject_pickle()
    for bp in bioprojects:
        try:
            print bp, get_org_name(bp) #, ', '.join([acc for acc in bp_dict.keys()
                                       # if bp_dict[acc] == bp])
        except:
            pass
def download_genome(acc):
    """download the genome of a given accession number"""
    print "Downloading genome %s" % acc
    filename = os.path.join(DATA_DIR, 'genomes', acc+'.gbk')
    net_handle = Entrez.efetch(db='nuccore', id=acc,
                               rettype='gbwithparts', retmode='text')
    out_handle = open(filename, 'w')
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()

def get_seq(acc, start, end):
    """Get [start,end] of the genome sequence"""
    gbk_file = os.path.join(DATA_DIR, 'genomes', acc+'.gbk')
    if not os.path.isfile(gbk_file):
        download_genome(acc)
    rec = SeqIO.read(gbk_file, 'genbank')
    return rec.seq[start:end]

def parse_homolog_seqs(blast_record):
    """Given a BLAST record, extract all hit sequences"""
    bp_cache = load_bioproject_pickle() # load genome acc - bioproject mappings
    subject_seqs = {}
    for aln in blast_record.alignments:
        hsps = sorted(aln.hsps, key=lambda x: x.expect) # sort hits by e-value
        hsp = hsps[0]
        if hsp.expect > EVAL_THRESHOLD:
            continue
        bp = get_bioproject(aln.accession, bp_cache)
        left_ext = hsp.query_start - 1
        right_ext = blast_record.query_length - hsp.query_end
        if hsp.sbjct_start < hsp.sbjct_end: # forward strand
            sbjct_start = max(1, hsp.sbjct_start - left_ext)
            sbjct_end = hsp.sbjct_end + right_ext
            seq = get_seq(aln.accession, sbjct_start-1, sbjct_end)
        else:
            sbjct_start = hsp.sbjct_start + left_ext
            sbjct_end = max(1, hsp.sbjct_end - right_ext)
            seq = get_seq(aln.accession, sbjct_end-1, sbjct_start)\
                .reverse_complement()
        if any((c for c in seq.tostring() if c not in "ACGTacgt")):
            print "ambiguous sequence, skipping %s" % bp
            pass
        elif len(seq) < 0.95 * blast_record.query_length:
            # if the whole sub-sequence cannot be retrieved (e.g. the
            # sub-sequence is at the end of the contig)
            print "not enough coverage, skipping %s" % bp
            pass
        else:
            subject_seqs[bp] = seq
    return subject_seqs

def prep_fasta_for_clustal(query, blast_record, bioprojects):
    """Given a dictionary of <bioproject, seq> and a list of bioproject
    identifiers, return the multiple sequence alignment of sequences that belong
    to one of the organisms in the bioproject list."""
    # generate fasta file for clustalw input
    subject_seqs = parse_homolog_seqs(blast_record)
    infile_name = blast_record.query.split()[0] + '.fasta'
    with open(os.path.join(DATA_DIR, 'alignments', infile_name), 'w') as f:
        f.write('>%s\n' % query.description.replace(' ', '_'))
        f.write('%s\n' % query.seq.tostring())
        for bp in bioprojects:
            if bp not in subject_seqs:
                continue
            f.write('>%s_%s\n' % (bp, get_org_name(bp).replace(' ', '_')))
            f.write('%s\n' % subject_seqs[bp])

def prep_all_fasta():
    """Make all fasta files for clustalw alignment"""
    blast_records = read_blast()
    universal_bps = get_species(blast_records)
    promoters = read_promoters()
    for promoter, blast_rec in zip(promoters, blast_records):
        print promoter.description
        prep_fasta_for_clustal(promoter, blast_rec, universal_bps)

def run_clustal(infile):
    cline = ClustalwCommandline("clustalw", infile=infile)
    print cline
    cline()

def run_clustal_all():
    cur_dir = os.getcwd()
    os.chdir(os.path.join(DATA_DIR, 'alignments'))
    for infile in glob.glob('*.fasta'):
        run_clustal(infile)
    os.chdir(cur_dir)

def read_alignments():
    alns = [AlignIO.read(aln_file, 'clustal')
            for aln_file in glob.glob(os.path.join(DATA_DIR, 'alignments', '*.aln'))]
    return alns

def find_eltor_seq(aln):
    # Find the V. cholareae el tor promoter in the alignment
    locus_tags = read_locus_tags()
    seqs = [seq for seq in aln.get_all_seqs()
            if any(lt in seq.description for lt in locus_tags)]
    assert len(seqs) == 1
    return seqs[0]

def frequencies(xs):
    length = 0
    counts = {}
    for x in xs:
        if not x in counts:
            counts[x] = 1
        else:
            counts[x] += 1
        length += 1
    length = float(length)
    return [count/length for count in counts.values()]

def h(ps):
    """compute entropy (in bits) of a probability distribution ps"""
    return -sum([p * safe_log2(p) for p in ps])

def dna_entropy(xs,correct=True):
    """compute entropy (in bits) of a DNA sanmple"""
    ps = frequencies(xs)
    correction = (3)/(2*log(2)*len(xs)) #Basharin 1959
    #print "correction:",correction
    return h(ps) + (correction if correct else 0)

def alignment_entropy(aln, correct=False):
    """Return the list of column-wise entropies for the alignment"""
    entropies = np.array([dna_entropy(aln.get_column(i), correct)
                          for i in xrange(aln.get_alignment_length())])
    return entropies

def load_sites(infile):
    with open(infile) as f:
        sites = f.read().strip().split()
    return sites

def build_motif(sites):
    instances = [Seq(site, IUPAC.unambiguous_dna) for site in sites]
    m = motifs.create(instances)
    m.pseudocounts = 0.8/4      # Nishida et al.
    # http://www.ncbi.nlm.nih.gov/pubmed/19106141
    return m

def score_seq(motif, seq):
    cleanseq = Seq(''.join([b for b in seq if b != '-']), alphabet=motif.alphabet)
    assert len(cleanseq) == len(motif), '%d %d' % (len(cleanseq), len(motif))
    return max(motif.pssm.calculate(cleanseq),
               motif.pssm.calculate(cleanseq.reverse_complement()))

def num_bases(seq):
    """Given a sequence with potentially gaps, return the number of ACTG
    bases"""
    return len([b for b in seq if b != '-'])

def score_promoter(motif, promoter):
    """PSSM search of the promoter region"""
    scores = {}
    for start_pos in xrange(len(promoter)):
        # If the site starts with a gap, it means that the sequence is scored in
        # one of the previous iterations, skip it. This doesn't apply to the
        # first iteration if the aligned seq starts with a gap.
        if scores and promoter[start_pos] == '-':
            continue
        # get len(motif) bps for pssm search (skip gaps)
        end_pos = start_pos
        while (end_pos < len(promoter)-1 and
               num_bases(promoter[start_pos:end_pos]) < len(motif)):
            end_pos += 1
        if num_bases(promoter[start_pos:end_pos]) < len(motif):
            break
        scores[start_pos] = score_seq(motif, promoter[start_pos:end_pos])
    return scores

def get_threshold(motif):
    """For a given motif return a threshold to be used for PSSM search
    in Vibrio cholera genome."""
    print 'Computing threshold for the motif.'
    dist = motif.pssm.distribution(precision=10**4,
                background={'A': 0.2625622082660725,
                            'C': 0.23743779173392748,
                            'G': 0.2625622082660725,
                            'T': 0.23743779173392748})
    #background={'A': 0.32, 'C': 0.18, 'G': 0.18, 'T': 0.32})
    alpha = -log(1-0.01)/350
    th = dist.threshold_fpr(alpha)
    print 'th', th
    return th

def promoter_pssm_search_test(motif, aln):
    """Get an alignment, find the eltor promoter and perform pssm search on
    it."""
    # compute the score distribution of the motif
    th = get_threshold(motif)
    print "Threshold:", th
    prom = find_eltor_seq(aln)
    scores = defaultdict(list)
    for seq in aln:
        tmp_scores = score_promoter(motif, seq)
        for pos, score in tmp_scores.items():
            if score < th:
                continue
            scores[pos].append(seq)
    return scores

def get_tf_file(tf):
    """Given the TF name, return the file name containing sites"""
    if tf == 'VpsT':
        return './matrices/VpsT_MEME_symm.txt'
    if tf == 'VpsR':
        return './matrices/VpsR_MEME_sym.txt'
    return "./matrices/%s.txt" % tf

TFS = ["AphA", "CRP", "HapR", "H-NS", "VpsR", "VpsT"]

def pssm_search_all(aln):
    """Given an alignment and a list of TFs of interest, perform PSSM search for
    all TFs."""
    tfs = TFS
    scores = {}
    for tf in tfs:
        print "PSSM search for %s" % tf
        sites = load_sites(get_tf_file(tf))
        motif = build_motif(sites)
        scores[tf] = promoter_pssm_search_test(motif, aln)
    return scores

def get_color_patterns():
    colors = "rgbcmyk"
    patterns = ('/', '\\', '.', 'o', 'x', '-')
    col_pats = zip(colors, patterns)
    color_pattern_dict = {}
    for tf, (col, pattern) in zip(TFS, col_pats):
        color_pattern_dict[tf] = (col, pattern)
    return color_pattern_dict

font = {'size'   : 18}
matplotlib.rc('font', **font)

def plot_legend_only(horizontal=False):
    """Only plot the legend"""
    color_patterns = get_color_patterns()
    fig = pylab.figure()
    if horizontal:
        figlegend = pylab.figure(figsize=(10, 1))
    else:
        figlegend = pylab.figure(figsize=(2, 5))
    ax = fig.add_subplot(111)
    for tf, (acolor, apattern) in color_patterns.items():
        rect = pylab.Rectangle((0, 0), 10, 10, alpha=0.6,
                               edgecolor=acolor, hatch=apattern,
                               facecolor='w',label=tf)
        ax.add_patch(rect)
    handles, labels = ax.get_legend_handles_labels()
    legend_dict = dict((l, h) for h, l in zip(handles, labels))
    sorted_keys = sorted(legend_dict.keys())
    if horizontal:
        figlegend.legend([legend_dict[k] for k in sorted_keys], sorted_keys,
                         loc='center', ncol=len(TFS))
    else:
        figlegend.legend([legend_dict[k] for k in sorted_keys], sorted_keys, loc='center')
    figlegend.savefig('./plots/legend.png', dpi=200)

def get_gene_positions(aln):
    locus_tag = find_eltor_seq(aln).description.split('_')[0]
    start, end = gene_positions(locus_tag)
    # where do these positions map in the alignment (with gaps)
    num_bases = 0               # count bases
    for (i, b) in enumerate(find_eltor_seq(aln).seq.tostring()):
        if b in 'ACTGactg':
            num_bases += 1
        if num_bases == start:
            correct_start = i
        if num_bases == end:
            correct_end = i
    return (correct_start, correct_end)

def plot(aln, scores, show=False):
    """Plot the alignment"""
    pylab.rcParams['figure.figsize'] = (15, 3)
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.subplots_adjust(hspace=0)

    # tf binding site layer
    ax1.set_xlim(0, aln.get_alignment_length())
    ax1.set_ylim(0, 105)
    ax2.yaxis.get_major_ticks()[0].label1.set_visible(False) # remove zero tick
    ax1.set_ylabel("% proms.")
    # promoter search
    #scores = pssm_search_all(aln)
    for tf in TFS:
        mlen = build_motif(load_sites(get_tf_file(tf))).length
        acolor, apattern = get_color_patterns()[tf]
        for start_pos in scores[tf]:
            conservation = len(scores[tf][start_pos]) / float(len(aln)) * 100
            rect = pylab.Rectangle((start_pos, 0), mlen, conservation, alpha=0.6,
                                   edgecolor=acolor, hatch=apattern, label=tf,
                                   facecolor='w',
                                   linewidth=2)
            ax1.add_patch(rect)

            # remove redundant legends
    handles, labels = ax1.get_legend_handles_labels()
    legend_dict = dict((l, h) for h, l in zip(handles, labels))
    sorted_keys = sorted(legend_dict.keys())
    #fig.legend([legend_dict[k] for k in sorted_keys], sorted_keys, loc='best')

    # plot start position of the gene of interest
    # and the end position of the upstream gene
    x, y = get_gene_positions(aln)
    ax1.axvline(x=x, color='r')
    ax2.axvline(x=x, color='r')
    ax1.axvline(x=y, color='g')
    ax2.axvline(x=y, color='g')

    # entropy layer
    ent = alignment_entropy(aln)
    ax2.plot(ent)
    ax2.set_ylim(0, 1.5)
    ax2.invert_yaxis()
    ax2.set_ylabel("H(x)")

    # ticks
    ax1.yaxis.set_ticks([0, 50, 100])
    ax2.yaxis.set_ticks([0, 0.4, 0.8, 1.2])

    desc = find_eltor_seq(aln).description
    locus_tag = desc.split('_')[0]
    gene = desc.split('_')[1]
    print locus_tag, gene
    plt.title("%s (%s)" % (locus_tag, gene), y=2)
    if show:
        plt.show()
    else:
        # save the figure
        plt.savefig('./plots/' + desc+'.png', frameon=False, dpi=200)
    return plt

def plot_all(scores, show=False):
    for aln,score in scores:
        plot(aln, score, show)

def list_of_isolates():
    alignments = read_alignments()
    seq_descriptions = [seq.description for seq in alignments[2]]
    bioprojects = [desc.split('_')[0] for desc in seq_descriptions]
    get_org_names(bioprojects)
    
def get_genome(acc):
    """Get genome record given the accession number."""
    gbk_file = os.path.join(DATA_DIR, 'genomes', acc+'.gbk')
    if not os.path.isfile(gbk_file):
        download_genome(acc)
    rec = SeqIO.read(gbk_file, 'genbank')
    return rec

def distance(loca, locb):
    """Distance between two locations."""
    assert loca[0] < loca[1] and locb[0] < locb[1]
    return max(loca[0], locb[0]) - min(loca[1], locb[1])

def get_promoter(name, chr_a, chr_b, genes_a, genes_b):
    """Given gene name (locus tag if gene name is not available), return the
    upstream region.

    genes_a and genes_b are list of (name, feature) tuples.
    """
    gene_names_a = [x for (x, _) in genes_a]
    gene_names_b = [x for (x, _) in genes_b]
    if name not in gene_names_a + gene_names_b:
        print 'WARNING |', 'gene', name, 'not found, skipping'
        return (-1, Seq(''))
    if name in gene_names_a:
        genes_to_look = genes_a
        chr_to_look = chr_a
    else:
        genes_to_look = genes_b
        chr_to_look = chr_b

    index = 0
    while genes_to_look[index][0] != name:
        index += 1
    goi = genes_to_look[index][1]
    if ((goi.strand == 1 and index == 0) or (goi.strand == -1 and index ==
                                             len(genes_to_look)-1)):
        print 'WARNING | first/last gene on the chromosome, skipping.', name
        return (-1, Seq(''))
    up = (genes_to_look[index-1][1] if goi.strand == 1 else
          genes_to_look[index+1][1])
    goi_loc = (goi.location.start.position, goi.location.end.position)
    up_loc = (up.location.start.position, up.location.end.position)
    dist = distance(goi_loc, up_loc)
    #print goi
    #print up
    if dist < 300:
        if goi.strand == 1:
            promoter_loc = (goi_loc[0]-300, goi_loc[0]+50)
        else:
            promoter_loc = (goi_loc[1]-50, goi_loc[1]+300)
    else:
        if goi.strand == 1:
            promoter_loc = (up_loc[1], goi_loc[0]+50)
        else:
            promoter_loc = (goi_loc[1]-50, up_loc[0])
    return (promoter_loc,
            chr_to_look.seq[promoter_loc[0]:promoter_loc[1]+1].tostring())


def get_genes(rec):
    """Given a genome record, return locus tags of its genes."""
    genes = [(f.qualifiers.get('gene', f.qualifiers['locus_tag'])[0], f)
             for f in [feature for feature in rec.features
                       if feature.type == 'gene']]
    genes.sort(key=lambda (name, f): f.location.start.position)
    return genes

def get_chra():
    return get_genome('NC_002505.1')

def get_chrb():
    return get_genome('NC_002506.1')

def chra_genes():
    return get_genes(get_chra())

def chrb_genes():
    return get_genes(get_chrb())

def read_rnaseq_table(tf):
    """Read the file containing RNA-Seq data and get the list of genes that have
    significantly regulated by VpsT or VpsR.
    """
    if tf == 'VpsT':
        motif = build_motif(load_sites(get_tf_file('VpsT')))
        df = pd.read_excel('Wt vs. VpsT-Rg1 removed subexperiment.xlsx', 1)
    else:
        motif = build_motif(load_sites(get_tf_file('VpsR')))
        df = pd.read_excel('Wt vs. VpsR-Rg1 removed subexperiment.1.xlsx', 0)

    threshold = get_threshold(motif)
    chra = get_chra()
    chrb = get_chrb()
    chra_genes = get_genes(chra)
    chrb_genes = get_genes(chrb)
    genes = df['Feature ID']
    promoters = {}
    for g in genes:
        loc, seq = get_promoter(g, chra, chrb, chra_genes, chrb_genes)
        pssm_search = score_promoter(motif, seq)
        for pos, score in pssm_search.items():
            if score > threshold:
                site = Seq(seq[pos:pos+motif.length], alphabet=motif.alphabet)
                rsite = site.reverse_complement()
                strand = ('+' if motif.pssm.calculate(site) >= motif.pssm.calculate(rsite)
                          else '-')
                print g, score, site if strand=='+' else rsite, strand,
                print '[%d, %d]' % (pos+loc[0]+1, pos+loc[0]+motif.length)
