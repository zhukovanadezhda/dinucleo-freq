from Bio.Seq import Seq
from Bio import Entrez, SeqIO
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import gzip



def count_nucleotides(seq):
    
    nucleotides = {"A":0, "C":0, "G":0, "T":0}
    
    for dna_forward_seq in tqdm(seq):
        for nucleotide in nucleotides:
                nucleotides[nucleotide] = dna_forward_seq.count(nucleotide)

    a = nucl_counts["A"] 
    t = nucl_counts["T"]
    c = nucl_counts["C"]
    g = nucl_counts["G"]

    nucl_counts["A"] += t
    nucl_counts["T"] += a
    nucl_counts["C"] += g
    nucl_counts["G"] += c 
    
    #total = sum(nucleotides.values())

    #for nucleotide in nucleotides:
    #    nucleotides[nucleotide] = round(nucleotides[nucleotide] / total, 4) 
    
    return nucleotides


def count_dinucleotide_freq_theo(nucl_freq):
    
    dinucleotides = ["TT", "GT", "CT", "AT", "TG", "GG", "CG", "AG", "TC", "GC", "CC", "AC", "TA", "GA", "CA", "AA"]
    theo_freq = {}
    
    for dinucleotide in dinucleotides:
        theo_freq[dinucleotide] = nucl_freq[dinucleotide[0]]*nucl_freq[dinucleotide[1]]

    return theo_freq


def count_dinucleotides(seq):
    
    di_counts = {}
    
    dinucleotides = ["TT", "GT", "CT", "AT", "TG", "GG", "CG", "AG", "TC", "GC", "CC", "AC", "TA", "GA", "CA", "AA"]
    
    for dinucleotide in dinucleotides:
        di_counts[dinucleotide] = 0
     
    for dna_forward_seq in seq:

        for i in tqdm(range(len(dna_forward_seq)-1)):
            di_forward = dna_forward_seq[i:i+2]
            di_reversed = di_forward[i+1]+di_forward[i+1]
            if di_forward in di_counts:
                di_counts[di_forward] += 1
            if di_reversed in di_counts:
                di_counts[di_reversed] += 1

    #total_di = sum(di_counts.values())

    #freq = {}
    #for di in sorted(di_counts.keys()):
    #    freq[di] = round(di_counts[di] / total_di, 4) 

    return di_counts

def obs_vs_theo_dinuc(gzipfastafile):
    
    print(f"Analysing {gzipfastafile.split('/')[-1]} genome \n")
    
    total_di = 0
    total_nu = 0
    di_count = {"TT":0, "GT":0, "CT":0, "AT":0, "TG":0, "GG":0, "CG":0, "AG":0, 
                "TC":0, "GC":0, "CC":0, "AC":0, "TA":0, "GA":0, "CA":0, "AA":0}
    nucl_count = {"A":0, "C":0, "G":0, "T":0}
  
    with gzip.open(gzipfastafile, "rt") as handle:
        for record in tqdm(SeqIO.parse(handle, "fasta")):
            dinuc = count_dinucleotides(record)
            for key in dinuc:
                di_count[key] += dinuc[key]
            nuc = count_nucleotides(record)
            for key in nuc:
                nucl_count[key] += nuc
            total_di += sum(di_count.values())
            total_nu += sum(nucl_count.values())
            
    di_freq = {"TT":0, "GT":0, "CT":0, "AT":0, "TG":0, "GG":0, "CG":0, "AG":0, 
               "TC":0, "GC":0, "CC":0, "AC":0, "TA":0, "GA":0, "CA":0, "AA":0}
    for di in sorted(di_count.keys()):
        di_freq[di] = round(di_count[di] / total_di, 4) 
        
    nucl_frec = {"A":0, "C":0, "G":0, "T":0} 
    for nucl in nucl_count:
        nucl_count[nucl] = round(nucl_count[nucl] / total_nu, 4) 
    
    di_freq_theo = count_dinucleotide_freq_theo(nucl_freq)

    
    #print(f'The genome consists of {len(genome)} chromosomes or scaffolds of total size of {round(n/1000000,3)} Mb\n')
 
    dinucleotides = []
    frequencies_obs = []
    frequencies_theo = []
    for key in di_freq:
        dinucleotides.append(str(key))
        frequencies_obs.append(di_freq[key])
        frequencies_theo.append(di_freq_theo[key])    

    bar_width = 0.45
    x_pos = np.arange(len(dinucleotides)) 

    plt.bar(x_pos - bar_width/2, frequencies_obs, width=bar_width, color='#00A3A6', zorder=2)
    plt.bar(x_pos + bar_width/2, frequencies_theo, width=bar_width, color='#ED6E6C',  zorder=2)

    plt.xlabel('Dinucleotide')
    plt.ylabel('Frequency')
    plt.title(f'Dinucleotide frequency in {fastafile[8:-6]} genome\n')
    plt.xticks(x_pos, dinucleotides, rotation=70) 


    obs_patch = mpatches.Patch(color='#00A3A6', label='observed')
    theo_patch = mpatches.Patch(color='#ED6E6C', label='theoretical')
    plt.legend(handles=[obs_patch, theo_patch])

    plt.grid(zorder=0, alpha=0.3)

    plt.show()
    
    return di_freq