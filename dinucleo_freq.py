from Bio.Seq import Seq
from Bio import Entrez, SeqIO
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np



def nucleotide_freq(genome):
    
    nucleotides = {"A":0, "C":0, "G":0, "T":0}
    
    for dna_forward_seq in tqdm(genome):
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
    
    total = sum(nucleotides.values())

    for nucleotide in nucleotides:
        nucleotides[nucleotide] = round(nucleotides[nucleotide] / total, 4) 
    
    return nucleotides


def dinucleotide_freq_theo(nucl_freq):
    
    dinucleotides = ["TT", "GT", "CT", "AT", "TG", "GG", "CG", "AG", "TC", "GC", "CC", "AC", "TA", "GA", "CA", "AA"]
    theo_freq = {}
    
    for dinucleotide in dinucleotides:
        theo_freq[dinucleotide] = nucl_freq[dinucleotide[0]]*nucl_freq[dinucleotide[1]]

    return theo_freq


def dinucleotide_freq(genome):
    
    di_counts = {}
    
    dinucleotides = ["TT", "GT", "CT", "AT", "TG", "GG", "CG", "AG", "TC", "GC", "CC", "AC", "TA", "GA", "CA", "AA"]
    
    for dinucleotide in dinucleotides:
        di_counts[dinucleotide] = 0
     
    for dna_forward_seq in genome:
    
        dna_forward_seq = Seq(dna_forward_seq)

        for i in tqdm(range(len(dna_forward_seq)-1)):
            di_forward = dna_forward_seq[i:i+2]
            di_reversed = di_forward.reverse_complement()
            if di_forward in di_counts:
                di_counts[di_forward] += 1
            if di_reversed in di_counts:
                di_counts[di_reversed] += 1

    total_di = sum(di_counts.values())

    freq = {}
    for di in sorted(di_counts.keys()):
        freq[di] = round(di_counts[di] / total_di, 4) 

    return freq

def obs_vs_theo_dinuc(fastafile):
    
    print(f'Analysing {fastafile[8:-6]} genome \n')
    
    n = 0
    genome = []
    for seq_record in SeqIO.parse(fastafile, "fasta"):
        genome.append(seq_record.seq)
        n += len(seq_record.seq)
    
    print(f'The genome consists of {len(genome)} chromosomes or scaffolds of total size of {round(n/1000000,3)} Mb\n')
        
    di_freq = dinucleotide_freq(genome)
    nucl_freq = nucleotide_freq(genome)
    di_freq_theo = dinucleotide_freq_theo(nucl_freq)
    
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