import gzip
import json
import requests

from Bio.Seq import Seq
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy.stats import ttest_rel
from scipy.stats import chi2_contingency
import tarfile
from tqdm import tqdm


def count_fasta_records(gzipfastafile):
    """Calculating the number of records in a .fa.gz file.

    Parameters
    ----------
    gzipfastafile : str
        .fa.gz file with genome records.

    Returns
    -------
    count : int
        The number of records.
    """
    count = 0
    with gzip.open(gzipfastafile, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            count += 1
            
    return count


def count_nucleotides(seq):
    """Calculating the number nucleotides of each type in the sequence.

    Parameters
    ----------
    seq : str
        Nucleotide sequence.

    Returns
    -------
    nucl_counts : dic
        The number of nucleotides in the sequence for each nucleotide.
    """
    
    nucl_counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    
    for nucleotide in nucl_counts:
        nucl_counts[nucleotide] = seq.count(nucleotide)
    
    a = nucl_counts["A"] 
    t = nucl_counts["T"]
    c = nucl_counts["C"]
    g = nucl_counts["G"]

    nucl_counts["A"] += t
    nucl_counts["T"] += a
    nucl_counts["C"] += g
    nucl_counts["G"] += c 
    
    return nucl_counts


def count_dinucleotide_freq_theo(nucl_freq):
    """Calculating the theoretical dinucleotide frequency for each type of dinucleotide.

    Parameters
    ----------
    nucl_freq : dic
        Nucleotide frequency for each type of nucleotide.

    Returns
    -------
    theo_freq : dic
        Theoretical dinucleotide frequency for each type of dinucleotide.
    
    """
    dinucleotides = ["TT", "GT", "CT", "AT", "TG", "GG", "CG", "AG", "TC", "GC", "CC", "AC", "TA", "GA", "CA", "AA"]
    theo_freq = {}
    
    for dinucleotide in dinucleotides:
        theo_freq[dinucleotide] = round(nucl_freq[dinucleotide[0]]*nucl_freq[dinucleotide[1]], 4)

    return theo_freq


def count_dinucleotides(seq):
    """Calculating the number dinucleotides of each type in the sequence.

    Parameters
    ----------
    seq : str
        Nucleotide sequence.

    Returns
    -------
    di_counts : dic
        The number of dinucleotides in the sequence for each nucleotide.
    """

    di_counts = {"TT": 0, "GT": 0, "CT": 0, "AT": 0, "TG": 0, "GG": 0, "CG": 0, "AG": 0, "TC": 0, "GC": 0, "CC": 0, "AC": 0, "TA": 0, "GA": 0, "CA": 0, "AA": 0}
    compliments = {"A": "T", "T": "A", "C": "G", "G": "C"}
    
    try:
        dna_forward_seq = seq.seq
    except:
        dna_forward_seq = seq

    for i in range(len(dna_forward_seq)-1):
        di_forward = dna_forward_seq[i:i+2]
        if di_forward[0] in compliments and di_forward[1] in compliments:
            di_reversed = compliments[di_forward[1]]+compliments[di_forward[0]]
            di_counts[di_forward] += 1
            di_counts[di_reversed] += 1

    return di_counts


def obs_vs_theo_dinuc(gzipfastafile):
    """Calculating the theoretical and observed dinucleotide frequency for each type of dinucleotide from a .fa.gz file.

    Parameters
    ----------
    gzipfastafile : str
        .fa.gz file with genome records.

    Returns
    -------
    di_freq : dic
        Observed dinucleotide frequency for each type of dinucleotide.
    di_freq_theo : dic 
        Theoretical dinucleotide frequency for each type of dinucleotide.
    """
    
    try:
        print(f"😼Analysing {gzipfastafile.split('.')[0].split('_')[0]+' '+gzipfastafile.split('.')[0].split('_')[1]} genome")
    except:
        print("😼Analysing " + str(gzipfastafile).split("'")[1][:-4].split("_")[0] + " " + str(gzipfastafile).split("'")[1][:-4].split("_")[1] + " genome")
    
    nb_rec = count_fasta_records(gzipfastafile)
    rest = nb_rec
    print(f"🙀Total number of records: {nb_rec}")    
    
    total_di = 0
    total_nu = 0
    di_count = {"TT": 0, "GT": 0, "CT": 0, "AT": 0, "TG": 0, "GG": 0, "CG": 0, "AG": 0, 
                "TC": 0, "GC": 0, "CC": 0, "AC": 0, "TA": 0, "GA": 0, "CA": 0, "AA": 0}
    nucl_count = {"A": 0, "C": 0, "G": 0, "T": 0}
    
    with gzip.open(gzipfastafile, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            rest -= 1
            if nb_rec < 20:
                print(f"😸Analysing record {nb_rec-rest}: {record.id}, sequence length: {len(record.seq)}, number of records rest: {rest}")
            else:
                if (nb_rec-rest-1) % round(nb_rec/20, 0) == 0:
                    print(f"😸Analysing record {nb_rec-rest}: {record.id}, sequence length: {len(record.seq)}, number of records rest: {rest}\n...")
             
            dinuc = count_dinucleotides(record)
            for key in dinuc:
                di_count[key] += dinuc[key]
            nuc = count_nucleotides(record)
            for key in nuc:
                nucl_count[key] += nuc[key]
                
            #if len(record) > 1000000:
            #    test_significance(record)
                
                
    total_di = sum(di_count.values())
    total_nu = sum(nucl_count.values())  
    
    di_freq = {"TT": 0, "GT": 0, "CT": 0, "AT": 0, "TG": 0, "GG": 0, "CG": 0, "AG": 0, 
                "TC": 0, "GC": 0, "CC": 0, "AC": 0, "TA": 0, "GA": 0, "CA": 0, "AA": 0}
    
    for di in sorted(di_count.keys()):
        di_freq[di] = round(di_count[di] / total_di, 4) 

    nucl_freq = {"A": 0, "C": 0, "G": 0, "T": 0} 
    for nucl in nucl_count:
        nucl_freq[nucl] = round(nucl_count[nucl] / total_nu, 4)  
    di_freq_theo = count_dinucleotide_freq_theo(nucl_freq)

    print(f"😻Completed")
    
    return di_freq, di_freq_theo


def test_significance(record):
    """Testing the significance of the difference between observed and theoretical dinucleotide frequencies.

    Parameters
    ----------
    record : record (biopython)
        A fasta record.
    """
    
    print(f"😸Sampling {record.id}")
    
    sequence = record.seq

    samples = []
    sample_index = random.sample(range(len(sequence)-1000), 1000)
    sample = ''
    for i in sample_index:
        samples.append(sequence[i:i+1000])
        
    di_freq_distribution = {"TT": [], "GT": [], "CT": [], "AT": [], "TG": [], "GG": [], "CG": [], "AG": [], 
                                "TC": [], "GC": [], "CC": [], "AC": [], "TA": [], "GA": [], "CA": [], "AA": []}
        
    di_freq_theo_distribution = {"TT": [], "GT": [], "CT": [], "AT": [], "TG": [], "GG": [], "CG": [], "AG": [], 
                                     "TC": [], "GC": [], "CC": [], "AC": [], "TA": [], "GA": [], "CA": [], "AA": []}
    print(f"😻Sampling completed")
    print(f"😸Analising samples")

    for i in tqdm(range(len(samples))):
        sample = samples[i]
        total_di = 0
        total_nu = 0
        di_count = {"TT": 0, "GT": 0, "CT": 0, "AT": 0, "TG": 0, "GG": 0, "CG": 0, "AG": 0, 
                    "TC": 0, "GC": 0, "CC": 0, "AC": 0, "TA": 0, "GA": 0, "CA": 0, "AA": 0}
        nucl_count = {"A": 0, "C": 0, "G": 0, "T": 0}
    
        dinuc = count_dinucleotides(sample)
        for key in dinuc:
            di_count[key] += dinuc[key]
        nuc = count_nucleotides(sample)
        for key in nuc:
            nucl_count[key] += nuc[key]
                
        total_di = sum(di_count.values())
        total_nu = sum(nucl_count.values())  

        di_freq = {"TT": 0, "GT": 0, "CT": 0, "AT": 0, "TG": 0, "GG": 0, "CG": 0, "AG": 0, 
                    "TC": 0, "GC": 0, "CC": 0, "AC": 0, "TA": 0, "GA": 0, "CA": 0, "AA": 0}

        for di in sorted(di_count.keys()):
            di_freq[di] = round(di_count[di] / total_di, 4) 

        nucl_freq = {"A": 0, "C": 0, "G": 0, "T": 0} 
        for nucl in nucl_count:
            nucl_freq[nucl] = round(nucl_count[nucl] / total_nu, 4)  
        di_freq_theo = count_dinucleotide_freq_theo(nucl_freq)
        for key in di_freq:
            di_freq_distribution[key].append(di_freq[key])
            di_freq_theo_distribution[key].append(di_freq_theo[key])
    
    for key in di_freq_distribution:
        
        observed = np.array(di_freq_distribution[key])
        theo = np.array(di_freq_theo_distribution[key]) 

        contingency_table = np.array([observed, theo])

        statistic, p_val, degrees_of_freedom, expected_values = chi2_contingency(contingency_table)


        print(f"Chi2 test results for {key}:")
        if p_val < 0.05:
            print(f"The difference is significant (p-value: {p_val})\n")
        else:
            print(f"The difference is not significant (p-value: {p_val})\n")

    
    