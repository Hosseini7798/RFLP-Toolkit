import sys
import os
import subprocess

try:
    import Bio 
except ModuleNotFoundError:
    print('The Biopython module is NOT installed. Wait for installing')

    # 👇️ optionally install module
    python = sys.executable
    subprocess.check_call(
        [python, '-m', 'pip', 'install', 'biopython'],
        stdout=subprocess.DEVNULL
    )
finally:
    from Bio import Restriction
    from Bio.Seq import Seq
    from Bio import Entrez, SeqIO

#--------------------------------------------------------------------------------
def read_txt(file_path):
    """
    Read DNA sequence submissions from a text file and parse the information.

    The file should contain submissions formatted with lines starting with '>'
    and four elements (Enzymes, Sequence, SNP_position, Alt_allele) for each submission.

    Parameters:
    - file_path (str): The path to the text file containing DNA sequence submissions.

    Returns:
    list: A list of dictionaries, each representing a DNA sequence submission.
          Each dictionary contains 'sub_id', 'Enzymes', 'Sequence', 'SNP_position', and 'Alt_allele' keys.
    """    
    data = []
    with open(file_path, 'r') as file:
        entry = {}
        for line in file:
            line = line.strip()
            if line.startswith("> Sub_id:"):
                if entry:
                    data.append(entry)
                    entry = {}
                entry = {'sub_id': line.split(': ')[1]}
            elif line.startswith("Enzymes:"):
                entry['enzymes'] = line.split(': ')[1].split(', ')
            elif line.startswith("Sequence:"):
                entry['sequence'] = line.split(': ')[1]
            elif line.startswith("SNP_position:"):
                entry['snp_pos'] = int(line.split(': ')[1])
            elif line.startswith("Alt_allele:"):
                entry['mut'] = line.split(': ')[1]
            elif line.startswith("Rs_number:"):
                entry['rs_num'] = line.split(': ')[1]
            elif line.startswith("Chromosome:"):
                entry['chrom'] = line.split(': ')[1]  
            else:
                continue
        data.append(entry)
    return data


#--------------------------------------------------------------------------------
def highlight_mutation(sequence, mutation_position, color='\033[91m', reset='\033[0m'):
    """
    Highlight the mutation position in a DNA sequence.

    Parameters:
    - sequence (str): The DNA sequence.
    - mutation_position (int): The position of the mutation.
    - color (str, optional): ANSI escape code for highlighting color. Default is red.
    - reset (str, optional): ANSI escape code to reset color. Default is to reset.

    Returns:
    str: The sequence with the mutation position highlighted.
    """
    if mutation_position < 0 or mutation_position >= len(sequence):
        raise ValueError("Invalid mutation position")

    return (
        sequence[:mutation_position]
        + f"{color}{sequence[mutation_position]}{reset}"
        + sequence[mutation_position + 1:])


#--------------------------------------------------------------------------------
def replace_char_at_index(string: str ,
                          index: int,
                          char: str):
    """
    Replace the character at the specified index in the given string.

    Parameters:
    - string (str): The input string.
    - index (int): The index at which the character should be replaced.
    - char (str): The character to substitute at the specified index.

    Returns:
    str: The modified string with the character replaced at the specified index.
    """
    return string[:index-1] + char + string[index:]


#--------------------------------------------------------------------------------
def insert_char_at_index(string: str ,
                          index: int,
                          char: str):
    
    return string[:index-1] + char + string[index-1:]


#--------------------------------------------------------------------------------
def get_snps_flank(rs_number, flank_length=40):
    Entrez.email = 'hosseini7798@gmail.com'
    snp_id = rs_number[2:]
    handle = Entrez.esummary(db="snp", id=snp_id)
    record = Entrez.read(handle)
    res = record['DocumentSummarySet']['DocumentSummary'][0]
    snp_pos = int(res['CHRPOS_SORT'])
    handle = Entrez.esearch(db="nucleotide", term=res['ACC'])
    record = Entrez.read(handle)
    ids = record['IdList']
    handle = Entrez.efetch(db="nucleotide", id=ids[0],
                           strand=1, seq_start=snp_pos-flank_length,
                           seq_stop=snp_pos+flank_length,
                           rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return str(record.seq), flank_length+1, res['SPDI'].split(':')[-1]


#--------------------------------------------------------------------------------
def get_loc_flank(chromosome, location, flank_length=40):
    chrom_id = {'1': 568815597,'2': 568815596,'3': 568815595,'4': 568815594,'5': 568815593,
                '6': 568815592,'7': 568815591,'8': 568815590,'9': 568815589,'10': 568815588,
                '11': 568815587,'12': 568815586,'13': 568815585,'14': 568815584,'15': 568815583,
                '16': 568815582,'17': 568815581,'18': 568815580,'19': 568815579,'20': 568815578,
                '21': 568815577,'22': 568815576,'X': 568815575,'Y': 568815574}
    handle = Entrez.efetch(db="nucleotide", id=chrom_id[str(chromosome)],
                           strand=1, seq_start=location-flank_length,
                           seq_stop=location+flank_length,
                           rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return str(record.seq), flank_length+1


#--------------------------------------------------------------------------------
def generate_possible_RFLP_seq_variants(sequence, mut, snp_pos, enzyme):
    """
    Generate a set of DNA sequence variants near a specified SNP position.

    This function takes a reference DNA sequence, SNP position, and enzyme name,
    and returns a set of new DNA sequences with variations introduced near the SNP position.

    Parameters:
    - sequence (str): The reference DNA sequence.
    - snp_pos (int): The position of the Single Nucleotide Polymorphism (SNP).
    - enzyme (str): The name of the enzyme (e.g., 'EcoRI') to determine cutting sites.

    Returns:
    set: A set of DNA sequences with variations near the specified SNP position.
    """
    enzyme = getattr(Restriction, enzyme)
    near_len = len(enzyme.site)
    new_seqs = set()
    for i in range(1, near_len):
        for j in ['A', 'T', 'C', 'G']:
            new_seqs.add(replace_char_at_index(sequence, snp_pos-i, j))
            new_seqs.add(replace_char_at_index(sequence, snp_pos+i, j))
    new_seqs.remove(sequence)
    rflp_seqs = set()
    for new_seq in new_seqs:
        if is_RFLP(new_seq, mut, snp_pos, str(enzyme)):
            rflp_seqs.add(new_seq)
    return list(rflp_seqs)


#--------------------------------------------------------------------------------
def is_RFLP(sequence, mut, snp_pos, enzyme, report=False):
    """
    Determine whether an enzyme can detect a difference between a reference DNA sequence and a mutant sequence.

    This function takes a reference DNA sequence, a mutation allele, SNP position, enzyme name,
    and an optional report flag. It checks if the enzyme can act on either the reference sequence
    or the mutant sequence, indicating a detectable difference.

    Parameters:
    - sequence (str): The reference DNA sequence.
    - mut (str): The mutation allele to replace at the SNP position.
    - snp_pos (int): The position of the Single Nucleotide Polymorphism (SNP).
    - enzyme (str): The name of the enzyme (e.g., 'EcoRI') to check for recognition sites.
    - report (bool, optional): If True, generate a report (not implemented yet).

    Returns:
    bool: True if the enzyme can detect a difference, False otherwise.
    """
    mut_sequence = replace_char_at_index(sequence, snp_pos, mut)
    ref_seq = Seq(sequence)
    mut_seq = Seq(mut_sequence)
    enzyme = getattr(Restriction, enzyme)
    ref_res = enzyme.search(ref_seq)
    mut_res = enzyme.search(mut_seq)

    if report:
        # Placeholder for generating a report, not implemented yet
        pass
    
    return len(ref_res) - len(mut_res)


#--------------------------------------------------------------------------------
def alignment(seq1, seq2):
    l2 = '|'*len(seq1)
    for idx, (n1, n2) in enumerate(zip(seq1, seq2)):
        if n1!=n2:
            l2 = replace_char_at_index(l2, idx+1, '-')
            break
    return [f'{idx+1}{n1}>{n2}', seq1, l2, seq2]


#--------------------------------------------------------------------------------
def snp_location(sequence, mut, snp_pos):
    """
    Display the reference and mutant DNA sequences with the position of a Single Nucleotide Polymorphism (SNP).

    This function takes a reference DNA sequence, a mutation allele, and the SNP position.
    It prints the reference and mutant sequences with an arrow indicating the position of the SNP.

    Parameters:
    - sequence (str): The reference DNA sequence.
    - mut (str): The mutation allele to replace at the SNP position.
    - snp_pos (int): The position of the Single Nucleotide Polymorphism (SNP).

    Returns:
    None
    """
    mut_sequence = replace_char_at_index(sequence, snp_pos, mut)
    l1 = sequence
    l3 = mut_sequence
    l2 = replace_char_at_index('|'*len(sequence), snp_pos, '-')
    return [f'{snp_pos}{sequence[snp_pos-1]}>{mut}', l1, l2, l3]


#--------------------------------------------------------------------------------
def restriction_sites(sequence, enzyme):
    """
    Display the restriction enzyme sites and cleaved DNA sequence.

    This function takes a DNA sequence, an enzyme name, and an optional flag for returning the results.
    It prints or returns a representation of the DNA sequence with highlighted restriction enzyme sites and cleaved positions.

    Parameters:
    - sequence (str): The DNA sequence.
    - enzyme (str): The name of the enzyme (e.g., 'EcoRI') to check for recognition sites.
    - returning (bool, optional): If True, return a list of strings representing the sequence. Default is False.

    Returns:
    None or list: If returning is True, returns a list of strings representing the sequence.
    If returning is False, prints the sequence representation.
    """
    enzyme = getattr(Restriction, enzyme)
    seq = Seq(sequence)
    rseq = seq.reverse_complement()
    seq_split = enzyme.catalyse(seq)
    rseq_split = enzyme.catalyse(rseq)
    seq_cut = enzyme.search(seq)
    rseq_cut = enzyme.search(rseq)
    rseq_cut = [len(seq)-i+2 for i in rseq_cut]
    
    l1 = ' '.join([str(i) for i in seq_split])
    l3 = ' '.join([str(i)[::-1] for i in rseq_split[::-1]])
    if enzyme.is_blunt():
        l2 = '|'*len(l1)
        for idx, i in enumerate(enzyme.search(seq)):
            l2 = replace_char_at_index(l2, i+idx, '⥯')
    elif enzyme.is_5overhang():
        l2 = '|'*(len(seq)-len(seq_cut)*(len(enzyme.ovhgseq)-1))
        symbol = replace_char_at_index('⌃ ⌄', 2, '-'*(len(enzyme.ovhgseq)-1))
        for idx, i in enumerate(seq_cut):
            l2 = replace_char_at_index(l2, i+idx, symbol)
    else:
        l2 = '|'*(len(seq)-len(rseq_cut)*(len(enzyme.ovhgseq)-1))
        symbol = replace_char_at_index('⌄ ⌃', 2, '-'*(len(enzyme.ovhgseq)-1))
        for idx, i in enumerate(rseq_cut):
            l2 = replace_char_at_index(l2, i+idx, symbol)
    
    return [l1, l2, l3]


#--------------------------------------------------------------------------------
def seq_variants_RFLP(sequence, mut, snp_pos, enzyme):
    res = {}
    res['ref_seq_RFLP'] = is_RFLP(sequence, mut, snp_pos, enzyme)
    if res['ref_seq_RFLP']:
        res['snp_res'] = snp_location(sequence, mut, snp_pos)
        if res['ref_seq_RFLP']==1:
            res['restriction_sites_res'] = restriction_sites(sequence, enzyme)
        else: 
            mut_sequence = replace_char_at_index(sequence, snp_pos, mut)
            res['restriction_sites_res'] = restriction_sites(mut_sequence, enzyme) 
    else:
        new_seqs = generate_possible_RFLP_seq_variants(sequence, mut, snp_pos, enzyme)
        if len(new_seqs) == 0:
            return res
        res['ref_new_align'] = alignment(sequence, new_seqs[0])
        res['ref_seq_RFLP'] = is_RFLP(new_seqs[0], mut, snp_pos, enzyme)
        res['snp_res'] = snp_location(new_seqs[0], mut, snp_pos)
        if res['ref_seq_RFLP']==1:
            res['restriction_sites_res'] = restriction_sites(new_seqs[0], enzyme)
        else: 
            mut_sequence = replace_char_at_index(new_seqs[0], snp_pos, mut)
            res['restriction_sites_res'] = restriction_sites(mut_sequence, enzyme)
    return res


#--------------------------------------------------------------------------------
def prepare_submission(sub):
    sub_id = sub['sub_id']
    enzymes = sub.get('enzymes', ['ALL'])
    if enzymes[0] == 'ALL':
        enzymes = [str(i) for i in Restriction.CommOnly if i.freq > 1024]
    rs_num = sub.get('rs_num')
    if rs_num:
        sequence, snp_pos, mut = get_snps_flank(rs_num, 40)
    else:
        sequence = sub.get('sequence')
        if sequence:
            if sequence.__contains__('/'):
                snp_pos = sequence.find('[')+1
                mut = sequence[snp_pos+2]
                sequence = sequence.replace('[', '',)
                sequence = sequence[:snp_pos] + sequence[snp_pos+3:]
            else:
                snp_pos = sub['snp_pos']
                mut = sub['mut']
        else:
            chrom = sub.get('chrom')
            snp_pos = sub.get('snp_pos')
            mut = sub['mut']
            sequence, snp_pos = get_loc_flank(chrom, snp_pos, flank_length=40)
    return {'sub_id':sub_id, 'enzymes':enzymes,
            'sequence':sequence,'snp_pos':snp_pos,
            'mut':mut}


#--------------------------------------------------------------------------------
def report_submission(sub):
    for key, value in sub.items():
        globals()[key] = value
    enzymes_results = {}
    for i in enzymes:
        enzymes_results[i] = seq_variants_RFLP(sequence, mut, snp_pos, i)
    already_RFLP = [k for k, v in enzymes_results.items() if len(v)==3]
    not_RFLP_at_all = [k for k, v in enzymes_results.items() if len(v)==1]
    RFLP_after_change = [k for k, v in enzymes_results.items() if len(v)>1]

    report = f"> sub_id: {sub_id}\n"
    report += f"Selected enzymes: {', '.join(enzymes)}\n"
    report += f"All enzymes: {len(enzymes)}\n"
    report += f"Already RFLP enzyme for this SNP: {len(already_RFLP)}\n"
    report += f"Cannot used for RFLP at all: {len(not_RFLP_at_all)}\n"
    report += f"RFLP enzymes after change: {len(RFLP_after_change)}\n\n"
    for k, v in enzymes_results.items():
        if len(v) == 1:
            continue

        report += f"\tEnzyme: {k}\n"
        if v['ref_seq_RFLP']==1:
            report += f"\tCutting sequence: Reference\n\n"
        elif v['ref_seq_RFLP']==-1:
            report += f"\tCutting sequence: SNP\n\n"
        s = 'ref'
        if len(v)==4:
            s = 'new'
            report += f"\tAlignment of reference and new sequence: {v['ref_new_align'][0]}\n"
            report += f"\tref_seq: 5' {v['ref_new_align'][1]} 3'\n"
            report += f"\t            {v['ref_new_align'][2]}   \n"
            report += f"\tnew_seq: 5' {v['ref_new_align'][3]} 3'\n\n"
        if s == 'ref':
            report += '\tThis enzyme already could be used in RFLP.\n'
        report += f"\tDisplaying SNP: {v['snp_res'][0]}\n"
        report += f"\t{s}_seq: 5' {v['snp_res'][1]} 3'\n"
        report += f"\t            {v['snp_res'][2]}   \n"
        report += f"\tsnp_seq: 5' {v['snp_res'][3]} 3'\n\n"

        report += f"\tDisplaying restriction enzyme site:\n"
        report += f"\t{s}_seq: 5' {v['restriction_sites_res'][0]} 3'\n"
        report += f"\t            {v['restriction_sites_res'][1]}   \n"
        report += f"\tsnp_seq: 3' {v['restriction_sites_res'][2]} 5'\n"
        report += f"{'='*(len(sequence)+25)}\n"
    return report

#--------------------------------------------------------------------------------
def processing():
    path_list = os.listdir('submissions/')
    for path in path_list:
        submissions = read_txt(f'submissions/{path}')
        prepared_submissions = [prepare_submission(sub) for sub in submissions]
        report = ''
        for sub in prepared_submissions:
            report += report_submission(sub)
            report += f"{'#'*150}\n"
        with open(f'results/{path}','w', encoding='utf-8') as f:
            f.write(report) 
            f.close()


#--------------------------------------------------------------------------------
def main():
    processing()
if __name__ == "__main__":
    main()
