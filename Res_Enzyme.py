import sys
import subprocess
import argparse

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
    Reads a text file containing submission data and stores each submission as a dictionary in a list.

    Parameters:
    - file_path (str): The path to the input text file.

    Returns:
    list: A list of dictionaries, where each dictionary represents a submission with the following keys:
        - 'sub_id' (str): Submission ID.
        - 'enzymes' (list): List of enzymes.
        - 'sequence' (str): DNA sequence.
        - 'snp_pos' (int): SNP position.
        - 'mut' (str): Alternate allele.
        - 'rs_num' (str): Rs number.
        - 'chrom' (str): Chromosome.

    Example:
    ```python
    file_path = "path/to/your/file.txt"
    submissions = read_txt(file_path)
    ```
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
def replace_char_at_index(string: str ,
                          index: int,
                          char: str):
    """
    Replaces a character at a specified index in the given string with a new character.

    Parameters:
    - string (str): The input string.
    - index (int): The index of the character to be replaced (1-based index).
    - char (str): The new character to replace the existing one.

    Returns:
    str: The updated string with the character replaced.

    Example:
    ```python
    input_str = "example"
    updated_str = replace_char_at_index(input_str, 3, "x")
    # Result: "exxmple"
    ```
    """
    return string[:index-1] + char + string[index:]


#--------------------------------------------------------------------------------
def get_snps_flank(rs_number, flank_length=40):
    """
    Retrieves the flank region sequence of a SNP based on its rs_number.

    Parameters:
    - rs_number (str): The rs_number of the SNP.
    - flank_length (int): The length of the flank region on either side of the SNP (default is 40).

    Returns:
    tuple: A tuple containing:
        - str: The flank region sequence.
        - int: The total length of the retrieved sequence (including both flanks).
        - str: The alternate allele introduced by the SNP.

    Example:
    ```python
    rs_number = "rs123456"
    flank_sequence, total_length, spdi = get_snps_flank(rs_number, flank_length=50)
    # Result: ("ATCG...TGC", 101, "T")
    ```
    
    Note:
    - Requires the Biopython library for sequence retrieval.
    - Ensure that the Entrez module from Biopython is configured before using this function.
    """
    snp_id = rs_number[2:]
    handle = Entrez.esummary(db="snp", id=snp_id)
    record = Entrez.read(handle)
    res = record['DocumentSummarySet']['DocumentSummary'][0]
    snp_pos = int(res['CHRPOS_SORT'])
    if len(res['GLOBAL_MAFS'])>0:
        mut = res['GLOBAL_MAFS'][0]['FREQ'][0]
    else:
        mut = res['SPDI'].split(',')[0][-1]
    handle = Entrez.esearch(db="nucleotide", term=res['ACC'])
    record = Entrez.read(handle)
    ids = record['IdList']
    handle = Entrez.efetch(db="nucleotide", id=ids[0],
                           strand=1, seq_start=snp_pos-flank_length,
                           seq_stop=snp_pos+flank_length,
                           rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return str(record.seq), flank_length+1, mut


#--------------------------------------------------------------------------------
def get_loc_flank(chromosome, location, flank_length=40):
    """
    Retrieves the flank region sequence of a nucleotide based on its chromosome and location.

    Parameters:
    - chromosome (str, int): The chromosome identifier (e.g., '1', 'X', 'Y').
    - location (int): The position of the nucleotide on the chromosome.
    - flank_length (int): The length of the flank region on either side of the nucleotide (default is 40).

    Returns:
    tuple: A tuple containing:
        - str: The flank region sequence.
        - int: The total length of the retrieved sequence (including both flanks).

    Example:
    ```python
    chromosome = 'X'
    nucleotide_location = 123456
    flank_sequence, total_length = get_loc_flank(chromosome, nucleotide_location, flank_length=50)
    # Result: ("ATCG...TGC", 101)
    ```
    
    Note:
    - Requires the Biopython library for sequence retrieval.
    - Ensure that the Entrez module from Biopython is configured before using this function.
    """
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
    Generates a list of candidate sequences with potential Restriction Fragment Length Polymorphism (RFLP) variants.

    Parameters:
    - sequence (str): The original DNA sequence containing the SNP.
    - mut (str): The alternate allele introduced by the SNP.
    - snp_pos (int): The position of the SNP on the sequence.
    - enzyme (str): The name of the restriction enzyme.

    Returns:
    list: A list of candidate sequences that may exhibit RFLP variations.

    Example:
    ```python
    original_sequence = "AAATTCCC" 
    alternate_allele = "C" 
    snp_position = 4       # --> "AAACTCCC"
    enzyme_name = "EcoRI"  # site: "GAATTC" --> Cannot differentiate between ref_seq and mut_seq
    rflp_variants = generate_possible_RFLP_seq_variants(original_sequence, alternate_allele, snp_position, enzyme_name)
    # Result: ['GAATTCCC'] # --> Now it can be detected by EcoRI
    ```

    Note:
    - Requires the Biopython library for enzyme-related operations.
    - Ensure that the Restriction module from Biopython is configured before using this function.
    - The enzyme parameter should be provided as a string (e.g., 'EcoRI').
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
def is_RFLP(sequence, mut, snp_pos, enzyme):
    """
    Determines whether a given restriction enzyme can catalyze either the reference or mutant sequence
    and makes a distinction between them based on Restriction Fragment Length Polymorphism (RFLP).

    Parameters:
    - sequence (str): The original DNA sequence.
    - mut (str): The alternate allele introduced by a SNP.
    - snp_pos (int): The position of the SNP on the sequence.
    - enzyme (str): The name of the restriction enzyme.

    Returns:
    int: Returns -1 if the SNP sequence is cleaved, 1 if the reference sequence is cleaved, and 0 if none of them are cleaved.

    Example:
    ```python
    original_sequence = "GAATTCCCC"
    alternate_allele = "G"
    snp_position = 4
    enzyme_name = "EcoRI"  # site: "GAATTC"
    cleavage_result = is_RFLP(original_sequence, alternate_allele, snp_position, enzyme_name)
    # Result: 1 # Because EcoRI can catalyze the original_sequence but cannot catalyze the mutant sequence
    ```

    Note:
    - Requires the Biopython library for enzyme-related operations.
    - Ensure that the Restriction module from Biopython is configured before using this function.
    - The enzyme parameter should be provided as a string (e.g., 'EcoRI').
    """  
    mut_sequence = replace_char_at_index(sequence, snp_pos, mut)
    ref_seq = Seq(sequence)
    mut_seq = Seq(mut_sequence)
    enzyme = getattr(Restriction, enzyme)
    ref_res = enzyme.search(ref_seq)
    mut_res = enzyme.search(mut_seq)
    return len(ref_res) - len(mut_res)


#--------------------------------------------------------------------------------
def alignment(seq1, seq2):
    """
    Creates a representation of the alignment between two sequences with a single SNP difference.

    Parameters:
    - seq1 (str): The first DNA sequence.
    - seq2 (str): The second DNA sequence.

    Returns:
    list: A list containing the alignment representation with one SNP difference.
        - The SNP representation with the position and substitution (e.g., "25A>T").
        - The first DNA sequence.
        - The alignment indicator (| for matching, - for SNP position).
        - The second DNA sequence.

    Example:
    ```python
    sequence1 = "AATTCGAGAATGTA"
    sequence2 = "AATTCGTGAATGTA"
    alignment_result = alignment(sequence1, sequence2)
    # Result: ["7A>T",
    #          "AATTCGAGAATGTA",
    #          "||||||-|||||||",
    #          "AATTCGTGAATGTA"]
    ```
    
    Note:
    - The function assumes that the sequences have a single SNP difference.
    - The alignment indicator uses '|' for matching positions and '-' for SNP position.
    """
    l2 = '|'*len(seq1)
    for idx, (n1, n2) in enumerate(zip(seq1, seq2)):
        if n1!=n2:
            l2 = replace_char_at_index(l2, idx+1, '-')
            break
    return [f'{idx+1}{n1}>{n2}', seq1, l2, seq2]


#--------------------------------------------------------------------------------
def snp_location(sequence, mut, snp_pos):
    """
    Displays the reference and mutant DNA sequences with the position of a Single Nucleotide Polymorphism (SNP).

    Parameters:
    - sequence (str): The original DNA sequence.
    - mut (str): The alternate allele introduced by the SNP.
    - snp_pos (int): The position of the SNP on the sequence.

    Returns:
    list: A list containing the representation of the SNP location.
        - The SNP representation with the position and substitution (e.g., "25A>T").
        - The reference DNA sequence.
        - The alignment indicator (| for matching, - for SNP position).
        - The mutant DNA sequence.

    Example:
    ```python
    original_sequence = "AATTCGAGAATGTA"
    alternate_allele = "T"
    snp_position = 7
    snp_repr = snp_location(original_sequence, alternate_allele, snp_position)
    # Result: ["7A>T",
    #          "AATTCGAGAATGTA",
    #          "||||||-|||||||",
    #          "AATTCGTGAATGTA"]
    ```

    Note:
    - The function assumes a Single Nucleotide Polymorphism (SNP) at the specified position.
    - The alignment indicator uses '|' for matching positions and '-' for the SNP position.
    """
    mut_sequence = replace_char_at_index(sequence, snp_pos, mut)
    l1 = sequence
    l3 = mut_sequence
    l2 = replace_char_at_index('|'*len(sequence), snp_pos, '-')
    return [f'{snp_pos}{sequence[snp_pos-1]}>{mut}', l1, l2, l3]


#--------------------------------------------------------------------------------
def restriction_sites(sequence, enzyme):
    """
    Displays the locations of restriction enzyme sites and the cleaved DNA sequence.

    Parameters:
    - sequence (str): The DNA sequence to be analyzed.
    - enzyme (str): The name of the restriction enzyme.

    Returns:
    list: A list containing the representation of restriction enzyme sites and cleaved sequence.
        - The cleaved DNA sequence with spaces separating the cleavage sites.
        - The alignment indicator with symbols representing cleavage sites.
        - The reverse complement of the cleaved DNA sequence with spaces separating the cleavage sites.

    Example:
    ```python
    dna_sequence = "AAGAATTCCC"
    enzyme_name = "EcoRI"
    restriction_repr = restriction_sites(dna_sequence, enzyme_name)
    # Result: ['AAG AATTCCC',
    #          '|||⌃---⌄|||',
    #          'TTCTTAA GGG']
    ```

    Note:
    - Requires the Biopython library for enzyme-related operations.
    - Ensure that the Restriction module from Biopython is configured before using this function.
    - The enzyme parameter should be provided as a string (e.g., 'EcoRI').
    - The alignment indicator with symbols representing cleavage sites ('⥯' for blunt cleavage, '⌃-⌄' for overhang sticky end).
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
    """
    Stores the results of RFLP-related analyses in a dictionary.

    Parameters:
    - sequence (str): The original DNA sequence.
    - mut (str): The alternate allele introduced by a SNP.
    - snp_pos (int): The position of the SNP on the sequence.
    - enzyme (str): The name of the restriction enzyme.

    Returns:
    dict: A dictionary containing the results of various analyses:
        - 'ref_seq_RFLP' (int): Result of is_RFLP for the reference sequence (-1, 0, 1).
        - 'snp_res' (list): Result of snp_location with the SNP representation and aligned sequences.
        - 'restriction_sites_res' (list): Result of restriction_sites with cleaved sequence representation.
        - 'ref_new_align' (list): Result of alignment for the reference and a new variant sequence.

    Example:
    ```python
    original_sequence = "GAATTCCCC"
    alternate_allele = "G"
    snp_position = 4
    enzyme_name = "EcoRI"  # site: "GAATTC"
    analysis_results = seq_variants_RFLP(original_sequence, alternate_allele, snp_position, enzyme_name)
    # Result: {'ref_seq_RFLP': 1,
               'snp_res': ["4A>G", "GAATTCCCC", '|||-|||||', "GAAGTCCCC"],
               'restriction_sites_res': ['AAG AATTCCC', '|||⌃---⌄|||','TTCTTAA GGG']}
    ```

    Note:
    - Requires the Biopython library for enzyme-related operations.
    - Ensure that the Restriction module from Biopython is configured before using this function.
    - The enzyme parameter should be provided as a string (e.g., 'EcoRI').
    - The 'ref_seq_RFLP' result indicates whether the reference sequence is cleaved (-1) or not (0).
    """
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
    """
    Prepares the submission data for further analysis.

    Parameters:
    - sub (dict): The submission data, typically obtained from the read_txt function.

    Returns:
    dict: A dictionary containing the prepared submission data with the following keys:
        - 'sub_id' (str): Submission identifier.
        - 'enzymes' (list): List of enzymes for analysis.
        - 'sequence' (str): DNA sequence for analysis.
        - 'snp_pos' (int): Position of the SNP on the sequence.
        - 'mut' (str): Alternate allele introduced by the SNP.

    Example:
    ```python
    submission_data = {'sub_id': 'sub_1', 'enzymes': ['EcoRI', 'HindIII'], 'rs_num': 'rs123456'}
    prepared_data = prepare_submission(submission_data)
    # Result: {'sub_id': 'sub_1', 'enzymes': ['EcoRI', 'HindIII'], 'sequence': 'ATCG...TGC', 'snp_pos': 25, 'mut': 'A'}
    ```

    Note:
    - The function extracts information from the provided submission data and organizes it for further analysis.
    - If 'rs_num' is provided, it retrieves the sequence, SNP position, and mutation from the NCBI SNP database.
    - If 'sequence' is provided, it handles cases where the SNP information is embedded in the sequence.
    - If 'chrom' and 'snp_pos' are provided, it retrieves the sequence and SNP position from the specified chromosome location.
    - If 'enzymes' is ['ALL'], it includes all enzymes with a cleavage frequency less than 1024.
    """
    sub_id = sub['sub_id']
    enzymes = sub.get('enzymes', ['ALL'])
    if enzymes[0] == 'ALL':
        enzymes = [str(i) for i in Restriction.CommOnly if i.freq > 1024]
    rs_num = sub.get('rs_num')
    if rs_num:
        sequence, snp_pos, mut = get_snps_flank(rs_num, 40)
        mut = sub.get('mut', mut)
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
    """
    Generates a report as a string for all obtained results.

    Parameters:
    - results (dict): Dictionary containing results obtained from various analyses.

    Returns:
    str: A formatted report string summarizing the analysis results.

    Example:
    ```python
    analysis_results = {'ref_seq_RFLP': -1, 'snp_res': ["4A>G", "GAATTCCCC", '|||-|||||', "GCAATCCCC"], ...}
    report_string = make_report(analysis_results)
    # Result: "Analysis Results:\nRef Seq RFLP: -1\nSNP Location: 4A>G\n..."
    ```

    Note:
    - The function takes a dictionary of results and formats them into a human-readable report.
    - The report includes information such as RFLP results, SNP location, restriction sites, etc.
    - Customize the function based on the specific results obtained from your analyses.
    """
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
def processing(input_file, output_file, email):
    """
    Processes DNA sequence submissions from an input file and generates a report.

    Parameters:
    - input_file (str): Path to the input file containing DNA sequence submissions in TXT format.
    - output_file (str): Path to the output file where the report will be saved.
    - email (str): Email address used for accessing NCBI Entrez services.

    Returns:
    None: The function processes submissions and writes the report to the specified output file.

    Example:
    ```python
    input_path = "submissions.txt"
    output_path = "report.txt"
    user_email = "user@example.com"
    processing(input_path, output_path, user_email)
    ```

    Note:
    - The function reads DNA sequence submissions from the input file, processes them, and generates a report.
    - Each submission is prepared using the prepare_submission function.
    - The report includes details on RFLP results, SNP location, restriction sites, etc.
    - The final report is written to the specified output file.
    """
    Entrez.email = email
    submissions = read_txt(input_file)
    prepared_submissions = [prepare_submission(sub) for sub in submissions]
    report = ''
    for sub in prepared_submissions:
        report += report_submission(sub)
        report += f"{'#'*150}\n"
    with open(output_file,'w', encoding='utf-8') as f:
        f.write(report) 
        f.close()


#--------------------------------------------------------------------------------
def main():
    """
    Main entry point for the RFLP enzyme analysis tool.

    Parses command-line arguments, processes DNA sequence submissions, and generates a report.

    Usage:
    ```
    python Res_Enzyme.py <submission_file> <result_file> [--email <your_email>]
    ```

    Parameters:
    - submission_file (str): Path to the input file containing DNA sequence submissions in TXT format.
    - result_file (str): Path to the output file where the report will be saved.
    - --email (str, optional): Your email address required for accessing NCBI Entrez services.
      If not provided, the default email 'A.N.Other@example.com' will be used.

    Returns:
    None: The function processes submissions and writes the report to the specified output file.

    Example:
    ```bash
    python Res_Enzyme.py submissions.txt report.txt --email user@example.com
    ```

    Note:
    - This tool analyzes DNA sequence submissions, identifies RFLP enzymes, and generates a report.
    - Ensure that the Biopython library and required dependencies are installed before running the tool.
    - The results are saved to the specified output file.
    """
    parser = argparse.ArgumentParser(description='This is a tool for finding RFLP enzymez.')
    parser.add_argument('submission', type=str, help='submission file(txt)')
    parser.add_argument('result', type=str, help='result file name(txt)')
    parser.add_argument('--email', type=str, help='Your email address', default='A.N.Other@example.com')
    args = parser.parse_args()

    input_file = args.submission 
    output_file = args.result
    email = args.email   
    processing(input_file, output_file, email)
    print(f"Processing complete. Results saved to {output_file}")

if __name__ == "__main__":
    main()
