#!/usr/bin/env python3
import sys, os, argparse, datetime, gzip

#
# Globals
#
DATE = datetime.datetime.now().strftime("%Y%m%d")
PROG = sys.argv[0].split('/')[-1]
DESC = 'Parse and merge the output from RepeatMasker.'

def parse_args():
    '''Command line arguments'''
    p = argparse.ArgumentParser(prog=PROG, description=DESC)
    p.add_argument('-c', '--cross-match', required=True,
                   help='(str) Path to the cross match (*.out) table from Repeat Masker.')
    p.add_argument('-d', '--divsum', required=True,
                   help='(str) Path to the divergence summary (*.divsum) table from Repeat Masker.')
    p.add_argument('-o', '--outdir', required=False, default='.',
                   help='(str) Path to output directory.')
    p.add_argument('-m', '--min-length', required=False, default=10, type=int,
                   help='(int) Minimum length of well characterized bases (wellCharLen) needed to keep a sequence  [default=10]')
    # Check input arguments
    args = p.parse_args()
    args.outdir = args.outdir.rstrip('/')
    if not os.path.exists(args.outdir):
        sys.exit(f"Error: '{args.outdir}' not found.")
    if not os.path.exists(args.cross_match):
        sys.exit(f"Error: '{args.cross_match}' not found.")
    if not os.path.exists(args.divsum):
        sys.exit(f"Error: '{args.divsum}' not found.")
    return args

def now():
    '''Print the current date and time.'''
    return f'{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'

class RepeatDivsum:
    '''Divergence summary statistics for a given repeat.'''
    def __init__(self, rep_class, rep_id, abs_len, well_char_len, kimura):
        self.r_class   = rep_class.split('/')[0]
        self.r_family  = rep_class
        self.r_id      = rep_id
        self.abs_len   = abs_len
        self.wchar_len = well_char_len
        self.kimura    = kimura/100
    def __str__(self):
        return f'{self.r_class} {self.r_familyq} {self.r_id} {self.abs_len} {self.wchar_len} {self.kimura:0.06g}'

class RepeatAnnot:
    '''Annotation metadata for a given repeat.'''
    def __init__(self, rep_name, rep_class, chromosome, start, end, sw_score):
        self.name     = rep_name
        self.r_class  = rep_class.split('/')[0]
        self.r_family = rep_class
        self.chrom    = chromosome
        self.start    = start
        self.end      = end
        self.score    = sw_score
    def __str__(self):
        return f'{self.name} {self.r_class} {self.r_family} {self.chrom} {self.start} {self.end} {self.score}'

def parse_divsum_table(divsum_f, min_len=10):
    '''Parse the divergence summary table.'''
    fh = open(divsum_f, 'r')
    if divsum_f.endswith('.gz'):
        fh = gzip.open(divsum_f, 'rt')
    print('\nParsing divergence summary table...')
    divsum = dict()
    records = 0
    discard = 0
    # Parse file
    for i, line in enumerate(fh):
        line = line.strip('\n')
        # Skip comments and empty lines
        if line.startswith('#'):
            continue
        if len(line) == 0:
            continue
        # Split into different fields.
        # The target fields have 5 entries
        fields = line.split('\t')
        if len(fields) != 5:
            continue
        # Skip headers or placeholders
        if fields[0] in {'Class', 'ARTEFACT', 'Simple_repeat', '-----'}:
            continue
        # Work with the desired records
        records += 1
        repeat_summary = RepeatDivsum(rep_class=fields[0],
                                       rep_id=fields[1],
                                       abs_len=int(fields[2]),
                                       well_char_len=int(fields[3]),
                                       kimura=float(fields[4]))
        # Remove entries if too small
        if repeat_summary.wchar_len < min_len:
            discard += 1
            continue
        # save for export
        divsum[repeat_summary.r_id] = repeat_summary
    print(f'    Extracted {records:,} records from divergence summary table file.\n    Discarded {discard:,} records.')
    return divsum

def parse_crossmatch_table(crossmatch_f):
    '''Parse the cross_match table output from RepeatMasker'''
    fh = open(crossmatch_f, 'r')
    if crossmatch_f.endswith('.gz'):
        fh = gzip.open(crossmatch_f, 'rt')
    print('Parsing cross match table...')
    cross_match = dict()
    records = 0
    kept = 0
    # Parse file
    for i, line in enumerate(fh):
        line = line.strip('\n')
        # Skip comments and empty lines
        if line.startswith('#'):
            continue
        if len(line) == 0:
            continue
        # Split into different fields, remove headers
        fields = line.split()
        if not fields[0].isnumeric():
            continue
        if len(fields) not in range(14,16):
            continue
        # Now, process the record lines
        records += 1
        sw_score = int(fields[0])
        chromosome = fields[4]
        start_bp = int(fields[5])
        end_bp = int(fields[6])
        re_name = fields[9]
        re_class = fields[10]
        annotation = RepeatAnnot(re_name, re_class, chromosome, start_bp, end_bp, sw_score)
        # Filter some of the unwanted classes
        if annotation.r_class in ('Simple_repeat', 'Low_complexity', 'Satellite'):
            continue
        # Add to the annotation dict
        kept += 1
        cross_match.setdefault(annotation.name, [])
        cross_match[annotation.name].append(annotation)
    print(f'    Read {records:,} records from the cross_match table file.\n    Retained {kept:,} records.')
    return cross_match

def merge_cross_divsum(cross_match, divsum, outdir='.'):
    '''Merge the cross_match and divsum datasets into a final output.'''
    print('\nMatching cross_match and divsum records...')
    fh = open(f'{outdir}/repeat_masked_merged.tsv', 'w')
    header = ['#Name', 'Class', 'Family', 'Chromosome', 'StartBP', 'EndBP', 'WellCharLen', 'Kimura','SwScore']
    header = '\t'.join(header)
    fh.write(f'{header}\n')
    matches = 0
    # Process each cross_match annotation
    for name in cross_match:
        name_matches = cross_match[name]
        for annotation in name_matches:
            assert isinstance(annotation, RepeatAnnot)
            # Find the divsum for the annotation
            divergence = divsum.get(annotation.name, None)
            if divergence is None:
                continue
            assert isinstance(divergence, RepeatDivsum)
            # Make sure the class and family match exactly
            if (annotation.r_class != divergence.r_class) or (annotation.r_family != divergence.r_family):
                continue
            # Process only fully matching records
            matches += 1
            row = f'{annotation.name}\t{annotation.r_class}\t{annotation.r_family}\t{annotation.chrom}\t{annotation.start}\t{annotation.end}\t{divergence.wchar_len}\t{divergence.kimura:0.6g}\t{annotation.score:0.6g}\n'
            fh.write(row)
    print(f'    Exported a total of {matches:,} matching records.')
    fh.close()

def main():
    print(f'{PROG} started on {now()}\n')
    args = parse_args()
    # Parse cross match table
    cross_match = parse_crossmatch_table(args.cross_match)
    # Parse divsum table
    divsum = parse_divsum_table(args.divsum, args.min_length)
    # Join cross_match and divsum
    merge_cross_divsum(cross_match, divsum, args.outdir)
    # Finish
    print(f'\n{PROG} finished on {now()}')


# Run Code
if __name__ == '__main__':
    main()
