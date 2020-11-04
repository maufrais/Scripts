'''
@author: corinne.maufrais@pasteur.fr
Bioinformatics and Biostatistics Hub - Institut Pasteur
Paris, France
'''
import sys
import argparse, textwrap



def get_indexof_type(genotype_specification_type_and_order, geno_type):
    """
    Return the corresponding index of one word in a regex expression: 'GT:AD:DP'
    word GT --> 0
    word AD --> 1
    """
    try:
        fid = genotype_specification_type_and_order.split(':').index(geno_type)
    except:
        fid = -1
    return fid


def parse_vcf_info(genotype_datas, genotype_specification_type_and_order='GT:AD:DP:GQ:PL'):
    """
    Parse VCF line.
    Return genotype, coverage, depth_per_allele_by_sample and snp_qual
    """

    genotype_datas = genotype_datas.split(':') 
    genotype_gt = genotype_datas[get_indexof_type(genotype_specification_type_and_order, 'GT')]
    if get_indexof_type(genotype_specification_type_and_order, 'AD') != -1:
        try:
            depth_per_allele_by_sample_ad = map(int,genotype_datas[get_indexof_type(genotype_specification_type_and_order, 'AD')].split(','))
        except:
            depth_per_allele_by_sample_ad = []
    else:
        depth_per_allele_by_sample_ad = []
    try:
        if get_indexof_type(genotype_specification_type_and_order, 'DP') != -1:
            coverage_dp = float(genotype_datas[get_indexof_type(genotype_specification_type_and_order, 'DP')])  
        else:
            coverage_dp = 0.0
    except:
        coverage_dp = 0.0
    
    if get_indexof_type(genotype_specification_type_and_order, 'GQ') != -1:
        snp_qual = int(genotype_datas[get_indexof_type(genotype_specification_type_and_order, 'GQ')])
    else:
        snp_qual = 0
    return genotype_gt, coverage_dp, depth_per_allele_by_sample_ad, snp_qual


def get_snps_combine(fhin, min_ABHom=98, SNPqual_GQ=80, cov_DP=10, type_alz='high'):
    """
    Parse combine vcf file
    Return a dict with all informations
    """
    line = fhin.readline()
    while '#CHROM' not in line:
        line = fhin.readline()
    header = line.split()
    line = fhin.readline()
    datas ={}
    nb_pass = 0
    while line:
        keep_position = True
        line_datas = {}
        fld = line.split()
        chromo = fld[0]
        if chromo not in datas:
            datas[chromo]= {}
        position = int(fld[1])
        REF = fld[3]
        ALT = fld[4]
        if ',' in ALT:
            ALT = ALT.split(',')[0]
        filtered_tag = fld[6]
        genotype_specification_type_and_order = fld[8]
        if filtered_tag != '.' and filtered_tag != 'PASS':
            nb_pass +=1
            keep_position = False
        else:
            for i in range(9, len(header), 1):
                strain_name = header[i].split('.')[0]
                genotype_datas = fld[i]
                if genotype_datas != '.':
                    genotype_gt, coverage_dp, depth_per_allele_by_sample_ad, snp_qual = parse_vcf_info(genotype_datas, genotype_specification_type_and_order) 
                    base = ALT
                    infos = [genotype_gt, coverage_dp, depth_per_allele_by_sample_ad, snp_qual]
                    
                    if coverage_dp < cov_DP:
                        keep_position = False
                    if keep_position and depth_per_allele_by_sample_ad:
                        ABHom = float(depth_per_allele_by_sample_ad[1]) * 100 / coverage_dp
                        if ABHom < min_ABHom:
                            keep_position = False
                    if keep_position and snp_qual < SNPqual_GQ:
                        keep_position = False
                        
                    if not keep_position and type_alz == 'all':
                        base = '-'
                        keep_position = True
                else:
                    infos = ['',0.0,[],0]
                    base = REF
                line_datas[strain_name] = {'base': base, 'info':infos}
        
        if keep_position:
            for strain_name in line_datas.keys():
                if strain_name not in datas[chromo].keys():
                    datas[chromo][strain_name] = {position:{'base': line_datas[strain_name]['base'], 'info':line_datas[strain_name]['info']}}
                else:
                    datas[chromo][strain_name][position] = {'base': line_datas[strain_name]['base'], 'info':line_datas[strain_name]['info']}
        
        line = fhin.readline()
    
    nTot = 0
    for chromo in datas.keys():
        st_names = datas[chromo].keys()
        if st_names:
            nTot += len(datas[chromo][st_names[0]])
        for strain_name in st_names:
            for pos in datas[chromo][strain_name].keys():
                base = datas[chromo][strain_name][pos]['base']
                info = datas[chromo][strain_name][pos]['info']
    return datas, nTot


def make_sequence(strain_data):
    positions = strain_data.keys()
    positions.sort()
    seq = ''
    for pos in positions:
        seq += strain_data[pos]['base']
    return seq

def rewrite_datas(datas):
    new_datas = {}
    chromosomes = datas.keys()
    chromosomes.sort()
    for chromo in chromosomes:
        st_names = datas[chromo].keys()
        for strain_name in st_names:
            if strain_name in new_datas.keys():
                new_datas[strain_name] += make_sequence(datas[chromo][strain_name])
            else:
                new_datas[strain_name] = make_sequence(datas[chromo][strain_name])
    return new_datas


    
def write_fasta_ind(fhout, name, sequence, nb_pb_in_line=100):
    """
    write individual fasta file
    """
    nb_car = len(sequence)

    print_f = '>%s\n' % (name)

    i = 0
    while i < nb_car:
        if i < nb_pb_in_line:
            print_f += sequence[i:i + nb_pb_in_line] + '\n'
        else:
            print_f += sequence[i:i + nb_pb_in_line] + '\n'
        i += nb_pb_in_line
    if print_f[-1] == '\n':
        print_f = print_f[:-1]
    print >>fhout, print_f    

def write_fasta(fhout, datas):
    """
    write fasta file
    """
    strains = datas.keys()
    strains.sort()
    for strain_name in strains:
        name = strain_name.split('-')[1]
        write_fasta_ind(fhout, name, datas[strain_name])


def write_phylip(fhout, datas):
    """
    write phylip file
    """
    nb_strains = len(datas.keys())
    nb_pb_in_line = 60
    nb_car = len(datas.values()[0])
    print_f = ' %s %s\n' % (str(nb_strains), str(nb_car))
    i = 0
    while i < nb_car:
        if i < nb_pb_in_line:
            for strain_name in datas.keys():
                name = strain_name.split('-')[1]
                print_f += '%s ' % name + ' ' * (10 - len(name))
                print_f += datas[strain_name][i:i + nb_pb_in_line] + '\n'
            print_f += '\n'
        else:
            for strain_name in datas.keys():
                name = strain_name.split('-')[1]
                print_f += '\t' + datas[strain_name][i:i + nb_pb_in_line] + '\n'
            print_f += '\n'
        i += nb_pb_in_line
    print >>fhout, print_f

def write_gtx(fhout, datas):
    """
    Write all SNPs by positions for all strains
    """
    chromosomes = datas.keys()
    chromosomes.sort()
    st_names = datas[chromosomes[0]].keys()
    st_names.sort()
    print >>fhout, 'Chr-Position\tref\t' + '\t'.join([strain_name for strain_name in st_names])

    for chromo in chromosomes:
        st_names = datas[chromo].keys()
        st_names.sort()
        if st_names:
            l_pos = datas[chromo][st_names[0]].keys()
            l_pos.sort()
        else:
            l_pos = []
        for pos in l_pos:
            to_print = "%s-%s" % (chromo, pos)
            for strain_name in st_names:
                to_print += "\t%s" % (datas[chromo][strain_name][pos]['base'])
            print >>fhout, to_print

def write_nb_gap(fhout, datas):
    """
    Write number of filtered SNPs by positions
    """
    chromosomes = datas.keys()
    chromosomes.sort()
    st_names = datas[chromosomes[0]].keys()
    st_names.sort()
    for chromo in chromosomes:
        st_names = datas[chromo].keys()
        st_names.sort()
        if st_names:
            l_pos = datas[chromo][st_names[0]].keys()
            l_pos.sort()
        else:
            l_pos = []
        
        for pos in l_pos:
            to_print = ''
            for strain_name in st_names:
                to_print += "%s" % (datas[chromo][strain_name][pos]['base'])
            nbGap = to_print.count('-')
            if nbGap != 0:
                print >>fhout,  "%s\t%s\t%s" %  (chromo, pos, to_print.count('-'))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='variant_filtarion.py ',
                                     description="Filtering variant calls based on certain criteria. Variant not matching any of the DP, GQ and ABHom conditions will be considered bad and filtered out", 
                                     formatter_class=argparse.RawTextHelpFormatter)
    general_options = parser.add_argument_group(title="Options", description=None)
    general_options.add_argument("-i", "--vcffile", dest="vcf_combine",
                                 help="combine vcf file",
                                 metavar="FILE",
                                 type=file,
                                 )

    general_options.add_argument("-a", "--ABHom",
                                 dest="min_ABHom", metavar="INT",
                                 help="Homozygous SNP threshold (default: 90)",
                                 type=int,
                                 default=90
                                 )
    general_options.add_argument("-c", "--DP",
                                 dest="cov_DP", metavar="INT",
                                 help="Coverage DP threshold. (default: 10)",
                                 type=int,
                                 default=10)
    general_options.add_argument("-m", "--matrix",
                                 metavar='FILE',
                                 action='store',
                                 dest='out_table',
                                 type=argparse.FileType('w'),
                                 help='Report the retained SNPs by position and by strains in file')
    general_options.add_argument("-n", "--nb",
                                 metavar='FILE',
                                 action='store',
                                 dest='outfh_nbGap',
                                 type=argparse.FileType('w'),
                                 
                                 help='With the "all" variant filtration type, report the number of replace filtered variants by a gap in file.')
    general_options.add_argument("-t", "--type",
                                 dest="type", metavar="STR",
                                 help=textwrap.dedent("""Variant filtration type: 
- high:  Requires that all variants at a site not be filtered to be kept.
- all:   Don't require that all variants at a site not be filtered to be kept. 
         Replace filtered variants by a gap."""),
                                 choices=['high', 'all'],
                                 required=True)
    general_options.add_argument("-g", "--GQ",
                                 dest="SNPqual_GQ", metavar="INT",
                                 help="Normalized SNPs quality threshold (default: 80)",
                                 type=int,
                                 default=80)
    general_options.add_argument("-f", "--fasta",
                                 metavar='FILE',
                                 action='store',
                                 dest='fasta_outname',
                                 type=argparse.FileType('w'),
                                 
                                 help='Report SNPs in fasta file')
    general_options.add_argument("-p", "--phylip",
                                 metavar='FILE',
                                 action='store',
                                 dest='phylip_outname',
                                 type=argparse.FileType('w'),
                                 
                                 help='Report SNPs in phylip file')
    args = parser.parse_args()

    datas,nTot = get_snps_combine(args.vcf_combine, min_ABHom=args.min_ABHom, SNPqual_GQ=args.SNPqual_GQ, cov_DP=args.cov_DP, type_alz=args.type)

    datas_new = rewrite_datas(datas)
    if args.fasta_outname:
        write_fasta(args.fasta_outname, datas_new)
    if args.phylip_outname:
        write_phylip(args.phylip_outname, datas_new)
    if args.out_table:
        write_gtx(args.out_table, datas)
    if args.outfh_nbGap and args.type == 'all':
        write_nb_gap(args.outfh_nbGap, datas)
