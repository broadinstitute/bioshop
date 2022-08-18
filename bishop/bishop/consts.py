hg38_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz'

vcf_paths = {
    '1000g': 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz',
    'dbsnp': 'https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz',
}

AS_Fields = ['AS_BaseQRankSum', 'AS_FS', 'AS_InbreedingCoeff', 'AS_MQ', 'AS_MQRankSum', 'AS_QD', 'AS_ReadPosRankSum', 'AS_SOR']

chrom_aliases = {
 'chr1': ('1', 'CM000663.2', 'NC_000001.11'),
 'chr2': ('2', 'CM000664.2', 'NC_000002.12'),
 'chr3': ('3', 'CM000665.2', 'NC_000003.12'),
 'chr4': ('4', 'CM000666.2', 'NC_000004.12'),
 'chr5': ('5', 'CM000667.2', 'NC_000005.10'),
 'chr6': ('6', 'CM000668.2', 'NC_000006.12'),
 'chr7': ('7', 'CM000669.2', 'NC_000007.14'),
 'chr8': ('8', 'CM000670.2', 'NC_000008.11'),
 'chr9': ('9', 'CM000671.2', 'NC_000009.12'),
 'chr10': ('10', 'CM000672.2', 'NC_000010.11'),
 'chr11': ('11', 'CM000673.2', 'NC_000011.10'),
 'chr12': ('12', 'CM000674.2', 'NC_000012.12'),
 'chr13': ('13', 'CM000675.2', 'NC_000013.11'),
 'chr14': ('14', 'CM000676.2', 'NC_000014.9'),
 'chr15': ('15', 'CM000677.2', 'NC_000015.10'),
 'chr16': ('16', 'CM000678.2', 'NC_000016.10'),
 'chr17': ('17', 'CM000679.2', 'NC_000017.11'),
 'chr18': ('18', 'CM000680.2', 'NC_000018.10'),
 'chr19': ('19', 'CM000681.2', 'NC_000019.10'),
 'chr20': ('20', 'CM000682.2', 'NC_000020.11'),
 'chr21': ('21', 'CM000683.2', 'NC_000021.9'),
 'chr22': ('22', 'CM000684.2', 'NC_000022.11'),
 'chrM': ('J01415.2', 'MT', 'NC_012920.1'),
 'chrX': ('CM000685.2', 'NC_000023.11', 'X'),
 'chrY': ('CM000686.2', 'NC_000024.10', 'Y')
}


alias_to_chrom = {}
for name in chrom_aliases:
    for alias in chrom_aliases[name]:
        assert alias not in alias_to_chrom
        alias_to_chrom[alias] = name

