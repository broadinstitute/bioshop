from tqdm import tqdm

def vcf_progress_bar(vcf=None):
    seqlen_map = dict(zip(vcf.seqnames, vcf.seqlens))
    ns = dict(
        chrom=None,
        pbar=None,
        last_pos=0
    )
    def update(chrom=None, pos=None):
        pbar = ns['pbar']
        if pbar and chrom is None:
            pbar.close()
            return
        last_chrom = ns['chrom']
        if last_chrom != chrom:
            seqlen = seqlen_map[chrom]
            if pbar:
                pbar.close()
            pbar = tqdm(
                total=seqlen,
                desc=chrom,
                unit="base",
                unit_scale=True,
                colour='green'
            )
            ns['chrom'] = chrom
            ns['pbar'] = pbar
            ns['last_pos'] = 0
        last_pos = ns['last_pos']
        pbar.update(pos - last_pos)
        ns['last_pos'] = pos
    return update


