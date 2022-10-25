import subprocess as proc
import multiprocessing as mp
import sys
import os

INPUT_INTERVALS = "~{sep=" " input_intervals}"
OUTPUT_VCF_DIR = "vcf_splits"
OUTPUT_BASENAME = "~{base_vcf}"

COMMAND = """echo gatk \
    --java-options "-Xms2000m -Xmx3000m" \
    SelectVariants \
    -V ~{input_vcf_file} \
    -R ~{ref_fasta} \
    -L {interval_file} \
    -O {output_vcf_file}
"""

def run_job(job_info):
    job_idx = job_info.pop("job_idx")
    job_info['output_vcf_file'] = f"{OUTPUT_BASENAME}.{job_idx:04d}.vcf.gz"
    cmd = COMMAND.format(**job_info)
    cmdobj = proc.run(cmd, shell=True)
    return cmdobj.returncode == 0

def main():
    intervals = INPUT_INTERVALS_DIR.split(' ')
    jobs = [dict(interval_file=ifp, job_idx=idx) for (idx, ifp) in enumerate(intervals)]
    n_jobs = len(jobs)
    n_procs = mp.cpu_count()
    msg = f'[[ Running {n_jobs} jobs across {n_procs} threads ]]'
    print(msg)
    with mp.Pool() as pool:
        return_codes = pool.map(run_job, jobs)
    n_failed = return_codes.count(False)
    if n_failed:
        msg = f'[[ ERROR :: {n_failed} out of {n_jobs} jobs failed ]]'
        print(msg)
        return_code = -1
    else:
        msg = f'[[ SUCCESS :: all {n_jobs} jobs successfully completed ]]'
        print(msg)
        return_code = 0
    return return_code

if __name__ == '__main__':
    os.makedirs(OUTPUT_VCF_DIR, exist_ok=True)
    rc = main()
    sys.exit(rc)
