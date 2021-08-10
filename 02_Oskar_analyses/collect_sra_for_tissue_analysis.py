import os, subprocess 

import pandas as pd 

from ftplib import FTP
from tqdm import tqdm 

def collect_urls(tsa_id):
    url_list = []
    prefix = tsa_id[:2]
    suffix = tsa_id[2:4]
    print('Getting data URLs...')
    ftp = FTP('ftp.ncbi.nih.gov')
    ftp.login()

    url_base = f'/sra/wgs_aux/{prefix}/{suffix}/{tsa_id}/'
    ftp.cwd(url_base)
    for i in range(2):
        i += 1 
        url = f'{url_base}{tsa_id}.{i}.fsa_nt.gz'
        url_list.append(url)
    return url_list

        
def get_tsa_data(tsa_id, folder_output):
    tsa = []
    urls = collect_urls(tsa_id)
    for url in urls:
        file_name = url.split('/')[-1]
        if not os.path.join(folder_output, file_name):
            print(f'Downloading rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov{url} {folder_output}')
            os.system(f'rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov{url} {folder_output}')
        else:
            print(f'{file_name} already downloaded')
            tsa.append(file_name)
    return tsa 


def build_tsa_input():
    folder_base = '/home/savvy/PROJECTS/Oskar_Evolution'
    folder_output = f'{folder_base}/Data/02_Oskar_analyses/2.12/aedes_aegypti/'
    tsa = get_tsa_data('GFNA01', folder_output)

    os.chdir('..')
    os.chdir(folder_output)
    transcript_output = 'GFNA01_transcripts.fasta.gz'
    if os.path.isfile(f'{folder_output}{transcript_output}'):
        print(f'{transcript_output} already generated')
    else: 
        print('Generating TSA master file')
        cmd = f'zcat {tsa[0]} {tsa[1]} > {transcript_output.replace(".gz", "")}'
        subprocess.run(cmd, shell=True)
        subprocess.run(f'gzip {transcript_output.replace(".gz", "")}', shell=True)

def prefetch_reads():
    folder_base = '/home/savvy/PROJECTS/Oskar_Evolution'
    folder_output = f'{folder_base}/Data/02_Oskar_analyses/2.12/aedes_aegypti/sra_files'
    sra_table = pd.read_table(f'{folder_output}/sra_info_table.txt', sep=',')
    prefetch = '/home/savvy/bin/sratoolkit.2.11.0-ubuntu64/bin/prefetch'
    for run in sra_table['Run'].values:
        cmd = f'{prefetch} --max-size 1000000000000 -p {run}'
        subprocess.run(cmd, shell=True)


def fasterq_dump_reads():
    folder_base = '/home/savvy/PROJECTS/Oskar_Evolution'
    folder_output = f'{folder_base}/Data/02_Oskar_analyses/2.12/aedes_aegypti/sra_files'
    sra_table = pd.read_table(f'{folder_output}/sra_info_table.txt', sep=',')
    folder_output = f'/media/savvy/DATA1/savvy/Oskar_Evolution/fastq/'
    fasterq_dump = '/home/savvy/bin/sratoolkit.2.11.0-ubuntu64/bin/fasterq-dump'
    print('Downloading fastq files')
    for run in tqdm(sra_table['Run'].values):
        if not os.path.isfile(os.path.join(folder_output, f'{run}.fastq')):
            cmd = f'{fasterq_dump} -e 12 -O {folder_output} -p {run}'
            subprocess.run(cmd, shell=True)


def run_kallisto_index():
    folder_base = '/home/savvy/PROJECTS/Oskar_Evolution'
    folder_output = f'{folder_base}/Data/02_Oskar_analyses/2.12/aedes_aegypti/'
    transcripts_gz = 'GFNA01_transcripts.fasta.gz'
    transcripts_idx = 'GFNA01_transcripts.idx'
    print(f'\nBuilding index for {transcripts_gz.replace(".gz", "")}')
    cmd = f'kallisto index -i {transcripts_idx} {transcripts_gz}'
    print(f'Run ${cmd} in {folder_output}')
    # subprocess.run(cmd, shell=True)
        

if __name__ == '__main__':
    # build_tsa_input()
    # run_kallisto_index() #### Take ~ 6 min
    fasterq_dump_reads()