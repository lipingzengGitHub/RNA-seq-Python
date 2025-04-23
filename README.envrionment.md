#Build conda envrionment

conda env create -f environment.yml

#Activate envrionment

conda activate rnaseq-pipeline

#Run python script
python RNA-seq-ReadsCount.py -i ./fastq -o ./output -g genome_index -a annotation.gtf --fastqc --align --count


