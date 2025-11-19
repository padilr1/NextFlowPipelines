#!/usr/bin/env bash
#SBATCH --time=0-11:59:00 # time (D-HH:MM)
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G
#SBATCH --account=def-majewski
#SBATCH --mail-user=reinnier.padilla@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END

# time was initially 5-12:59:00
# time was 0-48:59:00
# takes about 8 hours to align and generate mcool files
# about 2 hours to generate cooltools files
# total = 10-12 hours, 10 cores, 70G 
# mem was 88G
module load StdEnv/2020 gcc/9.3.0 python/3.10.2 scipy-stack fastqc/0.11.9 samtools/1.13 kentutils/401 fastp/0.23.4 bwa-mem2/2.2.1
source /lustre06/project/6007495/padilr1/venvs/HiC/bin/activate
pdir=/lustre06/project/6007495/padilr1/pipelines/hic
export PATH=$pdir/bin:"$PATH"

# odir="/lustre06/project/6007495/padilr1/projects/10T/Hi_C/checkout/Feb_2025_batch/downsampled"
# samp="D3KO_reprimed_d8_HiC_785999_S26_L008"
# gn="hg38"

NP=$SLURM_CPUS_PER_TASK
wdir=$odir/$samp
gdir=$pdir/genomes/$gn
mkdir -p $wdir

# #symlink source sequences
lnk () {
  if [[ "$1" =~ 'q'$ ]]; then
    ln -s $(realpath $1) ${2}.fq
    fastqc ${2}.fq
  elif [[ "$1" =~ 'gz'$ ]]; then
    ln -s $(realpath $1) ${2}.fq.gz
    fastqc ${2}.fq.gz
  else
    echo "$1 doesn't end with q or gz, please rename. Exiting..."
    exit 1
  fi
}

# # dump BAM to FASTQ in case needed
bam2fq () {
  bam=$1
  wdir=$2
  nprocs=$3
  samtools sort -n -@ $((nprocs-1)) -o $wdir/reads.bam $bam 
  samtools fastq -1 $wdir/read1.fq -2 $wdir/read2.fq -0 /dev/null -s /dev/null -n $wdir/reads.bam
  rm $wdir/reads.bam
  echo $wdir/read1.fq,$wdir/read2.fq
}

if [[ "$reads" == *.bam ]]; then
  reads=$(bam2fq $reads $wdir $nprocs)
fi
if [[ "$reads" == *,* ]]; then
  lnk ${reads#*,} $wdir/${samp}.R2 &
fi
lnk ${reads%,*} $wdir/${samp} &
wait

# trim
cd $wdir
if [[ -f "${samp}.fq.gz" ]]; then
  mv ${samp}.fq.gz R1.fq.gz
  fastp_cmd="-i R1.fq.gz"
else
  mv ${samp}.fq R1.fq
  fastp_cmd="-i R1.fq"
fi
fastp_cmd="$fastp_cmd -o ${samp}.R1.fq.gz"
if [[ -f "${samp}.R2.fq.gz" ]]; then
  mv ${samp}.R2.fq.gz R2.fq.gz
  fastp_cmd="$fastp_cmd -I R2.fq.gz"
else
  mv ${samp}.R2.fq R2.fq
  fastp_cmd="$fastp_cmd -I R2.fq"
fi
fastp_cmd="$fastp_cmd -O ${samp}.R2.fq.gz"
fastp ${fastp_cmd} \
  --thread $NP \
  --overrepresentation_analysis \
  --detect_adapter_for_pe \
  --json ${samp}.fastp.json \
  --html ${samp}.fastp.html
sed -i "s/-i R1.fq[.gz]\?[.gz]\?[.gz]\?/-i ${samp}/" ${samp}.fastp.json
rm -f R{1,2}.fq.gz R{1,2}.fq

# align
bwa-mem2 mem -t $NP -SP $gdir/bwa2/${gn}.fa ${samp}.R{1,2}.fq.gz \
  | samtools view -bh \
  > ${samp}.bam

# index bam
samtools sort ${samp}.bam > ${samp}_sorted.bam
samtools index ${samp}_sorted.bam
samtools flagstat ${samp}_sorted.bam > ${samp}_sorted.bam.flagstat



echo "Generated aligned files"

# bam -> pairs
echo "Starting to generate pairs"
samtools view -h ${samp}.bam | 
  pairtools parse --nproc-in $NP --nproc-out $NP --drop-sam --drop-seq -c $gdir/${gn}.chrom.sizes --assembly ${gn} --output-stats ${samp}.step1.stats -o ${samp}.pairs.gz
pairtools sort --nproc $NP -o ${samp}.sorted.pairs.gz ${samp}.pairs.gz
pairtools select 'True' --chrom-subset $gdir/${gn}.canon.sizes -o ${samp}.selected.pairs.gz ${samp}.sorted.pairs.gz
pairtools restrict -f $gdir/${gn}.Arima.bed -o ${samp}.restricted.pairs.gz ${samp}.selected.pairs.gz
pairtools dedup -p $NP --backend cython --mark-dups --output-stats ${samp}.dedup.stats -o ${samp}.dedup.pairs.gz ${samp}.restricted.pairs.gz
# pairs -> pairix
pairtools select --chrom-subset $gdir/${gn}.chroms '(rfrag1!=rfrag2) and ((pair_type=="UU" or pair_type=="RU" or pair_type=="UR"))' -o ${samp}.final.pairs.gz ${samp}.dedup.pairs.gz
pairix ${samp}.final.pairs.gz
# stats
pairtools stats ${samp}.final.pairs.gz -o ${samp}.final.stats
pairtools stats ${samp}.pairs.gz -o ${samp}.step1.stats
pairtools stats ${samp}.selected.pairs.gz -o ${samp}.selected.stats
pairtools stats ${samp}.restricted.pairs.gz -o ${samp}.restricted.stats
pairtools stats ${samp}.dedup.pairs.gz -o ${samp}.deduped.stats
echo "Generated pairs"

# run cooler
deactivate
module load StdEnv/2020 gcc/9.3.0 python/3.10.2 scipy-stack fastqc/0.11.9 samtools/1.13 kentutils/401
source /lustre06/project/6007495/padilr1/venvs/Anaconda/bin/activate
conda activate dchic

echo "start generating cool file"
cooler cload pairix -p $NP --assembly ${gn} $gdir/${gn}.canon.sizes:1000 ${samp}.final.pairs.gz 1000.cool

# matrix balancing
echo "start matrix balancing"
cooler balance 1000.cool -p $NP --blacklist $gdir/blacklist.bed 

# Coarsen the matrix
echo "generate mcool file"
cooler zoomify --nproc $NP --out ${samp}.mcool --resolutions 1000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000 --balance  --balance-args "--nproc $NP --max-iters 1000" 1000.cool
rm 1000.cool
echo "finished generating mcool file"

deactivate
module load StdEnv/2020 gcc/9.3.0 python/3.10.2 scipy-stack fastqc/0.11.9 samtools/1.13 kentutils/401 fastp/0.23.4 bwa-mem2/2.2.1
export PATH="/home/padilr1/.local/bin:$PATH"
export PATH="/lustre06/project/6007495/padilr1/venvs/computeCanada_venv/cooltools_venv/bin:$PATH"
source /lustre06/project/6007495/padilr1/venvs/computeCanada_venv/cooltools_venv/bin/activate

# contact decay
echo "calculating a contact probability decay curve for the generated mcool file"
mkdir -p exp
for r in 1000 5000 10000 25000 50000 100000 250000 500000 1000000 2500000; do 
  ~/.local/bin/cooltools expected-cis --smooth --aggregate-smoothed -o exp/${r}.cis.tsv -p $NP ${samp}.mcool::/resolutions/$r
  if [[ $r -ge 10000 ]]; then
    ~/.local/bin/cooltools expected-trans -o exp/${r}.trans.tsv -p $NP ${samp}.mcool::/resolutions/$r
  fi
done

# compartment score
echo "Perform eigendecomposition on the same matrix"
mkdir -p eig
for r in 10000 25000 50000 100000 250000 500000; do
  ~/.local/bin/cooltools eigs-cis --phasing-track $gdir/bin/${r}.gc -o eig/$r ${samp}.mcool::/resolutions/$r --bigwig
done

#boundary score
echo "compute insulation scores"
mkdir -p ins
for p in 10000_100000 25000_1000000; do
  r=${p%_*}
  w=${p#*_}
  ~/.local/bin/cooltools insulation --bigwig ${samp}.mcool::/resolutions/$r $w --output ins/${r} -p $NP
  mv ins/$r ins/${r}.${w}.tsv
done

#saddle plots
echo "compute saddle plot matrices"
mkdir -p sdl
for r in 10000 25000 50000 100000 250000; do
  grep -v chrY exp/${r}.cis.tsv > exp.tmp
  ~/.local/bin/cooltools saddle -o sdl/${r}.cis --view $gdir/${gn}.autoX.bed --qrange 0.02 0.98 ${samp}.mcool::/resolutions/$r eig/${r}.cis.vecs.tsv::E1 exp.tmp::balanced.avg
  if [[ $r -ge 50000 ]]; then
    grep -v chrY exp/${r}.trans.tsv > exp.tmp
    ~/.local/bin/cooltools saddle -o sdl/${r}.trans -t trans --view $gdir/${gn}.autoX.bed --qrange 0.02 0.98 ${samp}.mcool::/resolutions/$r eig/${r}.cis.vecs.tsv::E1 exp.tmp::balanced.avg
  fi
done
rm exp.tmp

# loop score
echo "compute dot matrices for loop pile-up"
mkdir -p dot
for r in 5000 10000 25000; do
  ~/.local/bin/cooltools dots -o dot/$r -p $NP ${samp}.mcool::/resolutions/$r exp/${r}.cis.tsv::balanced.avg 
done

# log bin
deactivate
conda deactivate
module load StdEnv/2023
source /lustre06/project/6007495/padilr1/venvs/Anaconda/bin/activate
conda activate final_cooltools
export PATH="/home/padilr1/.local/bin:$PATH"
export PATH="/lustre06/project/6007495/padilr1/venvs/Anaconda/envs/final_cooltools/bin:$PATH"

for r in 1000 5000 10000 25000 50000 100000 ; do
  /cvmfs/soft.computecanada.ca/gentoo/2023/x86-64-v3/usr/bin/python /lustre06/project/6007495/padilr1/projects/10T/Hi_C/work/logbin/logbin_expected.py --resolution $r --bins-per-order-magnitude 20 exp/${r}.cis.tsv::balanced.avg exp/${r}.cis
done