// params
params.reads      = "data/*_{R1,R2}.fastq.gz"
params.ref        = "ref/genome.fa"
params.knownSites = "ref/known_sites.vcf.gz"
params.PICARD_JAR = 
params.SNPEFF_JAR = 
params.GATK_JAR   = 
params.VCF        = 
params.outdir     = "results"
params.interval   = null   // e.g. "chr20" or "chr20:1-500000" (optional)


process exportGATK {

}

// align, sort and index reads
process alignReads {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)
    path ref

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")

    script:
    """
    bwa mem $ref ${reads.join(' ')} | samtools sort -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam

    samtools flagstat ${sample_id}.sorted.bam.flagstat
    """
}

// mark duplicates 
process markDups {
    input:

    output:

    script:
    """
    java -Xmx6g -jar $PICARD_JAR AddOrReplaceReadGroups \
      I=$samp.bam \
      O=$samp.rg.bam \
      RGID=$samp \
      RGLB=$samp \
      RGPL=illumina \
      RGSM=$samp \
      RGPU=$samp
      CREATE_INDEX=true

    java -Xmx4G -jar ${PICARD_JAR} MarkDuplicates \
	-REMOVE_DUPLICATES false \
	-CREATE_INDEX true \
	-I $samp.rg.bam \
	-O $samp.rg.dup.bam \
	-M $samp.rg.dup.metrics
    """
}

// Perform base quality score recalibration with GATK

process reCalib {
    input:

    output:

    script:
    """
    gatk --java-options "-Xmx4G" \
	BaseRecalibrator \
       	-R $REF \
	-known-sites $VCF/dbSNP.vcf \
	-known-sites $VCF/Mills.indels.vcf \
	-known-sites $VCF/1KG.snps.vcf \
	-O $samp.rg.dup.recalibration_report.grp \
	-I $samp.rg.dup.bam
    gatk --java-options "-Xmx4G" \
	ApplyBQSR \
	-R $REF \
	--bqsr-recal-file $samp.rg.dup.recalibration_report.grp \
	-O $samp.rg.dup.recal.bam \
	-I $samp.rg.dup.bam
    """
}

// call variants

process callVariants {
    input:

    output:

    script:
    """
    gatk --java-options "-Xmx4G" \
	 HaplotypeCaller \
	 -R ${REF} \
	 -I $samp.rg.dup.recal.bam \
	 -ERC GVCF \
	 -O $samp.g.vcf.gz \
	 -L intervals.bed \
	 -ip 300 \
	 --native-pair-hmm-threads 1
    """
}

// GATK hard filtering
process gatkHardFiltering {
    input:

    output:

    script:
    """
    gatk --java-options "-Xmx4G" \
    VariantFiltration \
    -V HG002.g.vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O HG002.filt.g.vcf 2> >(grep -v undefined)
    """
}

// variant annotation with SnpEFF
process snpEFFannotations {
    input:

    output:

    script:
    """
    zless output.filtered.vcf.gz \
      	| sed '1,100s/ID=AD,Number=./ID=AD,Number=R/' \
	| vt decompose -s - \
	| vt normalize -r $REF - \
	| awk '$5 !~/\*/' > output.filtered.vt.vcf

    java -Xmx4G -jar $SNPEFF_JAR \
	-c $UTILS/snpEff/snpEff.config \
	-i vcf -o vcf -v GRCh37.75 \
	-canon \
	output.filtered.vt.vcf > output.snpeff.vcf

    bgzip output.snpeff.vcf
    tabix -p vcf output.snpeff.vcf.gz
    """
}

// run GEMINI 

process runGemini {
    input:

    output:

    script:
    """
    gemini load --cores 5 --tempdir $PWD -p samples.ped -t snpEff -v output.snpeff.vcf.gz output.snpeff.db

    gemini query --header --sample-filter "phenotype=2" --show-samples --in only -q "SELECT chrom, start, end, ref, alt, qual, gene, transcript, biotype, impact, cadd_scaled, is_coding, clinvar_sig, clinvar_origin, clinvar_gene_phenotype, clinvar_disease_name FROM variants WHERE (impact_severity='HIGH' OR impact_severity='MED') AND aaf_exac_all < 0.001" output.snpeff.db  > gemini_out.txt
    """
}

// run PLINK

process runPLINK {
    input:

    output:

    script:
    """
    bcftools filter -i'HWE=1 && F_MISSING<0.25' -o out/filt_candidate_EOPC_variants.vcf candidate_EOPC_variants.vcf

    plink2 --vcf out/filt_candidate_EOPC_variants.vcf --make-pgen --out out/temp_out

    plink2 --pfile out/temp_out --update-sex out/sex.txt --make-pgen --pheno out/pheno.txt --out out/withPhenoSex_candidate_EOPC_variants --mac 20

    plink2 \
    --pfile out/withPhenoSex_candidate_EOPC_variants \
    --glm \
    --covar out/covar.txt --covar-name sex age smoking_history bmi PC1 PC2 PC3 PC4 PC5 \
    --covar-variance-standardize \
    --out out/EOPC_assoc \
    --ci 0.95
    """
}