# IMPORT Python libraries
from os.path import join

# Import config file & parameters

#CONFIGFILE='configs/config_NB.yaml'
#CONFIGFILE='configs/config.telomeres.hg38.yaml'
#configfile: CONFIGFILE


# Import paths from config file
DATAPATH=config['datapath'].rstrip('/')
SRCPATH=config['srcpath']
#RESULTPATH=config['resultpath']

# Import sample information from config file
SAMPLE_TXT=config['sampletxt']
SAMPLE_FH=open(SAMPLE_TXT,  'r')
TMPSAMPLES=[SAMPLE[:-1] for SAMPLE in SAMPLE_FH]
SAMPLES=['\t'.join(SAMPLE.split('\t')[0:4]) for SAMPLE in TMPSAMPLES]

# Parse Input/Control information.
CONTROL=config['control']
CONTROL_CELLS=[SAMPLE.split('\t')[0] for SAMPLE in SAMPLES if (SAMPLE.split('\t')[2]==CONTROL and '#' not in SAMPLE)]
CONTROL_TREAT=[SAMPLE.split('\t')[1] for SAMPLE in SAMPLES if (SAMPLE.split('\t')[2]==CONTROL and '#' not in SAMPLE)]
CONTROL_CHIPS=[SAMPLE.split('\t')[2] for SAMPLE in SAMPLES if (SAMPLE.split('\t')[2]==CONTROL and '#' not in SAMPLE)]
CONTROL_REPS=[SAMPLE.split('\t')[3] for SAMPLE in SAMPLES if (SAMPLE.split('\t')[2]==CONTROL and '#' not in SAMPLE)]

CONTROL_SAMPLE=['_'.join(SAMPLE.split('\t')[0:3]) for SAMPLE in SAMPLES if (SAMPLE.split('\t')[2]==CONTROL and '#' not in SAMPLE)]

# Parse ChIP information.
IP_CELLS=[SAMPLE.split('\t')[0] for SAMPLE in SAMPLES if (SAMPLE.split('\t')[2]!=CONTROL and '#' not in SAMPLE)]
IP_TREAT=[SAMPLE.split('\t')[1] for SAMPLE in SAMPLES if (SAMPLE.split('\t')[2]!=CONTROL and '#' not in SAMPLE)]
IP_CHIPS=[SAMPLE.split('\t')[2] for SAMPLE in SAMPLES if (SAMPLE.split('\t')[2]!=CONTROL and '#' not in SAMPLE)]
IP_REPS=[SAMPLE.split('\t')[3] for SAMPLE in SAMPLES if (SAMPLE.split('\t')[2]!=CONTROL and '#' not in SAMPLE)]

IP_SAMPLE=['_'.join(SAMPLE.split('\t')[0:3]) for SAMPLE in SAMPLES if (SAMPLE.split('\t')[2]!=CONTROL and '#' not in SAMPLE)]


# merged control and chip lists
CELL  = IP_CELLS + CONTROL_CELLS
TREAT = IP_TREAT + CONTROL_TREAT
CHIP  = IP_CHIPS + CONTROL_CHIPS
SAMPLE = IP_SAMPLE + CONTROL_SAMPLE
REPS =  IP_REPS + CONTROL_REPS

# define some suffixes for use everywhere
FASTQ_SUFFIX        = config['fqsuffix']
BAM_SUFFIX          = '_mkdup_sorted.bam'
BAM_MERGED_SUFFIX   = '_mkdup_sorted_merged.bam'
PEAK_SUFFIX         = config['macs2peaks']
 

# Collect pipeline results.
rule all:
    input:
        ## ==== SPP using R
        expand(join(DATAPATH, '{cell}/{treat}/{chip}/QC/SPP/{cell}_{treat}_{chip}_sppPeakStats.txt'), zip,
                               cell=IP_CELLS , treat=IP_TREAT, chip=IP_CHIPS ),
        ##        
        ## ==== Peak calling        
        expand(expand(join(DATAPATH, '{cell}/{treat}/{chip}/macs2/{{peak}}',
               '{cell}_{treat}_{chip}_peaks.{{peak}}'), zip, cell=IP_CELLS, treat=IP_TREAT, chip=IP_CHIPS),peak=PEAK_SUFFIX),
        ##        
        ## ==== Signal track generation
        expand(expand(join(DATAPATH, '{cell}/{treat}/{chip}/bw/{cell}_{treat}_{chip}_{{scale}}_{{ratio}}.bw'),
               zip, cell=IP_CELLS, treat=IP_TREAT, chip=IP_CHIPS), scale=config['scale'], ratio=config['ratio']),
        expand(join(DATAPATH, '{cell}/{treat}/{chip}/bw/{cell}_{treat}_{chip}.bw'),
               zip, cell=IP_CELLS + CONTROL_CELLS, treat=IP_TREAT + CONTROL_TREAT, chip=IP_CHIPS + CONTROL_CHIPS),       
        ##
        ## ==== QC
        expand(expand(join(DATAPATH, '{cell}/{treat}/{chip}/QC/{cell}_{treat}_{chip}_{{scale}}_{{ratio}}_TSSheatmap.svg'),
               zip, cell=IP_CELLS, treat=IP_TREAT, chip=IP_CHIPS), scale=config['scale'], ratio='subtract'),
        expand(expand(join(DATAPATH, '{cell}/{treat}/{chip}/QC/{cell}_{treat}_{chip}_{{peak}}_macs2_chipPeakStats.txt' ),zip, cell=IP_CELLS , treat=IP_TREAT, chip=IP_CHIPS), 
               peak=PEAK_SUFFIX),      
        expand(join(DATAPATH, '{cell}/{treat}/{chip}/QC/{cell}_{treat}_{chip}_fingerprint.svg'), zip,
               cell=IP_CELLS, treat=IP_TREAT, chip=IP_CHIPS),
        ##
        ## ==== Fastqc
        expand(join(DATAPATH, '{cell}/{treat}/{chip}/QC/fastqc/{cell}_{treat}_{chip}_{reps}_fastqc.zip'), zip,
                cell=CELL, treat=TREAT, chip=CHIP, reps=REPS),
        ##
        ## ==== Alignment
        expand(join(DATAPATH, '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}' + BAM_MERGED_SUFFIX), zip,
               cell=CELL, treat=TREAT, chip=CHIP)
        




#================================================================================#
#                   ChIP-seq Pre-Processing - merged Exp                         #
#================================================================================#
### Perform SICER peak calling
rule sicer_callBroadPeaks:
    input:
        ip=join(DATAPATH, '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}' + BAM_MERGED_SUFFIX),
        control=join(DATAPATH, '{cell}/{treat}/' + CONTROL + '/bam/{cell}_{treat}_' +
                     CONTROL + BAM_MERGED_SUFFIX)
    output:
        ipbed=temp(join(DATAPATH, '{cell}/{treat}/{chip}/sicer/{cell}_{treat}_{chip}_chip.bed')),
        controlbed=temp(join(DATAPATH, '{cell}/{treat}/{chip}/sicer/{cell}_{treat}_{chip}_control.bed')),
        peaks=join(DATAPATH, '{cell}/{treat}/{chip}/sicer/{cell}_{treat}_{chip}_chip' +
                   '-W' + str(config['window']) + '-G' + str(config['gap']) +
                   '-islands-summary' + '-FDR' + str(config['sicerfdr']))
    params:
        #bedtools=config['bedtools'],
        #sicer=config['sicer'],
        outdir     = join(DATAPATH, '{cell}/{treat}/{chip}/sicer/'),
        genome     = config['genome'],
        redthres   = config['redundancythres'],
        window     = config['window'],
        fragsize   = config['fragsize'],
        genomefrac = config['genomefrac'],
        gap        = config['gap'],
        fdr        = config['sicerfdr']
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        # Convert IP BAM --> BED
        bedtools bamtobed -i {input.ip} > {output.ipbed}
        IPBED=$(basename {output.ipbed})
        # Convert CONTROL BAM --> BED
        bedtools bamtobed -i {input.control} > {output.controlbed}
        CONTROLBED=$(basename {output.controlbed})
        # Change into outdir.
        cd {params.outdir}
        # RUN SICER
        # Usage: ./SICER.sh [InputDir] [bed file] [control file]
        #  [OutputDir] [Species] [redundancy threshold] [window size (bp)]
        #  [fragment size] [effective genome fraction] [gap size (bp)] [FDR]
        SICER.sh {params.outdir} $IPBED $CONTROLBED \
                       {params.outdir} {params.genome} {params.redthres} \
                       {params.window} {params.fragsize} {params.genomefrac} \
                       {params.gap} {params.fdr}
        """

### Perform peak calling
rule macs2_callNarrowPeaks:
    input:
        ip=join(DATAPATH, '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}' + BAM_MERGED_SUFFIX),
        control=join(DATAPATH, '{cell}/{treat}/' + CONTROL + '/bam/{cell}_{treat}_' +
                     CONTROL + BAM_MERGED_SUFFIX)
    output:
        peaks=join(DATAPATH, '{cell}/{treat}/{chip}/macs2/narrowPeak/{cell}_{treat}_{chip}_peaks.narrowPeak'),
    params:
        outdir=join(DATAPATH, '{cell}/{treat}/{chip}/macs2/narrowPeak'),
        g=config['macs2species'],
        fdr=config['macs2fdr'],
        name='{cell}_{treat}_{chip}'
    conda:
        "envs/macs2.yaml"
    shell:
        """
        macs2 callpeak -t {input.ip} -c {input.control} -g {params.g} \
                                --outdir {params.outdir} -n {params.name} \
                                -q {params.fdr} --SPMR;
        """

### Perform broad peak calling
rule macs2_callBroadPeaks:
    input:
        ip=join(DATAPATH, '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}' + BAM_MERGED_SUFFIX),
        control=join(DATAPATH, '{cell}/{treat}/' + CONTROL + '/bam/{cell}_{treat}_' +
                     CONTROL + BAM_MERGED_SUFFIX)
    output:
        peaks=join(DATAPATH, '{cell}/{treat}/{chip}/macs2/broadPeak/{cell}_{treat}_{chip}_peaks.broadPeak')
    params:
        outdir=join(DATAPATH, '{cell}/{treat}/{chip}/macs2/broadPeak'),
        g=config['macs2species'],
        fdr=config['macs2fdr'],
        name='{cell}_{treat}_{chip}'
    conda:
        "envs/macs2.yaml"
    shell:
        """
        macs2 callpeak -t {input.ip} -c {input.control} -g {params.g} \
                                --outdir {params.outdir} -n {params.name} \
                                --broad --broad-cutoff {params.fdr};
        """

### Compute Signal Tracks
rule deeptools_bamCompare:
    input:
        ip=join(DATAPATH, '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}' + BAM_MERGED_SUFFIX),
        control=join(DATAPATH, '{cell}/{treat}/' + CONTROL + '/bam/{cell}_{treat}_' +
                     CONTROL + BAM_MERGED_SUFFIX)
    output:
        bw=join(DATAPATH, '{cell}/{treat}/{chip}/bw/{cell}_{treat}_{chip}_{scale}_{ratio}.bw')
    params:
        scale='{scale}',
        ratio='{ratio}',
        norm=config['norm'],
        bs=config['binsize'],
        ns=config['numbersamples'],
        effectiveGenomeSize=config['effectiveGenomeSize'],
        mqthres=config['mqthres']
    conda:
        "envs/deeptools.yaml"
    shell:
        """
        # Geneate BW from BAM
        bamCompare -b1 {input.ip} -b2 {input.control} -o {output.bw} \
                            --scaleFactorsMethod {params.scale} --operation {params.ratio} \
                -e 200 --effectiveGenomeSize  {params.effectiveGenomeSize} \
                --normalizeUsing {params.norm} \
                --minMappingQuality {params.mqthres} \
                --numberOfSamples {params.ns} \
                            -bs {params.bs} --ignoreDuplicates -p 4 -v
        """


rule compute_bw:
    input:
        bam=join(DATAPATH, '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}' + BAM_MERGED_SUFFIX)
    output:
        bw=join(DATAPATH, '{cell}/{treat}/{chip}/bw/{cell}_{treat}_{chip}.bw'),
        bg=temp(join(DATAPATH, '{cell}/{treat}/{chip}/bw/{cell}_{treat}_{chip}.bg'))
    params:
        chrinfo=config['chrinfo']
    conda:
        "envs/macs2.yaml"
    shell:
        """
        samtools view -F 4 {input.bam} |
        awk '{{OFS=\"\t\"}}{{print $3,$4,$4+1}}' |
        genomeCoverageBed -i - -g {params.chrinfo} -bg |bedtools sort -i - > {output.bg}
        # BedGraph --> BigWig
        bedGraphToBigWig {output.bg} {params.chrinfo} {output.bw}
        """
 


#================================================================================#
#                               ChIP-seq QC                                      #
#================================================================================#
### Spatial resolution QC via Heatmap around known TSS
rule plot_heatmapTSS:
    input:
        bw=join(DATAPATH, '{cell}/{treat}/{chip}/bw/{cell}_{treat}_{chip}_{scale}_{ratio}.bw'),
    output:
        txt        = join(DATAPATH, '{cell}/{treat}/{chip}/QC/{cell}_{treat}_{chip}_{scale}_{ratio}_TSSmatrix.txt'),
        svg        = join(DATAPATH, '{cell}/{treat}/{chip}/QC/{cell}_{treat}_{chip}_{scale}_{ratio}_TSSheatmap.svg'),
        svgprofile = join(DATAPATH, '{cell}/{treat}/{chip}/QC/{cell}_{treat}_{chip}_{scale}_{ratio}_TSSprofile.svg')
    params:
        #computematrix=config['computematrix'],
        tss=config['tss'],
        dist=config['tssdist'],
        #heatmapper=config['heatmapper'],
        color=config['heatmapcol'],
        #profiler=config['profiler'],
        threads=4
    conda:
        "envs/deeptools.yaml"
    shell:
        """
        # Compute signal matrix +/- 1.5kb around TSS.
        computeMatrix reference-point -S {input.bw} \
                                               -R {params.tss} -a {params.dist} \
                                               -b {params.dist} -out {output.txt} \
                                               -p {params.threads}

        # Plot computed signal matrix as heatmap.
        plotHeatmap -m {output.txt} -out {output.svg} --colorMap {params.color}

        # Plot computed signal matrix as profile.
        plotProfile -m {output.txt} -out {output.svgprofile}
        """

### Peak statistics - number of reads, FRiP, NSC/RSC
rule compute_peakStats:
    input:
        bam        = join(DATAPATH, '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}' + BAM_MERGED_SUFFIX),
        macs2peaks = join(DATAPATH, '{cell}/{treat}/{chip}/macs2/{peak}'+ '/{cell}_{treat}_{chip}_peaks.{peak}'),
        sicerpeaks = join(DATAPATH, '{cell}/{treat}/{chip}/sicer/{cell}_{treat}_{chip}_chip' +
                        '-W' + str(config['window']) + '-G' + str(config['gap']) +
                        '-islands-summary' + '-FDR' + str(config['sicerfdr'])),
        spptxt=join(DATAPATH, '{cell}/{treat}/{chip}/QC/SPP/{cell}_{treat}_{chip}_sppPeakStats.txt')
    output:
        txtmacs2         = join(DATAPATH, '{cell}/{treat}/{chip}/QC/{cell}_{treat}_{chip}_{peak}_macs2_chipPeakStats.txt')
    params:
        chrsizes = config['chrinfo'],
        chrorder = config['chrorder'],
        sortedmacs2 = temp(join(DATAPATH, '{cell}/{treat}/{chip}/QC/{cell}_{treat}_{chip}_{peak}_sorted_macs2.bed' )),
        sortedsicer = temp(join(DATAPATH, '{cell}/{treat}/{chip}/QC/{cell}_{treat}_{chip}_sorted_sicer.bed' )),
        txtsicer         = join(DATAPATH, '{cell}/{treat}/{chip}/QC/{cell}_{treat}_{chip}_sicer_chipPeakStats.txt' ),
        txt         = join(DATAPATH, '{cell}/{treat}/{chip}/QC/{cell}_{treat}_{chip}_encode_chipPeakStats.txt' ),
        flagstat    = join(DATAPATH, '{cell}/{treat}/{chip}/QC/{cell}_{treat}_{chip}_flagstats.txt')
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        ## ==== Compute total library size of IP
        LIB_SIZE=$(samtools view -c {input.bam})
        ## ==== MACS2 FRiP
        cat {input.macs2peaks} | cut -f1-3 |  bedtools sort -faidx {params.chrorder} -i - > {params.sortedmacs2}
        ##
        ## ==== Compute number of reads in peaks
        N_MACS2_PEAKREADS=$(bedtools coverage -a {params.sortedmacs2} -b {input.bam} -counts -sorted -g {params.chrsizes} |
                           awk 'BEGIN{{sum=0}}{{sum=sum+$NF}}END{{print sum}}')
        # Compute fraction of reads in peaks (FRiP) for MACS2 peaks
        MACS2_FRiP=$(echo $N_MACS2_PEAKREADS |
                     awk -v libsize=$LIB_SIZE '{{print $1/libsize}}')
        ### SICER FRiP
        cat {input.sicerpeaks} | cut -f1-3 | bedtools sort  -faidx {params.chrorder} -i - > {params.sortedsicer}
        # Compute number of reads in peaks
        N_SICER_PEAKREADS=$(bedtools coverage -a {params.sortedsicer} -b {input.bam} -counts -sorted  -g {params.chrsizes}  |
                           awk 'BEGIN{{sum=0}}{{sum=sum+$NF}}END{{print sum}}')
        # Compute fraction of reads in peaks (FRiP) for SICER peaks
        SICER_FRiP=$(echo $N_SICER_PEAKREADS |
                     awk -v libsize=$LIB_SIZE '{{print $1/libsize}}')
        ### Get NSC/RSC
        NSC=$(awk '{{print $(NF-2)}}' {input.spptxt})
        RSC=$(awk '{{print $(NF-1)}}' {input.spptxt})
        ### PCR Bottleneck Coefficient PBC = N1/Nd
        # N1 = number of genomic locations to which EXACTLY one unique mapping read maps
        # Nd = the number of genomic locations to which AT LEAST one unique mapping read maps
        PBC=$(samtools view {input.bam} | cut -f 3-4 | uniq -c |
             awk 'BEGIN{{n1=0; nd=0}}{{nd+=1; if($1 == 1) n1+=1}}END{{print n1/nd}}')
        # Write results to file
        printf 'LibrarySize\t'$LIB_SIZE'\n' >> {output.txtmacs2}
        printf 'LibrarySize\t'$LIB_SIZE'\n' >> {params.txtsicer}
        printf 'MACS2_FRiP\t'$MACS2_FRiP'\n' >> {output.txtmacs2}
        printf 'SICER_FRiP\t'$SICER_FRiP'\n' >> {params.txtsicer}
        printf 'NSC\t'$NSC'\n' > {params.txt}
        printf 'RSC\t'$RSC'\n' >> {params.txt}
        printf 'PBC\t'$PBC'\n' >> {params.txt}
        # Index BAM and flagstat
        samtools flagstat {input.bam} > {params.flagstat}
        """

### Phantompeakqualtools
### Normalized/Relative Strand Cross-correlation coefficient (NSC/RSC)
rule run_phantomPeakQualTools:
    input:
        ip=join(DATAPATH, '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}' + BAM_MERGED_SUFFIX),
        control=join(DATAPATH, '{cell}/{treat}/' + CONTROL + '/bam/{cell}_{treat}_' +
                     CONTROL + BAM_MERGED_SUFFIX)
    output:
        plot=join(DATAPATH, '{cell}/{treat}/{chip}/QC/SPP/{cell}_{treat}_{chip}_crossCorPlot.pdf'),
        txt=join(DATAPATH, '{cell}/{treat}/{chip}/QC/SPP/{cell}_{treat}_{chip}_sppPeakStats.txt')
    params:
        peakqualtools = SRCPATH + "/run_spp.R ",
        outdir=join(DATAPATH, '{cell}/{treat}/{chip}/QC/SPP'),
        threads=4
    conda:
        "envs/RplusPackages.yaml"
    shell:
        """
        # RUN SPP
         Rscript {params.peakqualtools} -c={input.ip} -i={input.control} \
                                                 -odir={params.outdir} \
                                                 -savp={output.plot} \
                                                 -out={output.txt} -p={params.threads}
        """

### Fingerprint plot
rule plot_bamFingerprint:
    input:
        ip      = join(DATAPATH, 
                       '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}' 
                       + BAM_MERGED_SUFFIX),
        control = join(DATAPATH, '{cell}/{treat}/' + CONTROL + 
                       '/bam/{cell}_{treat}_' +
                       CONTROL + BAM_MERGED_SUFFIX)
    output:
        svg     = join(DATAPATH, '{cell}/{treat}/{chip}/QC/{cell}_{treat}_{chip}_fingerprint.svg'),
        qc     = join(DATAPATH, '{cell}/{treat}/{chip}/QC/{cell}_{treat}_{chip}_fingerprint.qc.txt')
    params:
        #fingerprint=config['fingerprint'],
        title   = '{cell}_{treat}_{chip}_ChIP-seq',
        labels  = '{cell}_{treat}_{chip}_merged {cell}_{treat}_' + CONTROL,
        threads = 4
    conda:
        "envs/deeptools.yaml"
    shell:
        """
        plotFingerprint -b {input.ip} {input.control} \
                             --labels {params.labels} \
                             -T {params.title} -plot {output.svg} \
                             -p {params.threads} --outQualityMetrics {output.qc}
        """
#================================================================================#
#                            ChIP-seq ALIGNMENT                                  #
#================================================================================#
### Merge replicates.

def find_replicates(wilcards):
    # function to find technical replicates for merging  
    # using samtools merge
    replicates = []
    for i in range(len(CELL)):
        if CELL[i] == wilcards.cell and TREAT[i] == wilcards.treat and CHIP[i] == wilcards.chip:
            # add path based on wildcard informaqtion
            # this bam file is needed as input
            fp = DATAPATH + '/' + wilcards.cell + '/' + wilcards.treat + '/' + wilcards.chip + '/bam/' + SAMPLE[i] + '_'  + REPS[i] + BAM_SUFFIX
            replicates.append(fp)
    return replicates
    
rule merge_bam:
    input:
        bams = find_replicates
    output:
        tmpbam=temp(join(DATAPATH, '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}_rmdup_merged.bam')),
        mergebam=join(DATAPATH, '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}' + BAM_MERGED_SUFFIX)
    params:
        prefix=join(DATAPATH, '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}' + BAM_MERGED_SUFFIX.split('.')[0]),
        bamdir=join(DATAPATH, '{cell}/{treat}/{chip}/bam/')
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        pattern=" "
        # Get all bams.
        #BAMS=$(find -L {params.bamdir} -type f | grep bam$)
        #if [[ $(echo $BAMS | wc -w) != 1 ]]; then
        if [[ '{input.bams}' =~ $pattern ]]; then # check if input bams contains at least a single space, bc then its more than one file
            echo "multiple files"
            samtools merge {output.tmpbam} {input.bams}   # syntax is: samtools megrge outfile.bam infile.bam infile2.bam
            echo "merged"
            samtools sort  -@ 10  -o {output.mergebam} {output.tmpbam}
            echo "sorted"
            samtools index {output.mergebam} 
            echo "indexed"
        else
            touch {output.tmpbam}
            rm -f {output.mergebam}
            rm -f {output.mergebam}".bai"
            ln -srf {input.bams} {output.mergebam}
            ln -srf {input.bams}".bai" {output.mergebam}".bai"
        fi
        """

#Single end read alignment, duplicate rmv & sorting.
rule bowtie2_alignReads:
    input:
        fastq=join(DATAPATH, '{cell}/{treat}/{chip}/fastq_trimmed/{cell}_{treat}_{chip}_{reps}' + FASTQ_SUFFIX.split('.')[0] + '_trimmed.fq.gz')
    output:
        bam        = temp(join(DATAPATH, '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}_{reps}.bam')),
        tmpbam     = temp(join(DATAPATH, '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}_{reps}.tmp.bam')),
        tmpbamMate = temp(join(DATAPATH, '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}_{reps}.tmp.mate.bam')),
        tmpNameBam = temp(join(DATAPATH, '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}_{reps}.name.sort.tmp.bam')),
        mkdupbam  = join(DATAPATH, '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}_{reps}' + BAM_SUFFIX),
        log        = join(DATAPATH, '{cell}/{treat}/{chip}/QC/{cell}_{treat}_{chip}_{reps}_bowtie2.log'),
        flagstat   = join(DATAPATH, '{cell}/{treat}/{chip}/QC/{cell}_{treat}_{chip}_{reps}_flagstats.txt')
    params:
        index=config['index'],
        algnpara=config['alignparam'], 
        prefix=join(DATAPATH, '{cell}/{treat}/{chip}/bam/{cell}_{treat}_{chip}_{reps}' + BAM_SUFFIX.split('.')[0]),
        bamdir=join(DATAPATH, '{cell}/{treat}/{chip}/bam/')
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        # Align Reads with BOWTIE2 and convert to BAM
        bowtie2 -p 10 -x {params.index} -U {input.fastq} {params.algnpara} 2> {output.log} | samtools view -h -bS - > {output.tmpbam}
        
        # sort bam with name flag for fixmate
        samtools sort -n -@ 10 -o {output.tmpNameBam} {output.tmpbam}
        samtools fixmate -m {output.tmpNameBam} {output.tmpbamMate}
          
        # Sort BAM
        samtools sort -@ 10 -o {output.bam} {output.tmpbamMate}
        # Remove Duplicates from BAM
        samtools markdup {output.bam} {output.mkdupbam}
        
        # Index BAM and flagstat
        samtools index {output.mkdupbam};
        echo === {output.mkdupbam} > {output.flagstat}
        samtools flagstat  {output.mkdupbam} >> {output.flagstat}
        """

#Perform Adapter trimming with trim galore.
#SUFFIX extension _R1_val_1.fq.gz
rule run_trimGalore:
   input:
       fastq=join(DATAPATH, '{cell}/{treat}/{chip}/fastq/{cell}_{treat}_{chip}_{reps}' + FASTQ_SUFFIX),
   output:
       fastq=join(DATAPATH, '{cell}/{treat}/{chip}/fastq_trimmed/{cell}_{treat}_{chip}_{reps}' + FASTQ_SUFFIX.split('.')[0] + '_trimmed.fq.gz'),
   params:
       outdir=join(DATAPATH, '{cell}/{treat}/{chip}/fastq_trimmed'),
       trimsettings='--fastqc --gzip'
   conda:
        "envs/trimGalore.yaml"
   shell:
       """       echo $PATH;
       if [ ! -d {params.outdir} ]; then
           mkdir {params.outdir}
       fi

       # RUN trimgalore
       trim_galore  {input.fastq} -o {params.outdir} \
       {params.trimsettings}
       """

### Run FastQC
rule run_fastqc:
    input:
        fastq=join(DATAPATH, '{cell}/{treat}/{chip}/fastq/{cell}_{treat}_{chip}_{reps}' + FASTQ_SUFFIX),
    output:
        zip=join(DATAPATH, '{cell}/{treat}/{chip}/QC/fastqc/{cell}_{treat}_{chip}_{reps}_fastqc.zip')
    params:
        #fastqc=config['fastqc'],
        outdir=join(DATAPATH, '{cell}/{treat}/{chip}/QC/fastqc')
    conda:
        "envs/fastqc.yaml"
    shell:
        """
        if [ ! -d {params.outdir} ]; then
            mkdir {params.outdir}
        fi

        fastqc -o {params.outdir} --extract {input.fastq}
        """
