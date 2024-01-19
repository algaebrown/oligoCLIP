locals().update(config)

sample_labels = config['rbps']
libnames = config['libnames']
rbps = config['rbps']

HUMAN_RNA_NUCLEOTIDE='/tscc/projects/ps-yeolab4/seqdata/20200622_gencode_coords_hsher/GRCh38_latest_rna.fna'
N_READ_TO_SAMPLE=5*10**3

rule gather_trimming_stat:
    input:
        tr1=expand("QC/{libname}.Tr.metrics", libname = libnames)
    output:
        tr1="QC/cutadapt_stat.csv",
    params:
        run_time = "00:10:00",
        cores="1",
        error_out_file = "error_files/qctrim.txt",
        out_file = "stdout/qctrim",
        memory = 5000,
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/trimming_stat.py "{input.tr1:q}" {output.tr1}
        """

rule gather_fastqc_report:
    input:
        expand("fastqc/{sample_label}.umi.fqTrTr_fastqc/fastqc_data.txt", sample_label = sample_labels)
    output:
        basic='QC/fastQC_basic_summary.csv',
        passfail='QC/fastQC_passfail.csv'
    params:
        run_time = "00:40:00",
        cores = "1",
        memory = 10000,
        error_out_file = "error_files/fastqc_stat.txt",
        out_file = "stdout/fastqc",
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/fastqc_io.py -i "{input}" -p {output.passfail} -b {output.basic}
        """

rule gather_mapstat:
    input:
        expand("{libname}/bams/{sample_label}.Log.final.out", libname = libnames, sample_label = rbps),
    output:
        "QC/mapping_stats.csv"
    conda:
        "envs/metadensity.yaml"
    params:
        error_out_file = "error_files/mapstat",
        run_time = "00:40:00",
        cores = "1",
        memory = 10000,
        out_file = "stdout/mapstat",
    shell:
        """
        python {SCRIPT_PATH}/star_mapping_stat_io.py -i "{input}" -o {output}
        """

rule duplication_rate:
    input:
        dup=expand("{libname}/bams/{sample_label}.Aligned.sortedByCoord.out.bam", libname = libnames, sample_label = rbps),
        rmdup=expand("{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam",libname = libnames, sample_label = rbps),
        rmdup_bai=expand("{libname}/bams/{sample_label}.rmDup.Aligned.sortedByCoord.out.bam.bai",libname = libnames, sample_label = rbps)
    output:
        'QC/dup_level.csv'
    params:
        error_out_file = "error_files/dup",
        out_file = "stdout/dup_rate",
        run_time = "00:40:00",
        cores = "1",
        memory = 10000,
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/dup_level.py "{input.dup}" "{input.rmdup}" {output}
        """

#echo haha $(echo $(zcat multiplex_HEK293_3/fastqs/ultraplex_demux_QKI_Rev.fastq.gz | wc -l)/4 | bc)
rule count_demultiplex_ultraplex:
    input:
        fq1=expand("{libname}/fastqs/ultraplex_demux_{sample_label}_Rev.fastq.gz", libname = libnames, sample_label = rbps)
    output:
        'QC/demux_read_count.txt'
    params:
        error_out_file = "error_files/demux_count",
        out_file = "stdout/readcount",
        run_time = "00:40:00",
        cores = "1",
        memory = 10000,
    container: None
    shell:
        """
        touch {output}
        for f in {input.fq1} ; do \
        echo "$f $(echo $(zcat $f | wc -l)/4 | bc)" >> {output}; \
        done
        """


rule make_read_count_summary:
    input:
        feature_annotations = config['FEATURE_ANNOTATIONS'],
        counts = "counts/genome/megatables/{libname}.tsv.gz",
    output:
        region_summary =  "QC/read_count/{libname}.region.csv",
        type_summary =  "QC/read_count/{libname}.genetype.csv",
        name_summary =  "QC/read_count/{libname}.genename.csv",
        dist = "QC/read_count/{libname}.cosine_similarity.csv"
    params:
        error_out_file = "error_files/{libname}.read_count_summary.err",
        out_file = "stdout/{libname}.read_count_summary.out",
        run_time = "1:20:00",
        cores = 1,
        memory = 40000,
    run:
        import os
        print(output.region_summary)
        try:
            os.mkdir('QC/read_count')
        except Exception as e:
            print(e)
        import pandas as pd
        cnt = pd.read_csv(input.counts, sep = '\t')
        feature_annotations = pd.read_csv(input.feature_annotations, sep = '\t')

        df = pd.concat([feature_annotations, cnt], axis = 1)

        by_type = df.groupby(by = 'feature_type_top')[cnt.columns].sum()
        by_gene = df.groupby(by = 'gene_type_top')[cnt.columns].sum()
        by_name = df.groupby(by = 'gene_name')[cnt.columns].sum()

        by_type.to_csv(output.region_summary)
        by_gene.to_csv(output.type_summary)
        by_name.to_csv(output.name_summary)

        # distance
        from scipy.spatial.distance import pdist, squareform
        cov_filter = 10
        dist = squareform(pdist(cnt.loc[cnt.sum(axis = 1)>cov_filter].T, 'cosine'))

        dist_df = pd.DataFrame(1-dist, columns = cnt.columns, index = cnt.columns)
        dist_df.to_csv(output.dist)

##### debugging ######
rule what_is_read_wo_barcode:
    input:
        target=HUMAN_RNA_NUCLEOTIDE,
        target_db=HUMAN_RNA_NUCLEOTIDE + '.nog',
        query_fq_gz="{libname}/fastqs/ultraplex_demux_5bc_no_match_Fwd.fastq.gz"
    output:
        blast_result="QC/nobarcode_blast_output/{libname}.blast.tsv",
        fasta="QC/nobarcode_blast_output/{libname}.fasta"
    params:
        error_out_file = "error_files/demux_count",
        out_file = "stdout/readwobar",
        run_time = "00:40:00",
        cores = "1",
        memory = 20000,
        nlines = N_READ_TO_SAMPLE * 4
    conda:
        "envs/blast.yaml"
    shell:
        """
        set +o pipefail; 
        zcat {input.query_fq_gz} | head -n {params.nlines} | sed -n '1~4s/^@/>/p;2~4p' > {output.fasta}
        blastn -db {input.target} -query {output.fasta} -out {output.blast_result} -outfmt 6 -max_target_seqs 1 
        """

rule blast_unmapped_reads:
    input:
        target=HUMAN_RNA_NUCLEOTIDE,
        target_db=HUMAN_RNA_NUCLEOTIDE + '.nog',
        unmapped1_fq= "{libname}/bams/{sample_label}.Unmapped.out.mate1",
        unmapped2_fq= "{libname}/bams/{sample_label}.Unmapped.out.mate2"
    output:
        blast_result1="QC/unmapped_blast_output/{libname}.{sample_label}.1.blast.tsv",
        blast_result2="QC/unmapped_blast_output/{libname}.{sample_label}.2.blast.tsv",
        fasta1="QC/unmapped_blast_output/{libname}.{sample_label}.1.fasta",
        fasta2="QC/unmapped_blast_output/{libname}.{sample_label}.2.fasta"
    params:
        error_out_file = "error_files/demux_count",
        out_file = "stdout/blastunmap",
        run_time = "00:40:00",
        cores = "1",
        memory = 20000,
        nlines = N_READ_TO_SAMPLE * 4
    conda:
        "envs/blast.yaml"
    shell:
        """
        set +o pipefail; 
        samtools fasta {input.unmapped1_fq} | head -n {params.nlines} > {output.fasta1}
        samtools fasta {input.unmapped2_fq} | head -n {params.nlines} > {output.fasta2}
        blastn -db {input.target} -query {output.fasta1} -out {output.blast_result1} -outfmt 6 -max_target_seqs 1 
        blastn -db {input.target} -query {output.fasta2} -out {output.blast_result2} -outfmt 6 -max_target_seqs 1 
        """

rule blast_unmapped_reads_too_short:
    input:
        target=HUMAN_RNA_NUCLEOTIDE,
        target_db=HUMAN_RNA_NUCLEOTIDE + '.nog',
        bam= "{libname}/bams/{sample_label}.Aligned.sortedByCoord.out.bam",
    output:
        blast_result="QC/unmapped_blast_output/{libname}.{sample_label}.short.blast.tsv",
        fasta="QC/unmapped_blast_output/{libname}.{sample_label}.short.fasta",
    params:
        error_out_file = "error_files/demux_count",
        out_file = "stdout/blastunmap",
        run_time = "00:40:00",
        cores = "1",
        nlines = N_READ_TO_SAMPLE,
        memory = 20000,
    conda:
        "envs/blast.yaml"
    shell:
        """
        set +o pipefail; 
        samtools view -f 4 {input.bam} | grep uT:A:1 | head -n {params.nlines} | samtools fasta >  {output.fasta}
        blastn -db {input.target} -query {output.fasta} -out {output.blast_result} -outfmt 6 -max_target_seqs 1 
        """

rule summary_QC_statistics:
    ''' Gigantic summary table on where we lose our reads '''
    input:
        fastqc_initial = 'QC/fastQC_initial_basic_summary.csv',
        fastqc_post_processing = 'QC/fastQC_basic_summary.csv',
        cutadapt_stat = 'QC/cutadapt_stat.csv',
        mapping_stat = 'QC/mapping_stats.csv',
        dup_stat = 'QC/dup_level.csv',
        repeat_cnts = expand("counts/repeats/megatables/class/{libname}.tsv.gz", libname = libnames),
        genome_cnts = expand('QC/read_count/{libname}.region.csv', libname = libnames)
    output:
        "QC/summary.csv"
    params:
        error_out_file = "error_files/QC_summary",
        out_file = "stdout/QC_summary",
        run_time = "00:20:00",
        cores = 1,
        memory = 20000,
    run:
        import pandas as pd
        import re
        from pathlib import Path
        import json

        #### functions for ultraplex statistics extraction ####
        def extract_nread_with_barcode(s):
            ''' extract nread with barcode from ultraplex.log'''
            pattern = r"(\d+)\s+\((\d+\.\d+)%\)\s+reads\s+correctly\s+assigned"
            match = re.search(pattern, s)
            if match:
                num_reads = match.group(1)
                percentage = match.group(2)
                return num_reads, percentage
            else:
                print(s)
                return None, None
        def extract_input_reads(s):
            '''Demultiplexing complete! 21532666 reads processed in 1095.0 seconds'''
            pattern = r'(?<=Demultiplexing complete! )\d+(?= reads processed)'
            
            matches = re.findall(pattern, s)

            if matches:
                num_reads = int(matches[0])
                
                return num_reads
            else:
                print(s)
                return None

        def extract_quality_trimmed(s):
            pattern = r"(\d+)\s+\((\d+\.\d+)%\)\s+reads\s+quality\s+trimmed"
            match = re.search(pattern, s)
            if match:
                num_reads = match.group(1)
                percentage = match.group(2)
                return num_reads, percentage
            else:
                print(s)
                return None, None

        def get_ultraplex_stat(basedir):
            ultraplex_stat = []
            for file in basedir.glob('*/fastqs/ultraplex*.log'):
                libname = file.parent.parent.name
                with open(file) as f:
                    lines = f.readlines()
                    barcode_nread, barcode_perc = extract_nread_with_barcode(lines[-1])
                    ultraplex_input_nread = extract_input_reads([l for l in lines if 'Demultiplex' in l][0])
                ultraplex_stat.append([libname, barcode_nread, barcode_perc, ultraplex_input_nread])
            ultraplex_stat = pd.DataFrame(ultraplex_stat, columns = ['libname', 'nread_w_barcode', '%barcode', 'input_to_ultraplex']
                                ).drop_duplicates('libname').set_index('libname')
            return ultraplex_stat
        
        def get_fastp_stat(basedir):
            ''' extracts fastp statistics '''
            fastp_stats = []
            for fastp_stat_file in basedir.glob('QC/*umi.json'):
                fastp_stat = json.load(open(fastp_stat_file))
                fastp_stats.append([fastp_stat_file.name.split('.')[0],
                                fastp_stat['summary']['before_filtering']['total_reads'],
                                fastp_stat['summary']['after_filtering']['total_reads']]
                                )
            fastp_stats = pd.DataFrame(fastp_stats, columns = ['libname', 'before_umi_polyG', 'after_umi_polyG']).set_index('libname')
            fastp_stats['%polyG or too short']=100*(fastp_stats['before_umi_polyG']-fastp_stats['after_umi_polyG'])/fastp_stats['before_umi_polyG']
            return fastp_stats

        # read inputs:
        cutadapt = pd.read_csv('QC/cutadapt_stat.csv', index_col = 0)
        genome_stat = pd.read_csv('QC/mapping_stats.csv', index_col = 0)
        dup_df = pd.read_csv('QC/dup_level.csv', index_col = 0)
        ultraplex_stat = get_ultraplex_stat(Path('.'))
        fastp_stats = get_fastp_stat(Path('.'))

        # generate mapping summary
        genome_stat['libname']=genome_stat['STAR Log filename'].str.split('/', expand = True)[0]
        mapping_summary = pd.concat([genome_stat.groupby(by = 'libname')['Uniquely mapped reads number'].sum(),
                genome_stat.groupby(by = 'libname')['Number of reads mapped to multiple loci'].sum(),
                                    genome_stat.groupby(by = 'libname')['Number of input reads'].sum(),
                ], axis = 1
                )
        mapping_summary['total_mapped']=mapping_summary[['Uniquely mapped reads number', 'Number of reads mapped to multiple loci']].sum(axis = 1)
        mapping_summary['%Unique map']=100*mapping_summary['Uniquely mapped reads number']/mapping_summary['Number of input reads']
        mapping_summary['%Multi map']=100*mapping_summary['Number of reads mapped to multiple loci']/mapping_summary['Number of input reads']
        mapping_summary['% mapped']=100*mapping_summary['total_mapped']/mapping_summary['Number of input reads']

        # generate dedup summary: is the sum of read 1 and read2
        dup_df['libname']=dup_df['dup_bam'].str.split('/', expand = True)[0]
        dup_summary_df=dup_df.groupby(by = ['libname'])['before_dedup', 'after_dedup'].sum()
        dup_summary_df['% Unique frag']=100*dup_summary_df['after_dedup']/dup_summary_df['before_dedup']

        # generate repeat counts summary
        repeat_cnt_stat = []
        for repeat_cnt_file in list(Path('counts/repeats/megatables/class').glob('*tsv.gz')):
            repeat_cnt = pd.read_csv(repeat_cnt_file, sep = '\t', index_col = 0)
            n_read_in_repeat = repeat_cnt.sum().sum()
            n_read_in_rRNA = repeat_cnt.loc['rRNA'].sum()
            
            repeat_cnt_stat.append([repeat_cnt_file.name.split('.')[0],
                                    n_read_in_repeat,
                                    n_read_in_rRNA]
                                )
        repeat_cnt_stat = pd.DataFrame(repeat_cnt_stat, columns = ['libname', 'n_read_in_repeat', 'n_read_in_rRNA']).set_index('libname')

        # generate genome counts summary
        reads_in_window = []
        for region_cnt in list(Path('QC/read_count').glob('*region.csv')):
            region_read_count = pd.read_csv(region_cnt, index_col = 0)
            reads_in_window.append([region_cnt.name.split('.')[0], region_read_count.sum().sum()])
        reads_in_window = pd.DataFrame(reads_in_window, columns = ['libname', 'nread_in_genome_window']).set_index('libname')

        # make summary
        cutadapt.index = [i.split('.')[0] for i in cutadapt.index]
        try:
            summary = pd.concat([cutadapt[['Total reads processed', '% Reads that were too short']],
                    fastp_stats,
                    ultraplex_stat,

                    mapping_summary,
                    dup_summary_df,

                    reads_in_window,
                                repeat_cnt_stat],
                    axis = 1)
            is_paired = False
        except:
            summary = pd.concat([cutadapt[['Total read pairs processed', '% Pairs that were too short']],
                    fastp_stats,
                    ultraplex_stat,

                    mapping_summary,
                    dup_summary_df,

                    reads_in_window,
                                repeat_cnt_stat],
                    axis = 1)
            is_paired = True

        if not is_paired:
            summary['%_in_genome_window']=100*summary['nread_in_genome_window']/summary['after_dedup']
            summary['%_in_repeat']=100*summary['n_read_in_repeat']/summary['after_dedup']
            summary['%_in_rRNA']=100*summary['n_read_in_rRNA']/summary['after_dedup']
        else:
            # counting is only counting 1 read, but samtools counts both
            summary['%_in_genome_window']=100*2*summary['nread_in_genome_window']/summary['after_dedup']
            summary['%_in_repeat']=100*2*summary['n_read_in_repeat']/summary['after_dedup']
            summary['%_in_rRNA']=100*2*summary['n_read_in_rRNA']/summary['after_dedup']

        summary.T.to_csv(output[0])