"""
Metagenomics Pipeline
Автор: Вы
Описание: Полный пайплайн для метагеномного анализа
"""

import os
from glob import glob

# ============================================
# КОНФИГУРАЦИЯ
# ============================================
configfile: "config.yaml"

# Создаём необходимые директории
os.makedirs(config["results_dir"], exist_ok=True)
os.makedirs(os.path.join(config["results_dir"], "trimmed"), exist_ok=True)
os.makedirs(os.path.join(config["results_dir"], "assembly"), exist_ok=True)
os.makedirs(os.path.join(config["results_dir"], "alignment"), exist_ok=True)
os.makedirs(os.path.join(config["results_dir"], "clusters"), exist_ok=True)
os.makedirs(os.path.join(config["results_dir"], "genes"), exist_ok=True)
os.makedirs(os.path.join(config["results_dir"], "viruses"), exist_ok=True)


# ============================================
# ПРАВИЛО ALL - ЧТО ЗАПУСКАЕМ
# ============================================
rule all:
    input:
        viruses = os.path.join(config["results_dir"], "viruses/predictions.txt"),
        genes = os.path.join(config["results_dir"], "genes/proteins.faa"),
        clusters = os.path.join(config["results_dir"], "clusters/clusters.tsv"),
        assembly = os.path.join(config["results_dir"], "assembly/contigs.fa")
    run:
        print("\n" + "="*60)
        print("ПАЙПЛАЙН УСПЕШНО ЗАВЕРШЁН")
        print("="*60)
        print(f"📁 Контиги:      {input.assembly}")
        print(f"📁 Кластеры:     {input.clusters}")
        print(f"📁 Белки:        {input.genes}")
        print(f"📁 Предсказания: {input.viruses}")
        print("="*60 + "\n")


# ============================================
# ШАГ 1: ОБРЕЗКА РИДОВ
# ============================================
rule trim_reads:
    """
    Обрезка ридов по качеству с помощью пользовательской функции.
    """
    input:
        r1 = os.path.join(config["raw_reads_dir"], "{sample}_R1.fastq"),
        r2 = os.path.join(config["raw_reads_dir"], "{sample}_R2.fastq")
    output:
        r1 = os.path.join(config["results_dir"], "trimmed/{sample}_R1.trimmed.fastq"),
        r2 = os.path.join(config["results_dir"], "trimmed/{sample}_R2.trimmed.fastq")
    run:
        from trim_reads import trim_reads
        
        print(f"Обрезка ридов для образца {wildcards.sample}...")
        
        trim_reads(
            input.r1, 
            output.r1, 
            min_quality=config.get("min_quality", 20),
            min_length=config.get("min_read_length", 50)
        )
        
        trim_reads(
            input.r2, 
            output.r2, 
            min_quality=config.get("min_quality", 20),
            min_length=config.get("min_read_length", 50)
        )


# ============================================
# ШАГ 2: СБОРКА КОНТИГОВ
# ============================================
rule assemble:
    """
    Сборка контигов из обрезанных ридов.
    """
    input:
        r1 = expand(os.path.join(config["results_dir"], "trimmed/{sample}_R1.trimmed.fastq"), 
                   sample=config["samples"]),
        r2 = expand(os.path.join(config["results_dir"], "trimmed/{sample}_R2.trimmed.fastq"), 
                   sample=config["samples"])
    output:
        contigs = os.path.join(config["results_dir"], "assembly/contigs.fa")
    params:
        min_len = config.get("min_contig_len", 1000),
        kmer = config.get("kmer_size", 31)
    run:
        from assembly import assemble_contigs
        
        print("Сборка контигов...")
        
        # Для сборки используем все риды (объединяем)
        assemble_contigs(
            forward_reads=input.r1[0],  # берём первый образец
            reverse_reads=input.r2[0],   # для простоты
            output_fasta=output.contigs,
            min_contig_len=params.min_len,
            kmer_size=params.kmer
        )


# ============================================
# ШАГ 3: ВЫРАВНИВАНИЕ
# ============================================
rule align:
    """
    Выравнивание исходных ридов на собранные контиги.
    """
    input:
        contigs = rules.assemble.output.contigs,
        r1 = expand(os.path.join(config["results_dir"], "trimmed/{sample}_R1.trimmed.fastq"), 
                   sample=config["samples"]),
        r2 = expand(os.path.join(config["results_dir"], "trimmed/{sample}_R2.trimmed.fastq"), 
                   sample=config["samples"])
    output:
        bam = os.path.join(config["results_dir"], "alignment/aligned.bam")
    params:
        threads = config.get("threads", 4)
    run:
        from alignment import align_reads_to_contigs
        
        print("Выравнивание ридов на контиги...")
        
        align_reads_to_contigs(
            contigs_fasta=input.contigs,
            forward_reads=input.r1[0],
            reverse_reads=input.r2[0],
            output_bam=output.bam,
            threads=params.threads
        )


# ============================================
# ШАГ 4: КЛАСТЕРИЗАЦИЯ (БИННИНГ)
# ============================================
rule binning:
    """
    Кластеризация контигов в бины с помощью VAMB-подобного подхода.
    """
    input:
        contigs = rules.assemble.output.contigs,
        bam = rules.align.output.bam
    output:
        clusters = os.path.join(config["results_dir"], "clusters/clusters.tsv")
    params:
        kmer = config.get("kmer_size", 4),
        latent = config.get("latent_dim", 32),
        eps = config.get("dbscan_eps", 0.5)
    run:
        from binning import bin_contigs
        
        print("Кластеризация контигов...")
        
        bin_contigs(
            contigs_fasta=input.contigs,
            bam_file=input.bam,
            output_clusters=output.clusters,
            kmer_size=params.kmer,
            latent_dim=params.latent,
            dbscan_eps=params.eps
        )


# ============================================
# ШАГ 5: ПРЕДСКАЗАНИЕ ГЕНОВ
# ============================================
rule predict_genes:
    """
    Предсказание генов в контигах.
    """
    input:
        contigs = rules.assemble.output.contigs
    output:
        proteins = os.path.join(config["results_dir"], "genes/proteins.faa"),
        gff = os.path.join(config["results_dir"], "genes/genes.gff")
    params:
        min_gene = config.get("min_gene_length", 90)
    run:
        from gene_prediction import predict_genes
        
        print("Предсказание генов...")
        
        predict_genes(
            contigs_fasta=input.contigs,
            output_proteins=output.proteins,
            output_gff=output.gff,
            min_gene_length=params.min_gene
        )


# ============================================
# ШАГ 6: ОБНАРУЖЕНИЕ ВИРУСОВ
# ============================================
rule detect_viruses:
    """
    Поиск вирусных контигов.
    """
    input:
        contigs = rules.assemble.output.contigs
    output:
        predictions = os.path.join(config["results_dir"], "viruses/predictions.txt"),
        viral_fasta = os.path.join(config["results_dir"], "viruses/viral_contigs.fa")
    params:
        threshold = config.get("viral_threshold", 0.7)
    run:
        from virus_detection import detect_viruses
        
        print("Обнаружение вирусных последовательностей...")
        
        detect_viruses(
            contigs_fasta=input.contigs,
            output_predictions=output.predictions,
            threshold=params.threshold
        )


# ============================================
# ДОПОЛНИТЕЛЬНЫЕ ПРАВИЛА
# ============================================

rule clean:
    """
    Очистка временных файлов.
    """
    shell:
        """
        rm -rf {config[results_dir]}/tmp
        rm -f {config[results_dir]}/alignment/*.sam
        find {config[results_dir]} -name "*.tmp" -delete
        echo "Очистка завершена"
        """


rule report:
    """
    Генерация отчёта.
    """
    input:
        clusters = rules.binning.output.clusters,
        proteins = rules.predict_genes.output.proteins,
        viruses = rules.detect_viruses.output.predictions
    output:
        report = os.path.join(config["results_dir"], "report.txt")
    run:
        with open(output.report, 'w') as f:
            f.write("=" * 60 + "\n")
            f.write("ОТЧЁТ МЕТАГЕНОМНОГО ПАЙПЛАЙНА\n")
            f.write("=" * 60 + "\n\n")
            
            f.write(f"Образцы: {', '.join(config['samples'])}\n")
            f.write(f"Результаты: {config['results_dir']}\n\n")
            
            # Считаем кластеры
            n_clusters = 0
            if os.path.exists(input.clusters):
                with open(input.clusters) as cf:
                    n_clusters = sum(1 for _ in cf) - 1  # минус заголовок
            f.write(f"Количество кластеров: {n_clusters}\n")
            
            # Считаем гены
            n_genes = 0
            if os.path.exists(input.proteins):
                with open(input.proteins) as pf:
                    n_genes = sum(1 for line in pf if line.startswith('>'))
            f.write(f"Количество предсказанных генов: {n_genes}\n")
            
            # Считаем вирусы
            n_viral = 0
            if os.path.exists(input.viruses):
                with open(input.viruses) as vf:
                    for line in vf:
                        if 'VIRAL' in line:
                            n_viral += 1
            f.write(f"Вирусные контиги: {n_viral}\n")
            
        print(f"Отчёт сохранён в {output.report}")
