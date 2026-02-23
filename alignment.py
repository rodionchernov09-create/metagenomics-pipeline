#!/usr/bin/env python3

import os
import subprocess
import logging
from typing import List, Tuple

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def build_index(contigs_fasta: str, index_prefix: str) -> None:
    logger.info(f"Построение индекса для {contigs_fasta}...")
    
    cmd = ["bwa", "index", "-p", index_prefix, contigs_fasta]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info("Индекс успешно построен")
    except subprocess.CalledProcessError as e:
        logger.error(f"Ошибка при построении индекса: {e.stderr}")
        raise
    except FileNotFoundError:
        logger.error("bwa не найден. Установите bwa или проверьте окружение")
        raise


def align_reads(index_prefix: str, forward_reads: str, reverse_reads: str, output_sam: str, threads: int = 4) -> None:
    logger.info(f"Выравнивание ридов: {forward_reads}, {reverse_reads}...")
    
    cmd = [
        "bwa", "mem",
        "-t", str(threads),
        index_prefix,
        forward_reads,
        reverse_reads
    ]
    
    try:
        with open(output_sam, 'w') as sam_file:
            result = subprocess.run(cmd, stdout=sam_file, stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                logger.error(f"Ошибка bwa mem: {result.stderr}")
                raise subprocess.CalledProcessError(result.returncode, cmd)
        logger.info(f"Выравнивание завершено, SAM файл: {output_sam}")
    except FileNotFoundError:
        logger.error("bwa не найден. Установите bwa или проверьте окружение")
        raise


def sam_to_bam(sam_file: str, bam_file: str, threads: int = 4) -> None:
    logger.info(f"Конвертация SAM в BAM: {sam_file} -> {bam_file}")
    
    cmd_view = ["samtools", "view", "-@", str(threads), "-bS", sam_file]
    cmd_sort = ["samtools", "sort", "-@", str(threads), "-o", bam_file]
    
    try:
        with subprocess.Popen(cmd_view, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as view_proc:
            with subprocess.Popen(cmd_sort, stdin=view_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as sort_proc:
                _, stderr = sort_proc.communicate()
                if sort_proc.returncode != 0:
                    logger.error(f"Ошибка samtools sort: {stderr.decode()}")
                    raise subprocess.CalledProcessError(sort_proc.returncode, cmd_sort)
        logger.info(f"BAM файл создан: {bam_file}")
    except FileNotFoundError:
        logger.error("samtools не найден. Установите samtools или проверьте окружение")
        raise


def index_bam(bam_file: str) -> None:
    logger.info(f"Индексация BAM файла: {bam_file}")
    
    cmd = ["samtools", "index", bam_file]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info("BAM индекс создан")
    except subprocess.CalledProcessError as e:
        logger.error(f"Ошибка индексации BAM: {e.stderr}")
        raise
    except FileNotFoundError:
        logger.error("samtools не найден. Установите samtools или проверьте окружение")
        raise


def align_reads_to_contigs(
    contigs_fasta: str,
    forward_reads: str,
    reverse_reads: str,
    output_bam: str,
    threads: int = 4,
    keep_sam: bool = False
) -> None:
    logger.info("=" * 60)
    logger.info("НАЧАЛО ВЫРАВНИВАНИЯ РИДОВ НА КОНТИГИ")
    logger.info("=" * 60)
    
    work_dir = os.path.dirname(output_bam)
    if work_dir:
        os.makedirs(work_dir, exist_ok=True)
    
    index_prefix = contigs_fasta + ".bwa"
    sam_file = output_bam.replace('.bam', '.sam')
    
    build_index(contigs_fasta, index_prefix)
    
    align_reads(index_prefix, forward_reads, reverse_reads, sam_file, threads)
    
    sam_to_bam(sam_file, output_bam, threads)
    
    index_bam(output_bam)
    
    if not keep_sam:
        os.remove(sam_file)
        logger.info(f"Промежуточный SAM файл удалён: {sam_file}")
    
    logger.info("=" * 60)
    logger.info("ВЫРАВНИВАНИЕ ЗАВЕРШЕНО")
    logger.info(f"Результат: {output_bam}")
    logger.info("=" * 60)


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 5:
        print("Использование: python alignment.py <contigs.fasta> <forward.fastq> <reverse.fastq> <output.bam> [threads]")
        sys.exit(1)
    
    contigs = sys.argv[1]
    forward = sys.argv[2]
    reverse = sys.argv[3]
    output = sys.argv[4]
    threads = int(sys.argv[5]) if len(sys.argv) > 5 else 4
    
    align_reads_to_contigs(contigs, forward, reverse, output, threads)
