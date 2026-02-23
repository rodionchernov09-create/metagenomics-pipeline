#!/usr/bin/env python3

import os
import random
import numpy as np
from typing import List, Dict
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def random_dna(length: int, gc_content: float = 0.5) -> str:
    """Генерирует случайную ДНК последовательность с заданным GC составом."""
    bases = ['A', 'T', 'G', 'C']
    probs = [(1 - gc_content) / 2, (1 - gc_content) / 2, gc_content / 2, gc_content / 2]
    
    return ''.join(random.choices(bases, weights=probs, k=length))


def create_viral_genome() -> str:
    """Создаёт искусственный вирусный геном."""
    length = random.randint(5000, 15000)
    
    # Вирусы часто имеют немного другой GC состав
    gc = random.uniform(0.3, 0.6)
    
    return random_dna(length, gc)


def create_bacterial_genome() -> str:
    """Создаёт искусственный бактериальный геном."""
    length = random.randint(50000, 200000)  # Для теста укороченные
    
    # Бактерии обычно имеют GC состав 0.25-0.75
    gc = random.uniform(0.3, 0.7)
    
    return random_dna(length, gc)


def introduce_mutations(genome: str, mutation_rate: float = 0.001) -> str:
    """Вводит случайные мутации в геном."""
    genome_list = list(genome)
    bases = ['A', 'T', 'G', 'C']
    
    for i in range(len(genome_list)):
        if random.random() < mutation_rate:
            current = genome_list[i]
            options = [b for b in bases if b != current]
            genome_list[i] = random.choice(options)
    
    return ''.join(genome_list)

def generate_reads(genome: str, genome_id: str, coverage: float = 10.0, 
                   read_length: int = 150, insert_size: int = 300) -> List[Dict]:
    """
    Генерирует парные риды из одного генома.
    """
    reads = []
    genome_len = len(genome)
    num_reads = int((genome_len * coverage) / read_length)
    
    # Качество (Phred+33)
    quality = ''.join(chr(33 + 40) for _ in range(read_length))
    
    for i in range(num_reads):
        # Случайная позиция начала
        start = random.randint(0, genome_len - read_length - insert_size)
        
        # Прямой рид
        seq_r1 = genome[start:start + read_length]
        
        # Обратный рид (реверс-комплемент) - ТОТ ЖЕ ГЕНОМ
        seq_r2 = genome[start + insert_size:start + insert_size + read_length]
        seq_r2 = reverse_complement(seq_r2)
        
        # ОДИНАКОВОЕ ИМЯ для пары (ОЧЕНЬ ВАЖНО!)
        read_name = f"@read_{genome_id}_{i:06d}"
        
        reads.append({
            'name': read_name + '/1',
            'seq': seq_r1,
            'qual': quality,
            'genome_id': genome_id,
            'mate': 'R1'
        })
        
        reads.append({
            'name': read_name + '/2',
            'seq': seq_r2,
            'qual': quality,
            'genome_id': genome_id,
            'mate': 'R2'
        })
    
    return reads

def reverse_complement(seq: str) -> str:
    """Возвращает обратно-комплементарную последовательность."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, base) for base in reversed(seq))


def write_fasta(sequences: Dict[str, str], filename: str):
    """Записывает последовательности в FASTA файл."""
    with open(filename, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')
    logger.info(f"Записано {len(sequences)} последовательностей в {filename}")


def write_fastq(reads: List[Dict], filename_r1: str, filename_r2: str):
    """Записывает риды в FASTQ файлы (R1 и R2 отдельно)."""
    # Сортируем по genome_id и индексу, чтобы R1 и R2 шли в одном порядке
    reads_r1 = [r for r in reads if r['mate'] == 'R1']
    reads_r2 = [r for r in reads if r['mate'] == 'R2']
    
    # Сортируем по имени (чтобы пары были на тех же позициях)
    reads_r1.sort(key=lambda x: x['name'])
    reads_r2.sort(key=lambda x: x['name'])
    
    with open(filename_r1, 'w') as f1, open(filename_r2, 'w') as f2:
        for r1, r2 in zip(reads_r1, reads_r2):
            # Проверка, что это действительно пара
            assert r1['name'].replace('/1', '') == r2['name'].replace('/2', ''), \
                   f"Имена не совпадают: {r1['name']} vs {r2['name']}"
            
            f1.write(f"{r1['name']}\n{r1['seq']}\n+\n{r1['qual']}\n")
            f2.write(f"{r2['name']}\n{r2['seq']}\n+\n{r2['qual']}\n")
    
    logger.info(f"Записано {len(reads_r1)} пар ридов в {filename_r1} и {filename_r2}")

def main():
    """Главная функция генерации тестовых данных."""
    logger.info("=" * 60)
    logger.info("ГЕНЕРАЦИЯ ТЕСТОВЫХ ДАННЫХ")
    logger.info("=" * 60)
    
    # Создаём структуру папок
    raw_dir = "../data/raw"
    os.makedirs(raw_dir, exist_ok=True)
    
    # Генерируем геномы
    genomes = {}
    
    # 2 бактерии
    for i in range(1, 3):
        genome = create_bacterial_genome()
        # Добавляем варианты (мутации)
        genome_var = introduce_mutations(genome, 0.001)
        genomes[f'bacteria_{i}'] = genome
        genomes[f'bacteria_{i}_var'] = genome_var
    
    # 2 вируса
    for i in range(1, 3):
        genome = create_viral_genome()
        genomes[f'virus_{i}'] = genome
    
    logger.info(f"Сгенерировано геномов: {len(genomes)}")
    
    # Сохраняем геномы (для справки)
    write_fasta(genomes, "../data/reference_genomes.fa")
    
    # Генерируем риды
    all_reads = []
    sample_name = "test_sample"
    
    for genome_id, genome in genomes.items():
        # Разное покрытие для разных геномов
        if 'virus' in genome_id:
            coverage = 5.0  # Вирусов меньше
        else:
            coverage = 20.0  # Бактерий больше
        
        reads = generate_reads(genome, genome_id, coverage=coverage)
        all_reads.extend(reads)
    
    # Перемешиваем риды (как в реальном эксперименте)
    random.shuffle(all_reads)
    
    # Записываем FASTQ
    write_fastq(
        all_reads,
        f"{raw_dir}/{sample_name}_R1.fastq",
        f"{raw_dir}/{sample_name}_R2.fastq"
    )
    
    # Создаём config.yaml для теста
    config_content = f"""# Тестовая конфигурация
samples: ["{sample_name}"]
raw_reads_dir: "../data/raw"
results_dir: "../results"
min_contig_len: 500
threads: 4
min_quality: 20
min_read_length: 50
kmer_size: 31
viral_threshold: 0.7
"""
    
    with open("config_test.yaml", 'w') as f:
        f.write(config_content)
    
    logger.info("=" * 60)
    logger.info("ГЕНЕРАЦИЯ ЗАВЕРШЕНА")
    logger.info(f"FASTQ файлы: {raw_dir}/{sample_name}_R1.fastq, {raw_dir}/{sample_name}_R2.fastq")
    logger.info(f"Референсные геномы: ../data/reference_genomes.fa")
    logger.info(f"Тестовый конфиг: config_test.yaml")
    logger.info("=" * 60)


if __name__ == "__main__":
    random.seed(42)  # Для воспроизводимости
    np.random.seed(42)
    main()
