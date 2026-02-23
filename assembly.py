#!/usr/bin/env python3

import itertools
from collections import defaultdict, Counter
from typing import Dict, List, Set, Tuple, Optional
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def read_fastq(filepath: str) -> List[Tuple[str, str]]:
    sequences = []
    open_func = open
    
    if filepath.endswith('.gz'):
        import gzip
        open_func = gzip.open
    
    try:
        with open_func(filepath, 'rt') as f:
            while True:
                header = f.readline().strip()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline()
                f.readline()
                
                if seq:
                    sequences.append((header, seq))
    except Exception as e:
        logger.error(f"Ошибка чтения файла {filepath}: {e}")
        raise
    
    logger.info(f"Прочитано {len(sequences)} последовательностей из {filepath}")
    return sequences


def build_debruijn_graph(sequences: List[str], k: int = 31) -> Dict[str, Set[str]]:
    if k % 2 == 0:
        logger.warning(f"Рекомендуется нечётный k (чётный {k} может вызвать проблемы с палиндромами)")
    
    graph = defaultdict(set)
    kmer_counts = Counter()
    
    for seq_idx, seq in enumerate(sequences):
        if len(seq) < k + 1:
            continue
            
        for i in range(len(seq) - k):
            kmer1 = seq[i:i + k]
            kmer2 = seq[i + 1:i + k + 1]
            
            if not all(c in 'ACGTacgt' for c in kmer1 + kmer2):
                continue
                
            kmer1 = kmer1.upper()
            kmer2 = kmer2.upper()
            
            graph[kmer1].add(kmer2)
            kmer_counts[kmer1] += 1
            kmer_counts[kmer2] += 1
    
    logger.info(f"Построен граф с {len(graph)} узлами и {sum(len(v) for v in graph.values())} рёбрами")
    return graph


def find_contigs(graph: Dict[str, Set[str]]) -> List[str]:
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    
    all_nodes = set(graph.keys())
    for node, neighbors in graph.items():
        out_degree[node] = len(neighbors)
        for neighbor in neighbors:
            in_degree[neighbor] += 1
            all_nodes.add(neighbor)
    
    start_nodes = []
    for node in all_nodes:
        if out_degree[node] > in_degree[node] or (out_degree[node] > 0 and in_degree[node] == 0):
            start_nodes.append(node)
    
    if not start_nodes:
        start_nodes = [node for node in all_nodes if out_degree[node] > 0]
    
    logger.info(f"Найдено {len(start_nodes)} стартовых узлов")
    
    contigs = []
    visited_edges = set()
    
    for start in start_nodes:
        if out_degree[start] == 0:
            continue
            
        path = [start]
        current = start
        
        while True:
            neighbors = [n for n in graph.get(current, []) 
                        if (current, n) not in visited_edges]
            
            if not neighbors:
                break
                
            if len(neighbors) == 1:
                next_node = neighbors[0]
                visited_edges.add((current, next_node))
                path.append(next_node)
                current = next_node
            else:
                break
        
        if len(path) > 1:
            contig = path[0]
            for node in path[1:]:
                contig += node[-1]
            contigs.append(contig)
    
    logger.info(f"Собрано {len(contigs)} контигов")
    return contigs


def filter_contigs(contigs: List[str], min_length: int = 1000) -> List[str]:
    filtered = [c for c in contigs if len(c) >= min_length]
    logger.info(f"После фильтрации по длине >= {min_length}: {len(filtered)} контигов")
    return filtered


def write_fasta(contigs: List[str], output_file: str, prefix: str = "contig"):
    with open(output_file, 'w') as f:
        for i, contig in enumerate(contigs, 1):
            f.write(f">{prefix}_{i} length={len(contig)}\n")
            for j in range(0, len(contig), 80):
                f.write(contig[j:j+80] + '\n')
    
    logger.info(f"Записано {len(contigs)} контигов в {output_file}")


def assemble_contigs(forward_reads: str, reverse_reads: str, output_fasta: str, 
                     min_contig_len: int = 1000, kmer_size: int = 31) -> None:
    logger.info("=" * 60)
    logger.info("НАЧАЛО СБОРКИ КОНТИГОВ")
    logger.info("=" * 60)
    
    logger.info(f"Чтение прямых ридов: {forward_reads}")
    forward_seqs = [seq for _, seq in read_fastq(forward_reads)]
    
    logger.info(f"Чтение обратных ридов: {reverse_reads}")
    reverse_seqs = [seq for _, seq in read_fastq(reverse_reads)]
    
    all_seqs = forward_seqs + reverse_seqs
    logger.info(f"Всего последовательностей: {len(all_seqs)}")
    
    if not all_seqs:
        logger.error("Нет последовательностей для сборки!")
        with open(output_fasta, 'w') as f:
            f.write("")
        return
    
    logger.info(f"Построение графа де Брёйна с k={kmer_size}...")
    graph = build_debruijn_graph(all_seqs, k=kmer_size)
    
    if not graph:
        logger.warning("Граф пуст! Проверьте входные данные.")
        with open(output_fasta, 'w') as f:
            f.write("")
        return
    
    logger.info("Поиск контигов в графе...")
    contigs = find_contigs(graph)
    
    if not contigs:
        logger.warning("Контиги не найдены!")
        with open(output_fasta, 'w') as f:
            f.write("")
        return
    
    logger.info(f"Фильтрация контигов (минимальная длина: {min_contig_len})...")
    filtered_contigs = filter_contigs(contigs, min_contig_len)
    
    logger.info(f"Запись результатов в {output_fasta}...")
    write_fasta(filtered_contigs, output_fasta)
    
    total_bases = sum(len(c) for c in filtered_contigs)
    logger.info("=" * 60)
    logger.info("ИТОГИ СБОРКИ:")
    logger.info(f"  Всего контигов: {len(filtered_contigs)}")
    logger.info(f"  Всего баз: {total_bases}")
    if filtered_contigs:
        logger.info(f"  N50: {calculate_n50(filtered_contigs)}")
        logger.info(f"  Макс. длина: {max(len(c) for c in filtered_contigs)}")
        logger.info(f"  Средняя длина: {total_bases // len(filtered_contigs)}")
    logger.info("=" * 60)


def calculate_n50(contigs: List[str]) -> int:
    if not contigs:
        return 0
        
    sorted_contigs = sorted(contigs, key=len, reverse=True)
    total_length = sum(len(c) for c in contigs)
    half_length = total_length // 2
    
    cumulative = 0
    for contig in sorted_contigs:
        cumulative += len(contig)
        if cumulative >= half_length:
            return len(contig)
    
    return 0


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 4:
        print("Использование: python assembly.py <forward.fastq> <reverse.fastq> <output.fasta> [min_length]")
        sys.exit(1)
    
    forward = sys.argv[1]
    reverse = sys.argv[2]
    output = sys.argv[3]
    min_len = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
    
    assemble_contigs(forward, reverse, output, min_len)
