#!/usr/bin/env python3

import os
import re
from collections import defaultdict
from typing import Dict, List, Tuple, Optional
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# Генетический код (стандартный)
GENETIC_CODE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

START_CODONS = {'ATG', 'GTG', 'TTG'}
STOP_CODONS = {'TAA', 'TAG', 'TGA'}


def read_fasta_contigs(contigs_fasta: str) -> Dict[str, str]:
    logger.info(f"Чтение контигов из {contigs_fasta}...")
    
    contigs = {}
    current_header = None
    current_seq = []
    
    with open(contigs_fasta, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    contigs[current_header] = ''.join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        if current_header:
            contigs[current_header] = ''.join(current_seq)
    
    logger.info(f"Загружено {len(contigs)} контигов")
    return contigs


def find_orfs_in_sequence(sequence: str, min_length: int = 90) -> List[Dict]:
    """
    Находит все открытые рамки считывания (ORF) в последовательности.
    
    Args:
        sequence: строка ДНК
        min_length: минимальная длина ORF в нуклеотидах
    
    Returns:
        Список ORF с информацией о позициях и рамке
    """
    orfs = []
    seq_len = len(sequence)
    
    for frame in range(3):
        i = frame
        while i < seq_len - 2:
            codon = sequence[i:i+3]
            
            if codon in START_CODONS:
                start_pos = i
                orf_seq = []
                
                j = i
                while j < seq_len - 2:
                    current_codon = sequence[j:j+3]
                    orf_seq.append(current_codon)
                    
                    if current_codon in STOP_CODONS and len(orf_seq) > 1:
                        end_pos = j + 3
                        orf_length = end_pos - start_pos
                        
                        if orf_length >= min_length:
                            orfs.append({
                                'start': start_pos,
                                'end': end_pos,
                                'frame': frame + 1,
                                'strand': '+',
                                'length': orf_length,
                                'sequence': sequence[start_pos:end_pos],
                                'codons': orf_seq
                            })
                        break
                    
                    j += 3
                
                i = j
            else:
                i += 3
    
    return orfs


def reverse_complement(seq: str) -> str:
    """Возвращает обратно-комплементарную последовательность."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, base) for base in reversed(seq))


def find_orfs_both_strands(sequence: str, min_length: int = 90) -> List[Dict]:
    """Находит ORF на обеих цепях ДНК."""
    orfs = find_orfs_in_sequence(sequence, min_length)
    
    rev_seq = reverse_complement(sequence)
    rev_orfs = find_orfs_in_sequence(rev_seq, min_length)
    
    for orf in rev_orfs:
        orf['strand'] = '-'
        orig_start = len(sequence) - orf['end']
        orig_end = len(sequence) - orf['start']
        orf['start'] = orig_start
        orf['end'] = orig_end
        orf['sequence'] = reverse_complement(orf['sequence'])
    
    orfs.extend(rev_orfs)
    
    orfs.sort(key=lambda x: x['start'])
    return orfs


def translate_orf(orf_seq: str) -> str:
    """Транслирует ORF в белковую последовательность."""
    protein = []
    
    for i in range(0, len(orf_seq) - 2, 3):
        codon = orf_seq[i:i+3]
        if len(codon) == 3:
            aa = GENETIC_CODE.get(codon, 'X')
            protein.append(aa)
    
    return ''.join(protein)


def score_orf(orf: Dict) -> float:
    """Оценивает качество ORF на основе старт-кодона и длины."""
    score = 0.0
    
    if orf['sequence'][:3] == 'ATG':
        score += 2.0
    elif orf['sequence'][:3] in {'GTG', 'TTG'}:
        score += 1.0
    
    score += orf['length'] / 100.0
    
    protein = translate_orf(orf['sequence'])
    if protein and protein[-1] == '*':
        score += 1.0
    
    return score


def select_best_orfs(orfs: List[Dict], max_overlap: float = 0.5) -> List[Dict]:
    """Выбирает неперекрывающиеся ORF с наивысшими оценками."""
    if not orfs:
        return []
    
    for i, orf in enumerate(orfs):
        orf['score'] = score_orf(orf)
        orf['id'] = i
    
    orfs.sort(key=lambda x: x['score'], reverse=True)
    
    selected = []
    used_positions = []
    
    for orf in orfs:
        overlapping = False
        
        for start, end in used_positions:
            overlap_start = max(orf['start'], start)
            overlap_end = min(orf['end'], end)
            
            if overlap_end > overlap_start:
                overlap_len = overlap_end - overlap_start
                orf_len = orf['end'] - orf['start']
                
                if overlap_len / orf_len > max_overlap:
                    overlapping = True
                    break
        
        if not overlapping:
            selected.append(orf)
            used_positions.append((orf['start'], orf['end']))
    
    selected.sort(key=lambda x: x['start'])
    return selected


def predict_genes_in_contig(contig_id: str, sequence: str, min_gene_length: int = 90) -> List[Dict]:
    """Предсказывает гены в одном контиге."""
    orfs = find_orfs_both_strands(sequence, min_gene_length)
    selected_orfs = select_best_orfs(orfs)
    
    genes = []
    for i, orf in enumerate(selected_orfs, 1):
        protein_seq = translate_orf(orf['sequence'])
        
        gene = {
            'gene_id': f"{contig_id}_gene_{i}",
            'contig': contig_id,
            'start': orf['start'] + 1,
            'end': orf['end'],
            'frame': orf['frame'],
            'strand': orf['strand'],
            'length': orf['length'],
            'nucleotide_seq': orf['sequence'],
            'protein_seq': protein_seq[:-1] if protein_seq.endswith('*') else protein_seq
        }
        genes.append(gene)
    
    return genes


def write_proteins_fasta(genes: List[Dict], output_file: str) -> None:
    """Записывает белковые последовательности в FASTA."""
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, 'w') as f:
        for gene in genes:
            f.write(f">{gene['gene_id']} | contig:{gene['contig']} "
                   f"strand:{gene['strand']} length:{gene['length']}bp\n")
            
            protein = gene['protein_seq']
            for i in range(0, len(protein), 60):
                f.write(protein[i:i+60] + '\n')
    
    logger.info(f"Записано {len(genes)} белков в {output_file}")


def write_genes_gff(genes: List[Dict], output_file: str) -> None:
    """Записывает информацию о генах в GFF формат."""
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, 'w') as f:
        f.write("##gff-version 3\n")
        
        for gene in genes:
            f.write(f"{gene['contig']}\tProdigal_v2.6.3\tCDS\t"
                   f"{gene['start']}\t{gene['end']}\t.\t{gene['strand']}\t0\t"
                   f"ID={gene['gene_id']};Name={gene['gene_id']}\n")
    
    logger.info(f"Записано {len(genes)} генов в {output_file}")


def predict_genes(
    contigs_fasta: str,
    output_proteins: str,
    output_gff: Optional[str] = None,
    min_gene_length: int = 90
) -> None:
    """
    Главная функция для предсказания генов.
    
    Args:
        contigs_fasta: путь к FASTA файлу с контигами
        output_proteins: путь к выходному FASTA файлу с белками
        output_gff: путь к выходному GFF файлу (опционально)
        min_gene_length: минимальная длина гена в нуклеотидах
    """
    logger.info("=" * 60)
    logger.info("НАЧАЛО ПРЕДСКАЗАНИЯ ГЕНОВ")
    logger.info("=" * 60)
    
    contigs = read_fasta_contigs(contigs_fasta)
    
    all_genes = []
    
    for contig_id, sequence in contigs.items():
        if len(sequence) < min_gene_length:
            continue
        
        genes = predict_genes_in_contig(contig_id, sequence, min_gene_length)
        all_genes.extend(genes)
        
        logger.info(f"Контиг {contig_id}: найдено {len(genes)} генов")
    
    logger.info(f"Всего предсказано {len(all_genes)} генов")
    
    write_proteins_fasta(all_genes, output_proteins)
    
    if output_gff:
        write_genes_gff(all_genes, output_gff)
    
    logger.info("=" * 60)
    logger.info("ПРЕДСКАЗАНИЕ ГЕНОВ ЗАВЕРШЕНО")
    logger.info(f"Белки: {output_proteins}")
    if output_gff:
        logger.info(f"GFF: {output_gff}")
    logger.info("=" * 60)


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 3:
        print("Использование: python gene_prediction.py <contigs.fasta> <output_proteins.faa> [output.gff]")
        sys.exit(1)
    
    contigs = sys.argv[1]
    proteins = sys.argv[2]
    gff = sys.argv[3] if len(sys.argv) > 3 else None
    
    predict_genes(contigs, proteins, gff)
