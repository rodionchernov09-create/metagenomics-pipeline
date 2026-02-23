#!/usr/bin/env python3

import os
import numpy as np
import pysam
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Set
import logging
import itertools
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


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
                current_seq.append(line)
        
        if current_header:
            contigs[current_header] = ''.join(current_seq)
    
    logger.info(f"Загружено {len(contigs)} контигов")
    return contigs


def extract_kmer_frequencies(contigs: Dict[str, str], k: int = 4) -> Dict[str, np.ndarray]:
    logger.info(f"Извлечение частот {k}-меров...")
    
    if not contigs:
        logger.warning("Нет контигов для извлечения частот")
        return {}
    
    if k > 6:
        logger.warning(f"k={k} слишком большой, использую k=4")
        k = 4
    
    import itertools
    bases = ['A', 'C', 'G', 'T']
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
    kmer_to_idx = {kmer: i for i, kmer in enumerate(kmers)}
    
    frequencies = {}
    
    for contig_id, sequence in contigs.items():
        seq = sequence.upper()
        kmer_counts = np.zeros(len(kmers))
        
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if all(c in 'ACGT' for c in kmer):
                if kmer in kmer_to_idx:
                    kmer_counts[kmer_to_idx[kmer]] += 1
        
        total = kmer_counts.sum()
        if total > 0:
            frequencies[contig_id] = kmer_counts / total
        else:
            frequencies[contig_id] = np.zeros(len(kmers))
    
    logger.info(f"Извлечены частоты для {len(frequencies)} контигов")
    return frequencies

def extract_coverage(bam_file: str, contigs: Dict[str, str]) -> Dict[str, float]:
    logger.info(f"Извлечение покрытия из {bam_file}...")
    
    coverage = {}
    
    if not contigs:
        logger.warning("Нет контигов, возвращаю пустое покрытие")
        return coverage
    
    try:
        with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
            for contig_id in contigs.keys():
                try:
                    coverages = []
                    for pileup in bam.pileup(contig_id):
                        coverages.append(pileup.n)
                    
                    if coverages:
                        coverage[contig_id] = np.mean(coverages)
                    else:
                        coverage[contig_id] = 0.0
                except ValueError:
                    coverage[contig_id] = 0.0
    except Exception as e:
        logger.warning(f"Ошибка при чтении BAM файла: {e}")
        for contig_id in contigs.keys():
            coverage[contig_id] = 0.0
    
    logger.info(f"Извлечено покрытие для {len(coverage)} контигов")
    return coverage

def create_feature_matrix(
    contigs: Dict[str, str],
    kmer_freqs: Dict[str, np.ndarray],
    coverage: Dict[str, float]
) -> Tuple[np.ndarray, List[str]]:
    logger.info("Создание матрицы признаков...")
    
    if not contigs:
        logger.warning("Нет контигов, возвращаю пустую матрицу")
        return np.array([]), []
    
    contig_ids = list(contigs.keys())
    
    if not kmer_freqs or not coverage:
        logger.warning("Нет признаков для создания матрицы")
        return np.array([]), contig_ids
    
    kmer_dim = len(next(iter(kmer_freqs.values())))
    
    X = np.zeros((len(contig_ids), kmer_dim + 1))
    
    for i, contig_id in enumerate(contig_ids):
        if contig_id in kmer_freqs:
            X[i, :kmer_dim] = kmer_freqs[contig_id]
        else:
            X[i, :kmer_dim] = np.zeros(kmer_dim)
        
        if contig_id in coverage:
            X[i, kmer_dim] = np.log1p(coverage[contig_id])
        else:
            X[i, kmer_dim] = 0.0
    
    logger.info(f"Матрица признаков: {X.shape}")
    return X, contig_ids

class VAE(nn.Module):
    def __init__(self, input_dim: int, hidden_dim: int = 256, latent_dim: int = 32):
        super(VAE, self).__init__()
        
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU()
        )
        
        self.mu_layer = nn.Linear(hidden_dim, latent_dim)
        self.logvar_layer = nn.Linear(hidden_dim, latent_dim)
        
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, input_dim),
            nn.Sigmoid()
        )
    
    def encode(self, x):
        h = self.encoder(x)
        return self.mu_layer(h), self.logvar_layer(h)
    
    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std
    
    def decode(self, z):
        return self.decoder(z)
    
    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recon_x = self.decode(z)
        return recon_x, mu, logvar


def train_vae(X: np.ndarray, latent_dim: int = 32, epochs: int = 50, batch_size: int = 64) -> VAE:
    logger.info("Обучение вариационного автоэнкодера...")
    
    # ПРОВЕРКА: убедимся, что данные в [0,1]
    if X.min() < 0 or X.max() > 1:
        logger.warning(f"Данные вне диапазона [0,1]: min={X.min():.3f}, max={X.max():.3f}. Нормализую...")
        X_min = X.min(axis=0)
        X_max = X.max(axis=0)
        X_range = X_max - X_min
        X_range[X_range == 0] = 1
        X = (X - X_min) / X_range
        logger.info(f"После нормализации: min={X.min():.3f}, max={X.max():.3f}")
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f"Используется устройство: {device}")
    
    X_tensor = torch.FloatTensor(X).to(device)
    dataset = TensorDataset(X_tensor)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
    
    input_dim = X.shape[1]
    model = VAE(input_dim, latent_dim=latent_dim).to(device)
    optimizer = optim.Adam(model.parameters(), lr=1e-3)
    
    def loss_function(recon_x, x, mu, logvar):
        BCE = nn.functional.binary_cross_entropy(recon_x, x, reduction='sum')
        KLD = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
        return BCE + KLD
    
    model.train()
    for epoch in range(epochs):
        train_loss = 0
        for batch in dataloader:
            batch = batch[0]
            optimizer.zero_grad()
            recon_batch, mu, logvar = model(batch)
            loss = loss_function(recon_batch, batch, mu, logvar)
            loss.backward()
            train_loss += loss.item()
            optimizer.step()
        
        if epoch % 10 == 0:
            avg_loss = train_loss / len(dataloader.dataset)
            logger.info(f"Epoch {epoch}: loss = {avg_loss:.4f}")
    
    logger.info("Обучение завершено")
    return model

def extract_latent_representations(model: VAE, X: np.ndarray) -> np.ndarray:
    logger.info("Извлечение латентных представлений...")
    
    device = next(model.parameters()).device
    X_tensor = torch.FloatTensor(X).to(device)
    
    model.eval()
    with torch.no_grad():
        mu, _ = model.encode(X_tensor)
    
    return mu.cpu().numpy()


def cluster_contigs(latent_reps: np.ndarray, eps: float = 0.5, min_samples: int = 5) -> np.ndarray:
    logger.info(f"Кластеризация DBSCAN с eps={eps}, min_samples={min_samples}...")
    
    scaler = StandardScaler()
    scaled_reps = scaler.fit_transform(latent_reps)
    
    pca = PCA(n_components=min(20, scaled_reps.shape[1]))
    pca_reps = pca.fit_transform(scaled_reps)
    
    clustering = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1)
    labels = clustering.fit_predict(pca_reps)
    
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = list(labels).count(-1)
    
    logger.info(f"Найдено {n_clusters} кластеров, {n_noise} шумовых точек")
    return labels


def write_clusters(contig_ids: List[str], labels: np.ndarray, output_tsv: str) -> None:
    logger.info(f"Запись кластеров в {output_tsv}...")
    
    os.makedirs(os.path.dirname(output_tsv), exist_ok=True)
    
    with open(output_tsv, 'w') as f:
        f.write("cluster_id\tcontig_id\n")
        
        for contig_id, label in zip(contig_ids, labels):
            if label != -1:
                f.write(f"cluster_{label}\t{contig_id}\n")
    
    logger.info(f"Записано {len([l for l in labels if l != -1])} контигов в кластеры")


def bin_contigs(
    contigs_fasta: str,
    bam_file: str,
    output_clusters: str,
    kmer_size: int = 4,
    latent_dim: int = 32,
    dbscan_eps: float = 0.5,
    dbscan_min_samples: int = 5,
    vae_epochs: int = 50
) -> None:
    logger.info("=" * 60)
    logger.info("НАЧАЛО КЛАСТЕРИЗАЦИИ КОНТИГОВ (VAMB)")
    logger.info("=" * 60)
    
    contigs = read_fasta_contigs(contigs_fasta)

    if not contigs:
        logger.warning("Нет контигов для кластеризации, создаю пустой выходной файл")
        os.makedirs(os.path.dirname(output_clusters), exist_ok=True)
        with open(output_clusters, 'w') as f:
            f.write("cluster_id\tcontig_id\n")
        logger.info("=" * 60)
        logger.info("КЛАСТЕРИЗАЦИЯ ЗАВЕРШЕНА (пустой результат)")
        logger.info("=" * 60)
        return
    
    kmer_freqs = extract_kmer_frequencies(contigs, k=kmer_size)
    
    coverage = extract_coverage(bam_file, contigs)
    
    X, contig_ids = create_feature_matrix(contigs, kmer_freqs, coverage)
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    vae_model = train_vae(X_scaled, latent_dim=latent_dim, epochs=vae_epochs)
    
    latent_reps = extract_latent_representations(vae_model, X_scaled)
    
    labels = cluster_contigs(latent_reps, eps=dbscan_eps, min_samples=dbscan_min_samples)
    
    write_clusters(contig_ids, labels, output_clusters)
    
    logger.info("=" * 60)
    logger.info("КЛАСТЕРИЗАЦИЯ ЗАВЕРШЕНА")
    logger.info(f"Результат: {output_clusters}")
    logger.info("=" * 60)


if __name__ == "__main__":
    import sys
    import itertools
    
    if len(sys.argv) < 4:
        print("Использование: python binning.py <contigs.fasta> <alignments.bam> <output_clusters.tsv>")
        sys.exit(1)
    
    contigs = sys.argv[1]
    bam = sys.argv[2]
    output = sys.argv[3]
    
    bin_contigs(contigs, bam, output)
