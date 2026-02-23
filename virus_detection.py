#!/usr/bin/env python3

import os
import numpy as np
from collections import defaultdict
from typing import Dict, List, Tuple, Optional
import logging
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader

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
                current_seq.append(line.upper())
        
        if current_header:
            contigs[current_header] = ''.join(current_seq)
    
    logger.info(f"Загружено {len(contigs)} контигов")
    return contigs


def encode_sequence(seq: str, max_len: int = 1000) -> np.ndarray:
    """
    Кодирует последовательность ДНК в one-hot представление.
    
    Args:
        seq: строка ДНК
        max_len: максимальная длина (обрезает или дополняет нулями)
    
    Returns:
        Массив shape (4, max_len) с one-hot кодированием
    """
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    encoding = np.zeros((4, max_len), dtype=np.float32)
    
    for i, base in enumerate(seq[:max_len]):
        if base in mapping:
            encoding[mapping[base], i] = 1.0
    
    return encoding


def compute_kmer_frequencies(seq: str, k: int = 4) -> np.ndarray:
    """
    Вычисляет частоты k-меров в последовательности.
    
    Args:
        seq: строка ДНК
        k: размер k-мера
    
    Returns:
        Массив частот для всех возможных k-меров
    """
    bases = ['A', 'C', 'G', 'T']
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
    kmer_to_idx = {kmer: i for i, kmer in enumerate(kmers)}
    
    frequencies = np.zeros(len(kmers), dtype=np.float32)
    total = 0
    
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if kmer in kmer_to_idx:
            frequencies[kmer_to_idx[kmer]] += 1
            total += 1
    
    if total > 0:
        frequencies /= total
    
    return frequencies


class ContigDataset(Dataset):
    """Dataset для загрузки контигов."""
    
    def __init__(self, contigs: Dict[str, str], max_len: int = 1000):
        self.contig_ids = list(contigs.keys())
        self.sequences = [contigs[cid] for cid in self.contig_ids]
        self.max_len = max_len
    
    def __len__(self):
        return len(self.contig_ids)
    
    def __getitem__(self, idx):
        seq = self.sequences[idx]
        encoding = encode_sequence(seq, self.max_len)
        return torch.FloatTensor(encoding)

class DeepVirFinder(nn.Module):
    """
    Нейросеть для обнаружения вирусных последовательностей.
    Архитектура: CNN + полносвязные слои.
    """
    
    def __init__(self, input_channels: int = 4, seq_len: int = 1000, use_kmer: bool = False):
        super(DeepVirFinder, self).__init__()
        self.use_kmer = use_kmer
        
        self.conv1 = nn.Conv1d(input_channels, 64, kernel_size=8, padding=4)
        self.conv2 = nn.Conv1d(64, 64, kernel_size=8, padding=4)
        self.conv3 = nn.Conv1d(64, 128, kernel_size=8, padding=4)
        self.conv4 = nn.Conv1d(128, 128, kernel_size=8, padding=4)
        
        self.pool = nn.MaxPool1d(2)
        self.dropout = nn.Dropout(0.3)
        
        # Вычисляем размер после свёрток
        with torch.no_grad():
            test_input = torch.randn(1, input_channels, seq_len)
            x = F.relu(self.conv1(test_input))
            x = self.pool(x)
            x = F.relu(self.conv2(x))
            x = self.pool(x)
            x = F.relu(self.conv3(x))
            x = self.pool(x)
            x = F.relu(self.conv4(x))
            x = self.pool(x)
            self.conv_out_size = x.view(1, -1).size(1)
        
        # Размер входа для fc1 зависит от того, используем ли мы kmer
        fc1_input_size = self.conv_out_size
        if use_kmer:
            fc1_input_size += 256
            self.kmer_fc = nn.Linear(256, 256)
        
        self.fc1 = nn.Linear(fc1_input_size, 256)
        self.fc2 = nn.Linear(256, 128)
        self.fc3 = nn.Linear(128, 1)
        
    def forward(self, x, kmer_features=None):
        batch_size = x.size(0)
        
        x = F.relu(self.conv1(x))
        x = self.pool(x)
        x = self.dropout(x)
        
        x = F.relu(self.conv2(x))
        x = self.pool(x)
        x = self.dropout(x)
        
        x = F.relu(self.conv3(x))
        x = self.pool(x)
        x = self.dropout(x)
        
        x = F.relu(self.conv4(x))
        x = self.pool(x)
        x = self.dropout(x)
        
        x = x.view(batch_size, -1)
        
        if self.use_kmer and kmer_features is not None:
            kmer_features = F.relu(self.kmer_fc(kmer_features))
            x = torch.cat([x, kmer_features], dim=1)
        
        x = F.relu(self.fc1(x))
        x = self.dropout(x)
        x = F.relu(self.fc2(x))
        x = self.dropout(x)
        x = self.fc3(x)
        
        return torch.sigmoid(x)

class ViralDetector:
    """Класс для обнаружения вирусных последовательностей."""
    
    def __init__(self, model_path: Optional[str] = None, device: str = 'cpu'):
        self.device = torch.device(device if torch.cuda.is_available() and device == 'cuda' else 'cpu')
        logger.info(f"Используется устройство: {self.device}")
        
        self.model = DeepVirFinder(use_kmer=False).to(self.device)
        
        if model_path and os.path.exists(model_path):
            self.load_model(model_path)
        else:
            logger.warning("Модель не загружена, будут использоваться случайные веса")
    
    def load_model(self, model_path: str):
        """Загружает предобученную модель."""
        logger.info(f"Загрузка модели из {model_path}...")
        self.model.load_state_dict(torch.load(model_path, map_location=self.device))
        self.model.eval()
    
    def predict(self, contigs: Dict[str, str], batch_size: int = 32) -> Dict[str, float]:
        """
        Предсказывает вероятность вирусности для каждого контига.
        
        Returns:
            Словарь: id контига → вероятность вирусности (0-1)
        """
        logger.info("Предсказание вирусных контигов...")
        
        dataset = ContigDataset(contigs)
        dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=False)
        
        predictions = {}
        self.model.eval()
        
        with torch.no_grad():
            for i, batch in enumerate(dataloader):
                batch = batch.to(self.device)
                
                outputs = self.model(batch)
                
                for j, prob in enumerate(outputs.cpu().numpy()):
                    idx = i * batch_size + j
                    if idx < len(contigs):
                        contig_id = list(contigs.keys())[idx]
                        predictions[contig_id] = float(prob[0])
        
        logger.info(f"Предсказания получены для {len(predictions)} контигов")
        return predictions


def train_viral_detector(
    train_contigs: Dict[str, str],
    train_labels: Dict[str, int],
    val_contigs: Optional[Dict[str, str]] = None,
    val_labels: Optional[Dict[str, int]] = None,
    epochs: int = 20,
    batch_size: int = 32,
    learning_rate: float = 0.001,
    model_output: Optional[str] = None
) -> DeepVirFinder:
    """
    Обучает детектор вирусов на размеченных данных.
    
    Args:
        train_contigs: словарь обучающих контигов
        train_labels: словарь меток (1 - вирус, 0 - не вирус)
        val_contigs: словарь валидационных контигов
        val_labels: словарь валидационных меток
        epochs: количество эпох
        batch_size: размер батча
        learning_rate: скорость обучения
        model_output: путь для сохранения модели
    
    Returns:
        Обученная модель
    """
    logger.info("=" * 60)
    logger.info("НАЧАЛО ОБУЧЕНИЯ МОДЕЛИ")
    logger.info("=" * 60)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f"Используется устройство: {device}")
    
    model = DeepVirFinder().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    criterion = nn.BCELoss()
    
    train_ids = list(train_contigs.keys())
    train_seqs = [train_contigs[cid] for cid in train_ids]
    train_encodings = torch.FloatTensor(np.array([encode_sequence(seq) for seq in train_seqs]))
    train_targets = torch.FloatTensor([train_labels[cid] for cid in train_ids]).view(-1, 1)
    
    train_dataset = TensorDataset(train_encodings, train_targets)
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    
    if val_contigs and val_labels:
        val_ids = list(val_contigs.keys())
        val_seqs = [val_contigs[cid] for cid in val_ids]
        val_encodings = torch.FloatTensor(np.array([encode_sequence(seq) for seq in val_seqs]))
        val_targets = torch.FloatTensor([val_labels[cid] for cid in val_ids]).view(-1, 1)
        val_dataset = TensorDataset(val_encodings, val_targets)
        val_loader = DataLoader(val_dataset, batch_size=batch_size)
    
    for epoch in range(epochs):
        model.train()
        train_loss = 0.0
        
        for batch_x, batch_y in train_loader:
            batch_x, batch_y = batch_x.to(device), batch_y.to(device)
            
            optimizer.zero_grad()
            outputs = model(batch_x)
            loss = criterion(outputs, batch_y)
            loss.backward()
            optimizer.step()
            
            train_loss += loss.item()
        
        avg_train_loss = train_loss / len(train_loader)
        
        if val_contigs and val_labels:
            model.eval()
            val_loss = 0.0
            correct = 0
            total = 0
            
            with torch.no_grad():
                for batch_x, batch_y in val_loader:
                    batch_x, batch_y = batch_x.to(device), batch_y.to(device)
                    outputs = model(batch_x)
                    loss = criterion(outputs, batch_y)
                    val_loss += loss.item()
                    
                    predicted = (outputs > 0.5).float()
                    total += batch_y.size(0)
                    correct += (predicted == batch_y).sum().item()
            
            avg_val_loss = val_loss / len(val_loader)
            accuracy = correct / total
            
            logger.info(f"Epoch {epoch+1}/{epochs} - "
                       f"Train Loss: {avg_train_loss:.4f}, "
                       f"Val Loss: {avg_val_loss:.4f}, "
                       f"Val Acc: {accuracy:.4f}")
        else:
            logger.info(f"Epoch {epoch+1}/{epochs} - Train Loss: {avg_train_loss:.4f}")
    
    if model_output:
        os.makedirs(os.path.dirname(model_output), exist_ok=True)
        torch.save(model.state_dict(), model_output)
        logger.info(f"Модель сохранена в {model_output}")
    
    logger.info("=" * 60)
    logger.info("ОБУЧЕНИЕ ЗАВЕРШЕНО")
    logger.info("=" * 60)
    
    return model


def detect_viruses(
    contigs_fasta: str,
    output_predictions: str,
    model_path: Optional[str] = None,
    threshold: float = 0.7
) -> None:
    """
    Главная функция для обнаружения вирусов.
    
    Args:
        contigs_fasta: путь к FASTA файлу с контигами
        output_predictions: путь к выходному файлу с предсказаниями
        model_path: путь к предобученной модели (если None, используются случайные веса)
        threshold: порог для классификации
    """
    logger.info("=" * 60)
    logger.info("НАЧАЛО ОБНАРУЖЕНИЯ ВИРУСОВ")
    logger.info("=" * 60)
    
    contigs = read_fasta_contigs(contigs_fasta)
    
    detector = ViralDetector(model_path)
    
    predictions = detector.predict(contigs)
    
    os.makedirs(os.path.dirname(output_predictions), exist_ok=True)
    
    viral_count = 0
    with open(output_predictions, 'w') as f:
        f.write("contig_id\tlength\tprobability\tprediction\n")
        
        for contig_id, prob in predictions.items():
            is_viral = prob >= threshold
            if is_viral:
                viral_count += 1
            
            f.write(f"{contig_id}\t{len(contigs[contig_id])}\t{prob:.6f}\t{'VIRAL' if is_viral else 'non-viral'}\n")
    
    logger.info(f"Предсказания записаны в {output_predictions}")
    logger.info(f"Найдено {viral_count} вирусных контигов (порог {threshold})")
    
    viral_contigs = [cid for cid, prob in predictions.items() if prob >= threshold]
    
    viral_fasta = os.path.join(os.path.dirname(output_predictions), 'viral_contigs.fa')
    with open(viral_fasta, 'w') as f:
        for cid in viral_contigs:
            f.write(f">{cid}\n")
            seq = contigs[cid]
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')
    
    logger.info(f"Вирусные контиги сохранены в {viral_fasta}")
    
    logger.info("=" * 60)
    logger.info("ОБНАРУЖЕНИЕ ВИРУСОВ ЗАВЕРШЕНО")
    logger.info("=" * 60)


if __name__ == "__main__":
    import sys
    import itertools
    
    if len(sys.argv) < 3:
        print("Использование: python virus_detection.py <contigs.fasta> <output_predictions.txt> [model_path] [threshold]")
        sys.exit(1)
    
    contigs = sys.argv[1]
    output = sys.argv[2]
    model = sys.argv[3] if len(sys.argv) > 3 else None
    thresh = float(sys.argv[4]) if len(sys.argv) > 4 else 0.7
    
    detect_viruses(contigs, output, model, thresh)
