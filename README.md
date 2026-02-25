# Metagenomics Pipeline

Полный пайплайн для метагеномного анализа, написанный с нуля на Python.

## Возможности
- Обрезка ридов по качеству
- Сборка контигов (граф де Брёйна)
- Выравнивание ридов на контиги (BWA)
- Кластеризация контигов (VAE + DBSCAN)
- Предсказание генов (Prodigal-подобный алгоритм)
- Обнаружение вирусных последовательностей (CNN)

## Установка

```bash
# Клонировать репозиторий
git clone https://github.com/rodionchernov09-create/metagenomics-pipeline.git
cd metagenomics-pipeline

# Создать окружение
conda create -n snakemake python=3.10
conda activate snakemake
pip install snakemake numpy scikit-learn torch pysam
conda install -c bioconda bwa samtools
