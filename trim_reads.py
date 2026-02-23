import gzip

def trim_reads(input_fastq, output_fastq, min_quality=20, min_length=50):
    """
    Обрезка ридов по качеству и длине.
    
    Args:
        input_fastq: путь к входному FASTQ файлу
        output_fastq: путь к выходному FASTQ файлу
        min_quality: минимальное качество (Phred)
        min_length: минимальная длина рида после обрезки
    """
    open_func = gzip.open if input_fastq.endswith('.gz') else open
    
    with open_func(input_fastq, 'rt') as f_in:
        with open_func(output_fastq, 'wt') as f_out:
            
            while True:
                header = f_in.readline().strip()
                if not header:
                    break
                seq = f_in.readline().strip()
                plus = f_in.readline().strip()
                qual = f_in.readline().strip()
                
                # Конвертируем качество в числа
                quality_scores = [ord(c) - 33 for c in qual]
                
                # Ищем точку обрезки с конца
                trim_pos = len(seq)
                for i in range(len(quality_scores) - 1, -1, -1):
                    if quality_scores[i] >= min_quality:
                        trim_pos = i + 1
                        break
                
                # Обрезаем
                trimmed_seq = seq[:trim_pos]
                trimmed_qual = qual[:trim_pos]
                
                # Фильтруем по длине
                if len(trimmed_seq) >= min_length:
                    f_out.write(header + '\n')
                    f_out.write(trimmed_seq + '\n')
                    f_out.write(plus + '\n')
                    f_out.write(trimmed_qual + '\n')
