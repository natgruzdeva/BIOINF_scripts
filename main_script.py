from helpers.bio_files import (
    derive_single_line_path,
    write_single_line_copy,
    collect_best_descriptions,
    save_lines,
)

from helpers.fastq import (
    is_read_ok,
    make_bounds,
    make_output_path,
    read_next_fastq_record,
    write_fastq_record,
)


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds=(0, 100),
    length_bounds=(0, 2 ** 32),
    quality_threshold=0,
) -> None:
    
    """
    Фильтрация ридов по GC-составу, длине и качеству рида.

    Аргументы:
        input_fastq: путь ко входному FASTQ-файлу.
        output_fastq: имя выходного файла внутри папки filtered.
        gc_bounds: кортеж (низ, верх) или одно число (верхняя граница) -
            интервал GC-состава в процентах, по умолчанию (0, 100).
        length_bounds: кортеж (низ, верх) или одно число (верхняя граница) -
            интервал длины рида, по умолчанию (0, 2**32).
        quality_threshold: минимальное среднее качество рида (phred33),
            по умолчанию 0; риды ниже порога отбрасываются.

    Возвращает:
        None.
    """
    output_path = make_output_path(output_fastq)

    gc_low, gc_high = make_bounds(gc_bounds, 0)
    len_low, len_high = make_bounds(length_bounds, 0)

    with open(input_fastq, "r", encoding="utf-8") as input_file, open(output_path, "w", encoding="utf-8") as output_file:
        record = read_next_fastq_record(input_file)
        while record is not None:
            if is_read_ok(record, gc_low, gc_high,
                          len_low, len_high, quality_threshold):
                write_fastq_record(output_file, record)
            record = read_next_fastq_record(input_file)

def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta=None) -> None:
    """
    Переводит многострочный fasta-формат,в однострочный и сохраняет результат в новый файл.

    В памяти одновременно живёт только одна запись. 

    Аргументы:
        input_fasta: путь до входного многострочного fasta файла.
        output_fasta: путь до выхода; по умолчанию None.

    Возвращает:
        None.
    """
    if output_fasta is None:
        output_fasta = derive_single_line_path(input_fasta)

    with open(input_fasta, "r", encoding="utf-8") as input_file, open(output_fasta, "w", encoding="utf-8") as output_file:
        write_single_line_copy(input_file, output_file)

def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Извлекает описание лучшего совпадения для каждого запроса BLAST,
    удаляет дубликаты и записывает результат одним столбцом в алфавитном
    порядке (без учёта регистра).
    
    Аргументы:
        input_file: путь до txt-файла с результатами BLAST.
        output_file: путь до выходного файла.

    Возвращает:
        None.
    """
    with open(input_file, "r", encoding="utf-8") as blast_file:
        descriptions = collect_best_descriptions(blast_file)

    with open(output_file, "w", encoding="utf-8") as result_file:
        save_lines(result_file, sorted(set(descriptions), key=str.lower))