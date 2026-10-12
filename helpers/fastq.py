import os

# Ввод-вывод одного рида и путь к выводу (папка filtered)

def read_next_fastq_record(file_handle):
    """
    Читает один рид fastq (а именно: 4 строки) из открытого файла.

    Аргументы:
        file_handle: открытый файл на чтение.

    Возвращает:
        Кортеж из четырёх строк (header, sequence, plus_line, quality_string)
        без переносов строки. Если файл закончился - None.
    """
    header = file_handle.readline()
    if not header:
        return None

    sequence = file_handle.readline()
    plus_line = file_handle.readline()
    quality_string = file_handle.readline()

    return (
        header.strip(),
        sequence.strip(),
        plus_line.strip(),
        quality_string.strip(),
    )

def write_fastq_record(file_handle, record) -> None:
    """
    Записывает один рид fastq (а именно: 4 строки) в виде кортежа в открытый файл.

    Аргументы:
        file_handle: открытый файл на запись.
        record: кортеж (header, sequence, plus_line, quality_string).

    Возвращает:
        None.
    """
    header, sequence, plus_line, quality_string = record
    file_handle.write(header + "\n")
    file_handle.write(sequence + "\n")
    file_handle.write(plus_line + "\n")
    file_handle.write(quality_string + "\n")

def make_output_path(output_fastq: str) -> str:
    """
    Создаёт папку filtered (если её ещё нет).

    Аргументы:
        output_fastq: имя выходного файла.

    Возвращает:
        Путь вида filtered/<output_fastq>.
    """
    os.makedirs("filtered", exist_ok=True)
    return os.path.join("filtered", output_fastq)

# Нормализация задаваемых границ

def make_bounds(bounds, default_low: float):
    """
    Приводит границы фильтра к единому виду (нижняя граница, верхняя).

    Если передан кортеж/список из двух чисел, то использует его.
    Если передано одно число, то считает его верхней границей, нижней задает default_low.

    Аргументы:
        bounds: кортеж (низ, верх) или одно число (т.е. верхняя граница).
        default_low: нижняя граница по умолчанию.

    Возвращает:
        Кортеж (низ, верх) из двух чисел.
    """
    if isinstance(bounds, (tuple, list)):
        return (float(bounds[0]), float(bounds[1]))
    return (float(default_low), float(bounds))

# Функции (в т.ч. вычислительные)

def is_length_ok(length: int, low: float, high: float) -> bool:
    """
    Проверяет, попадает ли длина рида в заданные границы (включительно).

    Аргументы:
        length: длина последовательности.
        low: нижняя граница интервала.
        high: верхняя граница интервала.

    Возвращает:
        True, если low <= length <= high, иначе False.
    """
    return low <= length <= high


def calculate_gc_content(sequence: str) -> float:
    """
    Подсчёт GC-состава в процентах.

    Аргументы:
        sequence: нуклеотидная последовательность (str).

    Возвращает:
        GC-состав в процентах. Для пустой строки - 0.
    """
    if not sequence:
        return 0.0
    upper = sequence.upper()
    gc_count = sum(1 for nucleotide in upper if nucleotide in "GC")
    return gc_count / len(upper) * 100.0


def is_gc_ok(gc_percent: float, low: float, high: float) -> bool:
    """
    Проверяет, попадает ли GC-состав в заданные границы (включительно).

    Аргументы:
        gc_percent: GC-состав в процентах.
        low: нижняя граница интервала.
        high: верхняя граница интервала.

    Возвращает:
        True, если low <= gc_percent <= high, иначе False.
    """
    return low <= gc_percent <= high


def calculate_average_quality(quality: str) -> float:
    """
    Подсчёт среднего качества рида (шкала phred33).

    Аргументы:
        quality: строка качества из FASTQ.

    Возвращает:
        Среднее арифметическое качество всех символов. Для пустой строки - 0.
    """
    if not quality:
        return 0.0
    total = sum(ord(symbol) - 33 for symbol in quality)
    return total / len(quality)


def is_quality_ok(average_quality: float, threshold: float) -> bool:
    """
    Проверяет, что среднее качество не ниже порогового значения.

    Аргументы:
        average_quality: среднее качество рида (phred+33).
        threshold: пороговое значение среднего качества.

    Возвращает:
        True, если average_quality >= threshold, иначе False.
    """
    return average_quality >= threshold


# КОмпилятор всех функций
def is_read_ok(record, gc_low: float, gc_high: float,
               len_low: float, len_high: float,
               threshold: float) -> bool:
    """
    Проверка рида по всем трём критериям сразу: gc-состав, длина и качество рида.

    Аргументы:
        record: кортеж (header, sequence, plus_line, quality_string).
        gc_low, gc_high: нижняя и верхняя границы GC-состава.
        len_low, len_high: нижняя и верхняя границы длины.
        threshold: минимально допустимое среднее качество.

    Возвращает:
        True, если рид проходит все три проверки, иначе False.
    """
    _, sequence, _, quality_string = record
    return (
        is_gc_ok(calculate_gc_content(sequence), gc_low, gc_high)
        and is_length_ok(len(sequence), len_low, len_high)
        and is_quality_ok(calculate_average_quality(quality_string), threshold)
    )
