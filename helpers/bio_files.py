import os

def write_single_line_copy(input_file, output_file) -> None:
    """
    Читает fasta file и переписывает его в однострочный формат для каждого fasta >.

    Строка, начинающаяся с '>' начинает новую запись, остальные
    строки (непустые) - части её последовательности. Конец записи детектируется в двух
    случаях: когда виден новый заголовок '>' и когда файл кончился.

    Аргументы:
        input_file: открытый файл на чтение.
        output_file: открытый файл на запись.

    Возвращает:
        None.
    """
    header = None
    parts = []
    for line in input_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:                  # Запись кончилась
                save_single_line_record(output_file, header, "".join(parts))
            header = line
            parts = []
        else:
            parts.append(line)
    if header is not None:                          # Файл кончился
        save_single_line_record(output_file, header, "".join(parts))

def derive_single_line_path(input_fasta: str) -> str:
    """
    Имя выходного файла по умолчанию: входное имя + суффикс '_single_line'.

    Аргументы:
        input_fasta: путь до входного FASTA-файла.

    Возвращает:
        Путь к выходному файлу рядом со входным.
    """
    folder = os.path.dirname(input_fasta)
    stem, extension = os.path.splitext(os.path.basename(input_fasta))
    return os.path.join(folder, stem + "_single_line" + extension)


def save_single_line_record(fasta_file, header: str, sequence: str) -> None:
    """
    Записывает одну fasta-запись двумя строками: заголовок и
    последовательность (в однострочном формате).

    Аргументы:
        fasta_file: открытый файл на запись.
        header: строка-заголовок (начинается с '>').
        sequence: последовательность одной строкой.

    Возвращает:
        None.
    """
    fasta_file.write(header + "\n")
    fasta_file.write(sequence + "\n")