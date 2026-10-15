# BIOINF_scripts

Пакет биоинформатических инструментов: фильтрация FASTQ по качеству/GC/длине,
конвертация многострочного FASTA в однострочный и извлечение лучших совпадений
из отчёта BLAST.

При написании использована только стандартная библиотека Python. Все вычисления и проверки
вынесены в `helpers/`, главный скрипт содержит только оркестраторы-утилиты.

## Структура пакета

```text
BIOINF_scripts/
├── helpers/
│   ├── fastq.py
│   └── bio_files.py
├── main_script.py
└── README.md
```

## Использование утилит

```python
from main_script import (
    filter_fastq,
    convert_multiline_fasta_to_oneline,
    parse_blast_output,
)
```
Параметры каждой утилиты описаны в соответствующих разделах.

## Утилиты

### filter_fastq

Фильтрует FASTQ-файл по GC-составу, длине и среднему качеству.

```python
filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds=(0, 100),
    length_bounds=(0, 2 ** 32),
    quality_threshold=0,
) 
```
**Параметры:**

- `gc_bounds` - кортеж `(низ, верх)` или одно число (верхняя граница, нижняя = 0); по умолчанию `(0, 100)`.
- `length_bounds` - кортеж или одно число; по умолчанию `(0, 2**32)`.
- `quality_threshold` - минимальное среднее качество (phred+33); по умолчанию `0`.

Интервалы включительны. 

### convert_multiline_fasta_to_oneline

Переводит многострочный формат FASTA в однострочный (каждая последовательность одной строкой).

```python
convert_multiline_fasta_to_oneline(
    input_fasta: str,
    output_fasta=None,
) 
```

`output_fasta` необязателен: без него создаётся `<имя>_single_line.<расширение>` рядом со входом. Обработка потоковая.

### parse_blast_output

Извлекает описания лучших хитов BLAST для каждого запроса, удаляет дубликаты и сохраняет отсортированным столбцом.

```python
parse_blast_output(
    input_file: str,
    output_file: str,
)
```

## Требования

Только стандартная библиотека Python 3.

## Спасибо, что дошли до конца <3
