import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from scipy.signal import argrelextrema

FILENAME = "reads_for_analysis.fastq"
FASTQ_FORMATS = [
    {"name": "Sanger", "min": 33, "max": 73, "format": "fastq-sanger"},
    {"name": "Solexa", "min": 59, "max": 104, "format": "fastq-solexa"},
    {"name": "Illumina 1.3+", "min": 64, "max": 104, "format": "fastq-illumina"},
    {"name": "Illumina 1.5+", "min": 66, "max": 104, "format": "fastq-illumina"},
    {"name": "Illumina 1.8+", "min": 33, "max": 74, "format": "fastq-sanger"}
]


def get_fastq_format():
    q_formats = FASTQ_FORMATS
    with open(FILENAME, "r") as file:
        read_next_line = False
        for line in file:
            if read_next_line:
                ascii_values = [ord(char) for char in line.strip("\r\n")]
                q_max = max(ascii_values)
                q_min = min(ascii_values)
                q_formats = get_suitable_formats(q_min, q_max, q_formats)
                if len(q_formats) <= 1:
                    break
                read_next_line = False
            if line[0] == "+":
                read_next_line = True
    if len(q_formats) == 0:
        return None
    else:
        return q_formats[0]


def get_suitable_formats(min_val, max_val, formats):
    new_formats = []
    for format_desc in formats:
        if format_desc["min"] <= min_val <= format_desc["max"] and format_desc["min"] <= max_val <= format_desc["max"]:
            new_formats.append(format_desc)
    return new_formats


if __name__ == '__main__':
    records = SeqIO.parse(FILENAME, get_fastq_format()["format"])
    values = []
    for record in records:
        seq = record.seq
        gc_count = seq.count("G") + seq.count("C")
        gc_ratio = gc_count / len(seq)
        values.append((record, gc_ratio))

    n, bins, _ = plt.hist([value[1] for value in values], bins=25, range=(0, 1))
    plt.title("G/C frequency")
    plt.xlabel("Probability")
    plt.ylabel("Count")
    plt.show()

    peaks = argrelextrema(n, np.greater)
    print(peaks)
