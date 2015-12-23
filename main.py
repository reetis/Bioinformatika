import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO


def peakIndexes(y, thres=0.3, min_dist=1):
    thres *= np.max(y) - np.min(y)

    # find the peaks by using the first order difference
    dy = np.diff(y)
    peaks = np.where((np.hstack([dy, 0.]) < 0.)
                     & (np.hstack([0., dy]) > 0.)
                     & (y > thres))[0]

    if peaks.size > 1 and min_dist > 1:
        highest = peaks[np.argsort(y[peaks])][::-1]
        rem = np.ones(y.size, dtype=bool)
        rem[peaks] = False

        for peak in highest:
            if not rem[peak]:
                sl = slice(max(0, peak - min_dist), peak + min_dist + 1)
                rem[sl] = True
                rem[peak] = False

        peaks = np.arange(y.size)[~rem]

    return peaks


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


def get_records():
    global values
    records = SeqIO.parse(FILENAME, get_fastq_format()["format"])
    values = []
    for record in records:
        seq = record.seq
        gc_count = seq.count("G") + seq.count("C")
        gc_ratio = gc_count / len(seq) * 100
        values.append((record, gc_ratio))
    return values


def get_peak_ranges(bins, bin_values):
    peaks = peakIndexes(bin_values, min_dist=15)
    ranges = map(lambda x: (bins[x], bins[x + 1]), peaks)
    return ranges


if __name__ == '__main__':
    values = get_records()

    n, bins, _ = plt.hist([value[1] for value in values], bins=100, range=(0, 100))
    plt.title("G/C frequency")
    plt.xlabel("Probability")
    plt.ylabel("Count")
    # plt.show()

    ranges = get_peak_ranges(bins, n)
    for r in ranges:
        print(r)
