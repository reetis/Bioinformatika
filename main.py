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


def smooth(x, window_len=7, window='blackman'):
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]  # print(len(s))
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='valid')
    return y[(window_len / 2 - 1):-(window_len / 2)]


if __name__ == '__main__':
    records = SeqIO.parse(FILENAME, get_fastq_format()["format"])
    values = []
    for record in records:
        seq = record.seq
        gc_count = seq.count("G") + seq.count("C")
        gc_ratio = gc_count / len(seq) * 100
        values.append((record, gc_ratio))

    plt.subplot(211)
    n, bins, _ = plt.hist([value[1] for value in values], bins=100, range=(0, 100))
    plt.title("G/C frequency")
    plt.xlabel("Probability")
    plt.ylabel("Count")
    plt.subplot(212)
    n_25, _, _ = plt.hist([value[1] for value in values], bins=25, range=(0, 100))
    plt.title("G/C frequency smoothed")
    plt.xlabel("Probability")
    plt.ylabel("Count")
    plt.show()

    coarse_peaks = argrelextrema(n_25, np.greater)
    print(coarse_peaks[0])
    peaks = [peak * 4 + np.argmax(n[peak * 4:peak * 4 + 4]) for peak in coarse_peaks[0]]
    ranges = map(lambda x: (bins[x], bins[x + 1]), peaks)
    print(peaks)
    for r in ranges:
        print(r)
