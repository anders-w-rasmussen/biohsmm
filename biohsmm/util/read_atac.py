import numpy as np
import pysam

obs_vec = list()
length = list()
start = list()
chrom = list()

def bam_to_obs(c, s, e, filename):

    '''
    :param chrom:
    :param start:
    :param end:
    :param filename:
    :return: data array (obs)
    '''

    obs_vec.append(reads_from_sam(filename, c, s, e))
    length.append(e - s)
    start.append(s)
    chrom.append(c)

    cuts, read_length = process_data(obs_vec, length, start)
    data = np.zeros([np.sum(length), 3])
    data[:, 0] = cuts[:, 0]
    data[:, 1:] = read_length

    return data


def process_data(data, length, start):

    # Create Observation Vectors
    # Create Cuts
    obs_vec = np.zeros([np.sum(length), 1])
    tracking = 0
    for i_ in np.arange(len(data)):
        for n_ in np.arange(data[i_].shape[0]):
            r_st = np.max([np.int(data[i_][n_, 0]) - start[i_], 0])
            r_end = np.min([np.int(data[i_][n_, 1]) - start[i_], length[i_] - 1])
            obs_vec[tracking + r_st, 0] += 1
            obs_vec[tracking + r_end, 0] += 1
        tracking += length[i_]

    # Create Read Length Distributions
    obs_reads = np.zeros([np.sum(length), 2])
    tracking = 0
    for i_ in np.arange(len(data)):
        for n_ in np.arange(data[i_].shape[0]):
            r_st = np.max([np.int(data[i_][n_, 0]) - start[i_], 0])
            r_end = np.min([np.int(data[i_][n_, 1]) - start[i_], length[i_] - 1])
            obs_reads[tracking + r_st:tracking + r_end, 0] += np.ones(r_end - r_st)
            obs_reads[tracking + r_end, 1] += 1
        tracking += length[i_]

    return obs_vec, obs_reads


def reads_from_sam(samfile_name, chr, window_start, window_end):

    reads_array = []
    sam_file = pysam.AlignmentFile(samfile_name)

    for read in sam_file.fetch(chr, window_start, window_end):
        if read.flag in [83, 99, 147, 163]:

            left_tn5_start = min(read.reference_start, read.next_reference_start) - 4
            right_tn5_end = left_tn5_start + abs(read.template_length) + 8
            reads_array.append([left_tn5_start, right_tn5_end])

    return np.array(reads_array)