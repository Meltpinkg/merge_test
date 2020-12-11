import numpy as np


'''
    node -> list(Record)
    return (center)start, end
'''
def cal_center(node):
    #return len(node) / 2
    r_start = 0
    r_end = 0
    start_list = []
    end_list = []
    for record in node:
        r_start += record.start
        r_end += record.end
        start_list.append(record.start)
        end_list.append(record.end)
    r_start = r_start / len(node)
    r_end = r_end / len(node)
    r_idx = 0
    start_dif = 0x3f3f3f3f
    for idx in range(len(node)):
        if abs(node[idx].start - r_start) < start_dif:
            start_dif = node[idx].start - r_start
            r_idx = idx
    return r_idx, cal_ci(start_list), cal_ci(end_list)


def cal_ci(input_list):
    pos = int(1.96 * np.std(input_list) / len(input_list) ** 0.5)
    return "-%d,%d"%(pos, pos)
