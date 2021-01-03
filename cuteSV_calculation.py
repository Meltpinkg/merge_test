import numpy as np


'''
    node -> dict(id -> Record)
    return (center)start, end
'''
def cal_center(node):
    #return len(node) / 2
    r_start = 0
    r_end = 0
    start_list = []
    end_list = []
    for record in node.values():
        r_start += record.start
        r_end += record.end
        start_list.append(record.start)
        end_list.append(record.end)
    r_start = r_start / len(node)
    r_end = r_end / len(node)
    r_id = 0
    start_dif = 0x3f3f3f3f
    for id in node:
        if abs(node[id].start - r_start) < start_dif:
            start_dif = node[id].start - r_start
            r_id = id
    return r_id, cal_ci(start_list), cal_ci(end_list)


def cal_ci(input_list):
    pos = int(1.96 * np.std(input_list) / len(input_list) ** 0.5)
    return "-%d,%d"%(pos, pos)
