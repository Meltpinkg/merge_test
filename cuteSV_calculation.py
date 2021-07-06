import numpy as np


'''
    node -> list(Record)
    return represent_idx, (center)start, end
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
    r_id = 0
    start_dif = 0x3f3f3f3f
    for id in range(len(node)):
        if abs(node[id].start - r_start) < start_dif:
            start_dif = abs(node[id].start - r_start)
            r_id = id
    return r_id, cal_ci(start_list), cal_ci(end_list)


def cal_ci(input_list):
    pos = int(1.96 * np.std(input_list) / len(input_list) ** 0.5)
    return "-%d,%d"%(pos, pos)

'''
    semi cluster for ins/del
    node -> dict{id -> list(Record)}
'''
def cal_can(node, threshold_gloab):
    candidates = list()
    idx = 0
    for id in node:
        for record in node[id]:
            candidates.append([id, record.start, record.end, record])
    read_tag2SortedList = sorted(candidates, key = lambda x:x[2])
    global_len = [i[2] for i in read_tag2SortedList]
    DISCRETE_THRESHOLD_LEN_CLUSTER_TEMP = threshold_gloab * np.mean(global_len)
    last_len = read_tag2SortedList[0][2]
    allele_collect = list()
    allele_collect.append([[read_tag2SortedList[0][0]], # id
                            [read_tag2SortedList[0][3]],  # record
                            [read_tag2SortedList[0][2]],  # len
                            []]) # support_num
    for i in read_tag2SortedList[1:]:
        if i[2] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_TEMP:
            allele_collect[-1][3].append(len(allele_collect[-1][0]))
            allele_collect.append([[],[],[],[]])

        allele_collect[-1][0].append(i[0])
        allele_collect[-1][1].append(i[3])
        allele_collect[-1][2].append(i[2])
        last_len = i[2]
    allele_collect[-1][3].append(len(allele_collect[-1][0]))
    allele_sort = sorted(allele_collect, key = lambda x:x[3])
    #print(allele_sort)

    vis_id = set()
    flag_id = set()
    ans = []
    for allele in allele_sort:
        while len(flag_id) != len(allele[0]):
            ans.append([])
            for i in range(len(allele[0])):
                '''print('!!')
                print(i)
                print(flag_id)
                print(allele[0][i])
                print(vis_id)'''
                if i not in flag_id:
                    if allele[0][i] not in vis_id:
                        flag_id.add(i)
                        vis_id.add(allele[0][i])
                        ans[-1].append(allele[1][i]) # record.source = id
            # next population SV
            vis_id = set()
        flag_id = set()
    #print(ans)
    return ans