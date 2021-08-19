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
    ans -> [list(Record)] -> ([clu1, clu2, clu3])
'''
def cal_can(node, threshold_gloab, flag):
    if flag == 1:
        print(node)
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
    if flag == 1:
        print(allele_sort)
        # [[[2, 1, 2, 1, 0], [r1, r2, r3, r4, r5], [48, 50, 58, 79, 89], [5]]]

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

'''
    semi cluster for ins/del
    node -> dict{id -> list(Record)}
    candidates -> [list(Record)] -> ([clu1, clu2, clu3, ...])
    clu -> list(Record)
'''
def cal_can_greedy(node, max_dist, max_inspro, debug):
    if debug:
        print('start cal can:')
        for id in node:
            print(str(id) + ': ')
            for r in node[id]:
                print(r.to_string())
    candidates = list()  # clusters
    can_len = 0
    can_idx = 0
    for id in node:
        if len(node[id]) > can_len:
            can_idx = id
    if len(node[can_idx]) == 1:
        candidates.append([])
        for id in node:
            candidates[0].append(node[id][0])
        return candidates
    # multi clusters
    for record in node[can_idx]:
        candidates.append([{can_idx}, record.start, record.end, [record], 1])
    for id in node:
        if id == can_idx:
            continue
        for record in node[id]:
            st_bias = [0 for i in range(len(candidates))]
            ed_bias = [0 for i in range(len(candidates))]
            merge_flag = False
            for i in range(len(candidates)):
                if id not in candidates[i][0]:
                    if abs(candidates[i][1] - record.start) < max_dist and min(candidates[i][2], record.end) / max(candidates[i][2], record.end) > max_inspro:
                        merge_flag = True
                    st_bias[i] = abs(candidates[i][1] - record.start)
                    ed_bias[i] = abs(candidates[i][2] - record.end)
            if merge_flag is False:
                candidates.append([{id}, record.start, record.end, [record], 1])
            st_mean = np.mean(st_bias)
            ed_mean = np.mean(ed_bias)
            if ed_mean == 0 or st_mean == 0:
                ratio = 1
            else:
                ratio = st_mean / ed_mean
            min_bias = 0x3f3f3f3f
            min_bias_idx = 0
            for i in range(len(st_bias)):
                if st_bias[i] / ratio + ed_bias[i] < min_bias:
                    min_bias = st_bias[i] / ratio + ed_bias[i]
                    min_bias_idx = i
            candidates[min_bias_idx][0].add(id)
            candidates[min_bias_idx][1] = (candidates[min_bias_idx][4] * candidates[min_bias_idx][1] + record.start) / (candidates[min_bias_idx][4] + 1)
            candidates[min_bias_idx][2] = (candidates[min_bias_idx][4] * candidates[min_bias_idx][2] + record.end) / (candidates[min_bias_idx][4] + 1)
            candidates[min_bias_idx][3].append(record)
            candidates[min_bias_idx][4] += 1
    # re-merge candidates
    ans_can = list()
    len_can = len(candidates)
    vis_id = set()
    for i in range(len_can):
        if i in vis_id:
            continue
        vis_id.add(i)
        id_temp = candidates[i][0]
        st_temp = candidates[i][1]
        ed_temp = candidates[i][2]
        record_temp = candidates[i][3]
        for j in range(i, len_can, 1):
            if j not in vis_id:
                if len(id_temp & candidates[j][0]) == 0 and abs(candidates[i][1] - candidates[j][1]) < max_dist and min(candidates[i][2], candidates[j][2]) / max(candidates[i][2], candidates[j][2]) > max_inspro:
                    id_temp = id_temp | candidates[j][0]
                    st_temp = (st_temp + candidates[j][1]) / 2
                    ed_temp = (ed_temp + candidates[j][2]) / 2
                    record_temp = record_temp + candidates[j][3]
                    vis_id.add(j)
        ans_can.append(record_temp)
    if debug:
        for i in range(len(ans_can)):
            print(str(i + 1) + ':')
            for r in ans_can[i]:
                print(r.to_string())
    return ans_can
