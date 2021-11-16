import math
def check_two_sv(start1, end1, start2, end2):
    max_dist = 1000
    max_length_bias = 500
    max_inspro = 0.5
    #mean = (end1 + end2) / 2
    mean = max(end1, end2)
    #'''
    if abs(start1 - start2) < max_dist:
        if mean < 100 and min(end1, end2) / max(end1, end2) > 0.3:
            return True
        if 100 <= mean < 500 and abs(end1 - end2) < 300:
            return True
        if 500 <= mean < 1000 and abs(end1 - end2) < 500:
            return True
        if mean >= 1000 and abs(end1 - end2) < 1000:
            return True
    elif abs(start1 - start2) < max_dist * 1.5:
        if min(end1, end2) / max(end1, end2) > 0.9:
            return True
    '''
    if abs(start1 - start2) < max_dist and min(end1, end2) / max(end1, end2) > max_inspro:
        return True
    '''
    return False

def bfs(k, linky, dbx, dby, dis, n):
    py = 0
    yy = 0
    linky[py] = k
    slack = [0x3f3f3f3f for i in range(n + 1)]
    pre = [0 for i in range(n + 1)]
    vy = [0 for i in range(n + 1)]
    while True:
        px = linky[py]
        delta = 0x3f3f3f3f
        vy[py] = 1
        for i in range(1, n + 1):
            if not vy[i]:
                if dbx[px] + dby[i] - dis[px][i] < slack[i]:
                    slack[i] = dbx[px] + dby[i] - dis[px][i]
                    pre[i] = py
                if slack[i] < delta:
                    delta = slack[i]
                    yy = i
        for i in range(n):
            if vy[i]:
                dbx[linky[i]] -= delta
                dby[i] += delta
            else:
                slack[i] -= delta
        py = yy
        if linky[py] == 0:
            break
    while py:
        linky[py] = linky[pre[py]]
        py = pre[py]


def dfs(x, n, visx, visy, lx, ly, dis, linker, slack):
    visx[x] = 1
    for y in range(n):
        if visy[y]:
            continue
        tmp = lx[x] + ly[y] - dis[x][y]
        #if tmp == 0:
        if tmp < 1e-6:
            visy[y] = 1
            if linker[y] == -1 or dfs(linker[y], n, visx, visy, lx, ly, dis, linker, slack):
                linker[y] = x
                return 1
        elif slack[y] > tmp:
            slack[y] = tmp
    return 0


def km(n, dis):
    linker = [-1 for i in range(n)]
    ly = [0 for i in range(n)]
    lx = [-0x3f3f3f3f for i in range(n)]
    for i in range(n):
        for j in range(n):
            if dis[i][j] > lx[i]:
                lx[i] = dis[i][j]
    #print(lx)
    for x in range(n):
        slack = [0x3f3f3f3f for i in range(n)]
        while True:
            visx = [0 for i in range(n)]
            visy = [0 for i in range(n)]
            if dfs(x, n, visx, visy, lx, ly, dis, linker, slack):
                break
            d = 0x3f3f3f3f
            for i in range(n):
                if (not visy[i]) and d > slack[i]:
                    d = slack[i]
            for i in range(n):
                if visx[i]:
                    lx[i] -= d
            for i in range(n):
                if visy[i]:
                    ly[i] += d
                else:
                    slack[i] -= d
    for i in range(n):
        print('%d %d: %d'%(linker[i], i, dis[linker[i]][i]))


def main():
    boys = [(4371, 124)]
    girls = [(3395, 68), (3581, 340), (4702, 209)] # candidates
    #boys = [(363, 509)]
    #girls = [(887, 51), (399, 509), (1012, 75)]
    dis = []
    for boy in boys:
        dis.append([])
        for girl in girls:
            if check_two_sv(boy[0], boy[1], girl[0], girl[1]):
                #distance = - math.log(abs(boy[0] - girl[0] + 1)) - math.log(abs(boy[1] - girl[1] + 1))
                #distance = - abs(boy[0] - girl[0]) - abs(boy[1] - girl[1])
                #distance = -int( math.log(abs(boy[0] - girl[0]) + 1) * 100 + math.log(abs(boy[0] + boy[1] + girl[0] - girl[1]) + 1) * 1000)
                distance = -int( math.log(abs(boy[0] - girl[0]) + 1) * 100 + math.log(abs(max(boy[1], girl[1]) / min(boy[1], girl[1]) * 10)) * 1000 )  # 可以处理 205 DEL
                #print(math.log(abs(boy[0] - girl[0]) + 1))
                print(math.log(abs(boy[0] - girl[0]) + 1))
                print(math.log(abs(max(boy[1], girl[1]) / min(boy[1], girl[1]) * 10)))
                #print(math.log(abs(boy[1] / girl[1])))
                #distance = abs(candidate[1] - record.start) + abs(candidate[2] - record.end)
                #distance = -math.log(abs(candidate[1] - record.start) + 1) - math.log(abs(candidate[2] - record.end) + 1)
                #distance = -math.log(abs(candidate[2] - record.end) + 1)
            else:
                distance = -0x3f3f3f3f
            dis[-1].append(distance)
    n = len(girls)
    for i in range(len(boys), len(girls), 1):
        dis.append([-0x3f3f3f3f for j in range(n)])
    #dis = [[-6.56244409369372, -9.092907275084087], [-1.3862943611198906, -1061109567]]
    #dis = [[-6.5, -9.0], [-1.3, -1061109567.0]]
    print(dis)
    km(n, dis)


main()