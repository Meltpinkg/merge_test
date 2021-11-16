def calculate(f1, mdr, overlap):
    if f1 < 0.81 or mdr > 0.88 or overlap < 0.553:
        return -1
    


file = open('handle_del.txt')
cal = 0
with open('test_del.txt', 'r') as f:
    for line in f:
        seq = line.strip().split(',')
        if len(seq) == 1:
            seq = line.strip().split(' ')
            cal = 0
            seq.append(cal)
        else:
            f1 = seq[6]
            mdr = seq[8]
            overlap = seq[-1]

