class Record(object):
    def __init__(self, record, idx):
        self.source = idx
        self.name = record.id
        self.type = 'BND' if record.info['SVTYPE'] == 'TRA' else record.info['SVTYPE']
        self.start = parse_to_int(record.pos)
        if self.type not in ['INS','DEL','INV','DUP','BND']:
            if 'INS' in record.alts[0]:
                self.type = 'INS'
            elif 'DEL' in record.alts[0]:
                self.type = 'DEL'
        if ('INS' in self.type or 'DEL' in self.type) and 'SVLEN' in record.info:
            self.end = abs(parse_to_int(record.info['SVLEN']))
        elif 'END' in record.info:
            self.end = parse_to_int(record.info['END'])
        else:
            try:
                self.end = parse_to_int(record.stop)
            except:
                self.end = 0
        if self.type == 'BND':
            if 'CHR2' in record.info:
                chrom2 = record.info['CHR2']
            if 'END' in record.info:
                end = parse_to_int(record.info['END'])
            try:
                alt = str(record.alts[0])
                if alt[0] == ']':
                    tra_type = ']]N'
                    chrom2 = alt.split(':')[0][1:]
                    end = int(alt.split(':')[1][:-2])
                elif alt[0] == '[':
                    tra_type = '[[N'
                    chrom2 = alt.split(':')[0][1:]
                    end = int(alt.split(':')[1][:-2])
                else:
                    # print(type(alt))
                    if alt[1] == ']':
                        tra_type = 'N]]'
                        chrom2 = alt.split(':')[0][2:]
                        end = int(alt.split(':')[1][:-1])
                    else:
                        tra_type = 'N[['
                        chrom2 = alt.split(':')[0][2:]
                        end = int(alt.split(':')[1][:-1])
            except:
                print(record)
            self.chrom2 = chrom2
            self.end = end
            self.strand = tra_type
        '''
        if self.type == 'BND':
            tra_alt = str(record.alts[0])
            if tra_alt[0] == 'N':
                if tra_alt[1] == '[':
                    tra_type = 'A'
                else:
                    tra_type = 'B'
            elif tra_alt[0] == '[':
                tra_type = 'C'
            else:
                tra_type = 'D'
            if tra_alt[0] == 'N':
                tra_alt = tra_alt[2:-1]
            else:
                tra_alt = tra_alt[1:-2]
            self.chrom2 = tra_alt.split(':')[0]
            self.end = int(tra_alt.split(':')[1])
            self.strand = tra_type
        '''
        self.chrom1 = record.chrom
        self.strand = '.'
        if self.type != 'BND':
            self.chrom2 = record.chrom
            if 'STRAND' in record.info:
                self.strand = record.info['STRAND']
            elif 'STRANDS' in record.info:
                self.strand = record.info['STRANDS']
        if isinstance(self.strand, list) or isinstance(self.strand, tuple):
            self.strand = str(self.strand[0])
        self.ref = record.ref
        self.alt = record.alts  # tuple
        self.qual = record.qual  # NoneType
        if 'GT' in record.format:
            if record.samples[0]['GT'] == None or record.samples[0]['GT'][0] == None:
                self.gt = './.'
            else:
                self.gt = str(record.samples[0]['GT'][0]) + '/' + str(record.samples[0]['GT'][1])

    def to_string(self):
        return self.name + ', ' + self.type + ', start: ' + str(self.start) + ', end: ' + str(self.end) + ', strand: ' + self.strand


class ListNode(object):
    def __init__(self, id, record, pre=None, next=None):
        if record == None:
            self.variant_dict = dict()  # {id -> Record}
            self.represent = None
        else:
            self.variant_dict = dict()
            self.variant_dict[id] = [record]
            self.represent = record
        self.pre = pre
        self.next = next
    def add(self, id, record):
        if id not in self.variant_dict:
            self.variant_dict[id] = list()
        self.variant_dict[id].append(record)
        # self.variant_dict[id] = record
    def to_string(self):
        string = 'List size = ' + str(len(self.variant_dict)) + ', '
        for id in self.variant_dict:
            string += '~~~' + str(id) + ':'
            for xx in self.variant_dict[id]:
                string += xx.to_string() + '; '
        return string + '\t'


def parse_to_int(sth):
    if sth == None:
        return 0
    elif isinstance(sth, str):
        return int(sth)
    elif isinstance(sth, list):
        return parse_to_int(sth[0])
    elif isinstance(sth, tuple):
        return parse_to_int(sth[0])
    elif isinstance(sth, int):
        return sth
    else:
        return sth

'''
    cur, input -> Record
    max_dist [default = 1000]
    return True / False
'''
def check_is_same(cur, input, max_dist, max_inspro):
    #if input.type == cur.type and input.strand == cur.strand:
    if 1074800 < input.start < 1075400 and input.type == 'DEL':
        print('check')
        print(cur.to_string())
        print(input.to_string())
    if input.type == cur.type:
        if abs(input.start - cur.start) < max_dist:
            if input.type == 'INS' or input.type == 'DEL':
                if max_inspro < min(input.end, cur.end) / max(input.end, cur.end) <= 1.0:
                    return True
            else:
                if abs(input.end - cur.end) < max_dist:
                    return True
    if 1074800 < input.start < 1075400 and input.type == 'DEL':
        print('False')
    return False

'''
    head -> ListNode
    record -> Record
    向head后某处插入node(record)
'''
def insert_node(head, id, record):
    node = ListNode(id, record)
    if head.represent.start > record.start:
        node.next = head
        node.pre = head.pre
        head.pre.next = node
        head.pre = node
        return node
    while head.next != None:
        if head.next.represent.start > record.start:
            node.next = head.next
            node.pre = head
            head.next.pre = node
            head.next = node
            return node
        head = head.next
    head.next = node
    node.pre = head
    return node

'''
    head -> ListNode 插入的起始节点
    id -> sample_id
    record -> Record
    在head附近某处添加(add/insert) ListNode(id, record)
'''
def add_node(head, id, record, max_dist, max_inspro):
    if head.represent == None:
        node = ListNode(id, record)
        node.pre = head
        head.next = node
        return node
    cur = head
    while cur != None:
        if check_is_same(cur.represent, record, max_dist, max_inspro) and id not in cur.variant_dict:
            cur.add(id, record)
            return cur
        if cur.represent.start - record.start > 2 * max_dist:  # cannot merge
            break
        cur = cur.next
        #print('next')
    cur = head.pre
    while cur.represent != None:
        if check_is_same(cur.represent, record, max_dist, max_inspro) and id not in cur.variant_dict:
            cur.add(id, record)
            return cur
        if record.start - cur.represent.start > 2 * max_dist:
            break
        cur = cur.pre
        #print('pre')
    cur = cur.next
    return insert_node(cur, id, record)

def add_node_indel(head, id, record, max_dist, max_inspro):
    if head.represent == None:
        node = ListNode(id, record)
        node.pre = head
        head.next = node
        return node
    cur = head
    while cur != None:
        if check_is_same(cur.represent, record, max_dist, max_inspro):
            cur.add(id, record)
            return cur
        if cur.represent.start - record.start > 2 * max_dist:  # cannot merge
            break
        cur = cur.next
        #print('next')
    cur = head.pre
    while cur.represent != None:
        if check_is_same(cur.represent, record, max_dist, max_inspro):
            cur.add(id, record)
            return cur
        if record.start - cur.represent.start > 2 * max_dist:
            break
        cur = cur.pre
        #print('pre')
    cur = cur.next
    return insert_node(cur, id, record)

def print_list(head):
    print('Linked List:')
    while head != None:
        print(head.to_string())
        head = head.next


def parse_svtype(sv_type):
    if 'DEL' in sv_type:
        return 'DEL'
    if 'INS' in sv_type:
        return 'INS'
    if 'INV' in sv_type:
        return 'INV'
    if 'DUP' in sv_type:
        return 'DUP'
    if 'BND' in sv_type or 'TRA' in sv_type:
        return 'BND'

'''
if max(input.end, cur.end) <= 100:
    if abs(input.end - cur.end) < 70:
        return True
if max(input.end, cur.end) <= 300:
    if abs(input.end - cur.end) < 150:
        return True
if max(input.end, cur.end) <= 500:
    if abs(input.end - cur.end) < 300:
        return True
if max(input.end, cur.end) <= 1000:
    if abs(input.end - cur.end) < 600:
        return True
'''
'''
if max(input.end, cur.end) <= 100:
    if 0.3 < min(input.end, cur.end) / max(input.end, cur.end) <= 1.0:
        return True
if max(input.end, cur.end) <= 300:
    if 0.5 < min(input.end, cur.end) / max(input.end, cur.end) <= 1.0:
        return True
if max(input.end, cur.end) <= 1000:
    if 0.6 < min(input.end, cur.end) / max(input.end, cur.end) <= 1.0:
        return True
'''