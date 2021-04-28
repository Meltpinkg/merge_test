class Record(object):
    def __init__(self, record, idx):
        self.source = idx
        self.name = record.id
        self.type = parse_svtype(record.info['SVTYPE'])
        self.start = parse_to_int(record.pos)
        if ('INS' in record.info['SVTYPE'] or 'DEL' in record.info['SVTYPE']) and 'SVLEN' in record.info:
            self.end = abs(parse_to_int(record.info['SVLEN']))
        elif 'END' in record.info:
            self.end = parse_to_int(record.info['END'])
        else:
            try:
                self.end = parse_to_int(record.stop)
            except:
                self.end = 0
        if record.info['SVTYPE'] == 'BND' or record.info['SVTYPE'] == 'TRA':
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
        self.chrom1 = record.chrom
        if record.info['SVTYPE'] != 'TRA' and record.info['SVTYPE'] != 'BND':
            self.chrom2 = record.chrom
            if 'STRAND' in record.info:
                self.strand = record.info['STRAND']
            elif 'STRANDS' in record.info:
                self.strand = record.info['STRANDS']
            else:
                self.strand = '.'
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
        return self.name + ', start: ' + str(self.start) + ', end: ' + str(self.end) + ', strand: ' + self.strand


class ListNode(object):
    def __init__(self, id, record, pre=None, next=None):
        if record == None:
            self.variant_dict = dict()  # {id -> Record}
            self.represent = None
        else:
            self.variant_dict = dict()
            self.variant_dict[id] = record
            self.represent = record
        self.pre = pre
        self.next = next
    def add(self, id, record):
        if id in self.variant_dict:
            print('sample id already in variant dict')
        else:
            self.variant_dict[id] = record
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
    max_inspro [default = 0.7]
    return True / False
'''
def check_is_same(cur, input, max_dist, max_insratio, max_delratio):
    #if input.type == cur.type and input.strand == cur.strand:
    if input.type == cur.type:
        if abs(input.start - cur.start) < 1000:
            if input.type == 'INS':
                if max_insratio < min(input.end, cur.end) / max(input.end, cur.end) <= 1.0:
                    return True
            elif input.type == 'DEL':
                if max_delratio < min(input.end, cur.end) / max(input.end, cur.end) <= 1.0:
                    return True
            else:
                if abs(input.end - cur.end) < max_dist:
                    return True
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
def add_node(head, id, record, max_dist, max_insratio, max_delratio):
    if head.represent == None:
        node = ListNode(id, record)
        node.pre = head
        head.next = node
        return node
    cur = head
    while cur != None:
        if check_is_same(cur.represent, record, max_dist, max_insratio, max_delratio) and id not in cur.variant_dict:
            cur.add(id, record)
            return cur
        if cur.represent.start - record.start > 2 * max_dist:  # cannot merge
            break
        cur = cur.next
        #print('next')
    cur = head.pre
    while cur.represent != None:
        if check_is_same(cur.represent, record, max_dist, max_insratio, max_delratio) and id not in cur.variant_dict:
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


def init_linkedlist():
    vcf_filenames, vcfgz_filenames, chrom_set, chrom_cnt, contiginfo = pre_vcf(argv[0])

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