class Record(object):
    def __init__(self, record, idx):
        self.source = idx
        self.name = record.id
        self.type = record.info['SVTYPE']
        self.start = int(record.pos)
        if record.info['SVTYPE'] == 'INS' and 'SVLEN' in record.info:
            self.end = int(record.info['SVLEN'])
        else:
            try:
                self.end = int(record.info['END'])
            except:
                try:
                    self.end = record.stop
                except:
                    pass     
        if record.info['SVTYPE'] == 'BND' or record.info['SVTYPE'] == 'TRA':
            tra_alt = str(record.alts[0])  # 若格式不正确可能下标越界
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
            if 'STRANDS' in record.info:
                self.strand = record.info['STRANDS']
            else:
                self.strand = '.'
        self.ref = record.ref
        self.alt = record.alts  # tuple
        self.qual = record.qual  # NoneType

    def to_string(self):
        return self.name + ', start: ' + str(self.start) + ', end: ' + str(self.end) + ', strand: ' + self.strand


class ListNode(object):
    def __init__(self, id, record, pre=None, next=None):
        self.variant_list = [record]
        self.start = record.start
        self.end = record.end
        self.vis = set()
        self.vis.add(id)
        self.pre = pre
        self.next = next
    def add(self, id, record):
        self.variant_list.append(record)
        self.vis.add(id)
    def to_string(self):
        string = 'List size = ' + str(len(self.variant_list)) + ', '
        for rec in variant_list:
            string += rec.to_string() + '; '
        return string
        

'''
    cur, input -> Record
    max_dist [default = 1000]
    max_inspro [default = 0.7]
    return True / False
'''
def check_is_same(cur, input, id, vis, max_dist, max_inspro):
    if input.type == cur.type and input.strand == cur.strand and id not in vis:
        if abs(input.start - cur.start) < 1000:
            if input.type == 'INS' and max_inspro < min(input.end, cur.end) / max(input.end, cur.end) <= 1.0:
                return True
            elif input.type != 'INS' and abs(input.end - cur.end) < max_dist:
                return True
    return False


'''
    head, node -> ListNode
    向head后某处插入node
'''
def insert_node(head, node):
    if head.start > record.start:
        if head.pre == null:
            node.next = head
            head.pre = node
            return node
        else:
            node.next = head
            node.pre = head.pre
            head.pre.next = node
            head.pre = node
            return node
    while head.next != null:
        if head.next.start > record.start:
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
    if head == null:
        head = ListNode(id, record)
        return head
    cur = head
    while cur != null:
        if check_is_same(cur.variant_list[0], record, id, cur.vis, max_dist, max_inspro):
            cur.add(id, record)
            return cur
        if cur.start - record.start > 2 * max_dist:  # cannot merge
            break            
        cur = cur.next
    cur = head.pre
    while cur != null:
        if check_is_same(cur, record):
            cur.add(id, record)
            return cur
        if record.start - cur.start > 2 * max_dist:
            break
        cur = cur.pre
    return insert_node(cur, ListNode(id, record))


def print_list(head):
    print('Linked List:')
    while head != null:
        print(ListNode.to_string())



def init_linkedlist():
    vcf_filenames, vcfgz_filenames, chrom_set, chrom_cnt, contiginfo = pre_vcf(argv[0])


if __name__ == '__main__':
