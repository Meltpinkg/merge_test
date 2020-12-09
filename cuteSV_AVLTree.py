#from functools import total_ordering
import re

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


#@total_ordering
class TreeNode(object):
    def __init__(self, id, record, parent=None, left=None, right=None):
        self.variant_list = [record]
        self.start = self.variant_list[0].start
        self.end = self.variant_list[0].end
        self.vis = set()
        self.vis.add(id)
        self.left = left
        self.right = right
        self.parent = parent
        self.height = 0
    def add(self, id, record):
        self.variant_list.append(record)
        self.vis.add(id)
    def to_string(self):
        string = 'Node size = ' + str(len(self.variant_list)) + ', start = ' + str(self.start) + ', end = ' + str(self.end)
        if self.left is not None:
            string += ' left: ' +  self.left.to_string()
        if self.right is not None:
            string += ' right ' + self.right.to_string()
        return string


'''
    node(原), input(待插入) -> Record
    return True, 0 -> input与record被认为可以merge
    return True, 非0 -> start很近，去子树里找有无end近的
    return False, 任意 -> 不可以merge
'''
def check_node(id, input, node, vis):
    '''
    print('====checkans')
    print(node.to_string())
    print(input.to_string())
    print('====check ans finished')
    '''
    if input.type == node.type and input.strand == node.strand:
        if abs(input.start - node.start) < 1000:
            if id not in vis and input.type == 'INS' and 0.7 < min(input.end, node.end) / max(input.end, node.end) <= 1.0:
                return True, 0
            elif id not in vis and input.type != 'INS' and abs(input.end - node.end) < 1000:
                return True, 0
            else:  # 继续尝试左右子树
                return True, input.start - node.start
        return False, input.start - node.start
    else:  # 肯定不能add
        return False, input.start - node.start


class AVLTree(object):
    def __init__(self):
        self.root = None

    def insert(self, id, record):
        if not self.root:
            self.root = TreeNode(id, record)
        else:
            self.root = self._insert(id, record, self.root)
    def _insert(self, id, record, node):  # id -> str, record -> Record, node -> TreeNode
        if node is None:
            node = TreeNode(id, record)
            return node
        check_ans = ()
        check_ans = check_node(id, record, node.variant_list[0], node.vis)
        '''
        if record.start == 10866:
            print(record.to_string())
            print(node.variant_list[0].to_string())
            print(check_ans)
        '''
        if check_ans[0] == True and check_ans[1] == 0:  # 可以merge
            node.add(id, record)
        elif check_ans[0] == True and self.second_insert(id, record, node):  # 尝试左and右子树
            pass
        else:  # 一定不能merge
            if check_ans[1] < 0:
                node.left = self._insert(id, record, node.left)
                node.left.parent = node
                if (self.height(node.left) - self.height(node.right)) == 2:
                    if record.start - node.left.start < 0:
                        node = self.singleLeftRotate(node)
                    else:
                        node = self.doubleLeftRotate(node)
            elif check_ans[1] > 0:
                node.right = self._insert(id, record, node.right)
                node.right.parent = node
                if (self.height(node.right) - self.height(node.left)) == 2:
                    if record.start - node.right.start > 0:
                        #node.to_string()
                        node = self.singleRightRotate(node)
                    else:
                        node = self.doubleRightRotate(node)

        node.height = max(self.height(node.right), self.height(node.left)) + 1
        return node
    # return 0 -> record被合并完了，否则要在子树继续insert
    def second_insert(self, id, record, node):  # 查看record能否与node子树的某个点合并，如果能返回True
        cur_node = node
        #print('start second insert')
        while abs(cur_node.start - record.start) < 1000:
            if cur_node.left is not None:
                cur_node = cur_node.left
                while cur_node.right is not None:
                    cur_node = cur_node.right
            else:
                while cur_node != node and cur_node.parent.start > cur_node.start:
                    cur_node = cur_node.parent
                if cur_node == node:
                    break
                cur_node = cur_node.parent
            check_ans = check_node(id, record, cur_node.variant_list[0], node.vis)
            '''
            if record.start == 10866:
                print(cur_node.to_string())
                print(check_ans)
            '''
            if check_ans[1] == 0:
                cur_node.add(id, record)
                return True
        cur_node = node
        while abs(cur_node.start - record.start) < 1000:
            if cur_node.right is not None:
                cur_node = cur_node.right
                while cur_node.left is not None:
                    cur_node = cur_node.left
            else:
                while cur_node != node and cur_node.parent.start < cur_node.start:
                    cur_node = cur_node.parent
                if cur_node == node:
                    break
                cur_node = cur_node.parent
            check_ans = check_node(id, record, cur_node.variant_list[0], node.vis)
            if check_ans[1] == 0:
                cur_node.add(id, record)
                return True
        return False


    def delete(self, key):
        if self.root is None:
            raise KeyError('Error,empty tree')
        else:
            self.root = self._delete(key, self.root)
    def _delete(self, key, node):
        if node is None:
            raise KeyError('Error,key not in tree')
        elif key < node.data:
            node.left = self._delete(key, node.left)
            if (self.height(node.right) - self.height(node.left)) == 2:
                if self.height(node.right.right) >= self.height(node.right.left):
                    node = self.singleRightRotate(node)
                else:
                    node = self.doubleRightRotate(node)
            node.height = max(self.height(node.left), self.height(node.right)) + 1
        elif key > node.data:
            node.right = self._delete(key, node.right)
            if (self.height(node.left) - self.height(node.right)) == 2:
                if self.height(node.left.left) >= self.height(node.left.right):
                    node = self.singleLeftRotate(node)
                else:
                    node = self.doubleLeftRotate(node)
            node.height = max(self.height(node.left), self.height(node.right)) + 1
        elif node.left and node.right:
            if node.left.height <= node.right.height:
                minNode = self._findMin(node.right)
                node.key = minNode.key
                node.right = self._delete(node.key, node.right)  
            else:
                maxNode = self._findMax(node.left)
                node.key = maxNode.key
                node.left = self._delete(node.key, node.left)
            node.height = max(self.height(node.left), self.height(node.right)) + 1
        else:
            if node.right:
                node = node.right
            else:
                node = node.left
        return node
    def inorder(self, start_node, node_list, vis_list):
        if start_node is not None:
            self.inorder(start_node.left, node_list, vis_list) 
            node_list.append(start_node.variant_list)
            vis_list.append(start_node.vis)
            self.inorder(start_node.right, node_list, vis_list)

    def find(self, key):
        if not self.root:
            return None
        else:
            return self._find(key, self.root)
    def _find(self, key, node):
        if not node:
            return None
        elif key < node.data:
            return self._find(key, node.left)
        elif key > node.data:
            return self._find(key, node.right)
        else:
            return node
    def findMin(self):
        if self.root is None:
            return None
        else:
            return self._findMin(self.root)
    def _findMin(self, node):
        if node.left:
            return self._findMin(node.left)
        else:
            return node
    def findMax(self):
        if self.root is None:
            return None
        else:
            return self._findMax(self.root)
    def _findMax(self, node):
        if node.right:
            return self._findMax(node.right)
        else:
            return node
    def height(self, node):
        if node is None:
            return -1
        else:
            return node.height
    #在node节点的左孩子k1的左子树添加了新节点，左旋转
    def singleLeftRotate(self, node):
        k1 = node.left
        node.left = k1.right
        if node.left is not None:
            node.left.parent = node
        k1.right = node
        k1.right.parent = k1
        node.height = max(self.height(node.right), self.height(node.left)) + 1
        k1.height = max(self.height(k1.left), node.height) + 1
        return k1
    #在node节点的右孩子k1的右子树添加了新节点，右旋转
    def singleRightRotate(self, node):
        k1 = node.right
        node.right = k1.left
        if node.right is not None:
            node.right.parent = node
        k1.left = node
        k1.left.parent = k1
        node.height = max(self.height(node.right), self.height(node.left)) + 1
        k1.height = max(self.height(k1.right), node.height) + 1
        return k1
    #在node节点的左孩子的右子树添加了新节点，先左后右
    def doubleRightRotate(self, node):
        node.right = self.singleLeftRotate(node.right)
        node.right.parent = node
        return self.singleRightRotate(node)
    #在node节点的右孩子的左子树添加了新节点,先右后左
    def doubleLeftRotate(self, node):
        node.left = self.singleRightRotate(node.left)
        node.left.parent = node
        return self.singleLeftRotate(node)
