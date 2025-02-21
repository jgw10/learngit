import parameters as pam
import lattice as lat
import time
from itertools import combinations
import bisect


def get_hole_uid(hole):
    """
    将空穴信息以进制的规则转化为数字信息
    :param hole: hole = ('x', 'y', 'z', 'orb', 's'), orb轨道, s自旋
    :return: hole_uid, 空穴所对应的数字
    """
    x, y, z, orb, s = hole

    # 依次将x, y, z, orb, s的最大个数求出来, 并作为进制
    b_x = 2 * pam.Mc + 1
    b_y = 2 * pam.Mc + 1
    b_z = 2 * pam.layer_num - 1
    b_orb = pam.Norb
    b_s = 2

    # 将x, y, z, orb, s依次转为数字
    # 保证i_x, i_y, i_z一定是非负数
    i_x = x + pam.Mc
    i_y = y + pam.Mc
    i_z = z
    # 轨道和自旋转为数字
    i_orb = lat.orb_int[orb]
    i_s = lat.spin_int[s]

    # 按将这些数, 按进制规则转为一个大数, 从低位到高位分别是s, orb, z, y, x
    hole_uid = 0
    b_list = [b_s, b_orb, b_z, b_y, b_x]
    i_list = [i_s, i_orb, i_z, i_y, i_x]
    b = 1
    for idx in range(len(i_list)):
        i = i_list[idx]
        hole_uid += i * b
        b *= b_list[idx]
    assert hole == get_hole(hole_uid), 'check hole and get_hole are tuple'

    return hole_uid


def get_hole(hole_uid):
    """
    将空穴的数字信息转化为空穴信息
    :param hole_uid: 空穴的数字信息
    :return: hole = ('x', 'y', 'z', 'orb', 's'), orb轨道, s自旋
    """
    # 依次将x, y, z, orb, s的最大个数求出来, 并作为进制
    b_x = 2 * pam.Mc + 1
    b_y = 2 * pam.Mc + 1
    b_z = 2 * pam.layer_num - 1
    b_orb = pam.Norb
    b_s = 2

    # 将大数依据进制规则, 提取每个进制上的数
    # 依次提取i_s, i_orb; i_z, i_y, i_x
    i_s = hole_uid % b_s
    hole_uid //= b_s
    i_orb = hole_uid % b_orb
    hole_uid //= b_orb
    i_z = hole_uid % b_z
    hole_uid //= b_z
    i_y = hole_uid % b_y
    hole_uid //= b_y
    i_x = hole_uid % b_x

    # 将这些数转为对应的空穴信息 hole = (x, y, z, orb, s)
    # i_x, i_y减去pam.Mc
    x = i_x - pam.Mc
    y = i_y - pam.Mc
    z = i_z
    # 得到轨道和自旋信息
    orb = lat.int_orb[i_orb]
    s = lat.int_spin[i_s]
    hole = (x, y, z, orb, s)

    return hole


def count_inversion(state):
    """
    计算一个态, 需要经过多少次交换才能得到规范化的顺序
    :param state: state是一个嵌套元组, state = (hole1, hole2, ...)
    :return: inversion, 交换次数
    """
    # 将state中的每个空穴转为hole_uid
    uid_state = map(get_hole_uid, state)  # 注意map是惰性求解器
    uid_state = tuple(uid_state)

    inversion = 0
    hole_num = len(uid_state)
    for i in range(1, hole_num):
        behind_uid = uid_state[i]
        for front_uid in uid_state[:i]:
            if front_uid > behind_uid:
                inversion += 1
    return inversion


def make_state_canonical(state):
    """
    将态中的每个空穴按照hole_uid升序排列,
    :param state: state = (hole1, hole2, ...)
    :return: canonical_state and phase
    canonical_state, 按照hole_uid升序后的态
    phase, 每交换一次顺序相位乘以-1
    """
    inversion = count_inversion(state)
    phase = 1.0 if inversion % 2 == 0 else -1.0
    canonical_state = sorted(state, key=get_hole_uid)
    canonical_state = tuple(canonical_state)
    return canonical_state, phase


def get_state_type(state):
    """
    按照Ni, 层内O, 层间O每层的数量给态分类
    :param state: state = ((x1, y1, z1, orb1, s1), ...)
    :return:state_type
    """
    layer_num = pam.layer_num
    # 统计每一层Ni, 层内O, 层间O的数量
    Ni_num = [0] * layer_num
    O_num = [0] * layer_num
    Oap_num = [0] * layer_num
    for hole in state:
        x, y, z, orb, _ = hole
        # Ni
        if orb in pam.Ni_orbs:      # 这里尽量用轨道来判断, 方便以后修改lattice.py
            idx = int(z / 2)
            Ni_num[idx] += 1
        # 层内O
        if orb in pam.O_orbs:
            idx = int(z / 2)
            O_num[idx] += 1
        # 层间O
        if orb in pam.Oap_orbs:
            idx = int((z - 1) / 2)
            Oap_num[idx] += 1

    # 对称等价的数量分布
    Ni_num_sym = Ni_num[::-1]
    O_num_sym = O_num[::-1]
    Oap_num_sym = Oap_num[::-1]

    # 选择其中一种等价的态类型
    # 先判断Ni, 选择前面Ni数量多的, 再依次判断层内O和层间O, 选择前面少的
    if Ni_num_sym > Ni_num:
        Ni_num, O_num, Oap_num = Ni_num_sym, O_num_sym, Oap_num_sym
    elif Ni_num_sym == Ni_num:
        if [O_num_sym, Oap_num_sym] < [O_num, Oap_num_sym]:
            O_num, Oap_num = O_num_sym, Oap_num_sym

    # 根据每一层Ni, 层内O和层间O的空穴数量, 生成态的类型
    state_type = []
    for idx, Ni in enumerate(Ni_num):
        On = O_num[idx]

        # 先根据Ni和层内O的数量生成
        if Ni != 0 or On != 0:
            dL = []
            if Ni != 0:
                dL.append(f'd{10-Ni}')
            if On == 1:
                dL.append('L')
            elif On != 0:
                dL.append(f'L{On}')
            dL = ''.join(dL)
            state_type.append(dL)

        # 再根据层间O的数量生成
        if idx < len(Oap_num):
            Oap = Oap_num[idx]
            if Oap != 0:
                if Oap == 1:
                    Oap_type = 'O'
                else:
                    Oap_type = f'O{Oap}'
                state_type.append(Oap_type)
    state_type = '-'.join(state_type)

    return state_type


def get_atomic_energy(state, A, Upp, Uoo, ep, eo):
    """
    大致计算该态在原子极限下的能量(以d8为基态, 设d8的能量为0)
    :param state:state = ((x1, y1, z1, orb1, s1), ...)
    :param A: 描述d轨道相互作用的一个参数
    :param Upp: 描述p轨道相互作用的参数
    :param Uoo: 描述apz轨道(层间O的pz)相互作用的参数
    :param ep: p轨道上的在位能
    :param eo: apz轨道上的在位能
    :return:energy
    """
    energy = 0.
    Ni_num = [0] * pam.layer_num
    O_num = {}
    Oap_num = {}
    for x, y, z, orb, _ in state:
        # 计算p, apz轨道上的在位能
        if orb in pam.O_orbs:
            energy += ep
        elif orb in pam.Oap_orbs:
            energy += eo

        # 统计在相同位置上的个数
        # 并根据在Ni, 层内O, 层间O依次存入Ni_num, O_num, Oap_num
        if orb in pam.Ni_orbs:
            idx = int(z / 2)
            Ni_num[idx] += 1
        if orb in pam.O_orbs:
            if (x, y, z) in O_num:
                O_num[(x, y, z)] += 1
            else:
                O_num[(x, y, z)] = 1
        if orb in pam.Oap_orbs:
            if (x, y, z) in Oap_num:
                Oap_num[(x, y, z)] += 1
            else:
                Oap_num[(x, y, z)] = 1

    # 只保留大于1的值
    O_num = [num for num in O_num.values() if num > 1]
    Oap_num = [num for num in Oap_num.values() if num > 1]

    # dn的能量, 比d8高A/2 * abs(n - 8)
    for num in Ni_num:
        energy += A/2 * abs(num-2)
    for num in O_num:
        double_num = num * (num - 1) / 2
        energy += Upp * double_num
    for num in Oap_num:
        double_num = num * (num - 1) / 2
        energy += Uoo * double_num

    return energy


class VariationalSpace:
    def __init__(self):
        self.lookup_tbl = self.create_lookup_tbl()
        self.dim = len(self.lookup_tbl)
        print(f'VS.dim = {self.dim}')

    def create_lookup_tbl(self):
        """
        找出所有可能的态，并根据能量的大小，砍去一部分态后存储在列表中
        :return: lookup_tbl
        """
        t0 = time.time()
        Mc = pam.Mc
        layer_num = pam.layer_num
        hole_num = pam.hole_num
        # 生成单个空穴的所有可能组合
        hole_list = []
        for x in range(-Mc, Mc + 1):
            B = Mc - abs(x)
            for y in range(-B, B + 1):
                for z in range(2 * layer_num - 1):
                    orbs = lat.get_unit_cell_rep(x, y, z)
                    if orbs == ['NotOnSublattice']:
                        continue
                    for orb in orbs:
                        for s in ['up', 'dn']:
                            hole_list.append((x, y, z, orb, s))

        # 生成所有可能组合的态, 并分类型, 根据能量范围得到所需要的态
        max_energy = 100
        state_energy = []
        type_num = {}
        lookup_tbl = []
        for state in combinations(hole_list, hole_num):
            canonical_state, _ = make_state_canonical(state)        # 将态按一定顺序排列, 避免重复
            state_type = get_state_type(canonical_state)
            energy = get_atomic_energy(canonical_state, 6., 4., 4., 2.9, 3.24)
            if (state_type, energy) not in state_energy:
                state_energy.append((state_type, energy))
                type_num[(state_type, energy)] = 1
            else:
                type_num[(state_type, energy)] += 1
            if energy < max_energy:
                uid = self.get_state_uid(canonical_state)
                lookup_tbl.append(uid)

        lookup_tbl.sort()       # 一定要有这一步, 这会影响get_index函数
        # 输出所有的类型
        state_energy.sort(key=lambda st_e: st_e[1])
        with open('./data/state_energy', 'w') as file:
            for state_type, energy in state_energy:
                num = type_num[(state_type, energy)]
                file.write(f'{state_type}: energy = {energy}, num = {num}\n')
        t1 = time.time()
        print('VS cost time', t1 - t0)
        return lookup_tbl

    def get_state_uid(self, state):
        """
        将态信息转化为数字信息
        :param state: state = (hole1, hole2, ....)
        :return: uid, 态对应的数字
        """
        Mc = pam.Mc
        # 计算存储一个空穴信息需要多大的进制, 并记录在b_hole
        b_x = 2 * Mc + 1
        b_y = 2 * Mc + 1
        b_z = 2 * pam.layer_num - 1
        b_orb = pam.Norb
        b_hole = b_x * b_y * b_z * b_orb * 2

        # 将每个空穴数字, 按b_hole进制转成一个大数
        uid = 0
        for idx, hole in enumerate(state):
            i_hole = get_hole_uid(hole)
            uid += i_hole * (b_hole ** idx)
        assert state == self.get_state(uid), 'check state and get_state are tuple'

        return uid

    def get_state(self, uid):
        """
        将态的数字信息转为态信息
        :param uid: 态的数字信息
        :return: state = (hole1, hole2, ....)
        """
        Mc = pam.Mc
        # 计算存储一个空穴信息需要多大的进制, 并记录在b_hole
        b_x = 2 * Mc + 1
        b_y = 2 * Mc + 1
        b_z = 2 * pam.layer_num - 1
        b_orb = pam.Norb
        b_hole = b_x * b_y * b_z * b_orb * 2

        # 将大数uid按照b_hole进制, 提取每个进制上的数
        state = []
        while uid:
            hole_uid = uid % b_hole
            hole = get_hole(hole_uid)
            state.append(hole)
            uid //= b_hole
        state = tuple(state)

        return state

    def get_index(self, state):
        """
        根据lookup_tbl, 找到state对应的索引
        :param state: state = ((x1, y1, z1, orb1, s1)...)
        :return: index, state对应在lookup_tbl中的索引
        """
        uid = self.get_state_uid(state)
        index = bisect.bisect_left(self.lookup_tbl, uid)

        # 判断索引是否超出lookup_tbl
        if index < self.dim:
            if self.lookup_tbl[index] == uid:
                return index
            else:
                return None
        else:
            return None
