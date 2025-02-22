import time
import numpy as np
import scipy.sparse as sps
import parameters as pam
import lattice as lat
import variational_space as vs

directions_to_vecs = {'UR': (1, 1, 0), 'UL': (-1, 1, 0), 'DL': (-1, -1, 0), 'DR': (1, -1, 0),
                      'L': (-1, 0, 0), 'R': (1, 0, 0), 'U': (0, 1, 0), 'D': (0,-1,0),
                      'T': (0, 0, 1), 'B': (0, 0, -1),
                      'L2': (-2, 0, 0), 'R2': (2, 0, 0), 'U2': (0, 2, 0), 'D2': (0, -2, 0), 'T2': (0, 0, 2), 'B2': (0, 0, -2),
                      'pzL': (-1, 0, 1), 'pzR': (1, 0, 1), 'pzU': (0, 1, 1), 'pzD': (0, -1, 1),
                      'mzL': (-1, 0, -1), 'mzR': (1, 0, -1), 'mzU': (0, 1, -1), 'mzD': (0, -1, -1)}
tpp_nn_hop_dir = ['UR', 'UL', 'DL', 'DR']


def set_tpd_tpp(tpd, tpp):
    """
    设置通过tpd, tpp跳跃的轨道和对应的方向, 以及对应的值
    :param tpd: p, d轨道的跳跃值
    :param tpp: p, p轨道的跳跃值
    :return: tpd_nn_hop_dir, tpd_nn_hop_fac, tpp_nn_hop_fac
    """
    if pam.Norb == 5 or pam.Norb == 8:
        uu = 1.
        tpd_nn_hop_dir = {'d3z2r2': ['L', 'R', 'U', 'D'],
                          'dx2y2': ['L', 'R', 'U', 'D']}

        tpd_nn_hop_fac = {('d3z2r2', 'L', 'px'): -tpd/np.sqrt(3),
                          ('d3z2r2', 'R', 'px'): tpd/np.sqrt(3),
                          ('d3z2r2', 'U', 'py'): tpd/np.sqrt(3),
                          ('d3z2r2', 'D', 'py'): -tpd/np.sqrt(3),
                          ('dx2y2', 'L', 'px'): tpd,
                          ('dx2y2', 'R', 'px'): -tpd,
                          ('dx2y2', 'U', 'py'): tpd,
                          ('dx2y2', 'D', 'py'): -tpd}
        tpp_nn_hop_fac = {('UR', 'px', 'py'): -tpp,
                          ('UL', 'px', 'py'): tpp,
                          ('DL', 'px', 'py'): -tpp,
                          ('DR', 'px', 'py'): tpp}
    else:
        print('the Norb is error')

    return tpd_nn_hop_dir, tpd_nn_hop_fac, tpp_nn_hop_fac


def create_tpd_nn_matrix(VS, tpd_nn_hop_dir, tpd_nn_hop_fac):
    """
    创建Tpd哈密顿矩阵, 只用遍历d到p轨道的跳跃,
    而p轨道到d轨道的跳跃只需将行列交换, 值不变
    :param VS:类, 含有lookup_tbl(存储要计算的态), 函数get_state_uid, get_state, get_index
    :param tpd_nn_hop_dir: 跳跃轨道和方向
    :param tpd_nn_hop_fac: 跳跃的值
    :return:
    """
    t0 = time.time()
    dim = VS.dim
    # tpd_orbs = [orbital 1, ...], tpd_keys = (orbital 1, direction, orbital 2)
    tpd_orbs = tpd_nn_hop_dir.keys()
    tpd_keys = tpd_nn_hop_fac.keys()
    data = []
    row = []
    col = []

    # 遍历整个态空间
    for row_idx in range(dim):
        state = VS.get_state(VS.lookup_tbl[row_idx])
        hole_num = len(state)

        # 其中一个空穴跳跃, 其他空穴不动
        for hole_idx in range(hole_num):
            hole = state[hole_idx]
            x, y, z, orb, s = hole
            # 根据轨道决定是否跳跃, 求出跳跃后的坐标
            if orb in tpd_orbs:
                for direction in tpd_nn_hop_dir[orb]:
                    vx, vy, vz = directions_to_vecs[direction]
                    hop_x, hop_y, hop_z = x + vx, y + vy, z + vz
                    # 由跳跃后的坐标得出跳跃后的轨道, 自旋不变
                    hop_orbs = lat.get_unit_cell_rep(hop_x, hop_y, hop_z)
                    if hop_orbs == ['NotOnSublattice']:
                        continue
                    for hop_orb in hop_orbs:
                        orb12 = (orb, direction, hop_orb)
                        if orb12 in tpd_keys:

                            # 跳跃后的空穴
                            hop_hole = (hop_x, hop_y, hop_z, hop_orb, s)
                            if hop_hole not in state:       # 检验是否满足Pauli不相容原理
                                # 将其中的一个空穴换成是跳跃后的空穴
                                hop_state = list(state)
                                hop_state[hole_idx] = hop_hole
                                hop_state, ph = vs.make_state_canonical(hop_state)
                                col_idx = VS.get_index(hop_state)
                                if col_idx is not None:
                                    value = tpd_nn_hop_fac[orb12] * ph
                                    data.extend((value, value))
                                    row.extend((row_idx, col_idx))
                                    col.extend((col_idx, row_idx))

    out = sps.coo_matrix((data, (row, col)), shape=(dim, dim))
    t1 = time.time()
    print('Tpd cost time', t1-t0)

    return out
