import pandas as pd
import time
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg

import parameters as pam
import variational_space as vs


def get_ground_state(matrix, VS):
    """
    求解矩阵的本征值和本征向量, 并对求解结果进行整理
    :param matrix: 哈密顿矩阵
    :param VS: 要求解的所有态
    :return:
    """
    t0 = time.time()
    dim = VS.dim
    Neval = pam.Neval
    vals, vecs = sps.linalg.eigsh(matrix, k=Neval, which='SA')
    print('lowest eigenvalue of H from np.linalg.eigsh = ')
    print(vals)

    # 计算不同的本征值第一次出现的索引，并存储在degen_idx中
    # 对应的简并度即为degen_idx[i+1] - degen_idx[i]
    val_num = 1
    degen_idx = [0]
    for _ in range(val_num):
        for idx in range(degen_idx[-1] + 1, pam.Neval):
            if abs(vals[idx] - vals[degen_idx[-1]]) > 1e-4:
                degen_idx.append(idx)
                break

    for i in range(val_num):
        print(f'Degeneracy of {i}th state is {degen_idx[i + 1] - degen_idx[i]}')
        print('val = ', vals[degen_idx[i]])
        weight_average = np.average(abs(vecs[:, degen_idx[i]:degen_idx[i + 1]]) ** 2, axis=1)

        # 创建MultiIndex DataFrame, 类似excel格式
        data = {'state_type': [], 'state': [], 'weight': []}
        for istate in range(dim):
            weight = weight_average[istate]
            if weight < 1e-3:       # data只存储weight > 1e-3的数据
                continue
            state = VS.get_state(VS.lookup_tbl[istate])
            state_type = vs.get_state_type(state)
            data['state_type'].append(state_type)
            data['state'].append(state)
            data['weight'].append(weight)

        df = pd.DataFrame(data)
        # 计算state_type总的weight
        df['type_weight'] = df.groupby('state_type')['weight'].transform('sum')

        # 按type_weight降序排列
        df = df.sort_values(by=['type_weight', 'weight'], ascending=[False, False])

        # 先输出state_type: type_weight, 再层次化输出state和weight
        current_type = None
        for _, row in df.iterrows():
            if row['state_type'] != current_type:
                current_type = row['state_type']
                print(f'{current_type}: {row['type_weight']}')
            state = row['state']
            weight = row['weight']
            state_string = []
            for hole in state:
                x, y, z, orb, s = hole
                hole_string = f'({x}, {y}, {z}, {orb}, {s})'
                state_string.append(hole_string)

            # 将字符串列表分成四个一组
            chunks = [state_string[i: i+4] for i in range(0, len(state_string), 4)]
            # 每个组内用', '连接, 并把组与组之间用'\n\t'连接
            state_string = '\n\t'.join([', '.join(chunk) for chunk in chunks])
            print(f'\t{state_string}\n\tweight = {weight}\n')

    t1 = time.time()
    print('gs cost time', t1-t0)
