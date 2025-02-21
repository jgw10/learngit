import parameters as pam

# 将轨道, 自旋与数字一一对应, 用来生成每个态的数字uid
if pam.Norb == 5:
    orb_int = {'d3z2r2': 0,
               'dx2y2': 1,
               'px': 2,
               'py': 3,
               'apz': 4}
    int_orb = {value: key for key, value in orb_int.items()}
spin_int = {'up': 1, 'dn': 0}
int_spin = {value: key for key, value in spin_int.items()}


def get_unit_cell_rep(x, y, z):
    """
    确定需要计算的晶格, 根据坐标确定轨道
    :return:orbs
    """
    layer_num = pam.layer_num
    if z < 0 or z > 2*layer_num-2:
        print('z is error')
    z_Ni = range(0, 2 * layer_num, 2)
    z_Ni = tuple(z_Ni)
    if layer_num > 1:
        z_Oap = range(1, 2*layer_num-1, 2)
        z_Oap = tuple(z_Oap)
    else:
        z_Oap = ()
    if x == 0 and y == 0 and z in z_Ni:
        return pam.Ni_orbs
    elif x == 0 and y == 0 and z in z_Oap:
        return pam.Oap_orbs
    elif abs(x) % 2 == 1 and abs(y) % 2 == 0 and z in z_Ni:
        return pam.O1_orbs
    elif abs(x) % 2 == 0 and abs(y) % 2 == 1 and z in z_Ni:
        return pam.O2_orbs
    else:
        return ['NotOnSublattice']
