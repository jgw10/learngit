import os
import time
import parameters as pam
import variational_space as vs
import hamiltonian as ham
# import ground_state as gs


def compute_Aw_main(A, Uoo, Upp, ed, ep, eo, tpd, tpp, tdo, tpo):
    """
    计算La3Ni4O10的主程序
    :param A:
    :param Uoo:
    :param Upp:
    :param ed:
    :param ep:
    :param eo:
    :param tpd:
    :param tpp:
    :param tdo:
    :param tpo:
    :return:
    """
    tpd_nn_hop_dir, tpd_nn_hop_fac, tpp_nn_hop_fac = ham.set_tpd_tpp(tpd, tpp)
    Tpd = ham.create_tpd_nn_matrix(VS, tpd_nn_hop_dir, tpd_nn_hop_fac)
    gs.get_ground_state(Tpd, VS)


if __name__ == '__main__':
    t0 = time.time()
    # 创建data的文件夹
    os.makedirs('data', exist_ok=True)
    # 计算前, 清空文件夹所有文件的内容
    for filename in os.listdir('data'):
        file_path = os.path.join('data', filename)
        with open(file_path, 'w') as file:
            file.truncate(0)

    VS = vs.VariationalSpace()
    A = pam.A
    Uoo = pam.Uoos[0]
    Upp = pam.Upps[0]

    ed = pam.ed_list[4]
    ep = pam.ep_list[4]
    eo = pam.eo_list[4]

    tpd = pam.tpd_list[4]
    tpp = pam.tpp_list[4]
    tdo = pam.tdo_list[4]
    tpo = pam.tpd_list[4]

    compute_Aw_main(A, Uoo, Upp, ed, ep, eo, tpd, tpp, tdo, tpo)
    t1 = time.time()
    print('total cost time', t1-t0)
