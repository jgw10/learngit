1.主程序compute.py
2.检验程序是否正确:
(1).可以在parameters.py中修改层数layer_num和空穴数目hole_num, 当写完一个哈密顿矩阵时(比如Tpd), 可以令layer_num = 2, hole_num = 4, 5 or 6, 
这时相当于Ni2O9的4, 5和6空穴程序, 单独求解Tpd这一个矩阵时(gs.get_ground_state(参数1, ...), 参数1调成Tpd矩阵), 结果要相同
(2).检验在原子极限下(A, B, C参数不为0外, 其他参数均为0), 求解总的H矩阵对应的能谱是否正确

