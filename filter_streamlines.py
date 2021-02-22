from generate_data import *
import numpy as np

class filter_streamlines(object):
    def __init__(self):
        '''Initialize member variables.
                '''
        self.all_cor = None
        self.all_scale = None
        self.data_num = 0
        self.line_num = 0
        self.lines = None
        self.filter_line_bool = None
        self.rft = None
        self.idx = 0

    def __clear_all_data(self):
        '''Clear data for all class member variables.
        '''
        self.__init__()

    def __get_data(self, lines, all_cor):
        filter_cor_id = []
        filter_line_id = self.filter_line_bool.copy()
        for i in range(lines.shape[0]):
            if i in filter_line_id:
                filter_cor_id.append(lines[i])

    def __get_rft_f(self, id):
        # 获取原ReasonForTermination下的信息经过筛选后的信息
        # id：被保留的流线id，一维
        f_rft = []
        pre_rft = self.rft.copy()
        f_rft = pre_rft[id]
        return np.array(f_rft)



    def __save_fs(self, vtk_file_path, output_file_path):
        filter_line_id = np.array(np.where(self.filter_line_bool == 1)).reshape(-1) # 被保留的流线id
        filter_lines = self.lines[filter_line_id].copy()      # 流线信息，其中每个数组表示该流线包括的点id
        line_size = len(filter_line_id)                  # 剩下的流线数
        new_count = 0     # 新的流线计数
        old_count = 0     # 文件中的流线比新的多出来的部分计数
        rft = self.__get_rft_f(filter_line_id)                # 剩下的流线代表的RFT，一维
        flag_LRS = np.zeros([3])
        with open(vtk_file_path) as vf:
            with open(output_file_path, 'w+') as ff:
                for w in vf.readlines():
                    if 'POINT_DATA' in w:
                        ff.write(w)
                        old_count -= 1
                        continue
                    if new_count < 0:
                        new_count = 0
                        # 循环完毕后，重置flag
                        flag_LRS = np.zeros([3])

                    if 'LINES' in w:
                        # 写入 LINES 流线数 流线数+所有流线包括的点数
                        flag_LRS[0] = 1
                        info = w.split(' ')[:]
                        old_count = len(self.lines) - line_size
                        info[1] = str(line_size)
                        total_data = 0
                        for i in filter_lines:
                            total_data += i.size
                        info[2] = str(total_data + line_size)
                        new_count = line_size
                        ff.write(info[0] + " " + info[1] + " " + info[2] + '\n')
                        continue
                    if 'CELL_DATA' in w:
                        info = w.split(' ')[:]
                        info[1] = line_size
                        ff.write(str(info[0]) + " " + str(info[1]) +  '\n')
                        continue
                    if 'ReasonForTermination' in w:
                        # 写入 ReasonForTermination 信息
                        flag_LRS[1] = 1
                        info = w.split(' ')[:]
                        old_count = (len(self.lines) - line_size) // 9 + 1
                        info[2] = str(line_size)
                        new_count = line_size
                        ff.write(info[0] + " " + info[1] + " " + info[2] + " " + info[3])
                        continue
                    if 'SeedIds' in w:
                        # 写入SeedIds信息
                        flag_LRS[2] = 1
                        info = w.split(' ')[:]
                        old_count = (len(self.lines) - line_size) // 9 + 1
                        info[2] = str(line_size)
                        new_count = line_size
                        ff.write(info[0] + " " + info[1] + " " + info[2] + " " + info[3])
                        continue
                    if new_count == 0 and old_count > 0:
                        old_count -= 1
                        continue
                    if new_count <= 0:
                        # 循环完毕后，重置flag
                        new_count = 0
                        flag_LRS = np.zeros([3])
                        ff.write(w)
                        continue
                    else:
                        if flag_LRS[0] == 1:
                            # LINES下面的行
                            this_line = filter_lines[line_size - new_count]
                            data_num = len(this_line)
                            write_into_file = ''
                            write_into_file += str(data_num)
                            for i in range(data_num):
                                write_into_file = write_into_file + ' ' + str(this_line[i])
                            ff.write(write_into_file + '\n')
                            new_count -= 1
                        elif flag_LRS[1] == 1:
                            # ReasonForTermination
                            write_into_file = ''
                            c_low = line_size - new_count
                            c_high = c_low + 9
                            if c_high > line_size:
                                c_high = line_size
                            for i in range(c_low, c_high):
                                write_into_file += str(rft[i]) + ' '
                            ff.write(write_into_file + '\n')
                            new_count -= 9
                        elif flag_LRS[2] == 1:
                            # SeedIds
                            len2 = line_size // 2
                            write_into_file = ''
                            c_low = line_size - new_count
                            c_high = c_low + 9
                            if c_high > line_size:
                                c_high = line_size
                            for i in range(c_low, c_high):
                                write_into_file += str(self.idx % len2) + ' '
                                self.idx += 1
                            ff.write(write_into_file + '\n')

                            new_count -= 9



    def __judge_is_bi(self, fb):
        # 因为流线是从一个点向两个方向延伸的（双向积分），所以判断是否流线总数是单数
        len2 = len(fb) // 2
        for i in range(len2):
            if fb[i] != fb[i + len2]:
                fb[i] = 1
                fb[i + len2] = 1
        return fb

    def f_s(self, vtk_file_path, result_path):
        # 筛选流线
        obj = generate_data()
        self.lines, self.all_cor = obj.get_xyz(vtk_file_path)
        self.rft = obj.ret_rft()
        filter_line_bool = np.load(result_path)    # 大小为流线数，只有0和1，1代表该流线被留下
        # filter_line_bool = np.random.randint(0, 2, (1052))
        self.filter_line_bool = filter_line_bool
        self.__save_fs(vtk_file_path, 'tec030_query_1.vtk') #suggest save your final file in 'name of flow field'_'name of experiment'_'number of streamlines you want to visualize'.vtk

if __name__ == '__main__':
        file_path = 'tec030.vtk' #file path of the source flow field in vtk format
        r_path = 'tec030_query_1.npy' #file path of the serial numbers of the streamlines that you want to visualize
        obj = filter_streamlines()
        obj.f_s(file_path, r_path)