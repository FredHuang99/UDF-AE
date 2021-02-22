import numpy as np

class generate_data(object):
    def __init__(self):
        '''Initialize member variables.
                '''
        self.reader = None
        self.all_cor = None
        self.all_vec = None
        self.all_scale = None
        self.scale_name = None
        self.data_num = 0
        self.line_num = 0
        self.lines = None
        self.max_cor, self.min_cor = [], []  # 最大最小坐标
        self.center = []
        self.sprinkle_type = 0
        self.streamlines_max_length = -1
        self.xyz_flag = 0
        self.rft = None

    def __clear_all_data(self):
        '''Clear data for all class member variables.
        '''
        self.__init__()

    def __save_data_info(self, vtk_format_filename):
        '''
        遍历得到所有坐标值和最大最小坐标值，以及得到所有vector和标量
        对于航天的变量，三维："t","p","g",    "d", "e","ox","oy","oz","om","M"
                          温度 压力 标量浓度 密度 能量 涡量                马数
                      二维："vis","den","e","oz","Ma"
                            粘性   密度 能量 涡量  马数
        '''
        with open(vtk_format_filename) as f:
            file_lines = f.readlines()
            count = 0
            all_data = []
            all_scale = []
            scale_name = []
            line_count = 0
            r_count = 0
            lines = []
            rft = []   # ReasonForTermination
            r_l = []

            for line in file_lines:
                # 去除该行末尾的\n和前后空格，再以空格分割
                data = line.strip('\n').strip().split(' ')
                # point数
                if data[0] == 'POINTS':
                    count, n = int(data[1])//3 + int(data[1]) % 3, int(data[1])
                    self.data_num = n
                    continue
                if data[0] == 'LINES':
                    line_count = int(data[1])
                    self.line_num = line_count
                    continue
                if 'ReasonForTermination' in data:
                    r_count = self.line_num
                    continue
                if count > 0:
                    count -= 1
                    for j in range(len(data) // 3):
                        data_float = [float(data[j*3 + i]) for i in range(3)]
                        all_data.append(data_float[:])
                if line_count > 0:
                    line_count -= 1
                    p_id = [int(data[i]) for i in range(1, int(data[0]) + 1)]
                    lines.append(np.array(p_id))
                if r_count>0:
                    r_count -= 9
                    for i in range(len(data)):
                        r_l.append(int(data[i]))
                    if r_count < 0:
                        r_count = -r_count
                        for _ in range(r_count):
                            r_l.append(None)
                        r_count = 0
                rft = np.array(r_l.copy())


        self.all_cor = np.array(all_data[:n])
        self.lines = np.array(lines)
        max_cor = np.max(self.all_cor, axis=0)
        min_cor = np.min(self.all_cor, axis=0)
        self.max_cor = max_cor
        self.min_cor = min_cor
        self.rft = np.array(rft)

    def ret_rft(self):
        return self.rft

    def get_xyz(self, vtk_file_path):
        self.__save_data_info(vtk_file_path)
        return self.lines, self.all_cor
