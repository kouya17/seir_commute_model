import csv
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import population
import json


class SEIRModel:
    def __init__(self, m, a, b, g, s, e, i, r, t):
        # print('SEIRModel init')
        # print('m=', m, 'a=', a, 'b=', b, 'g=', g)
        # print('S=', s, 'E=', e, 'I=', i, 'R=', r)
        self.m = m
        self.a = a
        self.b = b
        self.g = g
        self.s = s
        self.e = e
        self.i = i
        self.r = r
        self.t = t
        
    def iterate_rk4(self, dt):
        # 今回扱う4式を配列にまとめる
        equation_list = [self.ds_dt, self.de_dt, self.di_dt, self.dr_dt]
        # ルンゲクッタ係数の初期化
        k1_list = [0] * len(equation_list)
        k2_list = [0] * len(equation_list)
        k3_list = [0] * len(equation_list)
        k4_list = [0] * len(equation_list)
        
        # 1段目をまわす
        for index, eq in enumerate(equation_list):
            k1_list[index] = eq(self.t, self.s, self.e, self.i, self.r)
        # 2段目をまわす
        for index, eq in enumerate(equation_list):
            k2_list[index] = eq(self.t + (dt / 2), self.s + (dt / 2) * k1_list[0], self.e + (dt / 2) * k1_list[1],
                                self.i + (dt / 2) * k1_list[2], self.r + (dt / 2) * k1_list[3])
        # 3段目をまわす
        for index, eq in enumerate(equation_list):
            k3_list[index] = eq(self.t + (dt / 2), self.s + (dt / 2) * k2_list[0], self.e + (dt / 2) * k2_list[1],
                                self.i + (dt / 2) * k2_list[2], self.r + (dt / 2) * k2_list[3])
        # 4段目をまわす
        for index, eq in enumerate(equation_list):
            k4_list[index] = eq(self.t + dt, self.s + dt * k3_list[0], self.e + dt * k3_list[1],
                                self.i + dt * k3_list[2], self.r + dt * k3_list[3])
        
        # 値を更新する
        self.s = self.s + (dt / 6) * (k1_list[0] + 2 * k2_list[0] + 2 * k3_list[0] + k4_list[0])
        self.e = self.e + (dt / 6) * (k1_list[1] + 2 * k2_list[1] + 2 * k3_list[1] + k4_list[1])
        self.i = self.i + (dt / 6) * (k1_list[2] + 2 * k2_list[2] + 2 * k3_list[2] + k4_list[2])
        self.r = self.r + (dt / 6) * (k1_list[3] + 2 * k2_list[3] + 2 * k3_list[3] + k4_list[3])
        self.t = self.t + dt
        
        # 呼び出しもとにも更新後の値を返却する
        return self.s, self.e, self.i, self.r
    
    def ds_dt(self, t, s, e, i, r):
        return self.m * (e + i + r) - self.b * s * i
    
    def de_dt(self, t, s, e, i, r):
        return self.b * s * i - (self.m + self.a) * e
    
    def di_dt(self, t, s, e, i, r):
        return self.a * e - (self.m + self.g) * i
    
    def dr_dt(self, t, s, e, i, r):
        return self.g * i - self.m * r


class InjectionSEIRModel(SEIRModel):
    def __init__(self, m, a, b, g, s, e, i, r, t):
        super().__init__(m, a, b, g, s, e, i, r, t)
        self.__b_i_inject = 0
    
    def set_b_i(self, b_i_inject):
        self.__b_i_inject = b_i_inject
    
    def ds_dt(self, t, s, e, i, r):
        return self.m * (e + i + r) - self.__b_i_inject * s
    
    def de_dt(self, t, s, e, i, r):
        return self.__b_i_inject * s - (self.m + self.a) * e


class DayModel:
    def __init__(self, seir_models, x, r_0, infect_t, t):
        self.__seir_models = [0] * len(seir_models)
        for i in range(len(seir_models)):
            self.__seir_models[i] = InjectionSEIRModel(seir_models[i].m, seir_models[i].a, seir_models[i].b, seir_models[i].g,
                                                       seir_models[i].s, seir_models[i].e, seir_models[i].i, seir_models[i].r,
                                                       seir_models[i].t)
        
        self.__x = x
        self.__r_0 = r_0
        self.__infect_t = infect_t
        
        self.__move_seir_models = [[0 for j in range(len(seir_models))] for i in range(len(seir_models))]
        w = [[[0 for k in range(4)] for j in range(len(seir_models))] for i in range(len(seir_models))]
        for i in range(len(seir_models)):
            for j in range(len(seir_models)):
                n_i = seir_models[i].s + seir_models[i].e + seir_models[i].i + seir_models[i].r
                w[i][j][0] = x[i][j] * seir_models[i].s / n_i
                w[i][j][1] = x[i][j] * seir_models[i].e / n_i
                w[i][j][2] = x[i][j] * seir_models[i].i / n_i
                w[i][j][3] = x[i][j] * seir_models[i].r / n_i
                self.__move_seir_models[i][j] = InjectionSEIRModel(seir_models[j].m, seir_models[j].a, seir_models[j].b, seir_models[j].g,
                                                                   w[i][j][0], w[i][j][1], w[i][j][2], w[i][j][3], t)
                # 移動元の人数を減らす
                self.__seir_models[i].s = self.__seir_models[i].s - w[i][j][0]
                self.__seir_models[i].e = self.__seir_models[i].e - w[i][j][1]
                self.__seir_models[i].i = self.__seir_models[i].i - w[i][j][2]
                self.__seir_models[i].r = self.__seir_models[i].r - w[i][j][3]
    
    def iterate(self, dt):
        # 感染率計算用に都道府県別のI, Nを計算する
        n_pref_group = [0] * len(self.__seir_models)
        i_pref_group = [0] * len(self.__seir_models)
        for i in range(len(self.__seir_models)):
            n_pref_group[i] = 0
            i_pref_group[i] = 0
            for j in range(len(self.__seir_models)):
                n_pref_group[i] = n_pref_group[i] + self.__move_seir_models[j][i].s
                n_pref_group[i] = n_pref_group[i] + self.__move_seir_models[j][i].e
                n_pref_group[i] = n_pref_group[i] + self.__move_seir_models[j][i].i
                i_pref_group[i] = i_pref_group[i] + self.__move_seir_models[j][i].i
                n_pref_group[i] = n_pref_group[i] + self.__move_seir_models[j][i].r
            n_pref_group[i] = n_pref_group[i] + self.__seir_models[i].s
            n_pref_group[i] = n_pref_group[i] + self.__seir_models[i].e
            n_pref_group[i] = n_pref_group[i] + self.__seir_models[i].i
            i_pref_group[i] = i_pref_group[i] + self.__seir_models[i].i
            n_pref_group[i] = n_pref_group[i] + self.__seir_models[i].r
        
        data = [[0 for j in range(4)] for i in range(len(self.__seir_models))]
        for i in range(len(self.__seir_models)):
            # 都道府県の感染率を計算する
            infection_rate = self.__r_0 * i_pref_group[i] / self.__infect_t / n_pref_group[i]
            
            # 移動グループを時間発展させる
            for j in range(len(self.__seir_models)):
                self.__move_seir_models[j][i].set_b_i(infection_rate)
                s, e, i_, r = self.__move_seir_models[j][i].iterate_rk4(dt)
                data[i][0] = data[i][0] + s
                data[i][1] = data[i][1] + e
                data[i][2] = data[i][2] + i_
                data[i][3] = data[i][3] + r
            # 非移動グループを時間発展させる
            self.__seir_models[i].set_b_i(infection_rate)
            s, e, i_, r = self.__seir_models[i].iterate_rk4(dt)
            data[i][0] = data[i][0] + s
            data[i][1] = data[i][1] + e
            data[i][2] = data[i][2] + i_
            data[i][3] = data[i][3] + r
        
        return data
    
    # 移動グループを元の都道府県に帰してモデルを組み直す
    def reset_commuter(self):
        night_seir_model = [0] * len(self.__seir_models)
        for i in range(len(self.__seir_models)):
            return_people = [0] * 4
            for j in range(len(self.__seir_models)):
                return_people[0] = return_people[0] + self.__move_seir_models[i][j].s
                return_people[1] = return_people[1] + self.__move_seir_models[i][j].e
                return_people[2] = return_people[2] + self.__move_seir_models[i][j].i
                return_people[3] = return_people[3] + self.__move_seir_models[i][j].r
            night_seir_model[i] = SEIRModel(self.__seir_models[i].m, self.__seir_models[i].a, self.__seir_models[i].b, self.__seir_models[i].g,
                                            self.__seir_models[i].s + return_people[0], self.__seir_models[i].e + return_people[1],
                                            self.__seir_models[i].i + return_people[2], self.__seir_models[i].r + return_people[3],
                                            self.__seir_models[i].t)
        
        return night_seir_model


class NightModel:
    def __init__(self, seir_models):
        self.__seir_models = seir_models
    
    def iterate(self, dt):
        data = [[0 for j in range(4)] for i in range(len(self.__seir_models))]
        for index in range(len(self.__seir_models)):
            s, e, i, r = self.__seir_models[index].iterate_rk4(dt)
            data[index][0] = s
            data[index][1] = e
            data[index][2] = i
            data[index][3] = r
        return data
        

class SEIRCommuteModel:
    def __init__(self, seir_models, x, r_0, infect_t, t_initial_day):
        print('SEIRCommuteModel init')
        print('x: ', x)
        self.__seir_models = seir_models
        self.__x = x
        self.__r_0 = r_0
        self.__infect_t = infect_t
        self.__t = t_initial_day
        hour = (t_initial_day * 24) % 24
        if 8 <= hour < 17:
            self.__state_model = DayModel(seir_models, x, r_0, infect_t, t_initial_day)
        else:
            self.__state_model = NightModel(seir_models)
        print('t: ', self.__t)
    
    def iterate(self, dt):
        # 日中、夜間状態の切り替え
        pre_hour = (self.__t * 24) % 24
        next_hour = ((self.__t + dt) * 24) % 24
        state = '↓'
        if pre_hour <= 8 and next_hour > 8 and int(self.__t % 7) != 5 and int(self.__t % 7) != 6 \
           and type(self.__state_model).__name__ == 'NightModel':
            self.__state_model = DayModel(self.__seir_models, self.__x, self.__r_0, self.__infect_t, self.__t)
            state = 'day'
        elif pre_hour <= 17 and next_hour > 17 and type(self.__state_model).__name__ == 'DayModel':
            self.__seir_models = self.__state_model.reset_commuter()
            self.__state_model = NightModel(self.__seir_models)
            state = 'night'
        
        self.__t = self.__t + dt
        print('t[day]: ', '{:>6.2f}'.format(self.__t), 'state: ', state)
        return self.__state_model.iterate(dt)


class FixedOrderFormatter(ScalarFormatter):
    def __init__(self, order_of_mag=0, useOffset=True, useMathText=True):
        self._order_of_mag = order_of_mag
        ScalarFormatter.__init__(self, useOffset=useOffset, useMathText=useMathText)
    
    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self._order_of_mag
        

def update_data_list(time_list, s_list, e_list, i_list, r_list, time, data):
    time_list.append(time)
    for i in range(len(data)):
        s_list[i].append(data[i][0])
        e_list[i].append(data[i][1])
        i_list[i].append(data[i][2])
        r_list[i].append(data[i][3])


def make_output_info_file(pref_info_list, pref_path_list, out_file_path):
    output_info_dict = {}
    for index, pref_info in enumerate(pref_info_list):
        output_info_dict[pref_info[0]] = pref_path_list[index]
    with open(out_file_path, 'w') as f:
        json.dump(output_info_dict, f, indent=4, ensure_ascii=False)


def output_figure(output_file_path, time_list, s_list, e_list, i_list, r_list):
    # fig, ax = plt.subplots(1, len(s_list), sharey=True, figsize=(16, 4))
    fig, ax = plt.subplots(1, 7, sharey=True, figsize=(16, 4))
    # titles = ['tokyo', 'ibaraki', 'tochigi', 'gunma', 'saitama', 'chiba', 'kanagawa']
    titles = ['Hokkaido', 'Aomori', 'Iwate', 'Miyagi', 'Akita', 'Yamagata', 'Fukushima']
    ax[0].set_ylabel('number of people')
    ax[3].set_xlabel('days')
    for i in range(len(ax)):
        ax[i].plot(time_list, s_list[i], label="S")
        ax[i].plot(time_list, e_list[i], label="E")
        ax[i].plot(time_list, i_list[i], label="I")
        ax[i].plot(time_list, r_list[i], label="R")
        ax[i].set_title(titles[i])
        ax[i].grid(True, linestyle='-.')
        # ax[i].yaxis.set_major_formatter(FixedOrderFormatter(4, useMathText=True))
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.show()


if __name__ == "__main__":
    # パラメータの設定
    # 基本再生産数
    R_0 = 2
    # 平均潜伏期間[日]
    L = 5
    # 平均発症期間[日]
    I = 14
    # 時間刻み[日]
    DT = 0.1
    
    # 初期値の設定
    population_list = population.get_population_list('./input/population.csv')
    print(population_list)
    INITIAL_DATA = [[0 for j in range(len(population_list))] for i in range(len(population_list))]
    for i in range(len(population_list)):
        # 東京の場合
        if i == 12:
            INITIAL_DATA[i][0] = population_list[i][1] - 1
            INITIAL_DATA[i][2] = 1
        else:
            INITIAL_DATA[i][0] = population_list[i][1]
    '''
    INITIAL_DATA = [[13951635, 0, 1, 0],
                    [2866325, 0, 0, 0],
                    [1942313, 0, 0, 0],
                    [1938053, 0, 0, 0],
                    [7341794, 0, 0, 0],
                    [6280344, 0, 0, 0],
                    [9204965, 0, 0, 0]]
    '''
    # 全体人数
    N = [0] * len(INITIAL_DATA)
    for i in range(len(N)):
        N[i] = N[i] + INITIAL_DATA[i][0] + INITIAL_DATA[i][1] + INITIAL_DATA[i][2] + INITIAL_DATA[i][3]
    # 初期時刻[日]
    T_INIT = 0
    
    # SEIRモデルに合わせてパラメータを変形する
    # 出生率・死亡率は今回は無視する
    m = 0
    # 感染症の発生率
    a = 1 / L
    # 感染症への感染率
    b = [0] * len(INITIAL_DATA)
    for i in range(len(b)):
        b[i] = R_0 / I / N[i]
    # 感染からの回復率
    g = 1 / I
    
    # 何日分シミュレーションするか
    T_END = 1000
    # 結果出力ファイルパス
    FILE_PATH_NAME = './output/result'
    
    commute_list = population.get_commute_list('./input/commute.csv')
    seir_models_pref = [0] * len(INITIAL_DATA)
    x = [[0 if j == i else commute_list[i][j] for j in range(len(INITIAL_DATA))] for i in range(len(INITIAL_DATA))]
    for i in range(len(INITIAL_DATA)):
        for j in range(len(INITIAL_DATA)):
            x[i][j] = x[i][j] * 0.0001
    '''
    x = [[0, 7619, 2770, 2251, 140961, 82706, 238314],
         [67284, 0, 22098, 1166, 17807, 41734, 3748],
         [17301, 18175, 0, 23503, 12067, 1197, 1772],
         [13614, 1075, 16561, 0, 27904, 895, 1518],
         [936105, 14437, 10049, 29250, 0, 43074, 28111],
         [716882, 35427, 1145, 791, 41668, 0, 25966],
         [1068513, 2688, 1433, 1151, 14035, 14932, 0]]
    '''
    for i in range(len(seir_models_pref)):
        seir_models_pref[i] = SEIRModel(m, a, b[i], g, INITIAL_DATA[i][0], INITIAL_DATA[i][1], INITIAL_DATA[i][2], INITIAL_DATA[i][3], T_INIT)
    seir_commute_model = SEIRCommuteModel(seir_models_pref, x, R_0, I, T_INIT)
    
    now_time = T_INIT
    now_data = [[0 for j in range(4)] for i in range(len(INITIAL_DATA))]
    for i in range(len(INITIAL_DATA)):
        for j in range(4):
            now_data[i][j] = INITIAL_DATA[i][j]
    
    time_list = []
    s_list = [[] for i in range(len(INITIAL_DATA))]
    e_list = [[] for i in range(len(INITIAL_DATA))]
    i_list = [[] for i in range(len(INITIAL_DATA))]
    r_list = [[] for i in range(len(INITIAL_DATA))]
    
    pref_path_list = [FILE_PATH_NAME + str(i) + '.txt' for i in range(len(population_list))]
    make_output_info_file(population_list, pref_path_list, './output_info.txt')
    
    for pref_num, file_path in enumerate(pref_path_list):
        with open(file_path, mode='w', newline="") as f:
            writer = csv.writer(f)
        
            writer.writerow(['day', 'S', 'E', 'I', 'R'])
            writer.writerow([now_time, now_data[pref_num][0], now_data[pref_num][1],
                             now_data[pref_num][2], now_data[pref_num][3]])
            update_data_list(time_list, s_list, e_list, i_list, r_list,
                             now_time, now_data)
        
    while now_time < T_END:
        now_data = seir_commute_model.iterate(DT)
        now_time = now_time + DT
        now_time_int = round(now_time)
        # if abs(now_time - now_time_int) < 0.001:
        for pref_num, file_path in enumerate(pref_path_list):
            with open(file_path, mode='a', newline="") as f:
                writer = csv.writer(f)
                writer.writerow([now_time, now_data[pref_num][0], now_data[pref_num][1],
                                now_data[pref_num][2], now_data[pref_num][3]])
        update_data_list(time_list, s_list, e_list, i_list, r_list,
                         now_time, now_data)
    
    output_figure('', time_list, s_list, e_list, i_list, r_list)
