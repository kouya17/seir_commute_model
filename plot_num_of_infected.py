import matplotlib.pyplot as plt
import json
import pandas as pd
import numpy as np
from matplotlib import rcParams
from matplotlib.ticker import ScalarFormatter


def get_output_info_dict(info_file_path):
    with open(info_file_path) as f:
        return json.load(f)


class FixedOrderFormatter(ScalarFormatter):
    def __init__(self, order_of_mag=0, useOffset=True, useMathText=True):
        self._order_of_mag = order_of_mag
        ScalarFormatter.__init__(self, useOffset=useOffset, useMathText=useMathText)
    
    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self._order_of_mag


if __name__ == "__main__":
    OUTPUT_PATH = './output'
    COMMUTE_FOLDER_LIST = ['commute00', 'commute00001', 'commute02', 'commute10']
    RESULT_NAME = 'result'
    
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Yu Gothic', 'Meirio', 'Takao', 'IPAexGothic', 'IPAPGothic', 'VL PGothic', 'Noto Sans CJK JP']
    
    output_info_dict = get_output_info_dict('./output_info.txt')
    print('output_info_dict: ', output_info_dict)
    
    labels = []
    for k in output_info_dict.keys():
        labels.append(k)
    
    infected_list = []
    for series in range(len(COMMUTE_FOLDER_LIST)):
        df_list = []
        for i in range(len(output_info_dict)):
            df_list.append(pd.read_csv(OUTPUT_PATH + '/' + COMMUTE_FOLDER_LIST[series] +
                                       '/' + RESULT_NAME + str(i) + '.txt'))
        
        for d in df_list:
            d['EIR'] = d['E'] + d['I'] + d['R']
        
        infected = []
        for i in range(len(output_info_dict)):
            infected.append(df_list[i].iloc[-1]['EIR'])
        infected_list.append(infected)
    
    x = np.arange(len(labels))
    width = 0.2
    
    fig = plt.figure(figsize=(10.0, 6.0))
    ax = fig.add_subplot(111)
    rects1 = ax.bar(x - 1.5 * width, infected_list[0], width, label='通勤・通学なし')
    rects2 = ax.bar(x - 0.5 * width, infected_list[1], width, label='通勤・通学あり、人数加工なし')
    rects3 = ax.bar(x + 0.5 * width, infected_list[2], width, label='通勤・通学あり、人数20%')
    rects4 = ax.bar(x + 1.5 * width, infected_list[3], width, label='通勤・通学なし、人数0.01%')
    
    ax.set_ylabel('累計感染者数')
    ax.yaxis.set_major_formatter(FixedOrderFormatter(4, useMathText=True))
    ax.set_title('1000日目時点での累計感染者数')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    # x軸縦書き（90度回転）
    plt.xticks(rotation=90)
    ax.legend()
    
    fig.tight_layout()
    
    plt.show()
