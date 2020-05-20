import matplotlib.pyplot as plt
from japanmap import picture
import json
import pandas as pd


def get_output_info_dict(info_file_path):
    with open(info_file_path) as f:
        return json.load(f)


if __name__ == "__main__":
    output_info_dict = get_output_info_dict('./output_info.txt')
    print('output_info_dict: ', output_info_dict)
    
    df_list = []
    for v in output_info_dict.values():
        df_list.append(pd.read_csv(v))
    
    # 人口表示テスト
    '''
    columns = df_list[0].columns.values
    initial_list = []
    for d in df_list:
        initial_list.append(d.iloc[0].values.tolist())
    print(columns)
    print(initial_list)
    pref_list = []
    for k in output_info_dict.keys():
        pref_list.append(k)
    
    population_df = pd.DataFrame(initial_list, columns=columns, index=pref_list)
    cmap = plt.get_cmap('Reds')
    norm = plt.Normalize(vmin=population_df.S.min(), vmax=population_df.S.max())
    fcol = lambda x: '#' + bytes(cmap(norm(x), bytes=True)).hex()
    plt.colorbar(plt.cm.ScalarMappable(norm, cmap))
    print(population_df.S.apply(fcol))
    plt.imshow(picture(population_df.S.apply(fcol)))
    plt.show()
    '''
    
    for d in df_list:
        d['EI_ratio'] = (d['I'] + d['E']) / (d['S'] + d['E'] + d['I'] + d['R'])
    print(df_list[0])
    columns = df_list[0].columns.values
    pref_list = []
    for k in output_info_dict.keys():
        pref_list.append(k)
    num = 0
    for i in range(len(df_list[0])):
        print(i)
        if i % 100 != 0:
            continue
        data_list = []
        for d in df_list:
            data_list.append(d.iloc[i].values.tolist())
        df = pd.DataFrame(data_list, columns=columns, index=pref_list)
        fig = plt.figure()
        cmap = plt.get_cmap('Reds')
        norm = plt.Normalize(vmin=0, vmax=1)
        fcol = lambda x: '#' + bytes(cmap(norm(x), bytes=True)).hex()
        plt.colorbar(plt.cm.ScalarMappable(norm, cmap))
        print(df.EI_ratio.apply(fcol))
        plt.title('t = {:>4} [day]'.format(int(i * 0.1)))
        plt.imshow(picture(df.EI_ratio.apply(fcol)))
        fig.savefig('./figure/{:0=3}.png'.format(num))
        num = num + 1
