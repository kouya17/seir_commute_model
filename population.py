import csv
import re


def is_pref(str):
    return re.match(r'.*[都道府県]$', str)


def get_population_list(population_file_path):
    with open(population_file_path, encoding='utf-8') as input_file:
        reader = csv.reader(input_file)
        population_list = []
        for row in reader:
            if is_pref(row[6]):
                population_list.append([row[6], int(row[7])])
    return population_list


def get_commute_list(commute_file_path):
    with open(commute_file_path, encoding='utf-8') as input_file:
        reader = csv.reader(input_file)
        commute_list = [[0 for j in range(47)] for i in range(47)]
        rows = [row for row in reader]
        for i in range(47):
            for j in range(47):
                commute_list[i][j] = 0 if rows[i + 11][j + 8] == '-' else int(rows[i + 11][j + 8])
    return commute_list


if __name__ == "__main__":
    population_list = get_population_list('./population.csv')
    print(population_list)
    print(len(population_list))
    commute_list = get_commute_list('./commute.csv')
    print(commute_list)
