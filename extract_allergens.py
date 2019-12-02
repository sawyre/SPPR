import argparse

import pandas as pd
import numpy as np
import re
import codecs
import nltk
from nltk import word_tokenize
import time

# nltk.download()

#!/usr/bin/python
# -*- coding: utf-8 -*-

def levenstein_dist(string1, string2):
    # Заполнение крайнего левого столбца и первой строки числами от 0 до длины строк, так как расстояние
    # между нулевой подстрокой и другой строкой равно длине второй.
    lev_dist = [[(i + j) if i == 0 or j == 0 else 0 for j in range(len(string2) + 1)] for i in range(len(string1) + 1)]
    # Непосредственное вычисление расстояния Левинштейна для подстрок двух строк
    for i in range(1, len(string1) + 1):
        for j in range(1, len(string2) + 1):
            # Если крайние символы двух подстрок равны, то расстояние Левинштейна между двумя строками
            # равно расстоянию между строками без этих крайних символов
            if string1[i - 1] == string2[j - 1]:
                lev_dist[i][j] = lev_dist[i-1][j-1]
            else:
                # Если крайние символы не равны, то необходимо совершить как минимум одну праЭвку.
                # В зависимости от правки:
                # 1) Удаление символа из первой подстроки (или добавление символа во вторую)
                # 2) Удаление символа из второй подстроки (или добавление символа в первую)
                # 3) Замена символа в одной из строк
                # Могут потребоваться дополнительные изменения, из которых выбирается наименьшее
                lev_dist[i][j] = 1 + min(lev_dist[i-1][j], lev_dist[i][j-1], lev_dist[i-1][j-1])
    # Значение расположенное в последнем столбце и строке матрицы и является минимальным расстоянием Левенштейна
    # между двумя числами
    return lev_dist[len(string1)][len(string2)]

# В данном блоке удалено все кроме русских знаков и точек (скобки, тире, подчеркивания, палки, лишние пробелы)
# Также все переведено в нижний регистр (нужно для приведения всего текста в более менее стандартный вид)
def text_cleaner(text):
    text = text.replace('\n', ' ').replace('\r', ' ') # удаение переноса строк
    rules = [
        {r'[a-z0-9]': u' '},
        {r'/': u' '},
        {r'\\': u' '},
        {r'\{': u' '},
        {r'\}': u' '},
        {r'\(': u' '},
        {r'\)': u' '},
        {r'\[': u' '},
        {r'\]': u' '},
        {r'-': u' '},
        {r'_': u' '},
        {r',': u' '},  # До этого момента - удаление мусора
        {r' +': u' '}  # До этого момента - удаление лишних пробелов
    ]
    for rule in rules:
        for (k, v) in rule.items():
            regex = re.compile(k)
            text = regex.sub(v, text)
            text = text.rstrip()
    return text.lower().strip() # Перевод в нижний регистр и удаление пробелов в начале и в конце

parser = argparse.ArgumentParser()
parser.add_argument("anamnes", help="Path to anamnes file.")
parser.add_argument("data", help="Path to data file.")
parser.add_argument("dict", help="Path to dict file.")
parser.add_argument("predict", help="Path to predict file.")
parser.add_argument("mark", help="Path to mark file.")
args = parser.parse_args()
print(args.anamnes)
print(args.data)
print(args.dict)
print(args.predict)
print(args.mark)

with open(args.anamnes, 'r') as f:
    data = f.readlines();
    
print('количество данных', len(data), '\n')

find_string = 'Аллергологический анамнез'.lower()
acceptable_dist = 1
i_stop = 6000
i_start = 4000
count_alerg = 0

start = time.time()

with open(args.data, 'w') as f, open(args.predict, 'w') as pr:
    dct = codecs.open(args.dict, 'r', 'utf-8')
    allergens = dct.readlines();  # словарь с аллергенами
    
    for i in range(len(allergens)):
        allergens[i] = text_cleaner(allergens[i]) # отчистка строк от мусора
        print(allergens[i])
        # print(allergens[i])
    dct.close()
    
    for i in range(i_start, i_stop): # сами записи наблюдений у врача          
        data[i] = text_cleaner(data[i]) # отчистка строк от мусора
        
        if data[i].find(find_string) != -1:
            count_alerg += 1
        
        if len(data[i][data[i].find(find_string):]) < 5:
            pass
            #pr.write(str(i) + ' empty\n')
            #f.write(str(i) + '\n')
        elif data[i].find('без особенностей') != -1 or data[i].find('неотягощен') != -1\
            or data[i].find('спокойный') != -1  or data[i].find('не отягощен') == 0 or data[i].find('отр.') == 0: 
            pass
            #pr.write(str(i) + ' without_features\n')
            #f.write(str(i) + ' ' + data[i][data[i].find(find_string) + len(find_string) + 1:] + '\n')
        else:
            tmp = data[i][data[i].find(find_string) + len(find_string) + 1:]
            f.write(str(i) + ' ' + tmp + '\n')
            pr.write(str(i) + ' ')
            
            for allergen in allergens:
                if len(allergen.split(' ')) > 1:
                    if tmp.find(allergen) != -1:
                        pr.write(allergen + ' ')
                else:
                    for word in tmp.split(' '):
                        if levenstein_dist(word, allergen) <= acceptable_dist:
                            pr.write(allergen + ' ')
                            break
                            
            pr.write('\n')
    
print('Время работы:', time.time() - start, 'секунд')
print('Количество упоминаний', find_string + ':', count_alerg)

f = codecs.open(args.mark, 'r', 'utf-8')
marks = f.readlines();  # словарь с аллергенами

for i in range(len(marks)):
    marks[i] = [int(x) for x in marks[i].strip().split(' ')]
f.close()

fact = 0
tp = 0
fp = 0

for mark in marks:
    fact += mark[1]
    tp += mark[2]
    fp += mark[3]
    
print('Fact count allergens:', fact)
print('True possitive allergens:', tp)
print('False possitive allergens:', fp)
print('Recall:', tp / fact)
print('Precision:', tp / (tp + fp))

print('Completed')