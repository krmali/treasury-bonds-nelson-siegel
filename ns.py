from dataclasses import dataclass
from numpy.core.multiarray import result_type
from scipy.optimize import minimize
import numpy as np
import pandas as pd
import csv
from datetime import datetime
import os
# from nelson_siegel_svensson import NelsonSiegelSvenssonCurve
import matplotlib.pyplot as plt

@dataclass
class Bond:
    isin_code: str
    isin: str
    issue_date: datetime
    maturity_date: datetime
    date: datetime
    t: str
    month: int
    count: int
    avg_yield: float


arr = []

with open('Book2.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter = ',')
    line_count = 0
    ignore = False
    for row in csv_reader:
        if not ignore:
            ignore = True
            continue
        bond = Bond(row[0],row[1],datetime.strptime(row[2], "%M/%d/%Y"),
                datetime.strptime(row[3], "%M/%d/%Y"),
                datetime.strptime(row[4],"%M/%d/%y"),
                row[5],int(row[6]),
                int(row[7]),float(row[8]))
        arr.insert(line_count,bond)
        line_count = line_count + 1

arr_days = []
temp_arr = []
temp_date = arr[0].date
for b in arr:
    if b.date == temp_date:
        temp_arr.append(b)
    else:
        arr_days.append(temp_arr)
        temp_arr = []
        temp_date = b.date
        temp_arr.append(b)


def neslon_betas(yields):
    beta0 = np.nanmax(yields)
    beta1 = np.nanmax(yields) - np.nanmin(yields)
    beta2 = 2*np.nanmean(yields) - np.nanmax(yields) - np.nanmin(yields)
    beta3 = beta2
    lambda1 = 5
    lambda2 = 5
    return [beta0, beta1, beta2, beta3, lambda1, lambda2]


def NSS(m,vi):
    f = vi[0]+vi[1]*((1-np.exp(-m/vi[4]))/(m/vi[4])) + vi[2]*((1-np.exp(-m/vi[4]))/(m/vi[4])-np.exp(-m/vi[4])) + vi[3]*((1-np.exp(-m/vi[5]))/(m/vi[5])-np.exp(-m/vi[5]))
    return f

def calculate_yields(dict, bonds):
    for b in bonds:
        if not b.month in dict:
            dict[b.month] = {'yield': b.avg_yield, 'estimated': np.nan} 
        else:
            dict[b.month]['yield'] = b.avg_yield
    keys = list(dict.keys())
    keys.sort()
    new_dict = {}
    for k in keys:
        new_dict[k] = dict[k]
    return new_dict

def calculate_estimated(dict, vi):
    for i in dict.keys():
        ns = NSS(i,vi)
        dict[i]['estimated'] = ns
    return dict


# def month_yield_nss(month_dict, bonds):
#     for b in bonds: 
#         month_dict[b.month] = b.avg_yield
#     month_dict_keys = list(month_dict.keys())
#     month_dict_keys.sort()
#     month_dict_values = list(month_dict.values())
#     vi = neslon_betas(month_dict_values)
#     month_dict_values = []
#     nelsons = []
#     for i in month_dict_keys:
#         ns = NSS(i,vi)
#         nelsons.append(ns)
#         month_dict_values.append(month_dict[i])
#     tracer(month_dict_values,nelsons,month_dict_keys)

def tracer(date, type,x1,x2,y):
    plt.clf()
    plt.figure(figsize=(12.0, 7.0))
    p1 = plt.scatter(y,x1,marker='.', color='darkblue')
    p1 = plt.plot(y,x2,color='darkcyan')
    plt.title('yield curve for ' + date, color='darkblue')
    plt.xlabel('maturity')
    plt.ylabel('yield')
    plt.xticks(np.arange(0, 241,12))
    return plt.savefig(date+'/'+ type + '.png', dpi=100)

def generate_for_day(bonds, method, date):
    if not os.path.isdir(date):
        os.makedirs(date)
    global month_count_dict
    month_count_dict = {
            3:{'yield':np.nan, 'estimated':np.nan}
            ,6:{'yield':np.nan, 'estimated':np.nan}
            ,9:{'yield':np.nan, 'estimated':np.nan}
            ,12:{'yield':np.nan, 'estimated':np.nan}
            ,24:{'yield':np.nan, 'estimated':np.nan}
            ,36:{'yield':np.nan, 'estimated':np.nan}
            ,48:{'yield':np.nan, 'estimated':np.nan}
            ,60:{'yield':np.nan, 'estimated':np.nan}
            ,84:{'yield':np.nan, 'estimated':np.nan}
            ,120:{'yield':np.nan, 'estimated':np.nan}
            ,240:{'yield':np.nan, 'estimated':np.nan}
            } 

    month_count_dict = calculate_yields(month_count_dict, bonds)
    month_keys = list(month_count_dict.keys())
    yields = []
    for k in month_keys:
        yields.append(month_count_dict[k]['yield'])
    vi = neslon_betas(yields)
    month_count_dict = calculate_estimated(month_count_dict, vi)
    estimated = []
    for k in month_keys:
        estimated.append(month_count_dict[k]['estimated'])


    def error(vi):
        global month_count_dict 
        month_count_dict = calculate_estimated(month_count_dict, vi)
        estimated = []
        for k in month_keys:
            estimated.append(month_count_dict[k]['estimated'])
        index = 0
        sum = 0
        for i in yields:
            if i is not np.nan:
                sum = sum + (i - estimated[index])**2
            index = index + 1
        return sum

    err = error(vi)
    result_vi = minimize(error, vi, method=method)
    estimated = []
    for k in month_keys:
        estimated.append(month_count_dict[k]['estimated'])
    tracer(date, method, yields, estimated, month_keys)
    d = {'maturity': list(month_count_dict.keys())
            , 'yield': yields
            , 'estimated': estimated}

    data = pd.DataFrame(data=d)
    data.to_csv(date+'/' + method + '.csv', index=False)

    result_vi = list(result_vi['x'])
    d = {'beta0': result_vi[0]
            ,'beta1': result_vi[1]
            ,'beta2': result_vi[2]
            ,'beta3': result_vi[3]
            ,'lambda0': result_vi[4]
            ,'lambda1': result_vi[5]
            }
    data = pd.DataFrame([d])
    data.to_csv(date+'/' + method + '-parameters.csv', index=False)
    return month_count_dict


for i in range(0,len(arr_days)):
    num = i 
    generate_for_day(arr_days[num], 'Nelder-Mead' ,str(arr_days[num][0].date.date()))
    estimates = []
    res = generate_for_day(arr_days[num], 'SLSQP' ,str(arr_days[num][0].date.date()))
    estimates.append(res)
    # generate_for_day(arr_days[num], 'BFGS' ,str(arr_days[num][0].date.date()))

