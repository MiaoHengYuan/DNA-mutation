import pandas as pd
import numpy as np
import scipy
from matplotlib import pyplot as plt
from scipy import stats
import math
import scipy.optimize
import powerlaw


def hist_creat(df, iterate_no, c, distance_range='short', threshold=5000, length=50):
 y_axis_list = []
 x_axis_list = []

 x_axis_list_nnoise = []
 y_axis_list_nnoise = []

 inter_num = int((threshold / length) + 1)

 if distance_range == 'short':
  df_short = df.loc[:, 'distance_' + str(iterate_no) + 'chr' + c:'counts_' + str(iterate_no) + 'chr' + c][
   df['distance_' + str(iterate_no)+'chr'+c] <= threshold]
  for i in range(inter_num):
   x_axis_list.append((i + 1) * length)
   y_axis_list.append(sum(df_short['counts_' + str(iterate_no) + 'chr' + c][
                           (df_short['distance_' + str(iterate_no) + 'chr' + c] >= i * length) & (
                                    df_short['distance_' + str(iterate_no) + 'chr' + c] < (i + 1) * length)]))

   if (i + 1) * length <= 1200:
    x_axis_list_nnoise.append((i + 1) * length)
    y_axis_list_nnoise.append(sum(df_short['counts_' + str(iterate_no) + 'chr' + c][
                                   (df_short['distance_' + str(iterate_no) + 'chr' + c] >= i * length) & (
                                            df_short['distance_' + str(iterate_no) + 'chr' + c] < (i + 1) * length)]))

 if distance_range == 'long':
  df_long = df.loc[:, 'distance_' + str(iterate_no) + 'chr' + c:'counts_' + str(iterate_no) + 'chr' + c][
   df['distance_' + str(iterate_no) + 'chr' + c] > threshold]
  dis_len = df_long['distance_' + str(iterate_no) + 'chr' + c].max() - df_long['distance_' + str(iterate_no) + 'chr' + c].min()
  ite_num = np.ceil(dis_len / 15000)
  if pd.isna(ite_num):
   x_axis_list = [0]
   y_axis_list = [0]
   x_axis_list_nnoise = [0]
   y_axis_list_nnoise = [0]
  else:
   for i in range(int(ite_num)):
    x_axis_list.append((i + 1) * 15000)
    y_axis_list.append(sum(df_long['counts_' + str(iterate_no) + 'chr' + c][(df_long['distance_' + str(iterate_no) + 'chr' + c] >= i * 15000) & (
             df_long['distance_' + str(iterate_no) + 'chr' + c] < (i + 1) * 15000)]))

    if (i + 1) * 15000 <= 400000:
     x_axis_list_nnoise.append((i + 1) * 15000)
     y_axis_list_nnoise.append(sum(df_long['counts_' + str(iterate_no) + 'chr' + c][
                                    (df_long['distance_' + str(iterate_no) + 'chr' + c] >= i * 15000) & (
                                             df_long['distance_' + str(iterate_no) + 'chr' + c] < (i + 1) * 15000)]))

 return x_axis_list, y_axis_list, x_axis_list_nnoise, y_axis_list_nnoise

def plot_function(df, iterate_no, distance, threshold=5000, length=50):
 x_, y_, x_n, y_n = hist_creat(df, iterate_no, 'short', threshold, length)
 x_l, y_l, x_nl, y_nl = hist_creat(df, iterate_no, 'long', threshold, length)
 dist_exp = getattr(stats, 'expon')
 dist_powerlaw = getattr(stats, 'powerlaw')

 ####### for whole distance ######
 plt.figure(figsize=(15, 15))
 plt.subplot(311)
 plt.hist(np.log10(distance[iterate_no]), bins=100)
 # plt.title('histgram of interval'+str(iterate_no))

 y_raw = [x for x in distance[iterate_no] if (math.isnan(x) == False) & (x <= threshold)]
 y_raw_n = [x for x in distance[iterate_no] if (math.isnan(x) == False) & (x <= 1200)]
 y_rawl = [x for x in distance[iterate_no] if (math.isnan(x) == False) & (x > threshold)]
 y_rawl_n = [x for x in distance[iterate_no] if (math.isnan(x) == False) & (x > 1200)]

 ####### for short distance #######
 plt.subplot(323)
 plt.scatter(np.log10(x_), np.log10(y_))
 if len(y_raw) == 0:
  print("No data at this time sliding under the case of short distance with noise!")
 else:
  param_exp = dist_exp.fit(y_raw)
  D_exp, p_exp = stats.kstest(y_raw, 'expon', args=param_exp)
  param_law = dist_powerlaw.fit(y_raw)
  D_law, p_law = stats.kstest(y_raw, 'powerlaw', args=param_law)

  print('for short distance range:')
  print('exponential dist with D value equal', D_exp)
  print('exponential dist with p value equal', p_exp)
  print('powerlaw dist with D value equal', D_law)
  print('powerlaw dist with p value equal', p_law)
 print('---------------------------------------')

 plt.subplot(324)
 plt.scatter(np.log10(x_n), np.log10(y_n))
 if len(y_raw_n) == 0:
  print("No data at this time sliding under the case of short distance without noise!")
 else:
  param_exp = dist_exp.fit(y_raw_n)
  D_exp, p_exp = stats.kstest(y_raw_n, 'expon', args=param_exp)
  param_law = dist_powerlaw.fit(y_raw_n)
  D_law, p_law = stats.kstest(y_raw_n, 'powerlaw', args=param_law)

  print('delete noise for short distance range:')
  print('exponential dist with D value equal', D_exp)
  print('exponential dist with p value equal', p_exp)
  print('powerlaw dist with D value equal', D_law)
  print('powerlaw dist with p value equal', p_law)
 print('---------------------------------------')

 ####### for long distance #######
 plt.subplot(325)
 plt.scatter(np.log10(x_l), np.log10(y_l))
 if len(y_rawl) == 0:
  print('No data at this time sliding under the case of long distance with noise!')
 else:
  param_exp = dist_exp.fit(y_rawl)
  D_exp, p_exp = stats.kstest(y_rawl, 'expon', args=param_exp)
  param_law = dist_powerlaw.fit(y_rawl)
  D_law, p_law = stats.kstest(y_rawl, 'powerlaw', args=param_law)

  print('for long distance range:')
  print('exponential dist with D value equal', D_exp)
  print('exponential dist with p value equal', p_exp)
  print('powerlaw dist with D value equal', D_law)
  print('powerlaw dist with p value equal', p_law)
 print('---------------------------------------')

 plt.subplot(326)
 plt.scatter(np.log10(x_nl), np.log10(y_nl))
 if len(y_rawl_n)==0:
  print('No data at this time sliding under the case of long distance without noise!')
 else:
  param_exp = dist_exp.fit(y_rawl_n)
  D_exp, p_exp = stats.kstest(y_rawl_n, 'expon', args=param_exp)
  param_law = dist_powerlaw.fit(y_rawl_n)
  D_law, p_law = stats.kstest(y_rawl_n, 'powerlaw', args=param_law)

  print('delete noise for long distance range:')
  print('exponential dist with D value equal', D_exp)
  print('exponential dist with p value equal', p_exp)
  print('powerlaw dist with D value equal', D_law)
  print('powerlaw dist with p value equal', p_law)
 print('---------------------------------------')
 plt.show()

def lambda_sliding(df_1, chromsome, c_start, c_end, win_size, step_len, df_basket=pd.DataFrame(), autosomal=True):
 chromsome_list = df_1[chromsome].unique().tolist()

 chromsome_list = list(set(list(map(str, chromsome_list))))
 if 'MT' in chromsome_list:
  chromsome_list.remove('MT')

 if 'Y' in chromsome_list:
  chromsome_list.remove('Y')

 if autosomal == False:
  chromsome_list.remove('X')
 chr_dict = {k: [] for k in chromsome_list}

 total_len = df_1[c_end].max() - df_1[c_start].min() # may be too long for chromosome based solution
 start_loc = df_1[c_start].min()
 sliding_times = math.ceil((total_len - win_size) / step_len)
 sliding_times = max(sliding_times, 29)
 l = locals()
 distance_list_all = []

 for i in range(sliding_times):
  for c in chromsome_list:
   df = df_1[df_1[chromsome] == c]
   df = df.sort_values(by=[c_start]).reset_index()

   win_head = df[df[c_start] == df[c_end][df[c_start] >= (start_loc + (i * step_len))].min()].index
   win_tail = df[df[c_end] == df[c_end][df[c_end] < start_loc + win_size + i * step_len].max()].index

   if (len(win_head) == 0) or (len(win_tail) == 0):

    distance_list = []

   else:
    if len(win_head) != 1:
     win_head = pd.Index([win_head[0]])
    if len(win_tail) != 1:
     win_tail = pd.Index([win_tail[0]])

    distance_list =list(np.sort(df.loc[list(win_head)[0]:list(win_tail)[0], c_end].diff().dropna().astype(int)))

   chr_dict[c] = chr_dict[c] + [distance_list]

   if len(distance_list) == 0:

    l['win_df_' + str(i) + 'chr' + c] = pd.DataFrame(columns=["distance_" + str(i) + 'chr' + c, 'counts_' + str(i) +'chr' + c])
   else:
    l['win_df_' + str(i) + 'chr' + c] = pd.DataFrame(distance_list).dropna().astype(int)[0].value_counts()
    l['win_df_' + str(i) + 'chr' + c] = l['win_df_' + str(i) + 'chr' + c].reset_index()
    l['win_df_' + str(i) + 'chr' + c].rename(columns={'index': "distance_" + str(i) + 'chr' + c, 0: 'counts_' + str(i) +'chr' + c}, inplace=True)

   df_basket = pd.concat([df_basket, l['win_df_' + str(i) + 'chr' + c]], axis=1)

  distance_list_all = distance_list_all + [distance_list]

 return chr_dict, df_basket

def maxlikelihood(distance, iterate_no, value_range):
 likelihood_dict = {}
 for i in value_range:
  short_set = [x for x in distance[iterate_no] if (math.isnan(x) == False) & (x <= i)]
  long_set = [x for x in distance[iterate_no] if (math.isnan(x) == False) & (x > i)]

  dist_exp = getattr(stats, 'expon')
  dist_powerlaw = getattr(stats, 'powerlaw')

  ##### for short part #####
  if len(short_set)==0:
   indicator = 2
   probility_s=0
  else:
   param_exp = dist_exp.fit(short_set)
   D_exp, p_exp = stats.kstest(short_set, 'expon', args=param_exp)
   param_law = dist_powerlaw.fit(short_set)
   D_law, p_law = stats.kstest(short_set, 'powerlaw', args=param_law)
   indicator = list([D_exp,D_law]).index(min(D_exp, D_law))

   if indicator == 0:### expon
    prob_short = stats.expon.pdf(short_set, param_exp[0], param_exp[1])

   if indicator == 1:### powerlaw
    prob_short = stats.powerlaw.pdf(short_set, param_law[0], param_law[1], param_law[2])

   probility_s = sum(np.log(prob_short)+[0])

  ##### for long part ######
  if len(long_set)==0:
   indicator = 2
   probility_l = 0
  else:
   param_exp = dist_exp.fit(long_set)
   D_exp, p_exp = stats.kstest(long_set, 'expon', args=param_exp)
   param_law = dist_powerlaw.fit(long_set)
   D_law, p_law = stats.kstest(long_set, 'powerlaw', args=param_law)
   indicator = list([D_exp,D_law]).index(min(D_exp, D_law))

   if indicator == 0:### expon
    prob_long = stats.expon.pdf(long_set, param_exp[0], param_exp[1])
   if indicator == 1:### powerlaw
    prob_long = stats.powerlaw.pdf(long_set, param_law[0], param_law[1], param_law[2])

   probility_l = sum(np.log(prob_long)+[0])

  prob = probility_s + probility_l

  if prob == 0:
   likelihood_dict[str(i)] = np.nan
  else:
   likelihood_dict[str(i)] = prob
 threshold = max(likelihood_dict, key=likelihood_dict.get)
 return threshold, likelihood_dict

def expon_function(lambda_, x, b):
 return -lambda_ * np.exp(-lambda_ * np.array(x)) + b

def linear_function(k, x, b):
 y = k * x + b
 return y

def del_zero(x_n, y_n):
 zero_number = y_n.count(0)
 for j in range(zero_number):
  zero_loc = y_n.index(0)
  del y_n[zero_loc]
  del x_n[zero_loc]

 return x_n, y_n

def fit_function(df, distance,c, visualization_1=False, visualization_2=False, threshold_list=[]):
 sliding_time = int(len(distance))
 dist_exp = getattr(stats, 'expon')
 dist_powerlaw = getattr(stats, 'powerlaw')
 coefficient_dict = {}
 coefficient_dict_new = {}

 for i in range(sliding_time):
  threshold = threshold_list[i]
  length = threshold/100
  x_, y_, x_n, y_n = hist_creat(df, i, c, 'short', threshold, length)
  x_l, y_l, x_nl, y_nl = hist_creat(df, i, c, 'long', threshold, length)

  y_raw = [x for x in distance[i] if (math.isnan(x) == False) & (x <= threshold)]
  y_rawl = [x for x in distance[i] if (math.isnan(x) == False) & (x > threshold)]

  x_n, y_n = del_zero(x_n, y_n)
  x_nl, y_nl = del_zero(x_nl, y_nl)
  x_, y_ = del_zero(x_, y_)
  x_l, y_l = del_zero(x_l, y_l)

  ### for short distance ###
  short_dict = {}
  short_dict_new = {}
  if len(y_raw)==0:
   indicator = 2
  else:
   param_exp = dist_exp.fit(y_raw)
   D_exp, p_exp = stats.kstest(y_raw, 'expon', args=param_exp)
   param_law = dist_powerlaw.fit(y_raw)
   D_law, p_law = stats.kstest(y_raw, 'powerlaw', args=param_law)
   indicator = list([D_exp, D_law]).index(min(D_exp, D_law))

  if indicator == 0:  ### expon
   if len(y_) <=1:
    short_dict['expon'] = ([np.nan, np.nan], np.nan)
    short_dict_new['expon'] = ([np.nan, np.nan], np.nan)
   else:
    params_ = param_exp
    short_dict['expon'] = (params_, np.nan)
    short_dict_new['expon'] = (params_, np.nan)


  if indicator == 1: ### power-law
   if len(y_) <=1:
    short_dict['powerlaw'] = ([np.nan, np.nan], np.nan)
    short_dict_new['powerlaw'] = ([np.nan, np.nan], np.nan)
   else:
    params_ = param_law
    short_dict['powerlaw'] = (params_, np.nan)
    fit_short = powerlaw.Fit(y_raw)
    params_new = (fit_short.power_law.alpha, fit_short.power_law.sigma)
    short_dict_new['powerlaw'] = (params_new, np.nan)

  if indicator == 2:
   short_dict['empty'] = ([np.nan, np.nan], np.nan)
   short_dict_new['empty'] = ([np.nan, np.nan], np.nan)

  if visualization_1 == True:
   plt.figure()
   if indicator == 0:
    plt.plot(x_n, np.log10(y_n), '.', label="data")
    plt.plot(x_n, linear_function(params_[0], np.array(x_n), params_[1]), '--', label="fitted")
    # plt.title("Fitted Exponential Curve of sliding time_"+str(i)+'in short distance')

   if indicator == 1:
    plt.plot(np.log10(x_n), np.log10(y_n), '.', label="data")
    plt.plot(np.log10(x_n), linear_function(params_[0], np.array(np.log10(x_n)), params_[1]), '--', label="fitted")
    # plt.title("Fitted Powerlaw Curve of sliding time_"+str(i)+'in short distance')
   plt.show()

  ### for long distance ###
  long_dict = {}
  long_dict_new = {}
  if len(y_rawl)==0:
   indicator = 2
  else:
   param_exp = dist_exp.fit(y_rawl)
   D_exp, p_exp = stats.kstest(y_rawl, 'expon', args=param_exp)
   param_law = dist_powerlaw.fit(y_rawl)
   D_law, p_law = stats.kstest(y_rawl, 'powerlaw', args=param_law)
   indicator = list([D_exp, D_law]).index(min(D_exp, D_law))

  if indicator == 0:  ### expon
   if len(y_l) <=1:
    long_dict['expon'] = ([np.nan, np.nan], np.nan)
    long_dict_new['expon'] = ([np.nan, np.nan], np.nan)
   else:
    params_ = param_exp
    long_dict['expon'] = (params_, np.nan)
    long_dict_new['expon'] = (params_, np.nan)

  if indicator == 1:  ### powerlaw
   if len(y_l) <=1:
    long_dict['powerlaw'] = ([np.nan, np.nan], np.nan)
    long_dict_new['powerlaw'] = ([np.nan, np.nan], np.nan)
   else:
    params_ = param_law
    long_dict['powerlaw'] = (params_, np.nan)
    fit_long = powerlaw.Fit(y_rawl)
    params_new = (fit_long.power_law.alpha, fit_long.power_law.sigma)
    long_dict_new['powerlaw'] = (params_new, np.nan)

  if indicator == 2:
   long_dict['empty'] = ([np.nan, np.nan], np.nan)
   long_dict_new['empty'] = ([np.nan, np.nan], np.nan)

  if visualization_2 == True:
   plt.figure()
   if indicator == 0:
    plt.plot(x_nl, np.log10(y_nl), '.', label="data")
    plt.plot(x_nl, linear_function(params_[0], np.array(x_nl), params_[1]), '--', label="fitted")
    plt.title("Fitted Exponential Curve of sliding time_" + str(i) + 'in long distance')

   if indicator == 1:
    plt.plot(np.log10(x_nl), np.log10(y_nl), '.', label="data")
    plt.plot(np.log10(x_nl), linear_function(params_[0], np.array(np.log10(x_nl)), params_[1]), '--', label="fitted")
    plt.title("Fitted Powerlaw Curve of sliding time_" + str(i) + 'in long distance')
   plt.show()

  coefficient_dict[i] = (short_dict, long_dict)
  coefficient_dict_new[i] = (short_dict_new, long_dict_new)

 return coefficient_dict, coefficient_dict_new

def parameter_conclude_1(coefficient_dict):
 short_exp = []
 long_law = []
 short_law = []
 long_exp = []

 amount = len(coefficient_dict.keys())
 for i in range(amount):
  ### for short distance ###
  if len(coefficient_dict[i][0])!=0:
   if list(coefficient_dict[i][0].keys())[0] == 'expon':
    short_exp = short_exp + [(i, coefficient_dict[i][0]['expon'][0][0])]

   if list(coefficient_dict[i][0].keys())[0] == 'powerlaw':
    short_law = short_law + [(i, coefficient_dict[i][0]['powerlaw'][0][0])]

  ### for long distance ###
  if len(coefficient_dict[i][1])!=0:
   if list(coefficient_dict[i][1].keys())[0] == 'expon':
    long_exp = long_exp + [(i, coefficient_dict[i][1]['expon'][0][0])]

   if list(coefficient_dict[i][1].keys())[0] == 'powerlaw':
    long_law = long_law + [(i, coefficient_dict[i][1]['powerlaw'][0][0])]

 return short_exp, short_law, long_exp, long_law

def parameter_conclude_2(short_exp, short_law, long_exp, long_law, patient_id, c, path, plot_together=True, plot_withDistType=False,
                         plot_withDistance=False):
 if len(short_exp) == 0:
  x_e_corrd = []
  y_e_corrd = []

 if len(short_exp) != 0:
  x_e_corrd = list(zip(*short_exp))[0]
  y_e_corrd = list(zip(*short_exp))[1]

 if len(short_law) == 0:
  x_p_corrd = []
  y_p_corrd = []

 if len(short_law) != 0:
  x_p_corrd = list(zip(*short_law))[0]
  y_p_corrd = list(zip(*short_law))[1]

 if len(long_exp) == 0:
  x_le_corrd = []
  y_le_corrd = []

 if len(long_exp) != 0:
  x_le_corrd = list(zip(*long_exp))[0]
  y_le_corrd = list(zip(*long_exp))[1]

 if len(long_law) == 0:
  x_lp_corrd = []
  y_lp_corrd = []

 if len(long_law) != 0:
  x_lp_corrd = list(zip(*long_law))[0]
  y_lp_corrd = list(zip(*long_law))[1]

 if plot_together == True:
  plt.figure()
  plt.scatter(x_e_corrd, y_e_corrd, marker='o', c='#ff7f0e', label='exp in short')
  plt.scatter(x_p_corrd, y_p_corrd, marker='^', c='#ff7f0e', label='powerlaw in short')
  plt.scatter(x_le_corrd, y_le_corrd, marker="o", c='#2ca02c', label='exp in long')
  plt.scatter(x_lp_corrd, y_lp_corrd, marker="^", c='#2ca02c', label='powerlaw in long')

  plt.legend()
  plt.show()

 if plot_withDistType == True:
  plt.figure()
  plt.scatter(x_e_corrd, y_e_corrd, marker='o', c='#ff7f0e', label='exp in short')
  plt.scatter(x_le_corrd, y_le_corrd, marker="o", c='#2ca02c', label='exp in long')
  plt.legend()
  plt.show()

  plt.scatter(x_p_corrd, y_p_corrd, marker='^', c='#ff7f0e', label='powerlaw in short')
  plt.scatter(x_lp_corrd, y_lp_corrd, marker="^", c='#2ca02c', label='powerlaw in long')
  plt.legend()
  plt.show()

 if plot_withDistance == True:
  plt.figure()
  plt.scatter(x_e_corrd, y_e_corrd, marker='o', c='#ff7f0e', label='exp in short')
  plt.scatter(x_p_corrd, y_p_corrd, marker='^', c='#ff7f0e', label='powerlaw in short')
  plt.legend()
  plt.title("short_group{}".format(patient_id))
  plt.savefig(path+'/short_group_{}.jpg'.format(c))
  #plt.show()
  plt.close()

  plt.scatter(x_le_corrd, y_le_corrd, marker="o", c='#2ca02c', label='exp in long')
  plt.scatter(x_lp_corrd, y_lp_corrd, marker="^", c='#2ca02c', label='powerlaw in long')
  plt.legend()
  plt.title("long_group{}".format(patient_id))
  plt.savefig(path+'/long_group_{}.jpg'.format(c))
  #plt.show()
  plt.close()

def fit_function_KL(df,distance, visualization_1=False, visualization_2=False, threshold=5000, length=50):
 sliding_time = int(df.shape[1] / 2)
 dist_exp = getattr(stats, 'expon')
 dist_powerlaw = getattr(stats, 'powerlaw')
 coefficient_dict = {}

 for i in range(sliding_time):
  x_, y_, x_n, y_n = hist_creat(df, i, 'short', threshold, length)
  x_l, y_l, x_nl, y_nl = hist_creat(df, i, 'long', threshold, length)

  y_raw = [x for x in distance[i] if (math.isnan(x) == False) & (x <= threshold)]
  y_rawl = [x for x in distance[i] if (math.isnan(x) == False) & (x > threshold)]

  sind = df[df['distance_' + str(i)] <= threshold].index
  lind = df[df['distance_' + str(i)] > threshold].index

  zero_number = y_n.count(0)
  zero_number_n = y_nl.count(0)

  for j in range(zero_number):
   zero_loc = y_n.index(0)
   del y_n[zero_loc]
   del x_n[zero_loc]

  for k in range(zero_number_n):
   zero_loc = y_nl.index(0)
   del y_nl[zero_loc]
   del x_nl[zero_loc]

  ### for short distance ###
  short_dict = {}

  param_exp = dist_exp.fit(y_raw)
  #D_exp, p_exp = stats.kstest(y_raw, 'expon', args=param_exp)
  p = list(df.loc[sind,'counts_' + str(i)]/df.loc[sind,'counts_' + str(i)].sum())
  q = stats.expon.pdf(list(df.loc[sind,'distance_' + str(i)]), scale=param_exp[1], loc=param_exp[0])
  print('P real value of short expon is:', p)
  print("Q estimated value of short expon is:", q)

  KLDSE = stats.entropy(p, q)
  print('short expon KLD value is:', KLDSE)

  param_law = dist_powerlaw.fit(y_raw)
  #D_law, p_law = stats.kstest(y_raw, 'powerlaw', args=param_law)

  q = stats.powerlaw.pdf(list(df.loc[sind,'distance_' + str(i)]), a=param_law[0],scale=param_law[2], loc=param_law[1])
  print("Q estimate value of short powerlaw is:",q)

  KLDSP = stats.entropy(p, q)
  print('short powerlaw KLD values is:',KLDSP)

  indicator = list([KLDSE, KLDSP]).index(min(KLDSE, KLDSP))


  if indicator == 0:  ### expon
   params_, cv_ = scipy.optimize.curve_fit(linear_function, x_n, np.log10(y_n))
   squaredDiffs = np.square(np.log10(y_n) - linear_function(params_[0], np.array(x_n), params_[1]))
   squaredDiffsFromMean = np.square(np.log10(y_n) - np.mean(np.log10(y_n)))
   rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
   print(f"R² = {rSquared}")
   short_dict['expon'] = (params_, rSquared)

  if indicator == 1:
   params_, cv_ = scipy.optimize.curve_fit(linear_function, np.log10(x_n), np.log10(y_n))
   squaredDiffs = np.square(np.log10(y_n) - linear_function(params_[0], np.array(np.log10(x_n)), params_[1]))
   squaredDiffsFromMean = np.square(np.log10(y_n) - np.mean(np.log10(y_n)))
   rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
   print(f"R² = {rSquared}")
   short_dict['powerlaw'] = (params_, rSquared)

  if visualization_1 == True:
   plt.figure()
   if indicator == 0:
    plt.plot(x_n, np.log10(y_n), '.', label="data")
    plt.plot(x_n, linear_function(params_[0], np.array(x_n), params_[1]), '--', label="fitted")
    # plt.title("Fitted Exponential Curve of sliding time_"+str(i)+'in short distance')

   if indicator == 1:
    plt.plot(np.log10(x_n), np.log10(y_n), '.', label="data")
    plt.plot(np.log10(x_n), linear_function(params_[0], np.array(np.log10(x_n)), params_[1]), '--', label="fitted")
    # plt.title("Fitted Powerlaw Curve of sliding time_"+str(i)+'in short distance')
   plt.show()

  ### for long distance ###
  long_dict = {}

  param_exp = dist_exp.fit(y_rawl)
  #D_exp, p_exp = stats.kstest(y_rawl, 'expon', args=param_exp)
  pl = list(df.loc[lind,'counts_' + str(i)]/df.loc[lind,'counts_' + str(i)].sum())
  ql = stats.expon.pdf(list(df.loc[lind,'distance_' + str(i)]), scale=param_exp[1], loc=param_exp[0])
  print("P real value of long expon is:",pl)
  print("Q estimated value of long expon is:",ql)

  KLDLE = stats.entropy(pl, ql)
  print('long expon KLD value is:',KLDLE)

  param_law = dist_powerlaw.fit(y_rawl)
  #D_law, p_law = stats.kstest(y_rawl, 'powerlaw', args=param_law)

  ql = stats.powerlaw.pdf(list(df.loc[lind,'distance_' + str(i)]), a=param_law[0],scale=param_law[2], loc=param_law[1])
  print("Q real value of long powerlaw is:",ql)

  KLDLP = stats.entropy(pl, ql)
  print('long powerlaw KLD value is:', KLDLP)

  indicator = list([KLDLE, KLDLP]).index(min(KLDLE, KLDLP))

  if indicator == 0:  ### expon
   params_, cv_ = scipy.optimize.curve_fit(linear_function, x_nl, np.log10(y_nl))
   squaredDiffs = np.square(np.log10(y_nl) - linear_function(params_[0], np.array(x_nl), params_[1]))
   squaredDiffsFromMean = np.square(np.log10(y_nl) - np.mean(np.log10(y_nl)))
   rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
   print(f"R² = {rSquared}")
   long_dict['expon'] = (params_, rSquared)

  if indicator == 1:  ### powerlaw
   params_, cv_ = scipy.optimize.curve_fit(linear_function, np.log10(x_nl), np.log10(y_nl))
   squaredDiffs = np.square(np.log10(y_nl) - linear_function(params_[0], np.array(np.log10(x_nl)), params_[1]))
   squaredDiffsFromMean = np.square(np.log10(y_nl) - np.mean(np.log10(y_nl)))
   rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
   print(f"R² = {rSquared}")
   long_dict['powerlaw'] = (params_, rSquared)

  if visualization_2 == True:
   plt.figure()
   if indicator == 0:
    plt.plot(x_nl, np.log10(y_nl), '.', label="data")
    plt.plot(x_nl, linear_function(params_[0], np.array(x_nl), params_[1]), '--', label="fitted")
   # plt.title("Fitted Exponential Curve of sliding time_" + str(i) + 'in long distance')

   if indicator == 1:
    plt.plot(np.log10(x_nl), np.log10(y_nl), '.', label="data")
    plt.plot(np.log10(x_nl), linear_function(params_[0], np.array(np.log10(x_nl)), params_[1]), '--', label="fitted")
    # plt.title("Fitted Powerlaw Curve of sliding time_" + str(i) + 'in long distance')
   plt.show()

  coefficient_dict[i] = (short_dict, long_dict)

 return coefficient_dict



