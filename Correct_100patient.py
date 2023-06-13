import pandas as pd
import numpy as np
import scipy
from matplotlib import pyplot as plt
from scipy import stats
import math
import scipy.optimize
from Correct_Funtion_set import *
import tqdm
import sys
import os
import pickle


def get_results(df_breast, cancer_type, k=0, j=50):
    print('>>>>>>>>> code start!')
    patient_list = df_breast[0].value_counts()[k:j].index.tolist()
    #patient_list = ['PD4115a', 'PD4847a']
    counter = len(patient_list)
    if not os.path.exists('Correct_'+cancer_type):
        os.makedirs('Correct_'+cancer_type)

    PLVD = {}### create dictionary to save the predicted lambda values

    for patient_id in patient_list:
        print("now is the patient_{}!".format(patient_id))
        if not os.path.exists('Correct_' + cancer_type + '/{}'.format(patient_id)):
            os.makedirs('Correct_' + cancer_type + '/{}'.format(patient_id))
        #os.makedirs('Correct_' + cancer_type+'/{}'.format(patient_id))
        path = 'Correct_' + cancer_type+'/{}'.format(patient_id)
        path_2 = 'Correct_' + cancer_type

        df_1 = df_breast[(df_breast[0]==patient_id) & (df_breast[1]=='subs')]
        df_sort = df_1.sort_values(by=[3])
        distance, df_bucket = lambda_sliding(df_1,2,3,4,50000000,7000000,autosomal=True)
        threshold_list = {k: [] for k in distance.keys()}
        for c in tqdm.tqdm(distance.keys()):
            for i in range(0, 29):
                threshold_can, dict_can = maxlikelihood(distance[c], i, range(3000, 10000, 100))
                threshold_list[c] = threshold_list[c] + [int(threshold_can)]

        dic_chr = {}
        for c in distance.keys():
            filter_col = [s for s in list(df_bucket.columns) if s.endswith('chr' + c)]
            dic_chr['chr' + c] = fit_function(df_bucket[filter_col], distance[c], c, False, False, threshold_list[c])

        ### plot distribution exchange ###
        D = locals()
        for c in distance.keys():
            D['dict_chr' + c] = {}
            se_, sp_, le_, lp_ = parameter_conclude_1(dic_chr['chr' + c])
            parameter_conclude_2(se_, sp_, le_, lp_, patient_id, c, path, False, False, plot_withDistance=True)
            D['dict_chr' + c]['shrot_exp'] = se_
            D['dict_chr' + c]['shrot_pow'] = sp_
            D['dict_chr' + c]['long_exp'] = le_
            D['dict_chr' + c]['long_pow'] = lp_


        for c in distance.keys():
            ### short expon ###
            try:
                se_ = D['dict_chr' + c]['shrot_exp']
                plt.plot(list(zip(*se_))[0], -1 * np.array(list(zip(*se_))[1]))
                plt.scatter(list(zip(*se_))[0], -1 * np.array(list(zip(*se_))[1]), color='palevioletred')
                plt.title("short_expon_{}".format(patient_id))
                plt.savefig(path+'/shot_expon_{}.jpg'.format(c))
                plt.close()

            except IndexError:
                pass

            ### short powerlaw ###
            try:
                sp_ = D['dict_chr' + c]['shrot_pow']
                plt.plot(list(zip(*sp_))[0], -1 * np.array(list(zip(*sp_))[1]))
                plt.scatter(list(zip(*sp_))[0], -1 * np.array(list(zip(*sp_))[1]), color='palevioletred')
                plt.title("short_powerlaw_{}".format(patient_id))
                plt.savefig(path+'/shot_powerlaw_{}.jpg'.format(c))
                plt.close()
            except IndexError:
                pass

            ### long expon ###
            try:
                le_ = D['dict_chr' + c]['long_exp']
                plt.plot(list(zip(*le_))[0], -1 * np.array(list(zip(*le_))[1]))
                plt.scatter(list(zip(*le_))[0], -1 * np.array(list(zip(*le_))[1]), color='palevioletred')
                plt.title("long_expon_{}".format(patient_id))
                plt.savefig(path+'/long_expon_{}.jpg'.format(c))
                plt.close()
            except IndexError:
                pass

            ### long powerlaw ###
            try:
                lp_ = D['dict_chr' + c]['long_pow']
                plt.plot(list(zip(*lp_))[0], -1 * np.array(list(zip(*lp_))[1]))
                plt.scatter(list(zip(*lp_))[0], -1 * np.array(list(zip(*lp_))[1]), color='palevioletred')
                plt.title("long_powerlaw_{}".format(patient_id))
                plt.savefig(path+'/long_powerlaw_{}.jpg'.format(c))
                plt.close()
            except IndexError:
                pass

        ### add each chromosome together to plot:
        P = locals()
        pos_chr = {}
        max_length = 0
        for c in distance.keys():
            pos_chr[c] = {}
            pos_chr[c]['min'] = df_1[3][df_1[2] == c].min()
            pos_chr[c]['max'] = df_1[3][df_1[2] == c].max()
            pos_chr[c]['mlength'] = pos_chr[c]['max'] - pos_chr[c]['min']
            if pos_chr[c]['max'] > max_length:
                max_length = pos_chr[c]['max']

        df_empty = pd.DataFrame()
        for i in df_1[2].unique():
            df_kk = df_1[df_1[2] == i].sort_values(by=[3])
            df_kk['distance'] = df_1[df_1[2] == i].sort_values(by=[3])[3].diff()

            df_empty = pd.concat([df_empty, df_kk], axis=0)
        df_k = df_empty

        ### normalized each chromosome to make them plot in sequnence
        l = [str(a) for a in range(1, 23)] + ['X']
        l_chr = ['chr_' + x for x in l]
        chr_bound = []
        max_val = 0

        new_l = []
        for i in l:
            if i in distance.keys():
                new_l = new_l + [i]
        new_l_chr = ['chr_' + x for x in new_l]

        for c in new_l:
            val = df_k[3][df_k[2] == c]
            idx = df_k[df_k[2] == c].index
            df_k.loc[idx, 3] = val - pos_chr[c]['min'] + max_val
            max_val = pos_chr[c]['max'] - pos_chr[c]['min'] + max_val
            chr_bound = chr_bound + [max_val]

        ######## set each mutation type with particular color #######
        ### C>A
        x_bule = df_k[(df_k[6]=='C')&(df_k[5]=='A')][3]
        y_bule = df_k[(df_k[6]=='C')&(df_k[5]=='A')]['distance']

        x_bule_1 = df_k[(df_k[6]=='G')&(df_k[5]=='T')][3]
        y_bule_1 = df_k[(df_k[6]=='G')&(df_k[5]=='T')]['distance']
        ### C>G
        x_black = df_k[(df_k[6]=='C')&(df_k[5]=='G')][3]
        y_black = df_k[(df_k[6]=='C')&(df_k[5]=='G')]['distance']

        x_black_1 = df_k[(df_k[6]=='G')&(df_k[5]=='C')][3]
        y_black_1 = df_k[(df_k[6]=='G')&(df_k[5]=='C')]['distance']
        ### C>T
        x_red = df_k[(df_k[6]=='C')&(df_k[5]=='T')][3]
        y_red = df_k[(df_k[6]=='C')&(df_k[5]=='T')]['distance']

        x_red_1 = df_k[(df_k[6]=='G')&(df_k[5]=='A')][3]
        y_red_1 = df_k[(df_k[6]=='G')&(df_k[5]=='A')]['distance']
        ### T>A
        x_puple = df_k[(df_k[6]=='T')&(df_k[5]=='A')][3]
        y_puple = df_k[(df_k[6]=='T')&(df_k[5]=='A')]['distance']

        x_puple_1 = df_k[(df_k[6]=='A')&(df_k[5]=='T')][3]
        y_puple_1 = df_k[(df_k[6]=='A')&(df_k[5]=='T')]['distance']
        ### T>C
        x_yellow = df_k[(df_k[6]=='T')&(df_k[5]=='C')][3]
        y_yellow = df_k[(df_k[6]=='T')&(df_k[5]=='C')]['distance']

        x_yellow_1 = df_k[(df_k[6]=='A')&(df_k[5]=='G')][3]
        y_yellow_1 = df_k[(df_k[6]=='A')&(df_k[5]=='G')]['distance']
        ### T>G
        x_green = df_k[(df_k[6]=='T')&(df_k[5]=='G')][3]
        y_green = df_k[(df_k[6]=='T')&(df_k[5]=='G')]['distance']

        x_green_1 = df_k[(df_k[6]=='A')&(df_k[5]=='C')][3]
        y_green_1 = df_k[(df_k[6]=='A')&(df_k[5]=='C')]['distance']


        plt.figure(figsize=(15,6))
        plt.scatter(x_bule,y_bule,c='blue',label='C>A')
        plt.scatter(x_black,y_black,c='black',label='C>G')
        plt.scatter(x_red,y_red,c='red',label='C>T')
        plt.scatter(x_puple,y_puple,c='purple',label='T>A')
        plt.scatter(x_yellow,y_yellow,c='yellow',label='T>C')
        plt.scatter(x_green,y_green,c='lime',label='T>G')

        plt.scatter(x_bule_1,y_bule_1,c='royalblue',label='G>T')
        plt.scatter(x_black_1,y_black_1,c='dimgray',label='G>C')
        plt.scatter(x_red_1,y_red_1,c='firebrick',label='G>A')
        plt.scatter(x_puple_1,y_puple_1,c='fuchsia',label='A>T')
        plt.scatter(x_yellow_1,y_yellow_1,c='gold',label='A>G')
        plt.scatter(x_green_1,y_green_1,c='limegreen',label='A>C')


        plt.vlines(x=chr_bound, ymax=10 ** 6, ymin=0, linestyles='dashed', colors='lightgray')
        plt.yscale('log')
        plt.yticks([1, 10 ** 1, 10 ** 2, 10 ** 4, 10 ** 6])
        plt.xticks(chr_bound, new_l_chr, rotation=45)
        plt.xlabel('Genomic position')
        plt.legend(bbox_to_anchor=(0.5, -0.2), loc='center', ncol=12, edgecolor='w')
        plt.title("patient_{}".format(patient_id))
        plt.savefig(path_2+'/Kataegis_{}.jpg'.format(patient_id))
        plt.close()
        print("patient_{} is finished!".format(patient_id))
        print("---------- {} patients left ----------".format(counter-1))
        counter = counter-1

        PLVD[patient_id] = dic_chr
    f = open('/file_{}.pkl'.format(cancer_type), 'wb')
    pickle.dump(PLVD, f)
    f.close()


if __name__ == "__main__":

    #df_breast = pd.read_csv('/Users/mhy/Desktop/Project/AlexandrovEtAI/somatic_mutation_data/Lymphoma B-cell/Lymphoma B-cell_raw_mutations_data.txt',\
    #                    delimiter = "\t", header=None)
    #df_breast = pd.read_csv('/Users/mhy/Desktop/Project/AlexandrovEtAI/somatic_mutation_data/Lung Squamous/Lung Squamous_raw_mutations_data.txt', \
    #    delimiter="\t", header=None)

    #df_breast = pd.read_csv('/Users/mhy/Desktop/Project/AlexandrovEtAI/somatic_mutation_data/Melanoma/Melanoma_raw_mutations_data.txt', \
    #    delimiter="\t", header=None)

    #df_breast = pd.read_csv('/Users/mhy/Desktop/Project/AlexandrovEtAI/somatic_mutation_data/Colorectum/Colorectum_raw_mutations_data.txt',\
    #                    delimiter = "\t", header=None)

    df_breast = pd.read_csv('/Users/mhy/Desktop/Project/AlexandrovEtAI/somatic_mutation_data/Breast/Breast_raw_mutations_data.txt', \
                delimiter="\t", header=None)

    get_results(df_breast, 'Breast', 0, 50)