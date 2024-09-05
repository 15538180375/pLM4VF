import numpy as np
import pandas as pd
import math
import re,os,sys
import readFasta
import csv

def normalizePSSM(PSSM):
    #PSSM = PSSM[:, 1:21]
    PSSM = PSSM.astype(float)
    PSSM = np.array(PSSM)
    seq_cn = np.shape(PSSM)[0]
    PSSM_norm = [[0.0] * 20] * seq_cn
    PSSM_norm = np.array(PSSM_norm)
    mean_matrix = np.mean(PSSM, axis=1)
    std_matrix = np.std(PSSM, axis=1)

    for i in range(seq_cn):
        for j in range(20):
            if std_matrix[i] == 0.0:
                PSSM_norm[i][j] = PSSM[i][j] - mean_matrix[i]
            else:
                PSSM_norm[i][j] = (PSSM[i][j] - mean_matrix[i]) / std_matrix[i]
    return PSSM_norm


def handleMixed(PSSM, ALPHA):
    row1 = [0.0] * 20
    row2 = [0.0] * 20

    matrix_final = [[0.0] * 40] * 1
    row1 = np.array(row1)
    row2 = np.array(row2)
    matrix_final = np.array(matrix_final)

    PSSM_norm = normalizePSSM(PSSM)
    seq_cn = np.shape(PSSM)[0]
    for i in range(seq_cn):
        # print PSSM_norm[i]
        row1 = list(map(sum, zip(row1, PSSM_norm[i])))
    # print row1
    row1 = np.divide(row1, seq_cn)

    for j in range(20):
        for i in range(seq_cn - ALPHA):
            row2[j] += (PSSM_norm[i][j] - PSSM_norm[i + ALPHA][j]) * (PSSM_norm[i][j] - PSSM_norm[i + ALPHA][j])
    # print row2
    row2 = np.divide(row2, seq_cn - ALPHA)

    row = np.hstack((row1, row2))
    matrix_final[0] = row
    # print np.shape(matrix_final)
    return matrix_final


def pse_pssm(input_matrix, ALPHA):
    # print "start pse_pssm function"
    # ALPHA=1
    pse_pssm_matrix = handleMixed(input_matrix, ALPHA)
    # print "end pse_pssm function"
    return pse_pssm_matrix


def PSSM(fastas, pssmDir):

    if pssmDir == None:
        print('Error: please specify the directory of predicted protein disorder files by "--path" \n\n')
        return 0

    #AA = 'ARNDCQEGHILKMFPSTWYV'

    encodings = []
    PsePSSM_all = []   #存储所有序列的pssm_composition特征

    #给每个特征添加名字
    psepssm_header = ['#']
    for i in range(40):
        psepssm_header.append('pse_pssm' + str(i))
    PsePSSM_all.append(psepssm_header)

    for i in fastas:
        name, sequence = re.sub('\|', '', i[0]), i[1]
        code = [name]

        #if os.path.exists(pssmDir+'/'+name+'.pssm') == False:
            #print('Error: pssm prfile for protein' + name + 'does not exist.')
            #sys.exit(1)
        if os.path.exists(pssmDir + '/' + name + '.pssm') == False:
            psepssm1 = []
            psepssm1 = code + list(map(lambda x: 0, range(40)))
            PsePSSM_all.append(psepssm1)
            continue

        with open(pssmDir+'/'+name+'.pssm') as f:
            records = f.readlines()[3:-6]

        pssmMatrix = []
        for line in records:
            array = line.strip().split()
            pssmMatrix.append(array[2:22])
        PSSM = pd.DataFrame(pssmMatrix)

        psepssm = [] #暂时存储当前序列的psepssm特征，每循环一次更新一次
        psepssm = code + list(pse_pssm(PSSM,1).reshape(-1))
        PsePSSM_all.append(psepssm)

    return PsePSSM_all


def pssm_output(list,pssm_filepath):
    f = open(pssm_filepath,'w',newline='')
    writer = csv.writer(f)
    for i in list:
        writer.writerow(i)
    f.close()

