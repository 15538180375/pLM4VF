import sys, os
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)

def PSSM_composition(AA, proteinSeq, pssmMatrix):
    
    PSSM_composition = []
    pssm_composition = [[0.0 for m in range(20)] for n in range(20)]  
    for k in range(20):
        for i in range(len(proteinSeq)):
            if AA[k] == proteinSeq[i]:
                for j in range(20):
                    pssm_composition[k][j] = pssm_composition[k][j] + pssmMatrix[i][j] / len(pssmMatrix)
    for a in pssm_composition:
        for b in a:
            PSSM_composition.append(b)
    return PSSM_composition


def S_FPSSM(AA, proteinSeq, pssmMatrix):
    # S_FPSSM特征提取
    F_PSSM = pssmMatrix  # 转换F_PSSM矩阵
    for k in range(20):
        for i in range(len(proteinSeq)):
            if F_PSSM[i][k] <= 0.0:
                F_PSSM[i][k] = 0.0
            elif F_PSSM[i][k] >= 7.0:
                F_PSSM[i][k] = 7.0
    S_fpssm = []
    s_fpssm = [[0.0 for m in range(20)] for n in range(20)]
    for k in range(20):
        for i in range(len(proteinSeq)):
            if AA[k] == proteinSeq[i]:
                for j in range(20):
                    s_fpssm[k][j] = s_fpssm[k][j] + F_PSSM[i][j]
    for a in s_fpssm:
        for b in a:
            S_fpssm.append(b)
    return S_fpssm


def RPM_PSSM(AA, proteinSeq, pssmMatrix):
    # RPM_PSSM特征提取
    PPSSM = pssmMatrix  # 将原始PSSM矩阵转换为PPSSM
    RPM_PSSM = []
    rpm_PSSM = [[0.0 for m in range(20)] for n in range(20)]
    for k in range(20):
        for i in range(len(proteinSeq)):
            if PPSSM[i][k] <= 0.0:
                PPSSM[i][k] = 0
    for k in range(20):
        for i in range(len(proteinSeq)):
            if AA[k] == proteinSeq[i]:
                for j in range(20):
                    rpm_PSSM[k][j] = rpm_PSSM[k][j] + PPSSM[i][j] / len(proteinSeq)
    for a in rpm_PSSM:
        for b in a:
            RPM_PSSM.append(b)
    return RPM_PSSM


def PSSM(fastas, pssmDir):

    #if checkFasta.checkFasta(fastas) == False:
        #print('Error: for "PSSM" encoding, the input fasta sequences should be with equal length. \n\n')
        #return 0


    if pssmDir == None:
        print('Error: please specify the directory of predicted protein disorder files by "--path" \n\n')
        return 0

    AA = 'ARNDCQEGHILKMFPSTWYV'

    #encodings = []
    #header = ['#']
    #for p in range(1, len(fastas[0][1]) + 1):
        #for aa in AA:
            #header.append('Pos.' + str(p) + '.' + aa)
    #encodings.append(header)

    encodings = []

    PSSM_composition_all = []  #存储所有序列的pssm_composition特征
    S_FPSSM_all = []
    RPM_PSSM_all = []

    #给每个特征添加名字
    pssm_composition_header = ['#']
    s_fpssm_header = ['#']
    rpm_pssm_header = ['#']
    for i in range(400):
        pssm_composition_header.append('pssm_composition'+ str(i))
        s_fpssm_header.append('s_fpssm' + str(i))
        rpm_pssm_header.append('rpm_pssm' + str(i))
    PSSM_composition_all.append(pssm_composition_header)
    S_FPSSM_all.append(s_fpssm_header)
    RPM_PSSM_all.append(rpm_pssm_header)

    for i in fastas:

        name, sequence = re.sub('\|', '', i[0]), i[1]
        code = [name]

        if os.path.exists(pssmDir+'/'+name+'.pssm') == False:
            pssm_composition1 = []
            s_fpssm1 = []
            rpm_pssm1 = []
            pssm_composition1 = code + list(map(lambda x: 0, range(400)))
            s_fpssm1 = code + list(map(lambda x: 0, range(400)))
            rpm_pssm1 = code + list(map(lambda x: 0, range(400)))
            PSSM_composition_all.append(pssm_composition1)
            S_FPSSM_all.append(s_fpssm1)
            RPM_PSSM_all.append(rpm_pssm1)
            continue
            #print('Error: pssm prfile for protein' + name + 'does not exist.')
            #sys.exit(1)

        with open(pssmDir+'/'+name+'.pssm') as f:
            records = f.readlines()[3:-6]

        proteinSeq = ''
        pssmMatrix = []
        for line in records:
            array = line.strip().split()
            pssmMatrix.append(array[2:22])
            proteinSeq = proteinSeq + array[1]

        for i in range(len(pssmMatrix)):  #
            pssmMatrix[i] = [float(x) for x in pssmMatrix[i]]

        pssm_composition = []  
        s_fpssm = []
        rpm_pssm = []

        pssm_composition = code + PSSM_composition(AA,proteinSeq,pssmMatrix)
        s_fpssm = code + S_FPSSM(AA,proteinSeq,pssmMatrix)
        rpm_pssm = code + RPM_PSSM(AA,proteinSeq,pssmMatrix)
        #pssm_composition.extend(code)
        #pssm_composition.extend(PSSM_composition(AA,proteinSeq,pssmMatrix))
        PSSM_composition_all.append(pssm_composition)
        S_FPSSM_all.append(s_fpssm)
        RPM_PSSM_all.append(rpm_pssm)

    return PSSM_composition_all,S_FPSSM_all,RPM_PSSM_all


        

