from django.core.files import temp
from django.shortcuts import render
from django.conf import settings
from django.http import HttpResponse
from django.core.files.storage import FileSystemStorage
from django.http import FileResponse
from django.views.static import serve
import xlsxwriter
import pdfkit
import csv
import numpy
#import required libraries
import pandas as pd
import pyexcel
import xlrd
from matplotlib import pylab
from matplotlib import collections  as mc
from pylab import *
from pylev3 import Levenshtein
from matplotlib.ticker import PercentFormatter
from matplotlib import pyplot
import matplotlib.pyplot as plt
import PIL, PIL.Image
import uuid, os

try:
    from StringIO import BytesIO
except ImportError:
    from io import BytesIO

DIR = ""
def generateDirname():
    '''
    Generates an unique name for the temporary directory
    '''
    dirname = uuid.uuid4().hex
    return dirname


'''from google.colab import drive
drive.mount('/content/drive')'''

# Create your views here.

def welcome(request):
    return HttpResponse("welcome")

def ourResponse(request):
    return HttpResponse("OUR RESPONSE")

def takeInput(request):
    return render(request,'input.html')

def similarity(seq1, seq2):
    l1 , l2 = len(seq1), len(seq2)
    ldist = Levenshtein.wf(seq1, seq2)
    return (1 - ldist/max(l1, l2))*100

def performAlgo(request):
    global DIR
    DIR = "/tmp/temp-media/" + generateDirname()
    os.mkdir(DIR)

    myfile = request.FILES['document']
    print(myfile.name)
    fs = FileSystemStorage()
    '''fs.save(myfile.name, myfile)'''
    workbook = xlsxwriter.Workbook(DIR + '/new.xlsx')
    family = request.POST.get("input01")
    outpath = DIR + "/new.xlsx"
    df1 = pd.read_excel(myfile)
    df2 = df1
    for i in range((df1.shape[0] - 1)):
        A = df1.loc[i, "Sequence"]
        B = df1.loc[(i + 1), "Sequence"]
        percent_similarity = similarity(A, B)

        if (percent_similarity >= 90):
            df2 = df2.drop(df2[df2.Sequence == B].index)

    df2.to_excel(outpath, index=False)

    def df_gdomain_counter(df):
        df_count = df["ProteinID"].value_counts()
        return df_count

    def match(x, y, mm):
        mismatch = 0
        for i in range(len(x)):
            if (x[i] == 'X' or x[i] == y[i]):
                pass
            else:
                mismatch += 1
        if (mismatch <= mm):
            return True
        else:
            return False

    def H(protein_id, protein, x1, x2, x3, x4, mm1, mm2, mm3, mm4, min13, min34, min45, max13, max34, max45):
        pL1 = []
        pL2 = []
        pL3 = []
        pL4 = []
        L1 = []
        L2 = []
        L3 = []
        L4 = []
        for i in range(len(protein) - len(x1)):
            if (match(x1, protein[i:i + len(x1)], mm1) == True):
                #                       global L1
                pL1 = pL1 + [i]
                L1 = L1 + [protein[i:i + len(x1)]]
        # print "L1 = ", pL1,L1
        for j in range(len(protein) - len(x2)):
            if (match(x2, protein[j:j + len(x2)], mm2) == True):
                #                       global L2
                pL2 = pL2 + [j]
                L2 = L2 + [protein[j:j + len(x2)]]
        # print "L2 = ", pL2,L2
        for k in range(len(protein) - len(x3)):
            if (match(x3, protein[k:k + len(x3)], mm3) == True):
                #                       global L3
                pL3 = pL3 + [k]
                L3 = L3 + [protein[k:k + len(x3)]]
        # print "L3 = ", pL3,L3
        for l in range(len(protein) - len(x4)):
            if (match(x4, protein[l:l + len(x4)], mm4) == True):
                #                       global L3
                pL4 = pL4 + [l]
                L4 = L4 + [protein[l:l + len(x4)]]
        candidates = []
        for i in range(len(pL1)):
            for j in range(len(pL2)):
                for k in range(len(pL3)):
                    for l in range(len(pL4)):
                        if (min13 <= pL2[j] - pL1[i] <= max13 and min34 <= pL3[k] - pL2[j] <= max34 and min45 <=
                                pL4[l] - pL3[k] <= max45):
                            # if 80 <=pL2[j]-pL1[i]  <= 120 and 40 <=pL3[k]- pL2[j] <= 80 and 20 <=pL4[l]- pL3[k] <= 80
                            a = L1[i]
                            a_pos = pL1[i]
                            b = L2[j]
                            b_pos = pL2[j]
                            c = L3[k]
                            c_pos = pL3[k]
                            d = L4[l]
                            d_pos = pL4[l]

                            candidates.append((protein_id, a, a_pos, b, b_pos, c, c_pos, d, d_pos))

        return candidates

    def shuffler(word):
        word_to_scramble = list(word)
        numpy.random.shuffle(word_to_scramble)
        # O=seq= ''.join(seq_temp)
        new_word = ''.join(word_to_scramble)
        return new_word

    abc = []
    l1 = []

    inpath = DIR + "/new.xlsx"
    mismatch1 = int(request.POST.get("mismatch1"))
    mismatch2 = int(request.POST.get("mismatch2"))
    mismatch3 = int(request.POST.get("mismatch3"))
    mismatch4 = int(request.POST.get("mismatch4"))
    mismatch41 = mismatch4
    x1 = request.POST.get("x1")
    x2 = request.POST.get("x2")
    x3 = request.POST.get("x3")
    x4 = request.POST.get("x4")

    Min_G1_G3 = int(request.POST.get("Min_G1_G3"))
    Max_G1_G3 = int(request.POST.get("Max_G1_G3"))
    Min_G3_G4 = int(request.POST.get("Min_G3_G4"))
    Max_G3_G4 = int(request.POST.get("Max_G3_G4"))
    Min_G4_G5 = int(request.POST.get("Min_G4_G5"))
    Max_G4_G5 = int(request.POST.get("Max_G4_G5"))

    outpath = DIR + "/output_wo_bias.xlsx"
    workbook = xlsxwriter.Workbook(outpath)

    df1 = pd.read_excel(inpath)
    df2 = df1.set_index("Entry", drop=False)
    protein = df2.loc[:, "Sequence"]
    protein_id = df2.loc[:, "Entry"]
    protein_id
    for i in range(len(protein)):
        l = H(protein_id[i], protein[i], x1, x2, x3, x4, mismatch1, mismatch2, mismatch3, mismatch4, Min_G1_G3,
              Min_G3_G4, Min_G4_G5, Max_G1_G3, Max_G3_G4, Max_G4_G5)
        l1.append(l)
        abc = [item for sublist in l1 for item in sublist]
    gdomains = pd.DataFrame(abc,
                            columns=['ProteinID', 'G1-box', 'Position', 'G3-box', 'Position.1', 'G4-box', 'Position.2',
                                     'G5-box', 'Position.3'])
    gdomains = gdomains[gdomains['ProteinID'].astype(bool)]
    gdomains.head()
    gdomains.to_excel(outpath, index=False)

    abc = []
    l1 = []

    outpath = DIR + "/SA_nomismatch.xlsx"
    workbook = xlsxwriter.Workbook(outpath)

    str1 = "XXX"
    x41 = str1 + x4 + "X"
    mismatch41 = 0

    df1 = pd.read_excel(inpath)
    df2 = df1.set_index("Entry", drop=False)
    protein = df2.loc[:, "Sequence"]
    protein_id = df2.loc[:, "Entry"]
    protein_id
    for i in range(len(protein)):
        l = H(protein_id[i], protein[i], x1, x2, x3, x41, mismatch1, mismatch2, mismatch3, mismatch41, Min_G1_G3,
              Min_G3_G4, Min_G4_G5, Max_G1_G3, Max_G3_G4, Max_G4_G5)
        l1.append(l)
        abc = [item for sublist in l1 for item in sublist]
    gdomains = pd.DataFrame(abc,
                            columns=['ProteinID', 'G1-box', 'Position', 'G3-box', 'Position', 'G4-box', 'Position',
                                     'G5-box', 'Position'])
    gdomains = gdomains[gdomains['ProteinID'].astype(bool)]
    gdomains.head()
    gdomains.to_excel(outpath, index=False)

    abc = []
    l1 = []

    outpath = DIR + "/SA_mismatch.xlsx"
    workbook = xlsxwriter.Workbook(outpath)

    df1 = pd.read_excel(inpath)
    df2 = df1.set_index("Entry", drop=False)
    protein = df2.loc[:, "Sequence"]
    protein_id = df2.loc[:, "Entry"]
    protein_id
    for i in range(len(protein)):
        l = H(protein_id[i], protein[i], x1, x2, x3, x41, mismatch1, mismatch2, mismatch3, mismatch4, Min_G1_G3,
              Min_G3_G4, Min_G4_G5, Max_G1_G3, Max_G3_G4, Max_G4_G5)
        l1.append(l)
        abc = [item for sublist in l1 for item in sublist]
    gdomains = pd.DataFrame(abc,
                            columns=['ProteinID', 'G1-box', 'Position', 'G3-box', 'Position', 'G4-box', 'Position',
                                     'G5-box', 'Position'])
    gdomains = gdomains[gdomains['ProteinID'].astype(bool)]
    gdomains.head()
    gdomains.to_excel(outpath, index=False)

    abc = []
    l1 = []

    outpath = DIR + "/A_nomismatch.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/A_nomismatch.xlsx')

    y = x4[1:]
    z = y[:-1]
    x42 = str1 + z + str1
    df1 = pd.read_excel(inpath)
    df2 = df1.set_index("Entry", drop=False)
    protein = df2.loc[:, "Sequence"]
    protein_id = df2.loc[:, "Entry"]
    protein_id
    for i in range(len(protein)):
        l = H(protein_id[i], protein[i], x1, x2, x3, x42, mismatch1, mismatch2, mismatch3, mismatch41, Min_G1_G3,
              Min_G3_G4, Min_G4_G5, Max_G1_G3, Max_G3_G4, Max_G4_G5)
        l1.append(l)
        abc = [item for sublist in l1 for item in sublist]
    gdomains = pd.DataFrame(abc,
                            columns=['ProteinID', 'G1-box', 'Position', 'G3-box', 'Position', 'G4-box', 'Position',
                                     'G5-box', 'Position'])
    gdomains = gdomains[gdomains['ProteinID'].astype(bool)]
    gdomains.head()
    gdomains.to_excel(outpath, index=False)

    inpath_SA_mm = DIR + "/SA_mismatch.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/SA_mm_7mer_X_dict.xlsx')
    outpath1_SA_mm = DIR + "/SA_mm_7mer_X_dict.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/SA_mm_7mer_X_dict_count.xlsx')
    outpath2_SA_mm = DIR + "/SA_mm_7mer_X_dict_count.xlsx"

    def list_of_7mer_X(sevenmer):
        x_data = []
        for r1 in range(7):
            x = list(sevenmer)
            x[r1] = "X"
            x = ''.join(x)
            x_data.append(x)
        return x_data

    str2 = [["Rab", 470], ["Rac", 128], ["Ran", 29], ["Ras", 190], ["Roc", 19], ["Arf", 140], ["AlG1", 44],
            ["Era", 188], ["FeoB", 18], ["Hflx", 26], ["GB1", 116], ["EngB", 401], ["Dynamin", 115], ["IRG", 10],
            ["Obg", 659], ["Septin", 86], ["SRP", 99], ["Translational", 2869], ["tRme", 454], ["EngA", 424]]
    for i in str2:
        if (i[0] == family):
            total = i[1]

    data = pd.read_excel(inpath_SA_mm)
    unique_7mers = data['G5-box'].unique()
    temp = data

    id_set = {}

    for j in range(temp.shape[0]):
        seq = temp.loc[j, "G5-box"]
        ID = temp.loc[j, "ProteinID"]

        x_data = list_of_7mer_X(seq)
        for x in x_data:
            if (x not in id_set):
                id_set[x] = set()
                id_set[x].add(ID)
            else:
                id_set[x].add(ID)
    id_set.items()
    with open(outpath1_SA_mm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set.items()]

    with open(outpath2_SA_mm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
         id_set.items()]

    inpath_A_nomm = DIR + "/A_nomismatch.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/A_nomm_7mer_X_dict.xlsx')
    outpath1_A_nomm = DIR + "/A_nomm_7mer_X_dict.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/A_nomm_7mer_X_dict_count.xlsx')
    outpath2_A_nomm = DIR + "/A_nomm_7mer_X_dict_count.xlsx"

    data1 = pd.read_excel(inpath_A_nomm)
    unique_7mers = data1['G5-box'].unique()
    temp = data1

    id_set1 = {}

    for j in range(temp.shape[0]):
        seq = temp.loc[j, "G5-box"]
        ID = temp.loc[j, "ProteinID"]

        x_data = list_of_7mer_X(seq)
        for x in x_data:
            if (x not in id_set1):
                id_set1[x] = set()
                id_set1[x].add(ID)
            else:
                id_set1[x].add(ID)
    id_set1.items()
    with open(outpath1_A_nomm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set1.items()]

    with open(outpath2_A_nomm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
         id_set1.items()]

    inpath_SA_nomm = DIR + "/SA_nomismatch.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/SA_nomm_7mer_X_dict.xlsx')
    outpath1_SA_nomm = DIR + "/SA_nomm_7mer_X_dict.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/SA_nomm_7mer_X_dict_count.xlsx')
    outpath2_SA_nomm = DIR + "/SA_nomm_7mer_X_dict_count.xlsx"
    data2 = pd.read_excel(inpath_SA_nomm)
    unique_7mers = data2['G5-box'].unique()
    temp = data2

    id_set2 = {}

    for j in range(temp.shape[0]):
        seq = temp.loc[j, "G5-box"]
        ID = temp.loc[j, "ProteinID"]

        x_data = list_of_7mer_X(seq)
        for x in x_data:
            if (x not in id_set2):
                id_set2[x] = set()
                id_set2[x].add(ID)
            else:
                id_set2[x].add(ID)
    id_set2.items()
    with open(outpath1_SA_nomm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set2.items()]

    with open(outpath2_SA_nomm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
         id_set2.items()]

    workbook = xlsxwriter.Workbook(DIR + '/7mer_X_dict.xlsx')
    outpath1 = DIR + "/7mer_X_dict.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/7mer_X_count_dict.xlsx')
    outpath2 = DIR + "/7mer_X_count_dict.xlsx"

    SA_nomm = pd.read_excel(inpath_SA_nomm)
    A_nomm = pd.read_excel(inpath_A_nomm)
    SA_mm = pd.read_excel(inpath_SA_mm)

    table = [SA_nomm[['ProteinID', 'G5-box', 'Position']], A_nomm[['ProteinID', 'G5-box', 'Position']],
             SA_mm[['ProteinID', 'G5-box', 'Position']]]

    # to be used when SA with no mismatch doesn't give any result.
    # table= [A_nomm[['Entry', 'G5_box', 'Value']], SA_mm[['Entry', 'G5_box', 'Value']]]

    data3 = pd.concat(table)
    data3 = data3.reset_index(drop=True)

    unique_7mers = data3['G5-box'].unique()
    temp = data3

    id_set3 = {}

    for j in range(temp.shape[0]):
        seq = temp.loc[j, "G5-box"]
        ID = temp.loc[j, "ProteinID"]

        x_data = list_of_7mer_X(seq)
        for x in x_data:
            if (x not in id_set3):
                id_set3[x] = set()
                id_set3[x].add(ID)
            else:
                id_set3[x].add(ID)

    id_set3.items()
    with open(outpath1, 'w') as f:
        [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set3.items()]

    with open(outpath2, 'w') as f:
        [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
         id_set3.items()]

    def list_of_7mer_2X(sevenmer):
        x_data = []
        for r1 in range(7):
            for r2 in range(7):
                if (r1 != r2):
                    x = list(sevenmer)
                    x[r1] = "X"
                    x[r2] = "X"
                    x = ''.join(x)
                    x_data.append(x)
        return x_data

    workbook = xlsxwriter.Workbook(DIR + '/SA_mm_7mer_2X_dict.xlsx')
    outpath1_SA_mm = DIR + "/SA_mm_7mer_2X_dict.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/SA_mm_7mer_2X_dict_count.xlsx')
    outpath2_SA_mm = DIR + "/SA_mm_7mer_2X_dict_count.xlsx"

    data = pd.read_excel(inpath_SA_mm)
    unique_7mers = data['G5-box'].unique()
    temp = data

    id_set = {}

    for j in range(temp.shape[0]):
        seq = temp.loc[j, "G5-box"]
        ID = temp.loc[j, "ProteinID"]

        x_data = list_of_7mer_2X(seq)
        for x in x_data:
            if (x not in id_set):
                id_set[x] = set()
                id_set[x].add(ID)
            else:
                id_set[x].add(ID)
    id_set.items()
    with open(outpath1_SA_mm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set.items()]

    with open(outpath2_SA_mm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
         id_set.items()]

    workbook = xlsxwriter.Workbook(DIR + '/A_nomm_7mer_2X_dict.xlsx')
    outpath1_A_nomm = DIR + "/A_nomm_7mer_2X_dict.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/A_nomm_7mer_2X_dict_count.xlsx')
    outpath2_A_nomm = DIR + "/A_nomm_7mer_2X_dict_count.xlsx"

    data = pd.read_excel(inpath_A_nomm)
    unique_7mers = data['G5-box'].unique()
    temp = data

    id_set1 = {}

    for j in range(temp.shape[0]):
        seq = temp.loc[j, "G5-box"]
        ID = temp.loc[j, "ProteinID"]

        x_data = list_of_7mer_2X(seq)
        for x in x_data:
            if (x not in id_set1):
                id_set1[x] = set()
                id_set1[x].add(ID)
            else:
                id_set1[x].add(ID)
    id_set1.items()
    with open(outpath1_A_nomm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set1.items()]

    with open(outpath2_A_nomm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
         id_set1.items()]

    workbook = xlsxwriter.Workbook(DIR + '/SA_nomm_7mer_2X_dict.xlsx')
    outpath1_SA_nomm = DIR + "/SA_nomm_7mer_2X_dict.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/SA_nomm_7mer_2X_dict_count.xlsx')
    outpath2_SA_nomm = DIR + "/SA_nomm_7mer_2X_dict_count.xlsx"
    data = pd.read_excel(inpath_SA_nomm)
    unique_7mers = data['G5-box'].unique()
    temp = data

    id_set2 = {}

    for j in range(temp.shape[0]):
        seq = temp.loc[j, "G5-box"]
        ID = temp.loc[j, "ProteinID"]

        x_data = list_of_7mer_2X(seq)
        for x in x_data:
            if (x not in id_set2):
                id_set2[x] = set()
                id_set2[x].add(ID)
            else:
                id_set2[x].add(ID)
    id_set2.items()
    with open(outpath1_SA_nomm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set2.items()]

    with open(outpath2_SA_nomm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
         id_set2.items()]

    workbook = xlsxwriter.Workbook(DIR + '/7mer_2X_dict.xlsx')
    outpath1 = DIR + "/7mer_2X_dict.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/7mer_2X_count_dict.xlsx')
    outpath2 = DIR + "/7mer_2X_count_dict.xlsx"

    SA_nomm = pd.read_excel(inpath_SA_nomm)
    A_nomm = pd.read_excel(inpath_A_nomm)
    SA_mm = pd.read_excel(inpath_SA_mm)

    table = [SA_nomm[['ProteinID', 'G5-box', 'Position']], A_nomm[['ProteinID', 'G5-box', 'Position']],
             SA_mm[['ProteinID', 'G5-box', 'Position']]]

    # table= [A_nomm[['Entry', 'G5_box', 'Value']], SA_mm[['Entry', 'G5_box', 'Value']]]

    data = pd.concat(table)
    data = data.reset_index(drop=True)

    unique_7mers = data['G5-box'].unique()
    temp = data

    id_set3 = {}

    for j in range(temp.shape[0]):
        seq = temp.loc[j, "G5-box"]
        ID = temp.loc[j, "ProteinID"]

        x_data = list_of_7mer_2X(seq)
        for x in x_data:
            if (x not in id_set3):
                id_set3[x] = set()
                id_set3[x].add(ID)
            else:
                id_set3[x].add(ID)

    id_set3.items()
    with open(outpath1, 'w') as f:
        [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set3.items()]

    with open(outpath2, 'w') as f:
        [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
         id_set3.items()]

    def list_of_7mer_3X(sevenmer):
        x_data = []
        for r1 in range(7):
            for r2 in range(7):
                for r3 in range(7):
                    if (r1 != r2 and r1 != r3 and r2 != r3):
                        x = list(sevenmer)
                        x[r1] = "X"
                        x[r2] = "X"
                        x[r3] = "X"
                        x = ''.join(x)
                        x_data.append(x)
        return x_data

    workbook = xlsxwriter.Workbook(DIR + '/SA_mm_7mer_3X_dict.xlsx')
    outpath1_SA_mm = DIR + "/SA_mm_7mer_3X_dict.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/SA_mm_7mer_3X_dict_count.xlsx')
    outpath2_SA_mm = DIR + "/SA_mm_7mer_3X_dict_count.xlsx"

    data = pd.read_excel(inpath_SA_mm)
    unique_7mers = data['G5-box'].unique()
    temp = data

    id_set = {}

    for j in range(temp.shape[0]):
        seq = temp.loc[j, "G5-box"]
        ID = temp.loc[j, "ProteinID"]

        x_data = list_of_7mer_3X(seq)
        for x in x_data:
            if (x not in id_set):
                id_set[x] = set()
                id_set[x].add(ID)
            else:
                id_set[x].add(ID)
    id_set.items()
    with open(outpath1_SA_mm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set.items()]

    with open(outpath2_SA_mm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
         id_set.items()]

    workbook = xlsxwriter.Workbook(DIR + '/A_nomm_7mer_3X_dict.xlsx')
    outpath1_A_nomm = DIR + "/A_nomm_7mer_3X_dict.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/A_nomm_7mer_3X_dict_count.xlsx')
    outpath2_A_nomm = DIR + "/A_nomm_7mer_3X_dict_count.xlsx"

    data = pd.read_excel(inpath_A_nomm)
    unique_7mers = data['G5-box'].unique()
    temp = data

    id_set1 = {}

    for j in range(temp.shape[0]):
        seq = temp.loc[j, "G5-box"]
        ID = temp.loc[j, "ProteinID"]

        x_data = list_of_7mer_3X(seq)
        for x in x_data:
            if (x not in id_set1):
                id_set1[x] = set()
                id_set1[x].add(ID)
            else:
                id_set1[x].add(ID)
    id_set1.items()
    with open(outpath1_A_nomm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set1.items()]

    with open(outpath2_A_nomm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
         id_set1.items()]

    workbook = xlsxwriter.Workbook(DIR + '/SA_nomm_7mer_3X_dict.xlsx')
    outpath1_SA_nomm = DIR + "/SA_nomm_7mer_3X_dict.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/SA_nomm_7mer_3X_dict_count.xlsx')
    outpath2_SA_nomm = DIR + "/SA_nomm_7mer_3X_dict_count.xlsx"
    data = pd.read_excel(inpath_SA_nomm)
    unique_7mers = data['G5-box'].unique()
    temp = data

    id_set2 = {}

    for j in range(temp.shape[0]):
        seq = temp.loc[j, "G5-box"]
        ID = temp.loc[j, "ProteinID"]

        x_data = list_of_7mer_3X(seq)
        for x in x_data:
            if (x not in id_set2):
                id_set2[x] = set()
                id_set2[x].add(ID)
            else:
                id_set2[x].add(ID)
    id_set2.items()
    with open(outpath1_SA_nomm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set2.items()]

    with open(outpath2_SA_nomm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
         id_set2.items()]

    workbook = xlsxwriter.Workbook(DIR + '/7mer_3X_dict.xlsx')
    outpath1 = DIR + "/7mer_3X_dict.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/7mer_3X_count_dict.xlsx')
    outpath2 = DIR + "/7mer_3X_count_dict.xlsx"

    SA_nomm = pd.read_excel(inpath_SA_nomm)
    A_nomm = pd.read_excel(inpath_A_nomm)
    SA_mm = pd.read_excel(inpath_SA_mm)

    table = [SA_nomm[['ProteinID', 'G5-box', 'Position']], A_nomm[['ProteinID', 'G5-box', 'Position']],
             SA_mm[['ProteinID', 'G5-box', 'Position']]]

    # table= [A_nomm[['Entry', 'G5_box', 'Value']], SA_mm[['Entry', 'G5_box', 'Value']]]

    data = pd.concat(table)
    data = data.reset_index(drop=True)

    unique_7mers = data['G5-box'].unique()
    temp = data

    id_set3 = {}

    for j in range(temp.shape[0]):
        seq = temp.loc[j, "G5-box"]
        ID = temp.loc[j, "ProteinID"]

        x_data = list_of_7mer_3X(seq)
        for x in x_data:
            if (x not in id_set3):
                id_set3[x] = set()
                id_set3[x].add(ID)
            else:
                id_set3[x].add(ID)

    id_set3.items()
    with open(outpath1, 'w') as f:
        [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set3.items()]

    with open(outpath2, 'w') as f:
        [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
         id_set3.items()]

    def list_of_7mer_4X(sevenmer):
        x_data = []
        for r1 in range(7):
            for r2 in range(7):
                for r3 in range(7):
                    for r4 in range(7):
                        if (r1 != r2 and r1 != r3 and r1 != r4 and r2 != r3 and r2 != r4 and r3 != r4):
                            x = list(sevenmer)
                            x[r1] = "X"
                            x[r2] = "X"
                            x[r3] = "X"
                            x[r4] = "X"
                            x = ''.join(x)
                            x_data.append(x)
        return x_data

    workbook = xlsxwriter.Workbook(DIR + '/SA_mm_7mer_4X_dict.xlsx')
    outpath1_SA_mm = DIR + "/SA_mm_7mer_4X_dict.xlsx"
    workbook = xlsxwriter.Workbook(DIR + "/SA_mm_7mer_4X_dict_count.xlsx")
    outpath2_SA_mm = DIR + "/SA_mm_7mer_4X_dict_count.xlsx"

    data = pd.read_excel(inpath_SA_mm)
    unique_7mers = data['G5-box'].unique()
    temp = data

    id_set = {}

    for j in range(temp.shape[0]):
        seq = temp.loc[j, "G5-box"]
        ID = temp.loc[j, "ProteinID"]

        x_data = list_of_7mer_4X(seq)
        for x in x_data:
            if (x not in id_set):
                id_set[x] = set()
                id_set[x].add(ID)
            else:
                id_set[x].add(ID)
    id_set.items()
    with open(outpath1_SA_mm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set.items()]

    with open(outpath2_SA_mm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
         id_set.items()]

    workbook = xlsxwriter.Workbook(DIR + '/A_nomm_7mer_4X_dict.xlsx')
    outpath1_A_nomm = DIR + "/A_nomm_7mer_4X_dict.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/A_nomm_7mer_4X_dict_count.xlsx')
    outpath2_A_nomm = DIR + "/A_nomm_7mer_4X_dict_count.xlsx"

    data = pd.read_excel(inpath_A_nomm)
    unique_7mers = data['G5-box'].unique()
    temp = data

    id_set1 = {}

    for j in range(temp.shape[0]):
        seq = temp.loc[j, "G5-box"]
        ID = temp.loc[j, "ProteinID"]

        x_data = list_of_7mer_4X(seq)
        for x in x_data:
            if (x not in id_set1):
                id_set1[x] = set()
                id_set1[x].add(ID)
            else:
                id_set1[x].add(ID)
    id_set1.items()
    with open(outpath1_A_nomm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set1.items()]

    with open(outpath2_A_nomm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
         id_set1.items()]

    workbook = xlsxwriter.Workbook(DIR + '/SA_nomm_7mer_4X_dict.xlsx')
    outpath1_SA_nomm = DIR + "/SA_nomm_7mer_4X_dict.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/SA_nomm_7mer_4X_dict_count.xlsx')
    outpath2_SA_nomm = DIR + "/SA_nomm_7mer_4X_dict_count.xlsx"
    data = pd.read_excel(inpath_SA_nomm)
    unique_7mers = data['G5-box'].unique()
    temp = data

    id_set2 = {}

    for j in range(temp.shape[0]):
        seq = temp.loc[j, "G5-box"]
        ID = temp.loc[j, "ProteinID"]

        x_data = list_of_7mer_4X(seq)
        for x in x_data:
            if (x not in id_set2):
                id_set2[x] = set()
                id_set2[x].add(ID)
            else:
                id_set2[x].add(ID)
    id_set2.items()
    with open(outpath1_SA_nomm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set2.items()]

    with open(outpath2_SA_nomm, 'w') as f:
        [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
         id_set2.items()]

    workbook = xlsxwriter.Workbook(DIR + '/7mer_4X_dict.xlsx')
    outpath1 = DIR + "/7mer_4X_dict.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/7mer_4X_count_dict.xlsx')
    outpath2 = DIR + "/7mer_4X_count_dict.xlsx"

    SA_nomm = pd.read_excel(inpath_SA_nomm)
    A_nomm = pd.read_excel(inpath_A_nomm)
    SA_mm = pd.read_excel(inpath_SA_mm)

    table = [SA_nomm[['ProteinID', 'G5-box', 'Position']], A_nomm[['ProteinID', 'G5-box', 'Position']],
             SA_mm[['ProteinID', 'G5-box', 'Position']]]

    # table= [A_nomm[['Entry', 'G5_box', 'Value']], SA_mm[['Entry', 'G5_box', 'Value']]]

    data = pd.concat(table)
    data = data.reset_index(drop=True)

    unique_7mers = data['G5-box'].unique()
    temp = data

    id_set3 = {}

    for j in range(temp.shape[0]):
        seq = temp.loc[j, "G5-box"]
        ID = temp.loc[j, "ProteinID"]

        x_data = list_of_7mer_4X(seq)
        for x in x_data:
            if (x not in id_set3):
                id_set3[x] = set()
                id_set3[x].add(ID)
            else:
                id_set3[x].add(ID)

    id_set3.items()
    with open(outpath1, 'w') as f:
        [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set3.items()]

    with open(outpath2, 'w') as f:
        [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
         id_set3.items()]

    with open(outpath2, 'rU') as f:
        reader = csv.reader(f)
        next(reader)
        answer = max(int(column[1].replace(',', '').replace('[', '').replace(']', '')) for column in reader)

    for each in id_set3.items():
        if (len(each[1]) == answer):
            str3 = each[0]

    if (str3.startswith("X")):
        str3 = str3[1:]
        if (str3.endswith("X")):
            x43 = str3[:-1]
        elif (str3.startswith("X")):
            x43 = str3[1:]
        else:
            x43 = str3
    elif (str3.endswith("X")):
        str3 = str3[:-1]
        if (str3.startswith("X")):
            x43 = str3[1:]
        elif (str3.endswith("X")):
            x43 = str3[:-1]
        else:
            x43 = str3
    else:
        x43 = str3

    abc = []
    l1 = []

    inpath = DIR + "/new.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/output_new.xlsx')
    outpath = DIR + "/output_new.xlsx"

    df1 = pd.read_excel(inpath)
    df2 = df1.set_index("Entry", drop=False)
    protein = df2.loc[:, "Sequence"]
    protein_id = df2.loc[:, "Entry"]
    protein_id

    for i in range(len(protein)):
        l = H(protein_id[i], protein[i], x1, x2, x3, x43, mismatch1, mismatch2, mismatch3, mismatch4, Min_G1_G3,
              Min_G3_G4, Min_G4_G5, Max_G1_G3, Max_G3_G4, Max_G4_G5)
        l1.append(l)
        abc = [item for sublist in l1 for item in sublist]

    gdomains = pd.DataFrame(abc,
                            columns=['ProteinID', 'G1-box', 'Position', 'G3-box', 'Position.1', 'G4-box', 'Position.2',
                                     'G5-box', 'Position.3'])
    gdomains = gdomains[gdomains['ProteinID'].astype(bool)]

    gdomains.head()
    gdomains.to_excel(outpath, index=False)

    df = pd.read_excel(DIR + "/output_new.xlsx")
    df.to_html(DIR + "/output_new.html")

    df = pd.read_excel(DIR + "/output_wo_bias.xlsx")
    df.to_html(DIR + "/output_wo_bias.html")
    
    counter = {}

    inpath = DIR + "/new.xlsx"
    inpath_before = DIR + "/output_new.xlsx"
    workbook = xlsxwriter.Workbook(DIR + '/neg.xlsx')
    outpath = DIR + "/neg.xlsx"
    df1 = pd.read_excel(inpath)
    df2 = df1.set_index("Entry", drop=False)
    protein = df2.loc[:, "Sequence"]
    protein_id = df2.loc[:, "Entry"]
    protein_id

    for j in range(50):
        abc = []
        l1 = []
        for i in range(len(protein)):
            s = shuffler(protein[i])
            l = H(protein_id[i], s, x1, x2, x3, x43, mismatch1, mismatch2, mismatch3, mismatch4, Min_G1_G3,
                  Min_G3_G4, Min_G4_G5, Max_G1_G3, Max_G3_G4, Max_G4_G5)
            l1.append(l)
            abc = [item for sublist in l1 for item in sublist]
        gdomains = pd.DataFrame(abc, columns=['ProteinID', 'G1-box', 'Position', 'G3-box', 'Position', 'G4-box',
                                              'Position', 'G5-box', 'Position'])
        counter[j] = df_gdomain_counter(gdomains)

    outpath1 = DIR + "/after.csv"
    outpath2 = DIR + "/before.csv"

    xyz = pd.DataFrame(pd.concat(counter.values()))
    xyz['Protein'] = xyz.index
    xyz_mean_after = xyz.groupby(['Protein']).sum() / 50
    xyz_mean_after['Protein'] = xyz_mean_after.index
    xyz_mean_after.to_csv(outpath1)

    xyz_before = pd.read_excel(inpath_before)
    xyz_before = pd.DataFrame(df_gdomain_counter(xyz_before))
    xyz_before['Protein'] = xyz_before.index
    xyz_before.to_csv(outpath2)

    df = pd.read_csv(DIR + "/after.csv")
    df.to_html(DIR + "/after.html")

    df = pd.read_csv(DIR + "/before.csv")
    df.to_html(DIR + "/before.html")

    return render(request, 'view.html', {})

def downloadfile(request):
    with open(DIR + '/output_wo_bias.xlsx', 'rb') as xlsx:
        response = HttpResponse(xlsx.read())
        response['content_type'] = 'application/xlsx'
        response['Content-Disposition'] = 'attachment;filename=result.xlsx'
        return response

def downloadfilenew(request):
    with open(DIR + '/output_new.xlsx', 'rb') as xlsx:
        response = HttpResponse(xlsx.read())
        response['content_type'] = 'application/xlsx'
        response['Content-Disposition'] = 'attachment;filename=result_SMA.xlsx'
        return response

def negative1(request):
    with open(DIR + '/after.csv', 'rb') as csv:
        response = HttpResponse(csv.read())
        response['content_type'] = 'application/csv'
        response['Content-Disposition'] = 'attachment;filename=negative_before.csv'
        return response

def negative2(request):
    with open(DIR + '/before.csv', 'rb') as csv:
        response = HttpResponse(csv.read())
        response['content_type'] = 'application/csv'
        response['Content-Disposition'] = 'attachment;filename=negative_after.csv'
        return response

def negative3(request):
    with open(DIR + '/before.csv', 'r') as f:
        before = f.readlines()
    beforedata = [x.strip('\n').split(',') for x in before]
    beforedata = dict([(x[0], x[1]) for x in beforedata[1:]])

    with open(DIR + '/after.csv', 'r') as f:
        negcontrol = f.readlines()
    negdata = [x.strip('\n').split(',') for x in negcontrol]
    negdata = dict([(x[0], float(x[1])) for x in negdata[1:]])

    xvals, yvals = [], []
    for k in beforedata.items():
        x = k[1]
        if (k[0] in negdata.keys()):
            y = negdata[k[0]]
            xvals.append(x)
            yvals.append(y)
    xy = zip(xvals, yvals)
    xysorted = sorted(xy, key=lambda x: x[0])
    beforevals = [v[0] for v in xysorted]
    negvals = [v[1] for v in xysorted]
    xvals = range(1, len(beforevals) + 1)
    plt.figure(6)
    plot(xvals, beforevals, "b*", label="original")
    plot(xvals, negvals, "ro", label="neg. control")
    legend(loc="upper left")
    xlabel('Proteins')
    ylabel('Number of matches')
    title('Negative Control Plot')
    '''show()'''

    '''return HttpResponse("Success")'''

    buffer = BytesIO()
    canvas = pylab.get_current_fig_manager().canvas
    canvas.draw()
    pilImage = PIL.Image.frombytes("RGB", canvas.get_width_height(), canvas.tostring_rgb())
    pilImage.save(buffer, "PNG")
    pylab.close()

    # Send buffer in a http response the the browser with the mime type image/png set
    return HttpResponse(buffer.getvalue(), content_type="image/png")

def negative_control(request):
    return render(request, 'neg.html', {})

def neg_control(request):
    seq = request.POST.get("query")

    sheet_path = DIR + "/sheet.csv"
    data4 = pd.read_csv(sheet_path)
    data4.head()

    abc = []
    l1 = []

    for xx in range(27):
        protein_family = data4.loc[xx, "Protein"]

        def match(x, y, mm):
            mismatch = 0
            for i in range(len(x)):
                if (x[i] == 'X' or x[i] == y[i]):
                    pass
                else:
                    mismatch += 1
            if (mismatch <= mm):
                return True
            else:
                return False

        def H(protein_family, protein, x1, x2, x3, x4, mm1, mm2, mm3, mm4, min13, min34, min45, max13, max34,
              max45):
            pL1 = []
            pL2 = []
            pL3 = []
            pL4 = []
            L1 = []
            L2 = []
            L3 = []
            L4 = []
            for i in range(len(protein) - len(x1)):
                if (match(x1, protein[i:i + len(x1)], mm1) == True):
                    #                       global L1
                    pL1 = pL1 + [i]
                    L1 = L1 + [protein[i:i + len(x1)]]
            # print "L1 = ", pL1,L1
            for j in range(len(protein) - len(x2)):
                if (match(x2, protein[j:j + len(x2)], mm2) == True):
                    #                       global L2
                    pL2 = pL2 + [j]
                    L2 = L2 + [protein[j:j + len(x2)]]
            # print "L2 = ", pL2,L2
            for k in range(len(protein) - len(x3)):
                if (match(x3, protein[k:k + len(x3)], mm3) == True):
                    #                       global L3
                    pL3 = pL3 + [k]
                    L3 = L3 + [protein[k:k + len(x3)]]
            # print "L3 = ", pL3,L3
            for l in range(len(protein) - len(x4)):
                if (match(x4, protein[l:l + len(x4)], mm4) == True):
                    #                       global L3
                    pL4 = pL4 + [l]
                    L4 = L4 + [protein[l:l + len(x4)]]
            candidates = []
            for i in range(len(pL1)):
                for j in range(len(pL2)):
                    for k in range(len(pL3)):
                        for l in range(len(pL4)):
                            if (min13 <= pL2[j] - pL1[i] <= max13 and min34 <= pL3[k] - pL2[j] <= max34 and min45 <=
                                    pL4[l] - pL3[k] <= max45):
                                # if 80 <=pL2[j]-pL1[i]  <= 120 and 40 <=pL3[k]- pL2[j] <= 80 and 20 <=pL4[l]- pL3[k] <= 80
                                a = L1[i]
                                a_pos = pL1[i]
                                b = L2[j]
                                b_pos = pL2[j]
                                c = L3[k]
                                c_pos = pL3[k]
                                d = L4[l]
                                d_pos = pL4[l]
                                candidates.append((protein_family, a, a_pos, b, b_pos, c, c_pos, d, d_pos, b_pos-a_pos, c_pos-b_pos, d_pos-c_pos))

            return candidates

        inpath = seq
        workbook = xlsxwriter.Workbook(DIR + '/output_wo_bias.xlsx')
        outpath = DIR + "/output_wo_bias.xlsx"

        mismatch1 = data4.loc[xx, "Mismatch_G1"]
        mismatch2 = data4.loc[xx, "Mismatch_G3"]
        mismatch3 = data4.loc[xx, "Mismatch_G4"]
        mismatch4 = data4.loc[xx, "Mismatch_G5"]
        x1 = data4.loc[xx, "G1_box"]
        x2 = data4.loc[xx, "G3_box"]
        x3 = data4.loc[xx, "G4_box"]
        x4 = data4.loc[xx, "G5_box"]
        Min_G1_G3 = data4.loc[xx, "Min_G1_G3"]
        Max_G1_G3 = data4.loc[xx, "Max_G1_G3"]
        Min_G3_G4 = data4.loc[xx, "Min_G3_G4"]
        Max_G3_G4 = data4.loc[xx, "Max_G3_G4"]
        Min_G4_G5 = data4.loc[xx, "Min_G4_G5"]
        Max_G4_G5 = data4.loc[xx, "Max_G4_G5"]

        '''l = H(protein_family, inpath, x1, x2, x3, x4, mismatch1, mismatch2, mismatch3, mismatch4, Min_G1_G3,
              Min_G3_G4, Min_G4_G5, Max_G1_G3, Max_G3_G4, Max_G4_G5)
        l1.append(l)
        abc = [item for sublist in l1 for item in sublist]

        gdomains = pd.DataFrame(abc,
                                columns=['Protein', 'G1-box', 'Position', 'G3-box', 'Position', 'G4-box',
                                         'Position',
                                         'G5-box', 'Position'])
        gdomains = gdomains[gdomains['Protein'].astype(bool)]
        gdomains.head()
        gdomains.to_excel(outpath, index=False)

        abc = []
        l1 = []

        workbook = xlsxwriter.Workbook(DIR + '/SA_nomismatch.xlsx')
        outpath = DIR + "/SA_nomismatch.xlsx"

        str1 = "XXX"
        x41 = str1 + x4 + "X"
        mismatch41 = 0
        l = H(protein_family, inpath, x1, x2, x3, x41, mismatch1, mismatch2, mismatch3, mismatch41, Min_G1_G3,
              Min_G3_G4, Min_G4_G5, Max_G1_G3, Max_G3_G4, Max_G4_G5)
        l1.append(l)
        abc = [item for sublist in l1 for item in sublist]
        gdomains = pd.DataFrame(abc,
                                columns=['Protein', 'G1-box', 'Position', 'G3-box', 'Position', 'G4-box',
                                         'Position',
                                         'G5-box', 'Position'])
        gdomains = gdomains[gdomains['Protein'].astype(bool)]
        gdomains.head()
        gdomains.to_excel(outpath, index=False)

        abc = []
        l1 = []

        workbook = xlsxwriter.Workbook(DIR + '/SA_mismatch.xlsx')
        outpath = DIR + "/SA_mismatch.xlsx"
        l = H(protein_family, inpath, x1, x2, x3, x41, mismatch1, mismatch2, mismatch3, mismatch4, Min_G1_G3,
              Min_G3_G4, Min_G4_G5, Max_G1_G3, Max_G3_G4, Max_G4_G5)
        l1.append(l)
        abc = [item for sublist in l1 for item in sublist]
        gdomains = pd.DataFrame(abc,
                                columns=['Protein', 'G1-box', 'Position', 'G3-box', 'Position', 'G4-box',
                                         'Position',
                                         'G5-box', 'Position'])
        gdomains = gdomains[gdomains['Protein'].astype(bool)]
        gdomains.head()
        gdomains.to_excel(outpath, index=False)

        abc = []
        l1 = []

        workbook = xlsxwriter.Workbook(DIR + '/A_nomismatch.xlsx')
        outpath = DIR + "/A_nomismatch.xlsx"

        y = x4[1:]
        z = y[:-1]
        x42 = str1 + z + str1
        l = H(protein_family, inpath, x1, x2, x3, x42, mismatch1, mismatch2, mismatch3, mismatch41, Min_G1_G3,
              Min_G3_G4, Min_G4_G5, Max_G1_G3, Max_G3_G4, Max_G4_G5)
        l1.append(l)
        abc = [item for sublist in l1 for item in sublist]
        gdomains = pd.DataFrame(abc,
                                columns=['Protein', 'G1-box', 'Position', 'G3-box', 'Position', 'G4-box',
                                         'Position',
                                         'G5-box', 'Position'])
        gdomains = gdomains[gdomains['Protein'].astype(bool)]
        gdomains.head()
        gdomains.to_excel(outpath, index=False)

        inpath_SA_mm = DIR + "/SA_mismatch.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/SA_mm_7mer_X_dict.xlsx')
        outpath1_SA_mm = DIR + "/SA_mm_7mer_X_dict.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/SA_mm_7mer_X_dict_count.xlsx')
        outpath2_SA_mm = DIR + "/SA_mm_7mer_X_dict_count.xlsx"

        def list_of_7mer_X(sevenmer):
            x_data = []
            for r1 in range(7):
                x = list(sevenmer)
                x[r1] = "X"
                x = ''.join(x)
                x_data.append(x)
            return x_data

        total = data4.loc[xx, "total"]

        data = pd.read_excel(inpath_SA_mm)
        unique_7mers = data['G5-box'].unique()
        temp = data

        id_set = {}

        for j in range(temp.shape[0]):
            seq = temp.loc[j, "G5-box"]
            ID = temp.loc[j, "Protein"]

            x_data = list_of_7mer_X(seq)
            for x in x_data:
                if (x not in id_set):
                    id_set[x] = set()
                    id_set[x].add(ID)
                else:
                    id_set[x].add(ID)
        id_set.items()
        with open(outpath1_SA_mm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set.items()]

        with open(outpath2_SA_mm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
             id_set.items()]

        inpath_A_nomm = DIR + "/A_nomismatch.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/A_nomm_7mer_X_dict.xlsx')
        outpath1_A_nomm = DIR + "/A_nomm_7mer_X_dict.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/A_nomm_7mer_X_dict_count.xlsx')
        outpath2_A_nomm = DIR + "/A_nomm_7mer_X_dict_count.xlsx"

        data1 = pd.read_excel(inpath_A_nomm)
        unique_7mers = data1['G5-box'].unique()
        temp = data1

        id_set1 = {}

        for j in range(temp.shape[0]):
            seq = temp.loc[j, "G5-box"]
            ID = temp.loc[j, "Protein"]

            x_data = list_of_7mer_X(seq)
            for x in x_data:
                if (x not in id_set1):
                    id_set1[x] = set()
                    id_set1[x].add(ID)
                else:
                    id_set1[x].add(ID)
        id_set1.items()
        with open(outpath1_A_nomm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set1.items()]

        with open(outpath2_A_nomm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
             id_set1.items()]

        inpath_SA_nomm = DIR + "/SA_nomismatch.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/SA_nomm_7mer_X_dict.xlsx')
        outpath1_SA_nomm = DIR + "/SA_nomm_7mer_X_dict.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/SA_nomm_7mer_X_dict_count.xlsx')
        outpath2_SA_nomm = DIR + "/SA_nomm_7mer_X_dict_count.xlsx"
        data2 = pd.read_excel(inpath_SA_nomm)
        unique_7mers = data2['G5-box'].unique()
        temp = data2

        id_set2 = {}

        for j in range(temp.shape[0]):
            seq = temp.loc[j, "G5-box"]
            ID = temp.loc[j, "Protein"]

            x_data = list_of_7mer_X(seq)
            for x in x_data:
                if (x not in id_set2):
                    id_set2[x] = set()
                    id_set2[x].add(ID)
                else:
                    id_set2[x].add(ID)
        id_set2.items()
        with open(outpath1_SA_nomm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set2.items()]

        with open(outpath2_SA_nomm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
             id_set2.items()]

        workbook = xlsxwriter.Workbook(DIR + '/7mer_X_dict.xlsx')
        outpath1 = DIR + "/7mer_X_dict.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/7mer_X_count_dict.xlsx')
        outpath2 = DIR + "/7mer_X_count_dict.xlsx"

        SA_nomm = pd.read_excel(inpath_SA_nomm)
        A_nomm = pd.read_excel(inpath_A_nomm)
        SA_mm = pd.read_excel(inpath_SA_mm)

        table = [SA_nomm[['Protein', 'G5-box', 'Position']], A_nomm[['Protein', 'G5-box', 'Position']],
                 SA_mm[['Protein', 'G5-box', 'Position']]]

        # to be used when SA with no mismatch doesn't give any result.
        # table= [A_nomm[['Entry', 'G5_box', 'Value']], SA_mm[['Entry', 'G5_box', 'Value']]]

        data3 = pd.concat(table)
        data3 = data3.reset_index(drop=True)

        unique_7mers = data3['G5-box'].unique()
        temp = data3

        id_set3 = {}

        for j in range(temp.shape[0]):
            seq = temp.loc[j, "G5-box"]
            ID = temp.loc[j, "Protein"]

            x_data = list_of_7mer_X(seq)
            for x in x_data:
                if (x not in id_set3):
                    id_set3[x] = set()
                    id_set3[x].add(ID)
                else:
                    id_set3[x].add(ID)

        id_set3.items()
        with open(outpath1, 'w') as f:
            [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set3.items()]

        with open(outpath2, 'w') as f:
            [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
             id_set3.items()]

        def list_of_7mer_2X(sevenmer):
            x_data = []
            for r1 in range(7):
                for r2 in range(7):
                    if (r1 != r2):
                        x = list(sevenmer)
                        x[r1] = "X"
                        x[r2] = "X"
                        x = ''.join(x)
                        x_data.append(x)
            return x_data

        workbook = xlsxwriter.Workbook(DIR + '/SA_mm_7mer_2X_dict.xlsx')
        outpath1_SA_mm = DIR + "/SA_mm_7mer_2X_dict.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/SA_mm_7mer_2X_dict_count.xlsx')
        outpath2_SA_mm = DIR + "/SA_mm_7mer_2X_dict_count.xlsx"

        data = pd.read_excel(inpath_SA_mm)
        unique_7mers = data['G5-box'].unique()
        temp = data

        id_set = {}

        for j in range(temp.shape[0]):
            seq = temp.loc[j, "G5-box"]
            ID = temp.loc[j, "Protein"]

            x_data = list_of_7mer_2X(seq)
            for x in x_data:
                if (x not in id_set):
                    id_set[x] = set()
                    id_set[x].add(ID)
                else:
                    id_set[x].add(ID)
        id_set.items()
        with open(outpath1_SA_mm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set.items()]

        with open(outpath2_SA_mm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
             id_set.items()]

        workbook = xlsxwriter.Workbook(DIR + '/A_nomm_7mer_2X_dict.xlsx')
        outpath1_A_nomm = DIR + "/A_nomm_7mer_2X_dict.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/A_nomm_7mer_2X_dict_count.xlsx')
        outpath2_A_nomm = DIR + "/A_nomm_7mer_2X_dict_count.xlsx"

        data = pd.read_excel(inpath_A_nomm)
        unique_7mers = data['G5-box'].unique()
        temp = data

        id_set1 = {}

        for j in range(temp.shape[0]):
            seq = temp.loc[j, "G5-box"]
            ID = temp.loc[j, "Protein"]

            x_data = list_of_7mer_2X(seq)
            for x in x_data:
                if (x not in id_set1):
                    id_set1[x] = set()
                    id_set1[x].add(ID)
                else:
                    id_set1[x].add(ID)
        id_set1.items()
        with open(outpath1_A_nomm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set1.items()]

        with open(outpath2_A_nomm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
             id_set1.items()]

        workbook = xlsxwriter.Workbook(DIR + '/SA_nomm_7mer_2X_dict.xlsx')
        outpath1_SA_nomm = DIR + "/SA_nomm_7mer_2X_dict.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/SA_nomm_7mer_2X_dict_count.xlsx')
        outpath2_SA_nomm = DIR + "/SA_nomm_7mer_2X_dict_count.xlsx"
        data = pd.read_excel(inpath_SA_nomm)
        unique_7mers = data['G5-box'].unique()
        temp = data

        id_set2 = {}

        for j in range(temp.shape[0]):
            seq = temp.loc[j, "G5-box"]
            ID = temp.loc[j, "Protein"]

            x_data = list_of_7mer_2X(seq)
            for x in x_data:
                if (x not in id_set2):
                    id_set2[x] = set()
                    id_set2[x].add(ID)
                else:
                    id_set2[x].add(ID)
        id_set2.items()
        with open(outpath1_SA_nomm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set2.items()]

        with open(outpath2_SA_nomm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
             id_set2.items()]

        workbook = xlsxwriter.Workbook(DIR + '/7mer_2X_dict.xlsx')
        outpath1 = DIR + "/7mer_2X_dict.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/7mer_2X_count_dict.xlsx')
        outpath2 = DIR + "/7mer_2X_count_dict.xlsx"

        SA_nomm = pd.read_excel(inpath_SA_nomm)
        A_nomm = pd.read_excel(inpath_A_nomm)
        SA_mm = pd.read_excel(inpath_SA_mm)

        table = [SA_nomm[['Protein', 'G5-box', 'Position']], A_nomm[['Protein', 'G5-box', 'Position']],
                 SA_mm[['Protein', 'G5-box', 'Position']]]

        # table= [A_nomm[['Entry', 'G5_box', 'Value']], SA_mm[['Entry', 'G5_box', 'Value']]]

        data = pd.concat(table)
        data = data.reset_index(drop=True)

        unique_7mers = data['G5-box'].unique()
        temp = data

        id_set3 = {}

        for j in range(temp.shape[0]):
            seq = temp.loc[j, "G5-box"]
            ID = temp.loc[j, "Protein"]

            x_data = list_of_7mer_2X(seq)
            for x in x_data:
                if (x not in id_set3):
                    id_set3[x] = set()
                    id_set3[x].add(ID)
                else:
                    id_set3[x].add(ID)

        id_set3.items()
        with open(outpath1, 'w') as f:
            [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set3.items()]

        with open(outpath2, 'w') as f:
            [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
             id_set3.items()]

        def list_of_7mer_3X(sevenmer):
            x_data = []
            for r1 in range(7):
                for r2 in range(7):
                    for r3 in range(7):
                        if (r1 != r2 and r1 != r3 and r2 != r3):
                            x = list(sevenmer)
                            x[r1] = "X"
                            x[r2] = "X"
                            x[r3] = "X"
                            x = ''.join(x)
                            x_data.append(x)
            return x_data

        workbook = xlsxwriter.Workbook(DIR + '/SA_mm_7mer_3X_dict.xlsx')
        outpath1_SA_mm = DIR + "/SA_mm_7mer_3X_dict.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/SA_mm_7mer_3X_dict_count.xlsx')
        outpath2_SA_mm = DIR + "/SA_mm_7mer_3X_dict_count.xlsx"

        data = pd.read_excel(inpath_SA_mm)
        unique_7mers = data['G5-box'].unique()
        temp = data

        id_set = {}

        for j in range(temp.shape[0]):
            seq = temp.loc[j, "G5-box"]
            ID = temp.loc[j, "Protein"]

            x_data = list_of_7mer_3X(seq)
            for x in x_data:
                if (x not in id_set):
                    id_set[x] = set()
                    id_set[x].add(ID)
                else:
                    id_set[x].add(ID)
        id_set.items()
        with open(outpath1_SA_mm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set.items()]

        with open(outpath2_SA_mm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
             id_set.items()]

        workbook = xlsxwriter.Workbook(DIR + '/A_nomm_7mer_3X_dict.xlsx')
        outpath1_A_nomm = DIR + "/A_nomm_7mer_3X_dict.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/A_nomm_7mer_3X_dict_count.xlsx')
        outpath2_A_nomm = DIR + "/A_nomm_7mer_3X_dict_count.xlsx"

        data = pd.read_excel(inpath_A_nomm)
        unique_7mers = data['G5-box'].unique()
        temp = data

        id_set1 = {}

        for j in range(temp.shape[0]):
            seq = temp.loc[j, "G5-box"]
            ID = temp.loc[j, "Protein"]

            x_data = list_of_7mer_3X(seq)
            for x in x_data:
                if (x not in id_set1):
                    id_set1[x] = set()
                    id_set1[x].add(ID)
                else:
                    id_set1[x].add(ID)
        id_set1.items()
        with open(outpath1_A_nomm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set1.items()]

        with open(outpath2_A_nomm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
             id_set1.items()]

        workbook = xlsxwriter.Workbook(DIR + '/SA_nomm_7mer_3X_dict.xlsx')
        outpath1_SA_nomm = DIR + "/SA_nomm_7mer_3X_dict.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/SA_nomm_7mer_3X_dict_count.xlsx')
        outpath2_SA_nomm = DIR + "/SA_nomm_7mer_3X_dict_count.xlsx"
        data = pd.read_excel(inpath_SA_nomm)
        unique_7mers = data['G5-box'].unique()
        temp = data

        id_set2 = {}

        for j in range(temp.shape[0]):
            seq = temp.loc[j, "G5-box"]
            ID = temp.loc[j, "Protein"]

            x_data = list_of_7mer_3X(seq)
            for x in x_data:
                if (x not in id_set2):
                    id_set2[x] = set()
                    id_set2[x].add(ID)
                else:
                    id_set2[x].add(ID)
        id_set2.items()
        with open(outpath1_SA_nomm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set2.items()]

        with open(outpath2_SA_nomm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
             id_set2.items()]

        workbook = xlsxwriter.Workbook(DIR + '/7mer_3X_dict.xlsx')
        outpath1 = DIR + "/7mer_3X_dict.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/7mer_3X_count_dict.xlsx')
        outpath2 = DIR + "/7mer_3X_count_dict.xlsx"

        SA_nomm = pd.read_excel(inpath_SA_nomm)
        A_nomm = pd.read_excel(inpath_A_nomm)
        SA_mm = pd.read_excel(inpath_SA_mm)

        table = [SA_nomm[['Protein', 'G5-box', 'Position']], A_nomm[['Protein', 'G5-box', 'Position']],
                 SA_mm[['Protein', 'G5-box', 'Position']]]

        # table= [A_nomm[['Entry', 'G5_box', 'Value']], SA_mm[['Entry', 'G5_box', 'Value']]]

        data = pd.concat(table)
        data = data.reset_index(drop=True)

        unique_7mers = data['G5-box'].unique()
        temp = data

        id_set3 = {}

        for j in range(temp.shape[0]):
            seq = temp.loc[j, "G5-box"]
            ID = temp.loc[j, "Protein"]

            x_data = list_of_7mer_3X(seq)
            for x in x_data:
                if (x not in id_set3):
                    id_set3[x] = set()
                    id_set3[x].add(ID)
                else:
                    id_set3[x].add(ID)

        id_set3.items()
        with open(outpath1, 'w') as f:
            [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set3.items()]

        with open(outpath2, 'w') as f:
            [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
             id_set3.items()]

        def list_of_7mer_4X(sevenmer):
            x_data = []
            for r1 in range(7):
                for r2 in range(7):
                    for r3 in range(7):
                        for r4 in range(7):
                            if (r1 != r2 and r1 != r3 and r1 != r4 and r2 != r3 and r2 != r4 and r3 != r4):
                                x = list(sevenmer)
                                x[r1] = "X"
                                x[r2] = "X"
                                x[r3] = "X"
                                x[r4] = "X"
                                x = ''.join(x)
                                x_data.append(x)
            return x_data

        workbook = xlsxwriter.Workbook(DIR + '/SA_mm_7mer_4X_dict.xlsx')
        outpath1_SA_mm = DIR + "/SA_mm_7mer_4X_dict.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/SA_mm_7mer_4X_dict_count.xlsx')
        outpath2_SA_mm = DIR + "/SA_mm_7mer_4X_dict_count.xlsx"

        data = pd.read_excel(inpath_SA_mm)
        unique_7mers = data['G5-box'].unique()
        temp = data

        id_set = {}

        for j in range(temp.shape[0]):
            seq = temp.loc[j, "G5-box"]
            ID = temp.loc[j, "Protein"]

            x_data = list_of_7mer_4X(seq)
            for x in x_data:
                if (x not in id_set):
                    id_set[x] = set()
                    id_set[x].add(ID)
                else:
                    id_set[x].add(ID)
        id_set.items()
        with open(outpath1_SA_mm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set.items()]

        with open(outpath2_SA_mm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
             id_set.items()]

        workbook = xlsxwriter.Workbook(DIR + '/A_nomm_7mer_4X_dict.xlsx')
        outpath1_A_nomm = DIR + "/A_nomm_7mer_4X_dict.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/A_nomm_7mer_4X_dict_count.xlsx')
        outpath2_A_nomm = DIR + "/A_nomm_7mer_4X_dict_count.xlsx"

        data = pd.read_excel(inpath_A_nomm)
        unique_7mers = data['G5-box'].unique()
        temp = data

        id_set1 = {}

        for j in range(temp.shape[0]):
            seq = temp.loc[j, "G5-box"]
            ID = temp.loc[j, "Protein"]

            x_data = list_of_7mer_4X(seq)
            for x in x_data:
                if (x not in id_set1):
                    id_set1[x] = set()
                    id_set1[x].add(ID)
                else:
                    id_set1[x].add(ID)
        id_set1.items()
        with open(outpath1_A_nomm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set1.items()]

        with open(outpath2_A_nomm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
             id_set1.items()]

        workbook = xlsxwriter.Workbook(DIR + '/SA_nomm_7mer_4X_dict.xlsx')
        outpath1_SA_nomm = DIR + "/SA_nomm_7mer_4X_dict.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/SA_nomm_7mer_4X_dict_count.xlsx')
        outpath2_SA_nomm = DIR + "/SA_nomm_7mer_4X_dict_count.xlsx"
        data = pd.read_excel(inpath_SA_nomm)
        unique_7mers = data['G5-box'].unique()
        temp = data

        id_set2 = {}

        for j in range(temp.shape[0]):
            seq = temp.loc[j, "G5-box"]
            ID = temp.loc[j, "Protein"]

            x_data = list_of_7mer_4X(seq)
            for x in x_data:
                if (x not in id_set2):
                    id_set2[x] = set()
                    id_set2[x].add(ID)
                else:
                    id_set2[x].add(ID)
        id_set2.items()
        with open(outpath1_SA_nomm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set2.items()]

        with open(outpath2_SA_nomm, 'w') as f:
            [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
             id_set2.items()]

        workbook = xlsxwriter.Workbook(DIR + '/7mer_4X_dict.xlsx')
        outpath1 = DIR + "/7mer_4X_dict.xlsx"
        workbook = xlsxwriter.Workbook(DIR + '/7mer_4X_count_dict.xlsx')
        outpath2 = DIR + "/7mer_4X_count_dict.xlsx"

        SA_nomm = pd.read_excel(inpath_SA_nomm)
        A_nomm = pd.read_excel(inpath_A_nomm)
        SA_mm = pd.read_excel(inpath_SA_mm)

        table = [SA_nomm[['Protein', 'G5-box', 'Position']], A_nomm[['Protein', 'G5-box', 'Position']],
                 SA_mm[['Protein', 'G5-box', 'Position']]]

        # table= [A_nomm[['Entry', 'G5_box', 'Value']], SA_mm[['Entry', 'G5_box', 'Value']]]

        data = pd.concat(table)
        data = data.reset_index(drop=True)

        unique_7mers = data['G5-box'].unique()
        temp = data

        id_set3 = {}

        for j in range(temp.shape[0]):
            seq = temp.loc[j, "G5-box"]
            ID = temp.loc[j, "Protein"]

            x_data = list_of_7mer_4X(seq)
            for x in x_data:
                if (x not in id_set3):
                    id_set3[x] = set()
                    id_set3[x].add(ID)
                else:
                    id_set3[x].add(ID)

        id_set3.items()
        with open(outpath1, 'w') as f:
            [f.write('{0},{1}\n'.format(key, value)) for key, value in id_set3.items()]

        with open(outpath2, 'w') as f:
            [f.write('{0},{1}\n'.format(key, [len(value), round((100 * len(value) / total), 2)])) for key, value in
             id_set3.items()]
        if (len(id_set3) != 0) :

            with open(outpath2, 'rU') as f:
                reader = csv.reader(f)
                answer = max(int(column[1].replace(',', '').replace('[', '').replace(']', '')) for column in reader)

            for each in id_set3.items():
                if (len(each[1]) == answer):
                    str3 = each[0]

            if (str3.startswith("X")):
                str3 = str3[1:]
                if (str3.endswith("X")):
                    x43 = str3[:-1]
                elif (str3.startswith("X")):
                    x43 = str3[1:]
                else:
                    x43 = str3
            elif (str3.endswith("X")):
                str3 = str3[:-1]
                if (str3.startswith("X")):
                    x43 = str3[1:]
                elif (str3.endswith("X")):
                    x43 = str3[:-1]
                else:
                    x43 = str3
            else:
                x43 = str3'''

        x4 = data4.loc[xx, "G5_box_new"]

        workbook = xlsxwriter.Workbook(DIR + '/output_new.xlsx')
        outpath = DIR + "/output_new.xlsx"

        l = H(protein_family, inpath, x1, x2, x3, x4, mismatch1, mismatch2, mismatch3, mismatch4, Min_G1_G3,
                  Min_G3_G4, Min_G4_G5, Max_G1_G3, Max_G3_G4, Max_G4_G5)
        l1.append(l)
        abc = [item for sublist in l1 for item in sublist]

        gdomains = pd.DataFrame(abc,
                                    columns=['Protein', 'G1-box', 'Position', 'G3-box', 'Position', 'G4-box',
                                             'Position',
                                             'G5-box', 'Position', 'G1G3 spacing', 'G3G4 spacing', 'G4G5 spacing'])
        gdomains = gdomains[gdomains['Protein'].astype(bool)]

        gdomains.head()
        gdomains.to_excel(outpath, index=False)

    df = pd.read_excel(DIR + "/output_new.xlsx")
    df.to_html(DIR + "/output_new.html")

    return render(request, 'view1.html', {})

def downloadfile1(request):
    with open(DIR + '/output_new.xlsx', 'rb') as xlsx:
        response = HttpResponse(xlsx.read())
        response['content_type'] = 'application/xlsx'
        response['Content-Disposition'] = 'attachment;filename=result.xlsx'
        return response

def showimage(request):
    inpath_before = DIR + "/output_wo_bias.xlsx"
    inpath_after = DIR + "/output_new.xlsx"
    outpath_G1G3 = DIR + "/G1G3.png"
    outpath_G3G4 = DIR + "/G3G4.png"
    outpath_G4G5 = DIR + "/G4G5.png"
    temp_before = pd.read_excel(inpath_before)
    temp_after = pd.read_excel(inpath_after)
    temp_before_spacing = pd.DataFrame()
    temp_before_spacing['ProteinID'] = temp_before['ProteinID']
    temp_before_spacing['G1_G3'] = temp_before['Position.1'] - (temp_before['Position'] + 7)
    temp_before_spacing['G3_G4'] = temp_before['Position.2'] - (temp_before['Position.1'] + 4)
    temp_before_spacing['G4_G5'] = temp_before['Position.3'] - (temp_before['Position.2'] + 4)
    temp_after_spacing = pd.DataFrame()
    temp_after_spacing['ProteinID'] = temp_after['ProteinID']
    temp_after_spacing['G1_G3'] = temp_after['Position.1'] - (temp_after['Position'] + 7)
    temp_after_spacing['G3_G4'] = temp_after['Position.2'] - (temp_after['Position.1'] + 4)
    temp_after_spacing['G4_G5'] = temp_after['Position.3'] - (temp_after['Position.2'] + 4)
    x1 = temp_before_spacing['G1_G3']
    y1 = temp_after_spacing['G1_G3']

    x1_w = np.ones_like(x1) / float(len(x1))
    y1_w = np.ones_like(y1) / float(len(y1))
    plt.figure(1)
    plt.hist([x1, y1], weights=[x1_w, y1_w], label=['Before', 'After'])
    # plt.hist([x1, y1], density = True, label=['Before', 'After'])

    plt.legend(loc='upper right')

    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.style.use('ggplot')

    plt.xlabel('Spacing')
    plt.ylabel('Frequency (%)')

    plt.title("G1 − G3 Spacing")
    plt.savefig(outpath_G1G3)

    x2 = temp_before_spacing['G3_G4']
    y2 = temp_after_spacing['G3_G4']

    x2_w = np.ones_like(x2) / float(len(x2))
    y2_w = np.ones_like(y2) / float(len(y2))
    plt.figure(2)
    plt.hist([x2, y2], weights=[x2_w, y2_w], label=['Before', 'After'])
    # plt.hist([x2, y2], density = True, label=['Before', 'After'])

    plt.legend(loc='upper right')

    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.style.use('ggplot')

    plt.xlabel('Spacing')
    plt.ylabel('Frequency (%)')

    plt.title("G3 − G4 Spacing")
    plt.savefig(outpath_G3G4)

    x3 = temp_before_spacing['G4_G5']
    y3 = temp_after_spacing['G4_G5']

    x3_w = np.ones_like(x3) / float(len(x3))
    y3_w = np.ones_like(y3) / float(len(y3))
    plt.figure(3)
    plt.hist([x3, y3], weights=[x3_w, y3_w], label=['Before', 'After'])
    # plt.hist([x3, y3], density = True, label=['Before', 'After'])

    plt.legend(loc='upper right')

    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.style.use('ggplot')

    plt.xlabel('Spacing')
    plt.ylabel('Frequency (%)')

    plt.title("G4 − G5 Spacing")
    plt.savefig(outpath_G4G5)

    im3 = PIL.Image.open(outpath_G1G3)
    im4 = PIL.Image.open(outpath_G3G4)
    im5 = PIL.Image.open(outpath_G4G5)

    inpath_before = DIR + "/output_wo_bias.xlsx"
    outpath_Gbefore = DIR + "/Gbefore.png"
    outpath_Gafter = DIR + "/Gafter.png"

    in_file_before = pd.read_excel(inpath_before)

    G1_before = in_file_before[['ProteinID', 'G1-box', 'Position']]
    G3_before = in_file_before[['ProteinID', 'G3-box', 'Position.1']]
    G4_before = in_file_before[['ProteinID', 'G4-box', 'Position.2']]
    G5_before = in_file_before[['ProteinID', 'G5-box', 'Position.3']]

    G1_before.columns = ['ID', 'Sequence', 'value']
    G3_before.columns = ['ID', 'Sequence', 'value']
    G4_before.columns = ['ID', 'Sequence', 'value']
    G5_before.columns = ['ID', 'Sequence', 'value']

    G1_before['box'] = 'G1'
    G3_before['box'] = 'G3'
    G4_before['box'] = 'G4'
    G5_before['box'] = 'G5'

    temp_before = pd.concat([G1_before, G3_before, G4_before, G5_before])

    temp_before = temp_before.drop_duplicates()

    temp.agg_before = temp_before.groupby(["Sequence", "box"]).size().reset_index(name="Value").sort_values(by=['box']).reset_index()
    temp.agg_before['Serial'] = (temp.agg_before.index + 1)
    temp.agg_before
    plt.figure(4)
    groups_before = temp.agg_before.groupby("box")
    for name, group in groups_before:
        plt.plot(group["Serial"], group["Value"], marker="o", linestyle="", label=name)
    legend()
    style.use('ggplot')

    xlabel('Unique G-Boxes')
    ylabel('Frequency')

    title("No. of Unique G-Boxes Vs Frequency")
    savefig(outpath_Gbefore)

    inpath_after = DIR + "/output_new.xlsx"

    in_file_after = pd.read_excel(inpath_after)

    G1_after = in_file_after[['ProteinID', 'G1-box', 'Position']]
    G3_after = in_file_after[['ProteinID', 'G3-box', 'Position.1']]
    G4_after = in_file_after[['ProteinID', 'G4-box', 'Position.2']]
    G5_after = in_file_after[['ProteinID', 'G5-box', 'Position.3']]

    G1_after.columns = ['ID', 'Sequence', 'value']
    G3_after.columns = ['ID', 'Sequence', 'value']
    G4_after.columns = ['ID', 'Sequence', 'value']
    G5_after.columns = ['ID', 'Sequence', 'value']

    G1_after['box'] = 'G1'
    G3_after['box'] = 'G3'
    G4_after['box'] = 'G4'
    G5_after['box'] = 'G5'

    temp_after = pd.concat([G1_after, G3_after, G4_after, G5_after])

    temp_after = temp_after.drop_duplicates()

    temp.agg_after = temp_after.groupby(["Sequence", "box"]).size().reset_index(name="Value").sort_values(
        by=['box']).reset_index()
    temp.agg_after['Serial'] = (temp.agg_after.index + 1)
    temp.agg_after
    plt.figure(5)
    groups_after = temp.agg_after.groupby("box")
    for name, group in groups_after:
        plt.plot(group["Serial"], group["Value"], marker="o", linestyle="", label=name)
    legend()
    style.use('ggplot')

    xlabel('Unique G-Boxes')
    ylabel('Frequency')

    title("No. of Unique G-Boxes Vs Frequency")
    plt.savefig(outpath_Gafter)

    im1 = PIL.Image.open(outpath_Gbefore)
    im2 = PIL.Image.open(outpath_Gafter)

    def get_concat_h_multi_resize(im_list, resample=PIL.Image.BICUBIC):
        min_height = min(im.height for im in im_list)
        im_list_resize = [im.resize((int(im.width * min_height / im.height), min_height), resample=resample)
                          for im in im_list]
        total_width = sum(im.width for im in im_list_resize)
        dst = PIL.Image.new('RGB', (total_width, min_height))
        pos_x = 0
        for im in im_list_resize:
            dst.paste(im, (pos_x, 0))
            pos_x += im.width
        return dst

    def get_concat_v_multi_resize(im_list, resample=PIL.Image.BICUBIC):
        min_width = min(im.width for im in im_list)
        im_list_resize = [im.resize((min_width, int(im.height * min_width / im.width)), resample=resample)
                          for im in im_list]
        total_height = sum(im.height for im in im_list_resize)
        dst = PIL.Image.new('RGB', (min_width, total_height))
        pos_y = 0
        for im in im_list_resize:
            dst.paste(im, (0, pos_y))
            pos_y += im.height
        return dst

    def get_concat_tile_resize(im_list_2d, resample=PIL.Image.BICUBIC):
        im_list_v = [get_concat_h_multi_resize(im_list_h, resample=resample) for im_list_h in im_list_2d]
        return get_concat_v_multi_resize(im_list_v, resample=resample)

    get_concat_tile_resize([[im1, im2], [im3, im4, im5]]).save(DIR + '/concat_1.jpg')



    img = open(DIR + '/concat_1.jpg', 'rb')

    response = FileResponse(img)

    return response



def showimage1(request):
    with open(DIR + '/concat_1.jpg', 'rb') as jpg:
        response = HttpResponse(jpg.read())
        response['content_type'] = 'application/jpg'
        response['Content-Disposition'] = 'attachment;filename=plots.jpg'
        return response
