import pprint as p
import xlrd
import statistics as s
import matplotlib.pyplot as plt
from astropy.table import Table
import math
import numpy
from operator import itemgetter



LA_M = 11.34327273
LA_Dev = .5598118338
SR_M = 3.184423077
SR_Dev = .1565052649
SP_M = 3.199285714
SP_Dev = .089399927
SV_M = 29.47321429
SV_Dev = 2.551994393
MV_M = 36.02727273
MV_Dev = 3.313084503
B_M = 7.924528302
B_Dev = 4.820745486
Vol_adv = 100
L_adv =  4.103154365 * (2**.5)
BF_adv = 6.6
BF_max = 15
BF_adv2 = 4.2
W_per_H_adv = 2.706369635
W_per_H_dev = .1985065614


def get_data(file, sheet):
    book = xlrd.open_workbook(file)
    sheet = book.sheet_by_name(sheet)
    data = [[sheet.cell_value(r, c) for c in range(sheet.ncols)] for r in range(sheet.nrows)]
    return data

def identify_null(arr):
    players_del = []
    for row in range(len(arr)):
        for col in range(len(arr[row])):
            if arr[row][col] == '-':
                players_del.append(arr[row][0])
    return players_del

def remove_null(arr, null):

    for name in null:
        for row in range(len(arr) - 1):
            if arr[row][0] == name:
                arr.remove(arr[row])
    return arr

def check_last(arr):

    for col in range(len(arr[-1])):
        if arr[-1][col] == '-':
            arr.remove(arr[-1])

def feet_to_inch(str):

    str_arr = str.split("'")
    height_inch = (float(str_arr[0]) * 12) + float(str_arr[1])
    return height_inch

def int_data():
    data_meses, data_ats = get_data(r'/users/2020jswain/documents/comb.xlsx', 'mes'), get_data(r'/users/2020jswain/documents/comb.xlsx', 'ats')

    players_del_messes, players_del_ats = identify_null(data_meses), identify_null(data_ats)
    players_del = list(dict.fromkeys(players_del_messes + players_del_ats))

    data_meses_clean, data_ats_clean = remove_null(data_meses, players_del), remove_null(data_ats, players_del)
    check_last(data_meses_clean)
    check_last(data_ats_clean)

    players, pos, LA, SR, SP, MV, SV, B, BF, HL, HW, HNS, HS, R, W, WS, W_per_H = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
    ats, meses  = [players, pos, LA, SR, SP, SV, MV, B], [players, pos, BF, HL, HW, HNS, HS, R, W, WS, W_per_H ]


    for at in range(len(data_ats_clean[0])):
        for player in range(len(data_ats_clean) - 2):
            ats[at].append(data_ats_clean[player + 2][at])

    for mes in range(len(data_meses_clean[0])):
        for player in range(len(data_meses_clean) - 1):
            meses[mes].append(data_meses_clean[player + 1][mes])

    for i in range(round(len(players)/2)):
        del players[0]
    for i in range(round(len(pos)/2)):
        del pos[0]

    for i in range(len(BF)):
        BF[i] = BF[i] * 100

    for i in range(len(HS)):
        HS[i]= feet_to_inch(HS[i])

    for i in range(round(len(players))):
        W_per_H.append(W[i] / HW[i])

    return ats, meses

def find_angle(W_per_H_adv, W_per_H_dev, W, H):
    W_per_H = W/H
    dev = ( ((W_per_H_adv - W_per_H)/(W_per_H_dev)) )

    an_Ag = ((math.atan(dev) + (math.pi/2)) / (4))
    an_V = ((math.atan(dev) + (math.pi/2)) / (4))
    an_St = (45 - ((math.atan(dev) + (math.pi/2)) / (4)))
    an_S = ((math.atan(dev) + (math.pi/2)) / (4))
    return an_S, an_Ag, an_V, an_St

def body_fat(BF_p, BF_max):
    BF = (BF_max - BF_p)/4
    return BF

def atributes(LA, LA_M, SR, SR_M, SP, SP_M, SV, SV_M, MV, MV_M, B, B_M):
    Ag =  ( ((LA_M - LA)/(LA_Dev)) + ((SR_M - SR)/(SR_Dev)) ) /2
    V = ( ((SV_M - SV)/(SV_Dev)) + ( (MV_M - MV)/(MV_Dev)) ) /2  * -1
    St = ( ((B_M - B)/(B_Dev)) ) * -1
    S = ( ((SP_M - SP)/(SP_Dev)) )
    return Ag, V, St, S

def find_mag(Ag, V, St, S, Vol_adv, BF_adv):
    L_Ag = L_adv * (1.025 ** Ag)
    L_V = L_adv  * (1.025 ** V)
    L_S = L_adv * (1.025 ** S)
    L_St = L_adv * (1.025 ** St)
    return L_Ag, L_V, L_S, L_St

def find_vol(an_Ag, L_Ag, an_V, L_V ,an_S, L_S, an_St, L_St, BF):
    Ag_x = L_Ag * -math.cos(an_Ag)
    Ag_y = L_Ag * math.sin(an_Ag)

    V_x = L_V * math.cos(an_V)
    V_y = L_V * math.sin(an_V)

    S_x = L_S * -math.cos(an_S)
    S_y = L_S * -math.sin(an_S)

    St_x = L_St * math.cos(an_St)
    St_y = L_St * -math.sin(an_St)

    vol = .5 * BF * ( ((S_x - V_x) * (St_y - Ag_y)) + ((S_y - V_y) * (Ag_x - St_x)) )
    return vol

def run():

    ans, vols, players_arr = [], [], []
    ats, meses = int_data()

    for player in range(len(ats[0])):

        players = ats[0][player]
        pos = ats[1][player]
        LA = ats[2][player]
        SR = ats[3][player]
        SP = ats[4][player]
        SV = ats[5][player]
        MV = ats[6][player]
        B = ats[7][player]

        players  = meses[0][player]
        pos = meses[1][player]
        BF = meses[2][player]
        HL = meses[3][player]
        HW = meses[4][player]
        HNS = meses[5][player]
        HS = meses[6][player]
        R = meses[7][player]
        W = meses[8][player]
        WS = meses[9][player]
        W_per_H = meses[10][player]

        BF = body_fat(BF, BF_max)
        AT = atributes(LA, LA_M, SR, SR_M, SP, SP_M, SV, SV_M, MV, MV_M, B, B_M)
        mag = find_mag(AT[0], AT[1], AT[2], AT[3], Vol_adv, BF_adv)
        angle = find_angle(W_per_H_adv, W_per_H_dev, W, HS)
        vol = find_vol(angle[1], mag[0], angle[2], mag[1] ,angle[0], mag[2], angle[3], mag[3], BF)
        ans.append([players, vol/100])
    ans = sorted(ans, key=lambda x: -x[1])


    for row in ans:
        players_arr.append(row[0])
        vols. append(row[1])

    print(players_arr, vols)

    p.pprint(ans)




run()
