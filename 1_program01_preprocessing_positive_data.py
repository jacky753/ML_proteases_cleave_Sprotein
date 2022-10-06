# Spike proteinを開裂するproteasesの基質たちを学習して，新たなSpike protein を開裂するproteases を探す機械学習モデルを作成する．

import re
import numpy as np
import glob
import random
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from urllib.request import urlopen
from lxml import etree

# merops からcopy&pasteで作ったproteaseの基質たちの情報が入ったファイルの名前を配列に格納します．
filenames_posi = []
filenames_posi = glob.glob("trainingdata_positive/*")
print(filenames_posi)

#ここからはじまる．
import mysql.connector as mydb

# コネクションの作成
conn = mydb.connect(
    host='localhost',
    port='3306',　#これは文字列型ですが，場合によってはint型にする必要あり．つまり，クオーテーションを取るべきときもあるということ．
    user='XXXXXXX', #MySQLのユーザー名
    password='XXXXXXXX', #ユーザーパスワード
    database='meropsrefs01'
)
# DB操作用にカーソルを作成
cur = conn.cursor()

# MySQLでプロテアーゼの検索に使うMEROPS code辞書作成
positive_proteases = {"HAT":"S01.047", "DESC1":"S01.021", "HATL1":"S01.292", 
                      "TMPRSS2":"S01.247", "TMPRSS4":"S01.034", "TMPRSS13":"S01.087"}

#試しにHATをデータベースmeropsrerfs01で検索してみる
print(positive_proteases["HAT"])
merops_code = [positive_proteases["HAT"]]
cur.execute("SELECT uniprot_acc, p1 FROM cleavage WHERE code=(%s);", merops_code)
a = cur.fetchall()
print(a)

positive_proteases_txt = {"HAT":"S01-047(human_airway_trypsin-like_peptidase).txt", "DESC1":"S01-021(DESC1_peptidase).txt",  
                      "TMPRSS2":"S01-247(epitheliasin).txt", "TMPRSS4":"S01-034(transmembrane_peptidase-serine_4).txt"}

# Uniprot APIを用いてアミノ酸配列を取得するスクリプトを用意します．
#https://seiyakukenkyusya.com/programming/collecting-uniprot-information/
#上記のリンク先を参考にしました．

uid_eg = 'P78365'
def aaseq_from_uid(uid, substrate_turn, index):
    df = pd.DataFrame(np.arange(3).reshape(1, 3), columns=['uniprotKB_accession', 'function', 'sequence'], index=['protease'+str(substrate_turn)+'_substrate'+str(index)])

    for column_name in df:
        df[column_name] = df[column_name].astype(str)

    df['uniprotKB_accession'][0] = uid

    #display(df)    

    url = "https://www.uniprot.org/uniprot/" + uid + ".xml"
    f = urlopen(url)
    xml = f.read()
    root = etree.fromstring(xml)
    
    #以下のコードは下の説明を参照
    function = root.find('./entry/comment[@type="function"]', root.nsmap)
    if function==None:
        print("function was not detected.")
        pass
    else:
        df["function"][0] = function[0].text
        #print(function[0].text+"¥n")

    sequence = root.find('./entry/sequence', root.nsmap) 
    if sequence==None: 
        print("sequence was not detected.")
        pass 
    else: 
        df["sequence"][0] = sequence.text 
        #print(sequence.text+"¥n")

    #display(df) 
    #df0=pd.concat([df0, df], axis=0)
    #display(df0)
    display(df) 
    return df["sequence"][0]

# 取得したアミノ酸からP1, P1'を特定する．今回はP4部位を戻り値にした．
def get_cleave_point(fullaa, cleave_pattern):
    #cleave_point=-1
    p4=-1
    cleave_points = []

    print(cleave_pattern)
    cleave_pattern02 = re.sub("-", "[A-Z]{1}", cleave_pattern)
    print(cleave_pattern02)
    
    for w in range(len(fullaa)-8):
        if re.search(cleave_pattern02, fullaa[w:w+8]):
            cleave_point = w
            cleave_points.append(w)
            print(cleave_points)
            print("cleave_point: " + str(cleave_point))
            continue
    
    print(f'cleave_point_list: {cleave_points}')
    if cleave_points:
        p4 = cleave_points[0]
    if p4 == -1:
        print("there is no cleave_pattern.")
    else:
        pass
    return p4, cleave_points

#アミノ酸全長配列を取得できるかテスト
cp = 0
cpl = []
clpn="[A-Z]{1}[A-Z]{1}ARTLLA"
aaa = aaseq_from_uid("Q6ZS99", 0, 0)
for w in range(len(aaa)-8):
    if re.search(clpn, aaa[w:w+8]) != None:
        cp = w
        cpl.append(w)
        print(cpl)
        print("cleave_point: " + str(cp))
        continue
print(cp)
print(cpl)

# ファイルを各々読み込んで，Uniprot IDからアミノ酸配列を取得します．
# 注意!! Uniprot IDは正規表現で[A-Z]{1}[0-9A-Z]{5}のように表す．
#だけど，文字列"MERNUM"とかがヒットするのはどう対処するのか．

#開裂パターンをコピペした基質情報が並んだファイルから取得します．
def get_cleave_pattern(p4_column):
    seq = ""
    three_aa = ""
    for j in range(8):
        three_aa = str(df.iloc[i, p4_column+j])
        print(f'three_aa: {three_aa}')
        
        #print(one_aa) 
        if three_aa == "Ala":
            one_aa = "A"
        elif three_aa == "Cys":
            one_aa = "C"
        elif three_aa == "Asp":
            one_aa = "D"
        elif three_aa == "Glu":
            one_aa = "E"
        elif three_aa == "Phe":
            one_aa = "F"
        elif three_aa == "Gly":
            one_aa = "G"
        elif three_aa == "His":
            one_aa = "H"
        elif three_aa == "Ile":
            one_aa = "I"
        elif three_aa == "Lys":
            one_aa = "K"
        elif three_aa == "Leu":
            one_aa = "L"
        elif three_aa == "Met":
            one_aa = "M"
        elif three_aa == "Asn":
            one_aa = "N"
        elif three_aa == "Pro":
            one_aa = "P"
        elif three_aa == "Gln":
            one_aa = "Q"
        elif three_aa == "Arg":
            one_aa = "R"
        elif three_aa == "Ser":
            one_aa = "S"
        elif three_aa == "Thr":
            one_aa = "T"
        elif three_aa == "Val":
            one_aa = "V"
        elif three_aa == "Trp":
            one_aa = "W"
        elif three_aa == "Tyr":
            one_aa = "Y"
        else:
            print("This is not a major amino acid or \"- \".")
            print(f'This is {three_aa}.')
            one_aa = "-"
            pass 
        print(f'one_aa: {one_aa}')
        seq = seq + one_aa
    print("cleave_pattern: "+seq)
    return seq

#アミノ酸の電荷を数値に変換します．
def aa_charge(e, cleave_pattern):
    for j in range(8):   
        one_char_aa = cleave_pattern[j]
        if one_char_aa == "A":
            x_train[e][0][j][0] = charge["A"]
        elif one_char_aa == "C":
            x_train[e][0][j][0] = charge["C"]
        elif one_char_aa == "D":
            x_train[e][0][j][0] = charge["D"]
        elif one_char_aa == "E":
            x_train[e][0][j][0] = charge["E"]
        elif one_char_aa == "F":
            x_train[e][0][j][0] = charge["F"]
        elif one_char_aa == "G":
            x_train[e][0][j][0] = charge["G"]
        elif one_char_aa == "H":
            x_train[e][0][j][0] = charge["H"]
        elif one_char_aa == "I":
            x_train[e][0][j][0] = charge["I"]
        elif one_char_aa == "K":
            x_train[e][0][j][0] = charge["K"]
        elif one_char_aa == "L":
            x_train[e][0][j][0] = charge["L"]
        elif one_char_aa == "M":
            x_train[e][0][j][0] = charge["M"]
        elif one_char_aa == "N":
            x_train[e][0][j][0] = charge["N"]
        elif one_char_aa == "P":
            x_train[e][0][j][0] = charge["P"]
        elif one_char_aa == "Q":
            x_train[e][0][j][0] = charge["Q"]
        elif one_char_aa == "R":
            x_train[e][0][j][0] = charge["R"]
        elif one_char_aa == "S":
            x_train[e][0][j][0] = charge["S"]
        elif one_char_aa == "T":
            x_train[e][0][j][0] = charge["T"]
        elif one_char_aa == "V":
            x_train[e][0][j][0] = charge["V"]
        elif one_char_aa == "W":
            x_train[e][0][j][0] = charge["W"]
        elif one_char_aa == "Y":
            x_train[e][0][j][0] = charge["Y"]
        else:
            pass

#アミノ酸の疎水性を数値に変換します．
def aa_hypho(e, cleave_pattern):
    for j in range(8):
        one_char_aa = cleave_pattern[j]
        if one_char_aa == "A":
            x_train[e][1][j][0] = hydrophobicity["A"]
        elif one_char_aa == "C":
            x_train[e][1][j][0] = hydrophobicity["C"]
        elif one_char_aa == "D":
            x_train[e][1][j][0] = hydrophobicity["D"]
        elif one_char_aa == "E":
            x_train[e][1][j][0] = hydrophobicity["E"]
        elif one_char_aa == "F":
            x_train[e][1][j][0] = hydrophobicity["F"]
        elif one_char_aa == "G":
            x_train[e][1][j][0] = hydrophobicity["G"]
        elif one_char_aa == "H":
            x_train[e][1][j][0] = hydrophobicity["H"]
        elif one_char_aa == "I":
            x_train[e][1][j][0] = hydrophobicity["I"]
        elif one_char_aa == "K":
            x_train[e][1][j][0] = hydrophobicity["K"]
        elif one_char_aa == "L":
            x_train[e][1][j][0] = hydrophobicity["L"]
        elif one_char_aa == "M":
            x_train[e][1][j][0] = hydrophobicity["M"]
        elif one_char_aa == "N":
            x_train[e][1][j][0] = hydrophobicity["N"]
        elif one_char_aa == "P":
            x_train[e][1][j][0] = hydrophobicity["P"]
        elif one_char_aa == "Q":
            x_train[e][1][j][0] = hydrophobicity["Q"]
        elif one_char_aa == "R":
            x_train[e][1][j][0] = hydrophobicity["R"]
        elif one_char_aa == "S":
            x_train[e][1][j][0] = hydrophobicity["S"]
        elif one_char_aa == "T":
            x_train[e][1][j][0] = hydrophobicity["T"]
        elif one_char_aa == "V":
            x_train[e][1][j][0] = hydrophobicity["V"]
        elif one_char_aa == "W":
            x_train[e][1][j][0] = hydrophobicity["W"]
        elif one_char_aa == "Y":
            x_train[e][1][j][0] = hydrophobicity["Y"]
        else:
            pass    

#https://trade-and-develop.hatenablog.com/entry/2017/02/23/021119
cur = conn.cursor(buffered=True)

#adding secondary structure（アミノ酸2次構造をデータベースmeropsrefs01から取得する．）
def add_secondary_structure(e, uid, p4, cleave_pattern):
    uid = [uid]
    cur.execute("SELECT substrate_2d FROM substrate_2d where uniprot_acc=(%s);", uid)  
    ss = cur.fetchall()
    if not ss:
        print("substrate_2d is empty.")
        return 0    
    else:
        print(ss)
    ssl=list(ss[0][0])
    
    for j in range(8):
        one_char_aa = cleave_pattern[j]
        if p4-1 +j < 0:
            continue
        if p4-1+j > len(ssl)-1:
            return 0
        if ssl[p4-1 + j] == 'a':
            x_train[e][2][j][0]= 1
        elif ssl[p4-1 + j] =='b':
            x_train[e][3][j][0] = 1
        else:
            pass  
    return 1

#アミノ酸の特性値を読み込む．
#https://www.sigmaaldrich.com/JP/ja/technical-documents/technical-article/protein-biology/protein-structural-analysis/amino-acid-reference-chart
df_aap = pd.read_table('./amino_acids_properties/amino_acid_info_merck.tsv', index_col=0)
print(df_aap)

#アミノ酸の特性値をアミノ酸文字列から変換する．
def aa_properties(e, cleave_pattern):
    for j in range(8):
        one_char_aa = cleave_pattern[j]
        if one_char_aa == "A":
            x_train[e][4][j][0] = df_aap.at['Alanine', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Alanine', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Alanine', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Alanine', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Alanine', 'pl']
            x_train[e][9][j][0] = df_aap.at['Alanine', 'C']
            x_train[e][10][j][0] = df_aap.at['Alanine', 'H']
            x_train[e][11][j][0] = df_aap.at['Alanine', 'N']
            x_train[e][12][j][0] = df_aap.at['Alanine', 'O']
            x_train[e][13][j][0] = df_aap.at['Alanine', 'S']
        elif one_char_aa == "C":
            #x_train[e][12][j][0][2] = hydrophobicity["C"]
            x_train[e][4][j][0] = df_aap.at['Cysteine', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Cysteine', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Cysteine', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Cysteine', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Cysteine', 'pl']
            x_train[e][9][j][0] = df_aap.at['Cysteine', 'C']
            x_train[e][10][j][0] = df_aap.at['Cysteine', 'H']
            x_train[e][11][j][0] = df_aap.at['Cysteine', 'N']
            x_train[e][12][j][0] = df_aap.at['Cysteine', 'O']
            x_train[e][13][j][0] = df_aap.at['Cysteine', 'S']
        elif one_char_aa == "D":
            #x_train[e][15][j][0][2] = hydrophobicity["D"]
            x_train[e][4][j][0] = df_aap.at['Aspartic acid', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Aspartic acid', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Aspartic acid', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Aspartic acid', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Aspartic acid', 'pl']
            x_train[e][9][j][0] = df_aap.at['Aspartic acid', 'C']
            x_train[e][10][j][0] = df_aap.at['Aspartic acid', 'H']
            x_train[e][11][j][0] = df_aap.at['Aspartic acid', 'N']
            x_train[e][12][j][0] = df_aap.at['Aspartic acid', 'O']
            x_train[e][13][j][0] = df_aap.at['Aspartic acid', 'S']
        elif one_char_aa == "E":
            #x_train[e][16][j][0][2] = hydrophobicity["E"]
            x_train[e][4][j][0] = df_aap.at['Glutamic acid', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Glutamic acid', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Glutamic acid', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Glutamic acid', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Glutamic acid', 'pl']
            x_train[e][9][j][0] = df_aap.at['Glutamic acid', 'C']
            x_train[e][10][j][0] = df_aap.at['Glutamic acid', 'H']
            x_train[e][11][j][0] = df_aap.at['Glutamic acid', 'N']
            x_train[e][12][j][0] = df_aap.at['Glutamic acid', 'O']
            x_train[e][13][j][0] = df_aap.at['Glutamic acid', 'S']
        elif one_char_aa == "F":
            #x_train[e][7][j][0][2] = hydrophobicity["F"]
            x_train[e][4][j][0] = df_aap.at['Phenylalanine', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Phenylalanine', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Phenylalanine', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Phenylalanine', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Phenylalanine', 'pl']
            x_train[e][9][j][0] = df_aap.at['Phenylalanine', 'C']
            x_train[e][10][j][0] = df_aap.at['Phenylalanine', 'H']
            x_train[e][11][j][0] = df_aap.at['Phenylalanine', 'N']
            x_train[e][12][j][0] = df_aap.at['Phenylalanine', 'O']
            x_train[e][13][j][0] = df_aap.at['Phenylalanine', 'S']
        elif one_char_aa == "G":
            #x_train[e][0][j][0][2] = hydrophobicity["G"]
            x_train[e][4][j][0] = df_aap.at['Glycine', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Glycine', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Glycine', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Glycine', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Glycine', 'pl']
            x_train[e][9][j][0] = df_aap.at['Glycine', 'C']
            x_train[e][10][j][0] = df_aap.at['Glycine', 'H']
            x_train[e][11][j][0] = df_aap.at['Glycine', 'N']
            x_train[e][12][j][0] = df_aap.at['Glycine', 'O']
            x_train[e][13][j][0] = df_aap.at['Glycine', 'S']
        elif one_char_aa == "H":
            #x_train[e][19][j][0][2] = hydrophobicity["H"]
            x_train[e][4][j][0] = df_aap.at['Histidine', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Histidine', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Histidine', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Histidine', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Histidine', 'pl']
            x_train[e][9][j][0] = df_aap.at['Histidine', 'C']
            x_train[e][10][j][0] = df_aap.at['Histidine', 'H']
            x_train[e][11][j][0] = df_aap.at['Histidine', 'N']
            x_train[e][12][j][0] = df_aap.at['Histidine', 'O']
            x_train[e][13][j][0] = df_aap.at['Histidine', 'S']
        elif one_char_aa == "I":
            #x_train[e][5][j][0][2] = hydrophobicity["I"]
            x_train[e][4][j][0] = df_aap.at['Isoleucine', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Isoleucine', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Isoleucine', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Isoleucine', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Isoleucine', 'pl']
            x_train[e][9][j][0] = df_aap.at['Isoleucine', 'C']
            x_train[e][10][j][0] = df_aap.at['Isoleucine', 'H']
            x_train[e][11][j][0] = df_aap.at['Isoleucine', 'N']
            x_train[e][12][j][0] = df_aap.at['Isoleucine', 'O']
            x_train[e][13][j][0] = df_aap.at['Isoleucine', 'S']
        elif one_char_aa == "K":
            #x_train[e][17][j][0][2] = hydrophobicity["K"]
            x_train[e][4][j][0] = df_aap.at['Lysine', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Lysine', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Lysine', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Lysine', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Lysine', 'pl']
            x_train[e][9][j][0] = df_aap.at['Lysine', 'C']
            x_train[e][10][j][0] = df_aap.at['Lysine', 'H']
            x_train[e][11][j][0] = df_aap.at['Lysine', 'N']
            x_train[e][12][j][0] = df_aap.at['Lysine', 'O']
            x_train[e][13][j][0] = df_aap.at['Lysine', 'S']
        elif one_char_aa == "L":
            #x_train[e][4][j][0][2] = hydrophobicity["L"]
            x_train[e][4][j][0] = df_aap.at['Leucine', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Leucine', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Leucine', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Leucine', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Leucine', 'pl']
            x_train[e][9][j][0] = df_aap.at['Leucine', 'C']
            x_train[e][10][j][0] = df_aap.at['Leucine', 'H']
            x_train[e][11][j][0] = df_aap.at['Leucine', 'N']
            x_train[e][12][j][0] = df_aap.at['Leucine', 'O']
            x_train[e][13][j][0] = df_aap.at['Leucine', 'S']
        elif one_char_aa == "M":
            #x_train[e][6][j][0][2] = hydrophobicity["M"]
            x_train[e][4][j][0] = df_aap.at['Methionine', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Methionine', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Methionine', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Methionine', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Methionine', 'pl']
            x_train[e][9][j][0] = df_aap.at['Methionine', 'C']
            x_train[e][10][j][0] = df_aap.at['Methionine', 'H']
            x_train[e][11][j][0] = df_aap.at['Methionine', 'N']
            x_train[e][12][j][0] = df_aap.at['Methionine', 'O']
            x_train[e][13][j][0] = df_aap.at['Methionine', 'S']
        elif one_char_aa == "N":
            #x_train[e][13][j][0][2] = hydrophobicity["N"]
            x_train[e][4][j][0] = df_aap.at['Asparagine', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Asparagine', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Asparagine', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Asparagine', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Asparagine', 'pl']
            x_train[e][9][j][0] = df_aap.at['Asparagine', 'C']
            x_train[e][10][j][0] = df_aap.at['Asparagine', 'H']
            x_train[e][11][j][0] = df_aap.at['Asparagine', 'N']
            x_train[e][12][j][0] = df_aap.at['Asparagine', 'O']
            x_train[e][13][j][0] = df_aap.at['Asparagine', 'S']
        elif one_char_aa == "P":
            #x_train[e][1][j][0][2] = hydrophobicity["P"]
            x_train[e][4][j][0] = df_aap.at['Proline', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Proline', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Proline', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Proline', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Proline', 'pl']
            x_train[e][9][j][0] = df_aap.at['Proline', 'C']
            x_train[e][10][j][0] = df_aap.at['Proline', 'H']
            x_train[e][11][j][0] = df_aap.at['Proline', 'N']
            x_train[e][12][j][0] = df_aap.at['Proline', 'O']
            x_train[e][13][j][0] = df_aap.at['Proline', 'S']
        elif one_char_aa == "Q":
            #x_train[e][14][j][0][2] = hydrophobicity["Q"]
            x_train[e][4][j][0] = df_aap.at['Glutamine', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Glutamine', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Glutamine', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Glutamine', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Glutamine', 'pl']
            x_train[e][9][j][0] = df_aap.at['Glutamine', 'C']
            x_train[e][10][j][0] = df_aap.at['Glutamine', 'H']
            x_train[e][11][j][0] = df_aap.at['Glutamine', 'N']
            x_train[e][12][j][0] = df_aap.at['Glutamine', 'O']
            x_train[e][13][j][0] = df_aap.at['Glutamine', 'S']
        elif one_char_aa == "R":
            #x_train[e][18][j][0][2] = hydrophobicity["R"]
            x_train[e][4][j][0] = df_aap.at['Arginine', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Arginine', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Arginine', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Arginine', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Arginine', 'pl']
            x_train[e][9][j][0] = df_aap.at['Arginine', 'C']
            x_train[e][10][j][0] = df_aap.at['Arginine', 'H']
            x_train[e][11][j][0] = df_aap.at['Arginine', 'N']
            x_train[e][12][j][0] = df_aap.at['Arginine', 'O']
            x_train[e][13][j][0] = df_aap.at['Arginine', 'S']
        elif one_char_aa == "S":
            #x_train[e][10][j][0][2] = hydrophobicity["S"]
            x_train[e][4][j][0] = df_aap.at['Serine', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Serine', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Serine', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Serine', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Serine', 'pl']
            x_train[e][9][j][0] = df_aap.at['Serine', 'C']
            x_train[e][10][j][0] = df_aap.at['Serine', 'H']
            x_train[e][11][j][0] = df_aap.at['Serine', 'N']
            x_train[e][12][j][0] = df_aap.at['Serine', 'O']
            x_train[e][13][j][0] = df_aap.at['Serine', 'S']
        elif one_char_aa == "T":
            #x_train[e][11][j][0][2] = hydrophobicity["T"]
            x_train[e][4][j][0] = df_aap.at['Threonine', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Threonine', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Threonine', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Threonine', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Threonine', 'pl']
            x_train[e][9][j][0] = df_aap.at['Threonine', 'C']
            x_train[e][10][j][0] = df_aap.at['Threonine', 'H']
            x_train[e][11][j][0] = df_aap.at['Threonine', 'N']
            x_train[e][12][j][0] = df_aap.at['Threonine', 'O']
            x_train[e][13][j][0] = df_aap.at['Threonine', 'S']
        elif one_char_aa == "V":
            #x_train[e][3][j][0][2] = hydrophobicity["V"]
            x_train[e][4][j][0] = df_aap.at['Valine', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Valine', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Valine', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Valine', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Valine', 'pl']
            x_train[e][9][j][0] = df_aap.at['Valine', 'C']
            x_train[e][10][j][0] = df_aap.at['Valine', 'H']
            x_train[e][11][j][0] = df_aap.at['Valine', 'N']
            x_train[e][12][j][0] = df_aap.at['Valine', 'O']
            x_train[e][13][j][0] = df_aap.at['Valine', 'S']
        elif one_char_aa == "W":
            #x_train[e][9][j][0][2] = hydrophobicity["W"]
            x_train[e][4][j][0] = df_aap.at['Tryptophan', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Tryptophan', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Tryptophan', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Tryptophan', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Tryptophan', 'pl']
            x_train[e][9][j][0] = df_aap.at['Tryptophan', 'C']
            x_train[e][10][j][0] = df_aap.at['Tryptophan', 'H']
            x_train[e][11][j][0] = df_aap.at['Tryptophan', 'N']
            x_train[e][12][j][0] = df_aap.at['Tryptophan', 'O']
            x_train[e][13][j][0] = df_aap.at['Tryptophan', 'S']
        elif one_char_aa == "Y":
            #x_train[e][8][j][0][2] = hydrophobicity["Y"]
            x_train[e][4][j][0] = df_aap.at['Tyrosine', 'MolecularWeight']
            x_train[e][5][j][0] = df_aap.at['Tyrosine', 'ResidueWeight_without_water']
            x_train[e][6][j][0] = df_aap.at['Tyrosine', 'pKa']
            x_train[e][7][j][0] = df_aap.at['Tyrosine', 'pKb']
            x_train[e][8][j][0] = df_aap.at['Tyrosine', 'pl']
            x_train[e][9][j][0] = df_aap.at['Tyrosine', 'C']
            x_train[e][10][j][0] = df_aap.at['Tyrosine', 'H']
            x_train[e][11][j][0] = df_aap.at['Tyrosine', 'N']
            x_train[e][12][j][0] = df_aap.at['Tyrosine', 'O']
            x_train[e][13][j][0] = df_aap.at['Tyrosine', 'S']
        else:
            pass    

#uniprot id をコピペした基質情報ファイルから取得する．
def get_uid(m):
    for j in range(4):
        if re.match(r'^[A-Z]{1}[0-9A-Z]{5}$', str(df.iloc[i+m, j])):
            result = re.match(r'^[A-Z]{1}[0-9A-Z]{5}$', str(df.iloc[i+m, j]))
            uniprot_id = result.group(0)
            print(f'uniprot_id: {uniprot_id}')
            break
        else:
            uniprot_id = "null"
            print(f'uniprot_id: null')
    return uniprot_id

# 各アミノ酸について電荷と疎水性の値を定義
charge = {"A":0., "C":-0.0735876, "D":-0.9994991, "E":-0.9987427, "F":0., "G":0., "H":0.0593509, "I":0.,
          "K":0.9997489, "L":0., "M":0., "N":0., "P":0., "Q":0., "R":0.9999950, "S":0., "T":0., "V":0., "W":0., "Y":-0.0001995}
hydrophobicity = {"A":0.1630295, "C":-0.2554557, "D":-0.9794608, "E":-0.7458280, "F":0.9152760,
                  "G":-0.2554557, "H":-0.7869063, "I":0.8510911, "K":-0.8510911, "L":1.,
                  "M":0.7227214, "N":-1., "P":-0.0860077, "Q":-0.8305520, "R":-0.7021823,
                  "S":-0.4685494, "T":-0.3838254, "V":0.6816431, "W":0.8305520, "Y":0.5738126}

#main処理(creating positive data)
count = 0
protease_turn = 0
consecutive_uid = 0
uid_list = []
flag = 0
substrate_uniprot_id_list = []
check = 0
total_subs_count = 0
subs_sum = 0
positive_data_dict = {}

for key_name in positive_proteases_txt:
    #filename = "./trainingdata(positive)/" + filename
    print("==================================protease"+str(count)+"=======================================")
    filename = positive_proteases_txt[key_name]
    path = f'../trainingdata(positive)/{filename}'
    df = pd.read_csv(path, sep='\t')
    display(df)
    
    x_train = np.zeros([len(df), 14, 8, 1])
    print("len(x_train):"+str(len(x_train))+"\n")

    df_cleave_pattern=pd.DataFrame(
        data={'uniprot_id':[], 'cleave_pattern':[]}
    )
    display(df_cleave_pattern)
    for i in range(len(df)):
        print("-------------------protease "+str(protease_turn)+"substrate"+str(i)+"---------------------")
        total_subs_count = total_subs_count + 1
        print("total_subs_count: "+str(total_subs_count))
        secondory_structure = "null"
        

        
        #get Uniprot ID
        uniprot_id = get_uid(0)
        
        #get the position of p4 column
        for j in range(5):
            #print(f'j: {j}')
            if re.match(r'^(N)*(S)*(P)*$', str(df.iloc[i, j])):
                p4_column = j + 2
                print(f'p4_column: {p4_column}')
                break
        
        #func1
        cleave_pattern = ""
        cleave_pattern = get_cleave_pattern(p4_column)       

        uid = uniprot_id 
        uid_list = []
        uid = str(uid)
        uid_list.append(uid)
        substrate_uniprot_id_list.append(uid)


        print("uid:" + uid)
        if uid == "null":
            print("This substrate does not have Uniprot ID.")
            pass
        else:
            fullaa = aaseq_from_uid(uid, protease_turn, i)
            print(f'fullaa_len: {len(fullaa)}')
            
            

        p4, cleave_point_list = get_cleave_point(fullaa, cleave_pattern)
        if p4 == -1:
            print(f'Cleave point was not detected.')
            continue
        print(f'p4: {str(p4)}')

        mpc = search_double_pattern(fullaa, cleave_pattern)
        print("cleave_site_count: "+str(mpc))
        if mpc  > 1:
            flag = flag +1
        
        #func2
        #aa_onezero()
        
        #func3
        aa_charge(i, cleave_pattern)

        #func4
        aa_hypho(i, cleave_pattern)


        #func5
        if uid == "null":
            print("secondary_structure information is none, because uid is null.")
            pass
        else:
            secondary_str = add_secondary_structure(i, uid, p4, cleave_pattern)
            if secondary_str == 0:
                secondary_structure = "null"
            else:
                secondary_structure = "not_null"

        #func6
        aa_properties(i, cleave_pattern)

        positive_data_dict[check] = [total_subs_count, protease_turn, uid, cleave_pattern, secondory_structure, x_train[i]]
        check = check + 1
        print("check: "+str(check))
        print("total_subs_count: "+str(total_subs_count))     
        
        df_cleave_pattern.loc[f'{i}'] = [uid, cleave_pattern]
        display(df_cleave_pattern)
        
        positive_data_dict[check] = [str(total_subs_count), str(protease_turn), uid, cleave_pattern, secondory_structure, x_train[i]]
        
    filename = f'./trainingdata/cleave_pattern_one_letter_aa_{key_name}.csv'
    df_cleave_pattern.to_csv(filename)
    
    print("The number of substrates that had multiple same cleavage sites: "+str(flag))
    np.save('./saving_ndarray/x_train/positive/x_train_positive_protease'+str(protease_turn), x_train)
    protease_turn = protease_turn + 1
    
    if count == 0:
        x_train_positive_all=x_train
        
    else:
        x_train_positive_all = np.append(x_train_positive_all, x_train, 0)
    subs_sum = subs_sum + len(df)
    count = count + 1
    
#print("consecutive_uid: "+str(consecutive_uid))
print("subs_sum: "+str(subs_sum))
print("total_subs_count: "+str(total_subs_count))
print(f"check: {check}")
print("END")

print(positive_data_dict)

#main処理で作成したポジティイブデータの辞書配列を確認する．
import csv
field_name = ['total_subs_count', 'protease_turn', 'uniprot_id', 'cleave_pattern', 'secondory_structure', 'x_train[i]']
with open(r'./saving_ndarray/x_train/positive_data_dict_contents.csv', 'w', encoding='utf-8') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames = field_name)
    writer.writeheader()
    writer.writerows(positive_data_dic)

#ここまでの計算は約10分で終わる．

#save positive
np.save('./saving_ndarray/x_train/x_train_all_only_positive', x_train_positive_all)

#load positive
x_train_positive_all = np.load('./saving_ndarray/x_train/x_train_all_only_positive.npy')

print(len(x_train_positive_all))

#スパース行列を数える．
del_list = []
for i in range(len(x_train_positive_all)):
    #print(x_train_all[i])
    a = x_train_positive_all[i]
    #print(a)
    
    zero_row_count = 0
    for j in range(14):
        print(a[j])
        if np.all(a[j]==0):
            zero_row_count = zero_row_count + 1
    if zero_row_count == 14:
        del_list.append(i)
    print("=====================================================")
    

#x_train_all = np.delete(x_train_all, del_list, 0)
print(x_train_positive_all)
print(len(x_train_positive_all))
print(len(del_list))    
#https://engineeeer.com/check-numpy-array-all-zero/

#スパース行列を削除する．
del_list = []
for i in range(len(x_train_positive_all)):
    #print(x_train_all[i])
    a = x_train_positive_all[i]
    #print(a)
    
    zero_row_count = 0
    for j in range(14):
        print(a[j])
        if np.all(a[j]==0):
            zero_row_count = zero_row_count + 1
    if zero_row_count == 14:
        del_list.append(i)
    print("=====================================================")
    

x_train_positive_all = np.delete(x_train_positive_all, del_list, 0)
print(x_train_positive_all)
print(len(x_train_positive_all))
print(len(del_list))
#https://engineeeer.com/check-numpy-array-all-zero/

print(len(x_train_positive_all))

print(x_train_positive_all[100])
print(positive_data_dic[100][3])
five_cross = len(positive_data_dic)//5
print(five_cross)
rand_list_5cross = random.sample(range(0, len(positive_data_dic), 1), k=five_cross)
print(rand_list_5cross)
df1 = pd.DataFrame(
    data={'check': rand_list_5cross}
)
print(df1)

df1.to_csv("./saving_ndarray/x_train/x_train_positive_all_5crossval.csv")
print(len(rand_list_5cross))

# 検証データpositive_5crossvalを作る．
count = 0
for i in rand_list_5cross:
    if count == 0:
        x_train_positive_all_5crossval = [positive_data_dic[i][5]]
    else:
        x_train_positive_all_5crossval = np.append(x_train_positive_all_5crossval, [positive_data_dic[i][5]], 0)
    count = count + 1
    print(positive_data_dic[i])

print(len(x_train_positive_all_5crossval))
print(x_train_positive_all_5crossval[0])
positive_all_list = set(list(range(0, len(positive_data_dic)))) - set(rand_list_5cross)

print(positive_all_list)
print(len(positive_all_list))

dlist1= list(range(0, len(positive_data_dic)))
dlist2 = rand_list_5cross
dlist3 = [i for i in dlist1 if i not in dlist2]
print(dlist3)
print(len(dlist3))

# 訓練データpositive_allを作る．
count = 0
for i in positive_all_list:
    if count == 0:
        x_train_positive_all = [positive_data_dic[i][5]]
    else:
        x_train_positive_all = np.append(x_train_positive_all, [positive_data_dic[i][5]], 0)
    count = count + 1
    print(positive_data_dic[i])
print(x_train_positive_all)
print(len(x_train_positive_all))

# check what sequence to try create validataion data.
print(x_train_positive_all)

print(x_train_all_5crossval)
print(len(x_train_all_5crossval))

print(len(x_train_positive_all))
y_train_positive_all =  np.array([[1.0, 0.0]]*len(x_train_positive_all))
print(len(y_train_all))

y_train_positive_all_5crossval  = np.array([[1.0, 0.0]]*len(x_train_positive_all_5crossval))
print(len(y_train_all_5crossval))

print(substrate_uniprot_id_list)
print(len(substrate_uniprot_id_list))

ulist_mece = set(substrate_uniprot_id_list)
print(ulist_mece)
print(len(ulist_mece))

print(len(x_train_positive_all))
print(len(x_train_positive_all_5crossval))
print(len(y_train_positive_all))
print(len(y_train_positive_all_5crossval))
print(subs_sum)

#途中経過を保存する．
# saving numarray(ndarray)
np.save('./saving_ndarray/x_train/x_train_positive_all', x_train_positive_all)
np.save('./saving_ndarray/x_train/y_train_positive_all', y_train_positive_all)
np.save('./saving_ndarray/x_train/x_train_positive_all_ 5crossval', x_train_positive_all_5crossval)
np.save('./saving_ndarray/x_train/y_train_positive_all_5crossval', y_train_positive_all_5crossval)
print(len(x_train_positive_all))
print(len(x_train_all_5crossval))
print(len(y_train_all))
print(len(y_train_all_5crossval))
print(subs_sum)

#データをロードする．
#loding ndarray
print(type(np.load('../saving_ndarray/x_train/x_train_all_only_positive.npy')))
print(np.load('./saving_ndarray/x_train/x_train_all_only_positive.npy'))
x_train_positive_all = np.load('./saving_ndarray/x_train/x_train_all_only_positive.npy')
y_train_all = np.load('../saving_ndarray/x_train/y_train_all_only_positive.npy')
x_train_all_5crossval = np.load('./saving_ndarray/x_train/x_train_all_ 5crossval_only_positive.npy')
y_train_all_5crossval = np.load('./saving_ndarray/x_train/y_train_all_5crossval_only_positive.npy')
print(len(x_train_positive_all))
print(len(x_train_all_5crossval))
print(len(y_train_all))
print(len(y_train_all_5crossval))
print(subs_sum)
