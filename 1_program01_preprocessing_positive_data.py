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


