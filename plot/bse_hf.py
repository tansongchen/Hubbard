import matplotlib.pyplot as plt
import matplotlib.font_manager
myfont = matplotlib.font_manager.FontProperties(fname='/System/Library/Fonts/PingFang.ttc')
from Helper import get_data_from, get_colors_with
import math
from matplotlib import lines

with open('output/HFOrbitals.dat') as f:
    energyList = [int(float(x)*100000) for x in f.readlines()[0].strip().split()]

d = {}
temp = -100
for energy in energyList:
    if energy == temp:
        d[energy] = d[energy] + 1
    else:
        temp = energy
        d[energy] = 1

print(d)

energyList = [x/100000 for x in energyList]

padding = (energyList[-1] - energyList[0]) / 10 # 首先为能级图确定一个边距，以免最高能级和最低能级和边框相距过近
figure = plt.figure(figsize=(2,6))
ax = plt.gca()
ax.get_xaxis().set_visible(False) # x轴无实际意义，隐藏
ax.get_yaxis().set_label_text('Energy')
ax.set_ylim(energyList[0] - padding, energyList[-1] + padding) # 根据边距设定y轴范围
maxDeg = max(d.values())
width = 1/maxDeg
for key, value in d.items():
    start = (maxDeg - value) / 2
    for i in range(value):
        line = lines.Line2D(((start+i+0.3)*width, (start+i+0.7)*width), (key/100000, key/100000), c='black') # 逐个添加能级
        ax.add_line(line)
plt.show()
