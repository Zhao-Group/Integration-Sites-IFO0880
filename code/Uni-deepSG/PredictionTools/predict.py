from turtle import forward
import numpy as np
import torch 
import torch.nn as nn 
import torch.nn.functional as F
from torch.utils.data import Dataset
import torch.utils.data as Data
import data_loader as dl
import train as tt
import os
import matplotlib.pyplot as plt
import scipy.stats as ss
import pandas as pd

df = pd.read_csv('data/ifo0880/expr_efficiency.csv')

scores = []
for i in range(len(df)):
    sgrna=df['Guide Sequence'][i] + df['PAM'][i]
    state="type1"

    if state=="type1":
        input = sgrna+"gtttTagagctagaaatagcaagttAaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttt"[:36].upper()#add scar if type1
    elif state=="type2":
        input = sgrna+"gtttCagagctaTGCTGgaaaCAGCAtagcaagttGaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttt"[:36].upper()  #add scar if type2
    else:
        print("add scar after the 23bp. choose type2 or type1 for sgrna scarfford")
    input = dl.translateseq(input)
    #print(input)
    input =input.cuda()
    typ = torch.LongTensor([[0]]).cuda()

    net=torch.load("code/Uni-deepSG/PredictionTools/generate_data/model/100.pt")
    #print(net)
    
    net.eval()
    eff = net(input,typ, train = False)
    scores.append(eff.item())

df['Uni-deepSG_eff'] = scores

import scipy.stats
import seaborn as sns
import matplotlib.pyplot as plt

print(scipy.stats.pointbiserialr(list(df['Uni-deepSG_eff']), list(df['Efficiency (Binary)'])))

from scipy.stats import f_oneway

# Separate data into groups
groups = [df['Uni-deepSG_eff'][df['Efficiency'] == 'Good'],
          df['Uni-deepSG_eff'][df['Efficiency'] == 'Bad']]

# Perform one-way ANOVA
f_statistic, p_value = f_oneway(*groups)

print(f"F-statistic: {f_statistic}")
print(f"P-value: {p_value}")

ax = sns.violinplot(x="Efficiency", y="Uni-deepSG_eff", data=df, hue="Efficiency", palette={"Good": "lightblue", "Bad": "mistyrose"}, dodge=False)
ax = sns.swarmplot(x="Efficiency", y="Uni-deepSG_eff", data=df, color=".25")

ax.grid(False)

ax.set_ylabel('Predicted score')
ax.set_xlabel('gRNA editing efficiency')
ax.set_title('Uni-deepSG')

ylim = ax.get_ylim()
# Annotate the plot with p-value
plt.text(0.5, 0.9 * (ylim[1] - ylim[0]) +  ylim[0], f'p-value: {p_value:.3f}', ha='center', va='center', zorder=3)

ax.legend().remove()

plt.tight_layout()
plt.savefig('data/ifo0880/Uni-deepSG_vs_experimental_efficiency.png', dpi = 300)