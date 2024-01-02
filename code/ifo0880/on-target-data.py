#Modules
import numpy as np
import pandas as pd
import scipy.stats
import seaborn as sns
import matplotlib.pyplot as plt
from keras.models import Sequential, Model
from keras.layers.core import  Dropout, Activation, Flatten
from keras.regularizers import l1,l2,l1_l2
from keras.layers import Conv1D, MaxPooling1D, Dense, LSTM, Bidirectional, BatchNormalization, MaxPooling2D, AveragePooling1D, Input, Multiply, Add, UpSampling1D
import sys
import os

sys.path.append(os.path.abspath(os.path.join(sys.path[0], '..')))
print(sys.path[-1])

#Data
all_df = pd.read_csv('data/ifo0880/output_10_5.csv')
df = pd.read_csv('data/ifo0880/expr_efficiency.csv')

# Merge dataframes based on different columns
merged_df = pd.merge(df, all_df, on='Guide Sequence', how='left')[['ID_x', 'Guide Sequence', 'PAM_x', 'Efficiency', 'Efficiency (Binary)', 'Left HR', 'Right HR', 'On-target Score']]

print(scipy.stats.pointbiserialr(list(merged_df['On-target Score']), list(merged_df['Efficiency (Binary)'])))

from scipy.stats import f_oneway

# Separate data into groups
groups = [merged_df['On-target Score'][merged_df['Efficiency'] == 'Good'],
          merged_df['On-target Score'][merged_df['Efficiency'] == 'Bad']]

# Perform one-way ANOVA
f_statistic, p_value = f_oneway(*groups)

print(f"F-statistic: {f_statistic}")
print(f"P-value: {p_value}")

ax = sns.violinplot(x="Efficiency", y="On-target Score", data=merged_df, hue="Efficiency", palette={"Good": "lightblue", "Bad": "mistyrose"}, dodge=False)
ax = sns.swarmplot(x="Efficiency", y="On-target Score", data=merged_df, color=".25")

ax.grid(False)

ax.set_ylabel('Predicted score')
ax.set_xlabel('gRNA editing efficiency')
ax.set_title('Rule Set 2')

ylim = ax.get_ylim()
# Annotate the plot with p-value
plt.text(0.5, 0.9 * (ylim[1] - ylim[0]) +  ylim[0], f'p-value: {p_value:.3f}', ha='center', va='center', zorder=3)

plt.tight_layout()
plt.savefig('data/ifo0880/Rule Set 2_vs_experimental_efficiency.png', dpi = 300)
plt.clf()

def one_hot_encoding(lines):
    data_n = len(lines) 
    SEQ = np.zeros((data_n, len(lines[0]), 4), dtype=int)
    
    for l in range(0, data_n):
        seq = lines[l]
        for i in range(28):
            if seq[i] in "Aa":
                SEQ[l, i, 0] = 1
            elif seq[i] in "Cc":
                SEQ[l, i, 1] = 1
            elif seq[i] in "Gg":
                SEQ[l, i, 2] = 1
            elif seq[i] in "Tt":
                SEQ[l, i, 3] = 1

    return SEQ

def scores_guides_cas9(guides):

    seq = one_hot_encoding(guides)
    
    # model architecture
    SEQ = Input(shape=(28,4))
    conv_1 = Conv1D(activation="relu", padding="same", strides=1, filters=20, kernel_size=5, kernel_regularizer = l2(0.0001))(SEQ)
    bat_norm1 = BatchNormalization()(conv_1)
    pool = MaxPooling1D(pool_size=(2))(bat_norm1)
    conv_2 = Conv1D(activation="relu", padding="same", strides=1, filters=40, kernel_size=8, kernel_regularizer = l2(0.0001))(pool)
    bat_norm2 = BatchNormalization()(conv_2)
    pool_1 = AveragePooling1D(pool_size=(2))(bat_norm2)
    enc = pool_1
    dec_pool_1 =  UpSampling1D(size=2)(enc)
    dec_bat_norm2 = BatchNormalization()(dec_pool_1)
    dec_conv_2  = Conv1D(activation="relu", padding="same", strides=1, filters=40, kernel_size=8, kernel_initializer='glorot_uniform',kernel_regularizer = l2(0.0001))(dec_bat_norm2)
    dec_pool = UpSampling1D(size=2)(dec_conv_2)
    dec_conv_1 = Conv1D(activation="relu", padding="same", strides=1, filters=20, kernel_size=5, kernel_initializer='glorot_uniform',kernel_regularizer = l2(0.0001))(dec_pool)
    dec = Conv1D(activation="relu", padding="same", strides=1, filters=4, kernel_size=5, kernel_initializer='glorot_uniform',kernel_regularizer = l2(0.0001))(dec_pool)
    model_seq = Model(inputs = SEQ, outputs= dec) 
    flatten = Flatten()(enc)
    dropout_1 = Dropout(0.5)(flatten)
    dense_1 = Dense(80, activation='relu', kernel_initializer='glorot_uniform')(dropout_1)
    dropout_2 = Dropout(0.5)(dense_1)
    out = Dense(units=1,  activation="linear")(dropout_2) 
    model = Model(inputs = SEQ, outputs= out)
    model.load_weights(sys.path[-1] + "/model/cas9_seq_wtt.h5")
    pred_y = model.predict(seq)

    score = [-1*ele for ele in pred_y.flatten()]
    
    return score

on_target_seq = []
for i in range(len(merged_df)):
    on_target_seq.append(merged_df['Left HR'][i][-2:] + merged_df['Guide Sequence'][i] + merged_df['PAM_x'][i] + merged_df['Right HR'][i][0:3])
    
merged_df['On-target Score DeepGuide'] = scores_guides_cas9(on_target_seq)

print(scipy.stats.pointbiserialr(list(merged_df['On-target Score DeepGuide']), list(merged_df['Efficiency (Binary)'])))

from scipy.stats import f_oneway

# Separate data into groups
groups = [merged_df['On-target Score DeepGuide'][merged_df['Efficiency'] == 'Good'],
          merged_df['On-target Score DeepGuide'][merged_df['Efficiency'] == 'Bad']]

# Perform one-way ANOVA
f_statistic, p_value = f_oneway(*groups)

print(f"F-statistic: {f_statistic}")
print(f"P-value: {p_value}")

ax = sns.violinplot(x="Efficiency", y="On-target Score DeepGuide", data=merged_df, hue="Efficiency", palette={"Good": "lightblue", "Bad": "mistyrose"}, dodge=False)
ax = sns.swarmplot(x="Efficiency", y="On-target Score DeepGuide", data=merged_df, color=".25")

ax.grid(False)

ax.set_ylabel('Predicted score')
ax.set_xlabel('gRNA editing efficiency')
ax.set_title('DeepGuide')

ylim = ax.get_ylim()
# Annotate the plot with p-value
plt.text(0.5, 0.9 * (ylim[1] - ylim[0]) +  ylim[0], f'p-value: {p_value:.3f}', ha='center', va='center', zorder=3)

plt.tight_layout()
plt.savefig('data/ifo0880/DeepGuide_vs_experimental_efficiency.png', dpi = 300)
plt.clf()