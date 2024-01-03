

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns




# presets
hfont = {'fontname':'Helvetica'}


# read data
regdf = pd.read_csv('regdf_with_rxpkp.csv')


# process data
ordered_df1 = regdf[regdf['var']=='p']
ordered_df2 = regdf[regdf['var']=='v']
ordered_df1.reset_index(inplace=True, drop=True)
ordered_df2.reset_index(inplace=True, drop=True)
# log transformation
ordered_df1['rx_permember_log'] = np.log(ordered_df1.RXPKP)
ordered_df2['rx_permember_log'] = np.log(ordered_df2.RXPKP)



# plot
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(7,7), sharex=True, sharey=True)

# horizontal lines at 0.50
axs[0].axhline(y=0.5, color='lightgray')
axs[1].axhline(y=0.5, color='lightgray')

sns.regplot(y = ordered_df1['R2'], 
         x = ordered_df1['rx_permember_log'],
          marker = "o", 
         color='teal',
         ax=axs[0],
         ci=None)
# annotate
for r in range(0,ordered_df1.shape[0]):
     axs[0].text(ordered_df1.rx_permember_log[r]*1.01, 
              ordered_df1.R2[r], 
              ordered_df1.PRIMARYCOND[r], 
              horizontalalignment='left', size='medium', color='black',
              rotation=0, **hfont)

sns.regplot(y = ordered_df2['R2'], 
         x = ordered_df2['rx_permember_log'],
         marker = "o", 
         color = 'salmon', label='Visits per capita', 
         ax=axs[1], ci=None)
# annotate
for r in range(0,ordered_df2.shape[0]):
     axs[1].text(ordered_df2.rx_permember_log[r]*1.01, 
              ordered_df2.R2[r], 
              ordered_df2.PRIMARYCOND[r], 
              horizontalalignment='left', size='medium', color='black',
              rotation=0, **hfont)
# add panel annotations
axs[0].text(-0.8, 1.1, 'A)', 
            size=18, #weight='bold'
            )
# add panel annotations
axs[1].text(-0.8, 1.1, 'B)', 
            size=18, #weight='bold'
            )
axs[0].spines[['right', 'top']].set_visible(False)
axs[1].spines[['right', 'top']].set_visible(False)
axs[0].set_ylabel('$R^2$', fontsize=15, **hfont)
axs[1].set_ylabel('$R^2$', fontsize=15, **hfont)
axs[1].set_xlim((0., 5))
axs[1].set_xticks([ 0.69314718, 1.60943791, 2.30258509, 3.21887582, 3.91202301, 4.60517019])
axs[1].set_xticklabels([2, 5, 10, 25,50,100])
axs[0].set_xlabel('')
axs[1].set_xlabel('Prescriptions per 1000 people')

plt.tight_layout()
plt.savefig('alt_fig1F_2panels.pdf', dpi=300)
plt.savefig('alt_fig1F_2panels.tiff', dpi=300)

