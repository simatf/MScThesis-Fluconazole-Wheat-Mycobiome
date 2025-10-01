import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator
import itertools
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba

# -------------------------------------------

data = pd.read_csv('shannon_rem_dup_isolates.csv')

data.dropna(how= 'all', inplace= True)

data['Replicate'] = data['Replicate'].astype(int)

data = data[data['Putative species'] != 'Unknown']

# group and count isolates
richness = data.groupby(['Treatment', 'Replicate'])['Isolate'].nunique().reset_index(name='isolate_richness') # 'Medium'

ax = sns.boxplot(data= richness, x= 'Treatment', y= 'isolate_richness', hue= 'Treatment', palette={"Control": "skyblue", "Fluconazole": "mediumorchid"}, legend= False)
sns.stripplot(data= richness, x= 'Treatment', y= 'isolate_richness',
              jitter= True, color= 'black', alpha= 0.8)

# --- Add statistical annotation ---
pairs = [("Control", "Fluconazole")]  # what to compare
annot = Annotator(ax, pairs, data=richness, x='Treatment', y='isolate_richness')
annot.configure(test='Mann-Whitney', text_format='star', loc='outside', fontsize= 14)
annot.apply_and_annotate()

### configurations
plt.ylim(bottom= 0, top= 42) # top 16 when medium
#plt.legend(title= 'Treatment')
plt.ylabel('Number of Isolates (N)', fontsize= 16)
plt.xlabel('Treatment', fontsize= 16)
plt.xticks(ticks= [0, 1], labels= ['Control', 'Fluconazole'], fontsize= 14)
plt.yticks(ticks= plt.yticks()[0], fontsize= 14)
plt.tight_layout()

plt.show()