import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator
import itertools
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# -------------------------------------------

data = pd.read_csv('shannon_rem_dup_isolates.csv')

data.dropna(how= 'all', inplace= True)

data['Replicate'] = data['Replicate'].astype(int)

data = data[data['Putative species'] != 'Unknown']

# group and count isolates
richness = data.groupby(['Treatment', 'Replicate', 'Medium'])['Putative species'].nunique().reset_index(name='species_richness')
#richness = data.groupby(['Treatment', 'Replicate'])['Putative species'].nunique().reset_index(name='species_richness')

ax = sns.boxplot(data= richness, x= 'Treatment', y= 'species_richness', hue= 'Medium', dodge= True, gap= 0.05)
sns.stripplot(data= richness, x= 'Treatment', y= 'species_richness', hue= 'Medium',
             dodge= True, jitter= True, color= 'black', alpha= 0.8, legend= False)
############## 2. by replicate ##################
# ax = sns.boxplot(data= richness, x='Replicate', y='species_richness', hue='Treatment', palette={"Control": "skyblue", "Fluconazole": "mediumorchid"}, legend= True)
# ax.set_ylim(bottom= 0)

# pts = sns.stripplot(
#     data=richness,
#     x='Replicate', y='species_richness',
#     hue='Treatment', dodge=True, jitter=True,
#     palette={"Control": "black", "Fluconazole": "black"},  # temp color
#     alpha=0.9, ax=ax, legend=False)
#################################################

### statistics
# define comparisons: for each medium, compare Control vs Treated
#pairs = [('DMSO', 'Fluconazole')]

from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm

#model = ols('species_richness ~ Treatment * Medium', data=richness).fit()
#anova_table = anova_lm(model)
#print(anova_table)

# add annotations
#annotator = Annotator(ax, pairs, data= richness, x= 'Treatment', y= 'species_richness')
#annotator.configure(test= 't-test_ind', text_format= 'star')  # or test='Mann-Whitney'
#annotator.apply_and_annotate()
# # define comparisons: for each medium, compare Control vs Treated
media = richness['Medium'].unique()
pairs = [(('Control', m), ('Fluconazole', m)) for m in media]
################################# PREVIOUSLY HAD THIS BEFORE #######################################
# pvals = []
# for pair in pairs:
#     group1 = richness[(richness['Treatment'] == pair[0][0]) & (richness['Medium'] == pair[0][1])]['species_richness']
#     group2 = richness[(richness['Treatment'] == pair[1][0]) & (richness['Medium'] == pair[1][1])]['species_richness']
    
#     stat, p = ttest_ind(group1, group2)
#     pvals.append(p)

# # Step 3: Apply multiple testing correction (e.g., Bonferroni)
# reject, pvals_corrected, _, _ = multipletests(pvals, method='bonferroni')
################################# PREVIOUSLY HAD THIS BEFORE #######################################
from scipy.stats import ttest_ind, mannwhitneyu, shapiro, levene

pvals = []
tests_used = []

for pair in pairs:
    group1 = richness[(richness['Treatment'] == pair[0][0]) & (richness['Medium'] == pair[0][1])]['species_richness']
    group2 = richness[(richness['Treatment'] == pair[1][0]) & (richness['Medium'] == pair[1][1])]['species_richness']
    
    # Check assumptions for t-test
    p_normal1 = shapiro(group1)[1] if len(group1) >= 3 else 0  # Shapiro needs >=3
    p_normal2 = shapiro(group2)[1] if len(group2) >= 3 else 0
    p_var = levene(group1, group2)[1]
    
    stat, p = mannwhitneyu(group1, group2, alternative="two-sided")
    tests_used.append("mannwhitneyu")
    
    pvals.append(p)

# Apply Bonferroni (works regardless of test type)
reject, pvals_corrected, _, _ = multipletests(pvals, method='bonferroni')
# Step 4: Plot and manually annotate
#ax = sns.barplot(data=richness, x='Treatment', y='species_richness', hue='Medium')

########annotator = Annotator(ax, pairs, data=richness, x='Treatment', y='species_richness', hue='Medium')
#####annotator.set_pvalues_and_annotate(pvalues= pvals_corrected)
# # add annotations
# annotator = Annotator(ax, pairs, data= richness, x= 'Treatment', y= 'species_richness', hue= 'Medium')
# annotator.configure(test= 't-test_ind', text_format= 'star')  # or test='Mann-Whitney'
# annotator.apply_and_annotate()
#define comparisons: for each medium, compare Control vs Treated
# replicates = richness['Replicate'].unique()
# treatments = richness['Treatment'].unique()

# pairs = []

# for treatment in treatments:
#     reps_in_treatment = richness[richness['Treatment'] == treatment]['Replicate'].unique()
#     treatment_pairs = list(itertools.combinations(reps_in_treatment, 2))  # all pairwise combos

#     # Each pair is ((RepA, Treatment), (RepB, Treatment))
#     for pair in treatment_pairs:
#         pairs.append(((pair[0], treatment), (pair[1], treatment)))

# # add annotations
# annotator = Annotator(ax, pairs, data= richness, x= 'Replicate', y= 'species_richness', hue= 'Treatment')
# annotator.configure(test= 't-test_ind', text_format= 'star')  # or test='Mann-Whitney'
# annotator.apply_and_annotate()

### configurations
plt.legend(title= 'Treatment')
plt.ylabel('Species Richness', fontsize= 14)#plt.xlabel('Treatment', fontsize= 14)
#plt.xlabel('Replicate', fontsize= 14)
plt.xlabel('Treatment', fontsize= 14)
plt.xticks(ticks= [0, 1], labels= ['Control', 'Fluconazole'], fontsize= 12)
#plt.xticks(ticks= [0, 1, 2, 3], labels= [1,2,3,4], fontsize= 12)
plt.yticks(fontsize= 12)
plt.ylim(bottom= 0, top= 6)
plt.tight_layout()

#plt.savefig("00_NEW_boxplot_species_replicate.svg", format= 'svg')

plt.show()