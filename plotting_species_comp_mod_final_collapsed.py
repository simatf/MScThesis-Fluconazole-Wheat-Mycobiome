import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import shapiro, levene, ttest_ind, mannwhitneyu
from matplotlib.patches import Patch

# Load data
df = pd.read_csv('data.csv')

# Normalize species names (remove 'aff.')
df['Species_Grouped'] = df['Species'].str.replace(r'\baff\.\s*', '', regex=True)

# Filter only DMSO and 5Fluconazole375
dmso_df = df[df['Group'] == '0DMSO']
flu375_df = df[df['Group'] == '5Fluconazole375']

# Calculate DMSO mean per isolate
dmso_means = dmso_df.groupby('Isolate')['Value'].mean().rename('DMSO_Mean')

# Merge and compute kill %
flu375_df = flu375_df.merge(dmso_means, on='Isolate', how='left')
flu375_df['Kill_Percent'] = (1 - (flu375_df['Value'] / flu375_df['DMSO_Mean'])) * 100

# Keep species, kill percent, and source
species_data = flu375_df[['Species', 'Species_Grouped', 'Kill_Percent', 'Source']].copy()

# Group isolates per species
grouped = []
for species, sub_df in species_data.groupby('Species_Grouped'):
    combined_data = sub_df['Kill_Percent'].dropna().tolist()
    sources = sub_df['Source'].unique()
    
    if len(sources) == 1 and sources[0] == 'Control':
        group = 'Control'
    else:
        group = 'Fluconazole + Control'
    
    for val in combined_data:
        grouped.append({
            'Species_Grouped': species,
            'Kill_Percent': val,
            'Group': group
        })

plot_df = pd.DataFrame(grouped)

# Plot
plt.figure(figsize=(12, 6))

# Define palette mapping
palette_map = {'Control':'skyblue', 'Fluconazole + Control': 'salmon'}

group1 = plot_df[plot_df['Group'] == 'Control']['Kill_Percent'].dropna().values
group2 = plot_df[plot_df['Group'] == 'Fluconazole + Control']['Kill_Percent'].dropna().values

# Choose test
if len(group1) >= 3 and len(group2) >= 3:
    normal1 = shapiro(group1)[1] > 0.05
    normal2 = shapiro(group2)[1] > 0.05
else:
    normal1 = normal2 = False

equal_var = levene(group1, group2)[1] > 0.05

if normal1 and normal2:
    test_name = "t-test"
    stat, p_val = ttest_ind(group1, group2, equal_var=equal_var)
else:
    test_name = "Mann-Whitney U"
    stat, p_val = mannwhitneyu(group1, group2, alternative='two-sided')

# Plot
ax = sns.boxplot(
    x='Group', 
    y='Kill_Percent', 
    data=plot_df, 
    palette=palette_map,
    showfliers=False
)
sns.stripplot(
    x='Group',
    y='Kill_Percent',
    data=plot_df,
    color='black',
    jitter=True,
    size=4
)

# Add annotation above the two boxes
x1, x2 = 0, 1  # positions of the boxes on x-axis
y_max = plot_df['Kill_Percent'].max()
h = 3  # height offset for the line
ax.plot([x1, x1, x2, x2], [y_max+h, y_max+h+2, y_max+h+2, y_max+h], lw=1.2, color='black')
ax.text((x1+x2)/2, y_max+h+2, f"p = {p_val:.2f}", 
        ha='center', va='bottom', fontsize=18)
    
legend_elements = [
    Patch(facecolor='skyblue', edgecolor='black', label='Control exclusively'),
    Patch(facecolor='salmon', edgecolor='black', label='Fluconazole + Control')
]
#ax.legend(handles=legend_elements, loc='upper right', fontsize= 13)

plt.ylim(0, 105)
plt.ylabel('Relative Inhibition [%]', fontsize= 20)
plt.xlabel('')
plt.xticks([0, 1], labels= ['Control (exclusive)', 'Fluconazole + Control'], fontsize= 18)
plt.yticks(fontsize= 18)
plt.tight_layout()
plt.show()
