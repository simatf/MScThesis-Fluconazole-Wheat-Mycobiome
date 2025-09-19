import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import shapiro, levene, ttest_ind, mannwhitneyu
from matplotlib.patches import Patch

# Load data
df = pd.read_csv('data.csv')

# Normalize species names (remove 'aff.')
# df['Species_Grouped'] = df['Species'].str.replace(r'\baff\.\s*', '', regex=True)
df["Species_Grouped"] = df["Species"].str.replace(
    r"^(P.)\s+(?!aff\.)",  # match 'Penicillium ' not followed by 'aff.'
    r"\1 aff. ",                    # insert 'aff.' after genus
    regex=True
)

df.loc[df["Species"].str.strip().eq("P. bialowiezense"), "Species_Grouped"] = "P. bialowiezense"


# Filter only DMSO and 5Fluconazole375
dmso_df = df[df['Group'] == '0DMSO']
flu375_df = df[df['Group'] == '0Fluconazole3']

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
        group_color = 'blue'
    else:
        group_color = 'red'
    
    for val in combined_data:
        grouped.append({
            'Species_Grouped': species,
            'Kill_Percent': val,
            'Color': group_color
        })

plot_df = pd.DataFrame(grouped)

# Plot
plt.figure(figsize=(12, 6))

# Define palette mapping
palette_map = {
    species: 'skyblue' if color == 'blue' else 'salmon'
    for species, color in plot_df.groupby('Species_Grouped')['Color'].first().items()
}

ax = sns.boxplot(
    x='Species_Grouped', 
    y='Kill_Percent', 
    data=plot_df, 
    palette=palette_map,
    showfliers=False
)
sns.stripplot(
    x='Species_Grouped',
    y='Kill_Percent',
    data=plot_df,
    color='black',
    jitter=True,
    size=4
)

# # Stats: Control vs Fluconazole within species
# y_base = plot_df['Kill_Percent'].max() + 5
# pad = 5

# for species, sub_df in plot_df.groupby('Species_Grouped'):
#     # Re-merge source info for proper testing
#     srcs = species_data[species_data['Species_Grouped'] == species]
#     d1 = srcs[srcs['Source'] == 'Control']['Kill_Percent'].dropna().astype(float)
#     d2 = srcs[srcs['Source'] == 'Fluconazole']['Kill_Percent'].dropna().astype(float)
    
#     if len(d1) >= 2 and len(d2) >= 2:
#         normal1 = shapiro(d1)[1] > 0.05 if len(d1) >= 3 else False
#         normal2 = shapiro(d2)[1] > 0.05 if len(d2) >= 3 else False
#         equal_var = levene(d1, d2)[1] > 0.05
        
#         if normal1 and normal2:
#             stat, p = ttest_ind(d1, d2, equal_var=equal_var)
#         else:
#             stat, p = mannwhitneyu(d1, d2, alternative='two-sided')
        
#         xpos = list(plot_df['Species_Grouped'].unique()).index(species)
#         ax.text(xpos, y_base, f'p={p:.3g}', ha='center', fontsize=8)
    
legend_elements = [
    Patch(facecolor='skyblue', edgecolor='black', label='Control (exclusive)'),
    Patch(facecolor='salmon', edgecolor='black', label='Fluconazole + Control')
]
ax.legend(handles=legend_elements, loc='upper right')

plt.ylim(-20, 50)
plt.ylabel('Relative Inhibition [%]')
plt.xlabel('')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()
