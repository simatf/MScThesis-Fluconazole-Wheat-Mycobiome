import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, ttest_ind, levene, shapiro
import itertools
import numpy as np
from matplotlib.patches import Patch

# Load the data
file_path = 'data.csv'
df = pd.read_csv(file_path, sep=',')

df['Species_Grouped'] = df['Species'].str.replace(r'\baff\.\s*', '', regex=True)

# Filter only DMSO and 5Fluconazole375
dmso_df = df[df['Group'] == '0DMSO']
flu375_df = df[df['Group'] == '0Fluconazole3']

# Calculate DMSO mean per isolate
dmso_means = dmso_df.groupby('Isolate')['Value'].mean().rename('DMSO_Mean')

# Merge into 375 group and compute kill %
flu375_df = flu375_df.merge(dmso_means, on='Isolate', how='left')
flu375_df['Kill_Percent'] = (1 - (flu375_df['Value'] / flu375_df['DMSO_Mean'])) * 100

# Prepare isolate labels with species prefix for clarity
flu375_df['Label'] = flu375_df['Species'] + '\n' + flu375_df['Isolate'] + ' (' + flu375_df['Source'].str[0] + ')'

# Set plotting order (grouped by species)
label_order = flu375_df.sort_values(['Species_Grouped', 'Species', 'Isolate'])['Label'].unique().tolist()

# Plot
plt.figure(figsize=(16, 7))
ax = sns.boxplot(x='Label', y='Kill_Percent', data=flu375_df,
                 hue='Source', showfliers=False, order=label_order,
                 palette={'Control': 'skyblue', 'Fluconazole': 'mediumorchid'}, legend= True)

sns.stripplot(x='Label', y='Kill_Percent', data=flu375_df, hue='Source',
              dodge=True, jitter=True, color='black',
              order=label_order, legend= False)

# # Fix legend
# handles, labels = ax.get_legend_handles_labels()
# by_label = dict(zip(labels, handles))
# ax.legend(by_label.values(), by_label.keys(), title='Source', loc='upper right')

# Set y-axis limits
plt.ylim(-20, 70)
plt.ylabel('Relative Inhibition [%]')
plt.xlabel('')
#plt.title('Kill % at 5Fluconazole375 across all Isolates (Grouped by Species)')
plt.xticks(rotation=45, ha='right')

# Add statistical comparisons within each species
y_base = 50
pad = 5
species_list = flu375_df['Species_Grouped'].unique()

for species in species_list:
    species_df = flu375_df[flu375_df['Species_Grouped'] == species]
    isolates = species_df['Isolate'].dropna().unique()

    isolate_pairs = list(itertools.combinations(isolates, 2))

    for i, (iso1, iso2) in enumerate(isolate_pairs):
        # Extract kill % for both isolates
        d1 = species_df[species_df['Isolate'] == iso1]['Kill_Percent'].dropna()
        d2 = species_df[species_df['Isolate'] == iso2]['Kill_Percent'].dropna()

        if len(d1) >= 2 and len(d2) >= 2:
            # Test normality
            normal1 = shapiro(d1)[1] > 0.05 if len(d1) >= 3 else False
            normal2 = shapiro(d2)[1] > 0.05 if len(d2) >= 3 else False

            # Test variance
            equal_var = levene(d1, d2)[1] > 0.05

            # Choose test
            if normal1 and normal2:
                if equal_var:
                    stat, p = ttest_ind(d1, d2)
                    print('t-test')
                    test= 'T-test'
                else:
                    stat, p = ttest_ind(d1, d2, equal_var=False)
                    print('t-test non-equal')
            else:
                stat, p = mannwhitneyu(d1, d2, alternative='two-sided')
                test= 'Wilcoxon'

            # Annotate inside plot area
            label1 = species_df[species_df['Isolate'] == iso1]['Label'].iloc[0]
            label2 = species_df[species_df['Isolate'] == iso2]['Label'].iloc[0]
            if label1 in label_order and label2 in label_order:
                x1 = label_order.index(label1)
                x2 = label_order.index(label2)
                y = y_base + i * pad
                ax.plot([x1, x1, x2, x2], [y, y, y, y], lw=1.2, color='gray')
                ax.text((x1 + x2) / 2, y + 0.5, f'p={p:.3g} ({test})', ha='center', fontsize=8)

# Step: Add faint dotted vertical lines between species
species_per_label = flu375_df.set_index('Label')['Species_Grouped'].to_dict()

for i in range(1, len(label_order)):
    current_species = species_per_label[label_order[i]]
    previous_species = species_per_label[label_order[i - 1]]

    if current_species != previous_species:
        ax.axvline(i - 0.5, color='gray', linestyle='dotted', linewidth=0.8, alpha=0.6)

legend_elements = [
    Patch(facecolor='skyblue', edgecolor='black', label='Control'),
    Patch(facecolor='mediumorchid', edgecolor='black', label='Fluconazole')
]
ax.legend(handles=legend_elements, loc='upper right')

plt.tight_layout()
plt.show()