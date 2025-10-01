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
flu375_df = df[df['Group'] == '0Fluconazole3']

# Calculate DMSO mean per isolate
dmso_means = dmso_df.groupby('Isolate')['Value'].mean().rename('DMSO_Mean')

# Merge and compute kill %
flu375_df = flu375_df.merge(dmso_means, on='Isolate', how='left')
flu375_df['Kill_Percent'] = (1 - (flu375_df['Value'] / flu375_df['DMSO_Mean'])) * 100

# Keep species, kill percent, and source
species_data = flu375_df[['Species', 'Species_Grouped', 'Kill_Percent', 'Source']].copy()

# --- Build grouped dataset ---
grouped = []
for species, sub_df in species_data.groupby('Species_Grouped'):
    combined_data = sub_df['Kill_Percent'].dropna().tolist()
    sources = sub_df['Source'].unique()

    if len(sources) == 1 and sources[0] == 'Control':
        if species in ['Ascochyta phacae', 'Ramularia sp.']:
            group = 'Control (plant pathogens)'
        else:
            group = 'Control (non-plant pathogens)'
    else:
        group = 'Fluconazole + Control'
    
    for val in combined_data:
        grouped.append({
            'Species_Grouped': species,
            'Kill_Percent': val,
            'Group': group
        })

plot_df = pd.DataFrame(grouped)

# --- Define groups ---
groups = {
    "Control (plant pathogens)": plot_df[plot_df["Group"] == "Control (plant pathogens)"]["Kill_Percent"].dropna().values,
    "Control (non-plant pathogens)": plot_df[plot_df["Group"] == "Control (non-plant pathogens)"]["Kill_Percent"].dropna().values,
    "Fluconazole + Control": plot_df[plot_df["Group"] == "Fluconazole + Control"]["Kill_Percent"].dropna().values,
}

# --- Statistical test helper ---
def best_test(x, y):
    if len(x) >= 3 and len(y) >= 3:
        normal_x = shapiro(x)[1] > 0.05
        normal_y = shapiro(y)[1] > 0.05
    else:
        normal_x = normal_y = False
    
    equal_var = levene(x, y)[1] > 0.05

    if normal_x and normal_y:
        stat, p = ttest_ind(x, y, equal_var=equal_var)
        return "t-test", p
    else:
        stat, p = mannwhitneyu(x, y, alternative="two-sided")
        return "Mann-Whitney U", p

# Perform tests
test1, p1 = best_test(groups["Control (plant pathogens)"], groups["Fluconazole + Control"])
test2, p2 = best_test(groups["Control (non-plant pathogens)"], groups["Fluconazole + Control"])
print(f"{test1} (plant vs fluc): p={p1:.3g}")
print(f"{test2} (non-plant vs fluc): p={p2:.3g}")

# --- Plot ---
plt.figure(figsize=(12, 6))

palette_map = {
    "Control (plant pathogens)": "skyblue",
    "Control (non-plant pathogens)": "skyblue",
    "Fluconazole + Control": "salmon"
}
order_list = ["Control (plant pathogens)", "Control (non-plant pathogens)", "Fluconazole + Control"]

ax = sns.boxplot(
    x="Group", y="Kill_Percent", data=plot_df,
    palette=palette_map, showfliers=False, order=order_list
)
sns.stripplot(
    x="Group", y="Kill_Percent", data=plot_df,
    color="black", jitter=True, size=4, order=order_list
)

# --- Add statistical annotations ---
y_max = plot_df["Kill_Percent"].max()
h = 5  # vertical spacing

# # Plant vs Fluconazole
# x1, x2 = 0, 2  # positions of the boxes on x-axis
# y_max = plot_df['Kill_Percent'].max()
# h = -1  # height offset for the line
# ax.plot([x1, x1, x2, x2], [y_max+h+4, y_max+h+6, y_max+h+6, y_max+h+4], lw=1.2, color='black')
# ax.text((x1+x2)/2, y_max+h+6, f"p = {p1:.3g} (Wilcoxon test)", 
#         ha='center', va='bottom', fontsize=10)

# Non-plant vs Fluconazole
# Add annotation above the two boxes
x1, x2 = 1, 2  # positions of the boxes on x-axis
y_max = plot_df['Kill_Percent'].max()
h = -8  # height offset for the line
ax.plot([x1, x1, x2, x2], [y_max+h+4, y_max+h+6, y_max+h+6, y_max+h+4], lw=1.2, color='black')
ax.text((x1+x2)/2, y_max+h+6, f"p = {p2:.2g}", 
        ha='center', va='bottom', fontsize=14)

plt.ylim(-20, 60)
plt.ylabel('Relative Inhibition [%]', fontsize= 16)
plt.xlabel('')
plt.yticks(fontsize= 14)
plt.xticks(fontsize= 14)
plt.tight_layout()
plt.show()
