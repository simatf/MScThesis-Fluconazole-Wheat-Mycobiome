import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

#----------------------------------

# import
df = pd.read_csv('../list_all_mod.csv', usecols= [2, 10, 14]) # 5, 6 # 7, 8,

df['Putative species'] = df['Penicillium grouping'].where(df['Penicillium grouping'].notna(), df['Putative species'])

df = df[~df['Putative species'].str.lower().str.startswith("unknown")]

df["Putative species"] = df["Putative species"].str.replace(
    r"^Penicillium\s+(?:aff\.\s+)?",  # match "Penicillium " with optional "aff."
    "Penicillium aff. ",              # replace with standardized form
    regex=True
)

# sets
control_set = set(df[df['Treatment'] == 'Control']['Putative species'])
fungicide_set = set(df[df['Treatment'] == 'Fluconazole']['Putative species'])

### plot
plt.figure(figsize=(6,6))
out_1 = venn2(
    [control_set, fungicide_set],
)

# Set left circle color
out_1.get_patch_by_id('10').set_facecolor('skyblue')
# Set right circle color
out_1.get_patch_by_id('01').set_facecolor('mediumorchid')
# Set overlap color separately
out_1.get_patch_by_id('11').set_facecolor('salmon')

for label in out_1.set_labels:
    if label:
        label.set_text('')
for text in out_1.subset_labels:
    text.set_fontsize(13)

#plt.savefig('final/venn_final.svg')
plt.show()

