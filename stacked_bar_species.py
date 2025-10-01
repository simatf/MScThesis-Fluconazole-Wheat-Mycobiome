import pandas as pd
import matplotlib.pyplot as plt

#----------------------------------

# import
df = pd.read_csv('../list_all_mod.csv', usecols= [2, 10, 14]) # 5, 6 # 7, 8,

df['Putative species'] = df['Penicillium grouping'].where(df['Penicillium grouping'].notna(), df['Putative species'])

df = df[~df['Putative species'].str.lower().str.startswith("unknown")]

#df.drop(['Penicillium grouping'], inplace= True)
### STANDARDIZE
df["Putative species standard"] = df["Putative species"].str.replace(
    r"^Penicillium\s+(?:aff\.\s+)?",  # match "Penicillium " with optional "aff."
    "Penicillium aff. ",              # replace with standardized form
    regex=True
)

count_df = df.groupby(['Putative species standard', 'Treatment']).size().reset_index(name= 'Count')

count_df = count_df[count_df['Putative species standard'] != 'Unknown']

pivot_df = count_df.pivot(index='Treatment', columns='Putative species standard', values='Count').fillna(0)

relative_df = pivot_df.div(pivot_df.sum(axis=1), axis=0) * 100

# sort
sorted_taxa_piv = pivot_df.loc['Control'].sort_values(ascending=False).index.tolist()
pivot_df = pivot_df[sorted_taxa_piv]

sorted_taxa_rel = relative_df.loc['Control'].sort_values(ascending=False).index.tolist()
relative_df = relative_df[sorted_taxa_rel]

### plot
taxa = pivot_df.columns.tolist()
colors = plt.cm.get_cmap('tab20', 17).colors

ax = relative_df.plot(kind='bar', stacked=True, color= colors, edgecolor='black', figsize=(10, 6))
#ax = pivot_df.plot(kind='bar', stacked=True, color= colors, edgecolor='black', figsize=(10, 6))

####################
for container in ax.containers:
    labels = [f"{w:.0f}%" if w > 0 else "" for w in container.datavalues]
    ax.bar_label(container, labels=labels, label_type='center', fontsize=9)
####################
# for container in ax.containers:
#    labels = [f"{int(w)}" if w > 0 else "" for w in container.datavalues]
#    ax.bar_label(container, labels=labels, label_type='center', fontsize=9)
####################

plt.xlim(-0.25, 1.25)
plt.xlabel('Treatment', fontsize=14)
plt.ylabel('Relative Abundance [%]', fontsize=14)
#plt.ylabel('Isolate Count (N)', fontsize=14)

plt.xticks(ticks= [0, 1], labels= ['Control', 'Fluconazole'], rotation=0)

# Get handles and labels from the current plot
handles, labels = plt.gca().get_legend_handles_labels()

# Reverse them
plt.legend(handles[::-1], labels[::-1])

plt.tight_layout()
plt.gcf().set_size_inches(13, 8)
#plt.savefig("final/mycobiome_stacked_bar_relative_final.svg")
#plt.savefig("final/mycobiome_stacked_bar_absolute_counts_final.svg")
plt.show()