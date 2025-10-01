import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu  # or ttest_ind

df = pd.read_csv('Shannon/wide_count_df4.csv')

# rename
df['Putative species'] = df['Putative species'].str.replace(
    r"^(Penicillium)\s+(?!aff\.)",  # match 'Penicillium ' not followed by 'aff.'
    r"\1 aff. ",                    # insert 'aff.' after genus
    regex=True
)

# ----------------
# 1. Convert to relative abundance per replicate
# ----------------
numeric_df = df.set_index("Putative species")
col_sums = numeric_df.sum(axis=0)
rel_abund = numeric_df.div(col_sums, axis=1)

# ----------------
# 2. Reshape to long format
# ----------------
long_df = rel_abund.reset_index().melt(
    id_vars="Putative species",
    var_name="Sample",
    value_name="Relative abundance"
)
long_df["Treatment"] = long_df["Sample"].apply(
    lambda x: "Control" if "Control" in x else "Fluconazole"
)

# ----------------
# 3. Plot each species separately
# ----------------

for species in long_df["Putative species"].unique():
    sub = long_df[long_df["Putative species"] == species]

    plt.figure(figsize=(5, 4))
    
    palette = {"Control": "skyblue", "Fluconazole": "mediumorchid"}

    # Strip/box for Control (always plotted)
    sns.stripplot(
        data=sub[sub["Treatment"] == "Control"],
        x="Treatment", y="Relative abundance",
        size=8, jitter=True, color="black", alpha=0.8
    )
    sns.boxplot(
        data=sub[sub["Treatment"] == "Control"],
        x="Treatment", y="Relative abundance",
        palette=palette
    )

    # Extract Fluconazole values
    flu = sub[sub["Treatment"] == "Fluconazole"]["Relative abundance"]
    ctrl = sub[sub["Treatment"] == "Control"]["Relative abundance"]

    if (flu > 0).any():  
        # Plot Fluconazole if not all zero
        sns.stripplot(
            data=sub[sub["Treatment"] == "Fluconazole"],
            x="Treatment", y="Relative abundance",
            size=8, jitter=True, color="black", alpha=0.8
        )
        sns.boxplot(
            data=sub[sub["Treatment"] == "Fluconazole"],
            x="Treatment", y="Relative abundance",
            palette=palette
        )

        # Statistical test
        stat, p = mannwhitneyu(ctrl, flu)
        y_max = sub["Relative abundance"].max()
        plt.plot([0, 1], [y_max * 1.05, y_max * 1.05], color="black")
        plt.text(0.5, y_max*1.01, f"p = {p:.3f}", ha="center", fontsize= 18)

    else:
        sns.stripplot(
            data=sub[sub["Treatment"] == "Fluconazole"],
            x="Treatment", y="Relative abundance",
            size=8, jitter=True, color="black", alpha=0
        )
        # Reserve the space but replace plot with "NA"
        y_mid = sub["Relative abundance"].max() * 0.5  # halfway up plot
        plt.text(1, y_mid, "NA", ha="center", va="center", color="black", fontsize=18)
# for species in long_df["Putative species"].unique():
#     sub = long_df[long_df["Putative species"] == species]

#     # Create figure
#     plt.figure(figsize=(5, 4))
    
#     palette = {"Control": "skyblue", "Fluconazole": "mediumorchid"}

#     # Scatter points
#     sns.stripplot(data=sub, x="Treatment", y="Relative abundance", size=8, jitter=True, color= 'black', alpha= 0.8)
#     # Transparent boxplot overlay
#     sns.boxplot(data=sub, x="Treatment", y="Relative abundance", palette= palette)

#     # Statistical test
#     ctrl = sub[sub["Treatment"] == "Control"]["Relative abundance"]
#     flu = sub[sub["Treatment"] == "Fluconazole"]["Relative abundance"]
#     stat, p = mannwhitneyu(ctrl, flu)

#     # Annotation
#     y_max = sub["Relative abundance"].max()
#     plt.plot([0, 1], [0.78, 0.78], color="black")
#     plt.text(0.5, 0.8, f"p = {p:.3f}", ha="center")

    plt.title(species, fontsize= 20)
    plt.ylabel("Relative Abundance", fontsize= 18)
    plt.xlabel("Treatment", fontsize= 18)
    plt.ylim(bottom= -0.02)
    plt.xticks(ticks= [0, 1], labels= ['Control', 'Fluconazole'], fontsize= 16)
    plt.yticks(fontsize= 16)
    plt.gcf().set_size_inches(5, 12)
    plt.tight_layout()
    plt.savefig(f"NEW_box_each_species/NEW/000{species}.png", dpi=300)
    #plt.show()
