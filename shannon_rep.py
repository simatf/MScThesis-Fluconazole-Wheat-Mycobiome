import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

#-------------------------------------------------

df = pd.read_csv('shannon_rem_dup_no_isolates.csv')

# Remove "Unknown"
df = df[~df["Putative species"].str.lower().str.startswith("unknown")]

# Count occurrences per Putative species × Treatment × Replicate
counts_df = (
    df.groupby(["Putative species", "Treatment", "Replicate"])
      .size()
      .reset_index(name="Count")
)

# Pivot so columns are "Treatment_repX"
wide_df = counts_df.pivot_table(
    index="Putative species",
    columns=["Treatment", "Replicate"],
    values="Count",
    fill_value=0
)
wide_df.columns = [f"{treat}_rep{rep}" for treat, rep in wide_df.columns]
wide_df = wide_df.reset_index()

wide_df.to_csv('wide_count_df4.csv', index= False)

# Shannon diversity function
def shannon_index(counts):
    counts = np.array(counts)
    counts = counts[counts > 0]
    p = counts / counts.sum()
    return -np.sum(p * np.log(p))

# Compute Shannon index per replicate
shannon_per_rep = []
for col in wide_df.columns:
    if col != "Putative species":
        H = shannon_index(wide_df[col].values)
        treatment = col.split("_")[0]
        shannon_per_rep.append({"Treatment": treatment, "Replicate": col, "Shannon": H})

shannon_df = pd.DataFrame(shannon_per_rep)

# Separate by treatment
treatments = shannon_df["Treatment"].unique()
group1 = shannon_df.loc[shannon_df["Treatment"] == treatments[0], "Shannon"]
group2 = shannon_df.loc[shannon_df["Treatment"] == treatments[1], "Shannon"]

# Statistical test (two-sample t-test here)
t_stat, p_val = ttest_ind(group1, group2, equal_var=False)

# Plot
fig, ax = plt.subplots()
bp = ax.boxplot([group1, group2], labels=treatments, patch_artist=True)
colors = ['skyblue', 'mediumorchid']
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
for median in bp['medians']:
    median.set_color("black")
    median.set_linewidth(2)
# Show individual data points with jitter
for i, treatment in enumerate(treatments, start=1):
    y_vals = shannon_df.loc[shannon_df["Treatment"] == treatment, "Shannon"]
    x_vals = np.random.normal(i, 0.04, size=len(y_vals))  # jitter around box
    ax.scatter(x_vals, y_vals, color="black", zorder=3)
ax.set_ylabel("Shannon Diversity Index")
#ax.set_title("Shannon Diversity by Treatment")

# Annotate p-value
x1, x2 = 1, 2
y_max = max(shannon_df["Shannon"])
y, h = y_max + 0.05, 0.02
ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, color="black")

if p_val < 0.001:
    p_text = "p < 0.001"
else:
    p_text = f"p = {p_val:.3f}"

ax.text((x1 + x2) / 2, y + h + 0.001, p_text, ha="center", va="bottom", fontsize= 14)

plt.ylabel('Shannon Diversity Index', fontsize= 16)
plt.xlabel('Treatment', fontsize= 16)
plt.xticks(ticks= [1, 2], labels= ['Control', 'Fluconazole'], fontsize= 14)
plt.yticks(ticks= plt.yticks()[0], fontsize= 14)
plt.show()