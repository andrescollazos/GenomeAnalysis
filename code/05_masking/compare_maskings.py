import matplotlib.pyplot as plt

# ---- Data ----

categories = [
    "Retroelements",
    "DNA transposons",
    "Unclassified",
    "Simple repeats",
    "Low complexity",
    "Non-repetitive"
]

custom = [9.07, 2.19, 9.99, 0.71, 0.19, 77.86]
bryophyta = [1.55, 0.45, 0.00, 0.86, 0.25, 96.86]

colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99','#c2c2f0','#d9d9d9']

# ---- Pie Charts ----

fig, axes = plt.subplots(1, 2, figsize=(12, 6))

axes[0].pie(custom, colors=colors, autopct='%1.1f%%', startangle=90)
axes[0].set_title("Custom Library")

axes[1].pie(
    bryophyta,
    colors=colors,
    autopct=lambda p: f'{p:.1f}%' if p > 0.5 else '', # Only labels > 0.5%
    startangle=90
)
axes[1].set_title("Bryophyta Library")

# Legend to the right
fig.legend(categories, loc='center right', bbox_to_anchor=(1.15, 0.5), title="Repeat Categories")

plt.suptitle("Genome Composition by Repeat Type", fontsize=14)
plt.tight_layout()
plt.show()

# ---- Bar Plot ----

x = range(len(categories) - 1)  # Don't include non-repetitive
bar_width = 0.35

fig, ax = plt.subplots(figsize=(10, 6))

ax.bar([i - bar_width/2 for i in x], custom[:-1], width=bar_width, label='Custom', color='#4c72b0')
ax.bar([i + bar_width/2 for i in x], bryophyta[:-1], width=bar_width, label='Bryophyta', color='#55a868')

ax.set_xticks(x)
ax.set_xticklabels(categories[:-1], rotation=45, ha='right')
ax.set_ylabel('Percentage of Genome')
ax.set_title('Comparison of Masked Repeat Content')
ax.legend()
plt.tight_layout()
plt.show()
