from collections import defaultdict
import matplotlib.pyplot as plt

# Step 1: Load GO namespaces from go-basic.obo
# File obtained using: wget http://purl.obolibrary.org/obo/go/go-basic.obo -O go-basic.obo
go2namespace = {}
with open("go-basic.obo", "r") as f:
    current_id = None
    for line in f:
        if line.startswith("id: GO:"):
            current_id = line.strip().split()[-1]
        elif line.startswith("namespace:") and current_id:
            ns = line.strip().split()[-1]
            go2namespace[current_id] = ns
            current_id = None

# Step 2: Read GO raw counts (See line 23 in 08_functional_annotation/03_extract_statistics.sh)
namespace_counts = defaultdict(int)
with open("go_raw_counts.txt", "r") as f:
    for line in f:
        count, go_id = line.strip().split()
        count = int(count)
        if go_id in go2namespace:
            ns = go2namespace[go_id]
            namespace_counts[ns] += count
        else:
            namespace_counts["unknown"] += count

# Step 3: Plot pie chart
label_map = {
    "biological_process": "Biological Process",
    "molecular_function": "Molecular Function",
    "cellular_component": "Cellular Component",
    "unknown": "Unknown"
}

labels = [label_map[k] for k in namespace_counts]
sizes = [namespace_counts[k] for k in namespace_counts]

plt.figure(figsize=(7, 7))
plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140)
plt.title("GO Level 1 Term Distribution")
plt.axis('equal')
plt.tight_layout()
plt.show()
