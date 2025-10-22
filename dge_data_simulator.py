import pandas as pd
import numpy as np

# --- Configuration ---
NUM_GENES = 10000
NUM_SAMPLES = 12
OUTPUT_FILE = 'mock_expression_counts.csv'

samples_control = [f'Control_S{i+1}' for i in range(NUM_SAMPLES // 2)]
samples_treated = [f'Treated_S{i+1}' for i in range(NUM_SAMPLES // 2)]
sample_names = samples_control + samples_treated


condition_labels = ['Control'] * (NUM_SAMPLES // 2) + ['Treated'] * (NUM_SAMPLES // 2)

# ---  Generate Gene Names ---
gene_names = [f'Gene_{i+1}' for i in range(NUM_GENES)]

np.random.seed(42)
baseline_means = np.random.lognormal(mean=4, sigma=1.5, size=NUM_GENES)
count_data = []

for i in range(NUM_SAMPLES // 2):
    counts = np.random.poisson(baseline_means * 1.0, size=NUM_GENES)
    count_data.append(counts)

num_dge = 100
dge_indices = np.random.choice(NUM_GENES, num_dge, replace=False)
fold_changes = np.random.choice([3, 4, 5], num_dge) 

for i in range(NUM_SAMPLES // 2):
    treated_counts = []
    base_treated_means = baseline_means * 1.1 
    
    for j in range(NUM_GENES):
        mean_j = base_treated_means[j]
        
        if j in dge_indices:

            if np.random.rand() < 0.5:
                 mean_j *= fold_changes[np.where(dge_indices == j)[0][0]]
            else:
                 mean_j /= fold_changes[np.where(dge_indices == j)[0][0]]
        
        treated_counts.append(np.random.poisson(mean_j, size=1)[0])
        
    count_data.append(np.array(treated_counts))

count_matrix = np.stack(count_data, axis=1)

df_counts = pd.DataFrame(count_matrix, index=gene_names, columns=sample_names)

df_counts.to_csv(OUTPUT_FILE, index=True)

df_metadata = pd.DataFrame({
    'SampleName': sample_names,
    'Condition': condition_labels
})
df_metadata.set_index('SampleName', inplace=True)
df_metadata.to_csv('sample_metadata.csv', index=True)

print(f"Successfully generated DGE count matrix: {OUTPUT_FILE}")
print("Metadata file saved: sample_metadata.csv")
