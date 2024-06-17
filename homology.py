import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def gene_associations(sum_df, species, hist_folder):

    sum_df['Num_Associated_Species'] = sum_df['Associated_Species'].apply(lambda x: len(x.split('; ')))

    species_count_distribution = sum_df['Num_Associated_Species'].value_counts().sort_index()

    plt.figure(figsize=(10, 8))
    bars = species_count_distribution.plot(kind='bar', color='steelblue')
    plt.xlabel(f'Number of Associated Species (apart from H. {species.lower()}) ', labelpad=20, fontsize=16)
    plt.ylabel('Total Gene Count', labelpad=20, fontsize=16)
    plt.title(f'Number of Associations for H. {species.lower} genes', fontsize=16)
    plt.xticks(rotation = 0, ha ='right', fontsize = 14)
    plt.yticks(ha = "center", fontsize = 14)
    plt.tick_params(axis='x', pad=10)
    plt.tick_params(axis='y', pad= -10)
    plt.tight_layout()

    for bar in bars.patches:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, yval + 5, int(yval), ha='center', va='bottom')

    hist_file_path = os.path.join(hist_folder, f'H.{species}_gene_homology_distribution.png')
    plt.savefig(hist_file_path)
    plt.close()

def hist_similarity(gene_dict, species_name, hist_folder):

    species_counts = {}
    for gene_info in gene_dict.values():
        for associated_species in gene_info['Associated_Species']:
            if associated_species not in species_counts:
                species_counts[associated_species] = 0
            species_counts[associated_species] += 1
    
    sorted_species = sorted(species_counts.items(), key=lambda x: x[1], reverse=True)
    sorted_species = sorted_species[1:]
    
    plt.figure(figsize=(10, 8))
    species_names, counts = zip(*sorted_species)
    bars = plt.bar(species_names, counts, color = "firebrick")
    plt.xlabel('Associated Species', fontsize = 16, labelpad = 20)
    plt.ylabel(f'Number of Homologous Genes', labelpad = 15, fontsize=16)
    plt.title(f'Gene Homology Distribution for H. {species_name}', fontsize=16)
    plt.xticks(range(len(species_names)), ["H. " + str(sp) for sp in species_names], rotation = 45, ha='right', fontsize = 12, fontstyle='italic')  
  

    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, yval + 5, int(yval), ha='center', va='bottom')

    plt.tight_layout()

    hist_file_path = os.path.join(hist_folder, f'H.{species_name}_species_count_histogram.png')
    plt.savefig(hist_file_path)
    plt.close()

def wrap_labels(labels, max_species_per_line=3):
    wrapped_labels = []
    for label in labels:
        species_list = label.split(', ')
        wrapped_label = ''
        for i in range(0, len(species_list), max_species_per_line):
            if i != 0:
                wrapped_label += '\n'
            wrapped_label += ', '.join(species_list[i:i + max_species_per_line])
        wrapped_labels.append(wrapped_label)
    return wrapped_labels

def association_analysis(gene_dict, species_name, hist_folder, association_length):
    
    species_gene_count = {}

    for gene, info in gene_dict.items():
        associated_species = info['Associated_Species']
        if len(associated_species) == association_length:
            species_key = ', '.join(associated_species) 
            species_key = species_key.split(', ', 1)[1] 
            if species_key not in species_gene_count:  
                species_gene_count[species_key] = 0
            species_gene_count[species_key] += 1

    sorted_species_gene_count = sorted(species_gene_count.items(), key=lambda x: x[1], reverse=True)

    species = [item[0] for item in sorted_species_gene_count]
    gene_counts = [item[1] for item in sorted_species_gene_count]
    wrapped_labels = wrap_labels(species)


    plt.figure(figsize=(10, 8))
    bars = plt.bar(species, gene_counts, color='steelblue')
    plt.xlabel(f'Associated Species (apart from H. {species_name})', labelpad=20, fontsize=14)
    plt.ylabel('Gene Count', labelpad=20, fontsize=14)
    plt.title(f'{association_length}- Species Associations for H. {species_name} Genes', fontsize=16)

    plt.xticks(range(len(wrapped_labels)), wrapped_labels, rotation=45, ha='right', rotation_mode='anchor', fontstyle = 'italic')

    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, yval, int(yval), ha='center', va='bottom')  # Adjust y position

    plt.tight_layout()
    hist_file_path = os.path.join(hist_folder, f'H.{species_name}_{int(association_length)-1}_species_association_gene_count_histogram.png')
    plt.savefig(hist_file_path)
    plt.close()

def association_distribution(gene_dict, species_name, hist_folder):
    association_counts = {}

    for gene, info in gene_dict.items():
        associated_species = info['Associated_Species']
        assoc_length = len(associated_species)
        if assoc_length not in association_counts:
            association_counts[assoc_length] = 0
        association_counts[assoc_length] += 1

    association_lengths = sorted(association_counts.keys())
    gene_counts = [association_counts[length] for length in association_lengths]

    plt.figure(figsize=(10, 8))
    bars = plt.bar(association_lengths, gene_counts, color='firebrick')
    plt.xlabel(f'Number of Associated Species (apart from H. {species_name})', labelpad=20, fontsize=16)
    plt.ylabel('Gene Count', labelpad=20, fontsize=16)
    plt.title(f'Gene Distribution Across Association Types for {species_name}', fontsize=16)

    plt.xticks(association_lengths, labels=[f'{int(length) - 1}- species' for length in association_lengths], rotation=45, ha='right')

    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, yval - 0.5, int(yval), ha='center', va='bottom')  # Adjust y position

    plt.tight_layout()
    dist_file_path = os.path.join(hist_folder, f'H.{species_name}_gene_distribution_across_association_types.png')
    plt.savefig(dist_file_path)
    plt.close()

def count_gene_occurrences(file_path, hist_folder, output_file):
    xls = pd.ExcelFile(file_path)

    for sheet_name in xls.sheet_names:
        gene_counts = {}
        print(f"Processing sheet: {sheet_name}")
        if sheet_name == "clermontiae":
            sheet_df = pd.read_excel(xls, sheet_name=sheet_name, skiprows = 2, usecols=lambda x: x not in [0])
            sheet_df = sheet_df.iloc[:, 1:]
        else:
            sheet_df = pd.read_excel(xls, sheet_name=sheet_name, skiprows = 1, usecols=lambda x: x not in [0])
    
        for species in sheet_df.columns:
            species_counts = sheet_df[species].value_counts()
            for gene, count in species_counts.items():
                if gene not in gene_counts:
                    gene_counts[gene] = {'Total_Count': 1, 'Associated_Species': [sheet_name.lower()]}
                if species not in gene_counts[gene]:
                    gene_counts[gene][species] = 0
                gene_counts[gene][species] += count
                gene_counts[gene]['Total_Count'] += count
                if species not in gene_counts[gene]['Associated_Species']:
                    gene_counts[gene]['Associated_Species'].append(species.lower())
        
        summary_df = pd.DataFrame(gene_counts).T.fillna(0)
        summary_df['Associated_Species'] = summary_df['Associated_Species'].apply(lambda x: '; '.join(x))

        sheet_folder = os.path.join(hist_folder, sheet_name)
        if not os.path.exists(sheet_folder):
            os.makedirs(sheet_folder)

        association_distribution(gene_counts, sheet_name, sheet_folder)

        hist_similarity(gene_counts, sheet_name, sheet_folder)

        for length in range(2, 9):
            association_analysis(gene_counts, sheet_name, sheet_folder, length)


        print(f"H. {sheet_name}: Results successfully saved")

        
        with pd.ExcelWriter(output_file, mode='a', engine='openpyxl', if_sheet_exists='replace') as writer:
            summary_df.to_excel(writer, sheet_name=f'{sheet_name}_results', index=True)

    return #THE FUNCTION DOES NOT RETURN ANYTHING



