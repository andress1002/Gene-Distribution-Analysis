import os
import homology as homology
import pandas as pd

def main():

    hist_folder = r'C:\Users\andre\Desktop\MCF\MCF_Project\Gene_Homolgy\Histograms' 
    df_file = r'C:\Users\andre\Desktop\MCF\MCF_Project\Gene_Homolgy\Similar_Genes.xlsx'
    output_file = r'C:\Users\andre\Desktop\MCF\MCF_Project\Gene_Homolgy\Similarity_Results2.xlsx'

    homology.count_gene_occurrences(df_file, hist_folder, output_file)

if __name__ == "__main__":
    main()
