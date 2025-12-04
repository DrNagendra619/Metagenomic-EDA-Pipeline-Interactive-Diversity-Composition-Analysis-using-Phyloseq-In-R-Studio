# Metagenomic-EDA-Pipeline-Interactive-Diversity-Composition-Analysis-using-Phyloseq-In-R-Studio
Metagenomic EDA Pipeline Interactive Diversity  Composition Analysis using Phyloseq In R Studio
# Metagenomic EDA Pipeline: Interactive Diversity & Composition Analysis

This repository contains an automated R pipeline for Exploratory Data Analysis (EDA) of metagenomic/microbiome data. It utilizes the **Phyloseq** package to clean, normalize, and visualize microbiome census data, producing both static (publication-ready) and interactive (HTML) visualizations.

## 📌 Features

* **Automated Setup:** Automatically creates output directories and manages dependencies.
* **Speed Optimized:** Filters data to the top 2,000 abundant taxa to ensure rapid execution (<1 minute on standard laptops).
* **Interactive Visualizations:** Uses `plotly` and `htmlwidgets` to generate zoomable, hover-able plots.
* **Comprehensive Workflow:**
    * **Quality Control:** Rarefaction curves.
    * **Alpha Diversity:** Shannon & Observed richness indices.
    * **Beta Diversity:** PCoA using Bray-Curtis dissimilarity.
    * **Composition Analysis:** Stacked bar plots of Phyla and Heatmaps of top Genera.
* **Auto-Export:** Automatically saves all results (PNG images, HTML widgets, CSV tables) to a local directory.

## 🛠️ Prerequisites

You need **R** (and preferably **RStudio**) installed on your system. 

The script depends on the following R packages:
* `phyloseq` (Bioconductor)
* `ggplot2`
* `vegan`
* `plotly`
* `DT`
* `htmlwidgets`

## 🚀 Usage

1.  **Clone this repository** or download the `.R` script.
2.  Open the script in RStudio.
3.  **Configure Output Path:**
    By default, the script saves outputs to `D:/DOWNLOADS`. You can change this line in **Section 1**:
    ```r
    output_dir <- "path/to/your/output/folder"
    ```
4.  **Run the Script:** Select all code (Ctrl+A) and run (Ctrl+Enter).

## 📊 Pipeline Steps & Outputs

The script processes the built-in `GlobalPatterns` dataset (or your own data if modified) through the following steps:

### 1. Data Loading & Optimization
* Loads the dataset.
* Filters for the top 2,000 most abundant taxa to improve processing speed.
* Normalizes counts to relative abundance.

### 2. Quality Control (Rarefaction)
* **Output:** `01_Rarefaction_Curve.png`
* Checks if the sequencing depth is sufficient to capture the microbial diversity.

### 3. Alpha Diversity Analysis
* **Output:** `02_Alpha_Diversity.png` (Static) & `.html` (Interactive)
* Visualizes Species Richness (Observed) and Evenness (Shannon) across different sample types using boxplots.

### 4. Beta Diversity Analysis (PCoA)
* **Output:** `03_Beta_Diversity_PCoA.png` (Static) & `.html` (Interactive)
* Performs Principal Coordinate Analysis (PCoA) with Bray-Curtis dissimilarity to visualize how samples cluster based on their microbial composition.

### 5. Taxonomic Composition
* **Output:** `04_Taxonomy_Barplot.png` (Static) & `.html` (Interactive)
* Displays a stacked bar chart showing the relative abundance of the top 10 Phyla across samples.

### 6. Heatmap Analysis
* **Output:** `05_Heatmap.png` (Static) & `.html` (Interactive)
* Visualizes the abundance of the top 20 Genera across different sample types.

### 7. Data Export
* **Output:** * `06_GlobalPatterns_OTU.csv`: Raw OTU count table.
    * `06_Taxonomy_Table.html`: A searchable, interactive table of the taxonomy.

## 📝 Example Code Snippet

Here is how the pipeline handles the critical Alpha Diversity step:

```r
# Alpha Diversity Plotting
p_alpha <- plot_richness(physeq_filt, x = "SampleType", measures = c("Shannon", "Observed")) + 
  geom_boxplot(aes(fill = SampleType), alpha = 0.6) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Alpha Diversity (Top 2000 Taxa)")

# Save Interactive HTML
saveWidget(ggplotly(p_alpha), file.path(output_dir, "02_Alpha_Diversity_Interactive.html"))
