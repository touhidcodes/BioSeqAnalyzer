
# 🧬 BioSeqAnalyzer

Biological Sequence Analysis and Visualization with Python


## 📖 Table of Contents

- [About the Project](#about-the-project)
- [Features](#features)
- [Project Overview](#project_overview)
- [Library Uses](#library_uses)
- [Installation](#installation)
- [Running Script](#running_script)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## 📃 About the Project

BioSeqAnalyzer is a Python-based tool for analyzing and visualizing biological sequences from FASTA files. It extracts essential sequence properties, calculates GC content, nucleotide composition, and sequence length statistics, and generates insightful plots. Additionally, it compiles results into a detailed PDF report for easy interpretation.


## 📊 Project Overview

BioSeqAnalyzer is a bioinformatics tool designed to analyze DNA, RNA, and protein sequences efficiently. It provides functionalities for sequence manipulation, statistical analysis, and visualization to assist researchers and bioinformaticians in genomic studies.

The tool leverages NumPy, Pandas, Matplotlib, Seaborn, Biopython, and FPDF to process biological data, generate insights, and export structured reports.

### Use Cases: 
- **Genomics research:** Analyze DNA sequences to understand GC content and sequence properties.

- **Bioinformatics education:** A practical tool for learning sequence analysis concepts.

- **Comparative genomics:** Identify variations in GC content and sequence lengths across different samples.


## 🔥 Features  

1. **FASTA File Selection:** 
- Automatically detects available FASTA files.

2. **Sequence Analysis:**
- Calculates GC content (%)
- Determines nucleotide composition (A, T, G, C percentages)
- Computes sequence length statistics
- Generates reverse complement and RNA transcription

3. **Statistical Analysis:**  
- Mean, Median, Variance, and Standard Deviation for GC content and sequence lengths

4. **Data Visualization:**  
- GC Content Distribution
- Sequence Length Distribution
- GC Content vs. Sequence Length Scatter Plot
- Pair Plot of GC Content and Sequence Length.  

5. **PDF Report Generation:**  
- Includes statistical results and visual plots for better insights

 
## 📈 Library Uses

1. **NumPy (numpy):**
- Used for numerical computing, handling large arrays, and performing mathematical operations efficiently.

2. **Pandas (pandas):**
- Used for data manipulation and analysis, providing data structures like DataFrames and Series for handling tabular data.

3. **Matplotlib (matplotlib):**
- Used for data visualization, allowing the creation of plots, charts, and graphs to analyze data.

4. **Seaborn (seaborn):**
- A statistical data visualization library built on top of Matplotlib, used for making informative and attractive graphs.

5. **Biopython (biopython):**
- Used for computational biology and bioinformatics, including DNA sequence analysis, protein structure handling, and interaction with biological databases.

6. **FPDF (fpdf):**
- Used for generating PDF documents, allowing the creation of custom reports, tables, and formatted text in PDF format.  
## 📦 Installation

Before running the script, ensure you have Python installed along with the required libraries. You can install them using the following command:

```bash
pip install numpy pandas matplotlib seaborn biopython fpdf
```

## 🚀 Running Script

1. Clone the repository or download the BioSeqAnalyzer.py file and FASTA files.

```bash
git clone https://github.com/touhidcodes/BioSeqAnalyzer.git
cd BioSeqAnalyzer
```

2. Ensure the necessary libraries are installed (see [Installation](#installation) section).

3. Place your FASTA file in the project directory.

4. Run the script

```bash
python BioSeqAnalyzer.py
```
5. Select a FASTA file when prompted.

6. The script will analyze the sequences, generate plots, and create a PDF report in a folder named after the FASTA file.
## 🤝 Contributing

Contributions are what make the open-source community an amazing place! 

### Steps to contribute:
  - Fork the Project
  - Create a branch (`git checkout -b feature/AmazingFeature`)
  - Commit changes (`git commit -m 'Add some AmazingFeature'`)
  - Push the branch (`git push origin feature/AmazingFeature`)
  - Open a Pull Request

## 📬 Contact 

Touhidur Zaman - [@touhidcodes](https://www.linkedin.com/in/touhidur-zaman/) - touhidcodes@gmail.com  

Project Link: [BioSeqAnalyzer](https://github.com/touhidcodes/BioSeqAnalyzer)
