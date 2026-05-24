# 🦠 ApexPredict

**Pathotype Classification, MLST, Serogroup Detection, and Resistance Profiling for *APEC***

##  Overview

ApexPredict is a comprehensive bioinformatics tool for analyzing bacterial genome sequences (FASTA files) with a focus on **Avian Pathogenic *E. coli*** (APEC). The tool integrates multiple analyses:

- **MLST Typing**: Sequence type determination and high-risk ST identification
- **Serogroup Detection**: ECTyper-based O-antigen serogrouping
- **Virulence Profiling**: Detection of virulence factors
- **AMR Profiling**: Antimicrobial resistance gene detection
- **Machine Learning**: APEC pathotype prediction with confidence scoring

##  Features

- **One-click Analysis**: Upload FASTA, get comprehensive results
- **Real-time Processing**: MLST, serogroup, virulence, and AMR detection
- **Risk Assessment**: High-risk MLST and APEC-associated serogroup identification
- **Visual Analytics**: Interactive tables and metrics display
- **Export Results**: Download virulence profiles as CSV

## Databases & Tools

- **APEC virulence database** — sourced from APECTyper (https://github.com/JohnsonSingerLab/APECtyper), 
  originally described in: Johnson et al. (2022). Refining the definition of the avian pathogenic Escherichia coli (APEC) pathotype through inclusion of high-risk clonal groups. Poultry Science, 101(10), 102009. https://doi.org/10.1016/j.psj.2022.102009

- **ECTyper** — https://github.com/phac-nml/ectyper
- **MLST** — https://github.com/tseemann/mlst
- **Abricate** — https://github.com/tseemann/abricate

### Prerequisites

Install the required command-line bioinformatics tools to your system or Conda environment:

```bash
# Install MLST and Abricate via Bioconda
conda install -c bioconda mlst abricate

# Install ECTyper
conda install -c bioconda ectyper
# OR via pip if preferred: pip install ectyper

# Update Abricate's default databases (NCBI, VirulenceFinder, etc.)
abricate --update
```` `

### Local Deployment
````bash
# Clone the repository
git clone https://github.com/Aychakhr/ApexPredict.git
cd ApexPredict

# Install Python package dependencies
pip install -r requirements.txt

# Place your model files in the directory
# - apec_random_forest_model.joblib
# - model_features.joblib
# - apec_db/

# Run the app
streamlit run app.py