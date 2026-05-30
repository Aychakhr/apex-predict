FROM continuumio/miniconda3

# Install system dependencies
RUN apt-get update && apt-get install -y \
    perl \
    wget \
    curl \
    ncbi-blast+ \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install bioinformatics tools via conda
RUN conda install -c bioconda -c conda-forge \
    mlst \
    abricate \
    ectyper \
    -y

# Install Python dependencies
COPY requirements.txt .
RUN pip install -r requirements.txt

# Copy app files
WORKDIR /app
COPY . /app

# Make abricate databases available
RUN abricate --setupdb

EXPOSE 8501

CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]