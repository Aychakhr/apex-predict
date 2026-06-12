FROM continuumio/miniconda3:latest

# ── System dependencies ────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y \
    perl \
    wget \
    curl \
    git \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# ── Bioinformatics environment (apec_a equivalent) ─────────────────────────
RUN conda create -n apec_a -c bioconda -c conda-forge -y \
    python=3.6 \
    mlst \
    abricate \
    ectyper \
    blast \
    gsl=2.7 \
    mash

# Fix mash/gsl library mismatch (libgsl.so.27 -> libgsl.so.25)
RUN ln -sf /opt/conda/envs/apec_a/lib/libgsl.so.27 /opt/conda/envs/apec_a/lib/libgsl.so.25

# Patch ECTyper's get_valid_format bug (fastq exception swallows fasta check)
RUN python3 - << 'EOF'
path = "/opt/conda/envs/apec_a/lib/python3.6/site-packages/ectyper/genomeFunctions.py"
content = open(path).read()

old = """    for fm in ['fastq', 'fasta']:
        try:
            with open(file, "r") as handle:
                data = SeqIO.parse(handle, fm)
                if any(data):
                    if is_tarfile(file):
                        LOG.warning("Compressed file is not supported: {}".format(file))
                        return None
                    return fm
        except FileNotFoundError as err:
            LOG.warning("{0} is not found".format(file))
            return None
        except UnicodeDecodeError as err:
            LOG.warning("{0} is not a valid file".format(file))
            return None
        except:
            LOG.warning("{0} is an unexpected file".format(file))
            return None
    LOG.warning("{0} is not a fasta/fastq file".format(file))
    return None"""

new = """    for fm in ['fastq', 'fasta']:
        try:
            with open(file, "r") as handle:
                data = SeqIO.parse(handle, fm)
                try:
                    result = any(data)
                except Exception:
                    continue
                if result:
                    if is_tarfile(file):
                        LOG.warning("Compressed file is not supported: {}".format(file))
                        return None
                    return fm
        except FileNotFoundError as err:
            LOG.warning("{0} is not found".format(file))
            return None
        except UnicodeDecodeError as err:
            LOG.warning("{0} is not a valid file".format(file))
            return None
    LOG.warning("{0} is not a fasta/fastq file".format(file))
    return None"""

if old in content:
    content = content.replace(old, new, 1)
    open(path, 'w').write(content)
    print("ECTyper patched successfully")
else:
    print("WARNING: ECTyper patch target not found - may already be patched or version differs")
EOF

# ── ML / Web app environment  ─────────────────────────
RUN conda create -n apec_new -c conda-forge -y \
    python=3.10 \
    streamlit \
    pandas \
    scikit-learn \
    joblib \
    requests \
    numpy

# ── App setup ───────────────────────────────────────────────────────────────
WORKDIR /app
COPY . /app

# Point app.py to the bioinformatics env's bin directory
ENV BIOTOOLS_BIN=/opt/conda/envs/apec_a/bin

EXPOSE 8501

# Run streamlit using the apec_new environment
SHELL ["/bin/bash", "-c"]
CMD source activate apec_new && \
    streamlit run app.py --server.port=8501 --server.address=0.0.0.0
