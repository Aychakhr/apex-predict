import streamlit as st
import pandas as pd
import joblib
import subprocess
import os
import base64
import shutil
import re
import glob
import tempfile

# ENVIRONMENT BRIDGE
BIOTOOLS_BIN = os.environ.get(
    "BIOTOOLS_BIN",
    "/home/khera_aycha/miniconda3/envs/apec_a/bin"
)
os.environ["PATH"] = BIOTOOLS_BIN + os.pathsep + os.environ.get("PATH", "")

def get_base64_of_bin_file(bin_file):
    with open(bin_file, 'rb') as f:
        data = f.read()
    return base64.b64encode(data).decode()

def check_ectyper_available():
    try:
        result = subprocess.run("ectyper --help", shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                timeout=5)
        return result.returncode in (0, 1)
    except:
        return False

def detect_serogroup_with_ectyper(fasta_path):
    try:
        output_dir = "temp_ectyper_out"
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
            
        linux_tmp = tempfile.mkdtemp(prefix="/tmp/ectyper_in_")
        linux_fasta = os.path.join(linux_tmp, "sample.fna")
        shutil.copy2(fasta_path, linux_fasta)
        linux_output_dir = tempfile.mkdtemp(prefix="/tmp/ectyper_out_")
        shutil.rmtree(linux_output_dir)
        cmd = f"ectyper --input {linux_fasta} --output {linux_output_dir}"
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, universal_newlines=True,
                                timeout=40)

        output_path = None
        for candidate in [
            os.path.join(linux_output_dir, "output.csv"),
            os.path.join(linux_output_dir, "output.tsv"),
        ]:
            if os.path.exists(candidate) and os.path.getsize(candidate) > 0:
                output_path = candidate
                break

        if output_path is None:
            search_paths = (
                glob.glob(os.path.join(linux_output_dir, "**", "*.csv"), recursive=True) or
                glob.glob(os.path.join(linux_output_dir, "**", "*.tsv"), recursive=True)
            )
            if search_paths:
                output_path = search_paths[0]

        if output_path is not None and os.path.getsize(output_path) > 0:
            sep = '\t' if output_path.endswith('.tsv') else ','
            df_ec = pd.read_csv(output_path, sep=sep)

            if not df_ec.empty:
                serogroup = "Unknown"
                serotype = "Unknown"

                for col in ["O_prediction", "O.type", "O_type", "Serogroup"]:
                    if col in df_ec.columns:
                        serogroup = str(df_ec.iloc[0][col]).strip()
                        break

                if "Serotype" in df_ec.columns:
                    serotype = str(df_ec.iloc[0]["Serotype"]).strip()
                elif "O_prediction" in df_ec.columns and "H_prediction" in df_ec.columns:
                    o = str(df_ec.iloc[0]["O_prediction"]).strip()
                    h = str(df_ec.iloc[0]["H_prediction"]).strip()
                    serotype = f"{o}:{h}"

                if not serogroup or serogroup in ["-", "nan", "NA"]:
                    serogroup = "Unknown"
                if not serotype or serotype in ["-", "nan", "NA"]:
                    serotype = "Unknown"

                for d in [linux_output_dir, linux_tmp, output_dir]:
                    if os.path.exists(d): shutil.rmtree(d)
                return serogroup, serotype
                
        for d in [linux_output_dir, linux_tmp, output_dir]:
            if os.path.exists(d): shutil.rmtree(d)
        return "Detection Failed", "Detection Failed"
            
    except subprocess.TimeoutExpired:
        for d in ["linux_output_dir", "linux_tmp", "output_dir"]:
            p = locals().get(d)
            if p and os.path.exists(p): shutil.rmtree(p)
        return "Timeout", "Timeout"
    except Exception as e:
        for d in ["linux_output_dir", "linux_tmp", "output_dir"]:
            p = locals().get(d)
            if p and os.path.exists(p): shutil.rmtree(p)
        return "Error", f"{str(e)[:50]}"

st.set_page_config(page_title="ApexPredict", page_icon="🦠🔬", layout="wide")

HIGH_RISK_STS = ["117", "95", "101", "23", "355", "428", "131"]

APEC_SEROGROUPS = ["O1", "O2", "O6", "O7", "O8", "O15",
                   "O18", "O25", "O36", "O45", "O73", "O86",
                   "O88", "O115", "O117", "O119", "O153", "O161", "O166"]

bg_img_path = "bg.png"
if os.path.exists(bg_img_path):
    img_base64 = get_base64_of_bin_file(bg_img_path)
    bg_style = f'background-image: url("data:image/png;base64,{img_base64}");'
else:
    bg_style = ""

st.markdown(
    f"""
    <style>
    .stApp {{
        {bg_style}
        background-position: bottom right;
        background-repeat: no-repeat;
        background-size: 35%;
        background-color: #f4f6f8;
    }}
    .block-container {{
        padding-top: 3.5rem;
        padding-left: 3rem;
        padding-right: 3rem;
        padding-bottom: 2rem;
        max-width: 1200px;
    }}
    section.main > div {{
        background-color: transparent !important;
    }}
    header[data-testid="stHeader"] {{
        background: transparent !important;
    }}
    .main-title, div.main-title, p.main-title {{
        color: #064949 !important;
        font-size: 52px !important;
        font-weight: 800 !important;
        margin-top: 0.5rem !important;
        margin-bottom: 0 !important;
        line-height: 1.1 !important;
    }}
    .sub-title, div.sub-title, p.sub-title {{
        color: #064949 !important;
        font-size: 25px !important;
        margin-top: 6px !important;
        margin-bottom: 30px !important;
        opacity: 0.85;
    }}
    table {{
        width: 100%;
        border-collapse: collapse;
        color: #064949 !important;
        margin-top: 20px;
    }}
    th {{
        background-color: #064949 !important;
        color: white !important;
        text-align: left;
        padding: 12px;
        font-weight: 700;
    }}
    [data-testid="stTable"] th,
    [data-testid="stTable"] thead tr th,
    .stDataFrame thead tr th {{
        background-color: #064949 !important;
        color: white !important;
        font-weight: 700;
    }}
    td {{
        padding: 10px;
        border-bottom: 1px solid #c1f4da;
        color: #064949 !important;
    }}
    [data-testid="stTable"] table td,
    [data-testid="stTable"] table th,
    .stDataFrame td,
    .stDataFrame th {{
        color: #064949 !important;
    }}
    [data-testid="stMetricValue"] {{
        color: #064949 !important;
        font-size: 28px;
    }}
    [data-testid="stHorizontalBlock"] {{
        background-color: white !important;
        border-radius: 10px;
        padding: 1rem;
        box-shadow: 0 2px 8px rgba(0,0,0,0.07);
    }}
    .streamlit-expanderContent {{
        background-color: rgba(255,255,255,0.96) !important;
        padding: 1rem;
        border-radius: 10px;
    }}
    [data-testid="stTable"] {{
        background-color: white !important;
        border-radius: 10px;
        padding: 0.5rem;
    }}
    .streamlit-expanderHeader {{
        background-color: rgba(255,255,255,0.95) !important;
        border-radius: 8px;
        padding: 0.5rem;
    }}
    </style>
    """,
    unsafe_allow_html=True
)

st.markdown('<div class="main-title">🦠ApexPredict</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-title">Pathotype Classification, MLST, Serogroup Detection, and Resistance Profiling</div>', unsafe_allow_html=True)


@st.cache_resource
def load_analytical_assets():
    model = joblib.load("apec_random_forest_model.joblib")
    features = joblib.load("model_features.joblib")
    return model, features

model, model_features = load_analytical_assets()

ectyper_available = check_ectyper_available()
if not ectyper_available:
    st.sidebar.warning("⚠️ ECTyper not detected. Serogroup detection will be limited.")

uploaded_file = st.file_uploader("Upload Bacterial FASTA", type=["fna", "fasta"])

if uploaded_file is not None:
    temp_path = "temp_sample.fna"
    raw_bytes = uploaded_file.getbuffer().tobytes()
    raw_text = raw_bytes.decode("utf-8", errors="replace").replace("\r\n", "\n").replace("\r", "\n")
    clean_lines = []
    for line in raw_text.splitlines():
        line = line.strip()
        if line:
            clean_lines.append(line)
    with open(temp_path, "w", newline="\n") as f:
        f.write("\n".join(clean_lines) + "\n")
    
    status_indicator = st.empty()
    status_indicator.info("Running Genomic Analysis (MLST, Virulence, AMR, Serogroup) ...")

    mlst_st = "Unknown"
    try:
        mlst_cmd = f"mlst {temp_path}"
        mlst_proc = subprocess.run(mlst_cmd, shell=True, stdout=subprocess.PIPE,
                                   universal_newlines=True, timeout=30)
        if mlst_proc.stdout:
            mlst_data = mlst_proc.stdout.strip().split('\t')
            if len(mlst_data) >= 3:
                mlst_st = mlst_data[2]
    except Exception as e:
        st.warning(f"MLST failed: {e}")

    is_high_risk = mlst_st in HIGH_RISK_STS
    risk_label = "⚠️ HIGH RISK" if is_high_risk else "Standard"

    serogroup = "Unknown"
    serotype = "Unknown"
    serogroup_risk = "Unknown"
    
    if ectyper_available:
        try:
            status_indicator.info("Running Genomic Analysis ...")
            serogroup, serotype = detect_serogroup_with_ectyper(temp_path)
            
            if serogroup in APEC_SEROGROUPS:
                serogroup_risk = "⚠️ APEC-associated"
            elif serogroup not in ["Unknown", "Detection Failed", "Timeout", "Error"]:
                serogroup_risk = "ℹ️ Other serogroup"
            else:
                serogroup_risk = "❓ Undetermined"
        except Exception as e:
            serogroup = "Detection Error"
            serotype = "Detection Error"
            serogroup_risk = "Error"
    else:
        serogroup = "ECTyper not available"
        serotype = "Install ECTyper for detection"
        serogroup_risk = "Not available"

    detailed_results = []
    resistance_results = []
    
    current_dir = os.path.dirname(os.path.abspath(__file__))
    db_directory = os.path.join(current_dir, "apec_db")

    for database in ["ecoli_vf", "apec"]:
        if database == "ecoli_vf":
            cmd = f"abricate --db ecoli_vf {temp_path}"
        else:
            cmd = f"abricate --datadir {db_directory} --db apec {temp_path}"
        
        try:
            proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                                  universal_newlines=True, timeout=60)
            if proc.stdout:
                lines = proc.stdout.strip().split('\n')
                if len(lines) > 1:
                    for entry in lines[1:]:
                        cols = entry.split('\t')
                        if len(cols) >= 10:
                            clean_name = cols[5].replace('~~~', '|').replace('~', '|').split('|')[0].strip().upper()
                            detailed_results.append({
                                "Gene Name": clean_name,
                                "Contig": cols[1],
                                "Identity (%)": cols[9],
                                "Coverage (%)": cols[6]
                            })
        except Exception as e:
            st.warning(f"Virulence screening for {database} failed: {e}")

    try:
        amr_cmd = f"abricate --db ncbi {temp_path}"
        amr_proc = subprocess.run(amr_cmd, shell=True, stdout=subprocess.PIPE,
                                  universal_newlines=True, timeout=60)
        if amr_proc.stdout:
            lines = amr_proc.stdout.strip().split('\n')
            if len(lines) > 1:
                for entry in lines[1:]:
                    cols = entry.split('\t')
                    if len(cols) >= 14:
                        resistance_results.append({
                            "Gene": cols[5],
                            "Class": cols[13],
                            "Identity (%)": cols[9]
                        })
    except Exception as e:
        st.warning(f"AMR screening failed: {e}")

    status_indicator.empty()

    found_genes_list = [r['Gene Name'].strip().upper() for r in detailed_results]
    input_row = {feat: [1 if feat in found_genes_list else 0] for feat in model_features}
    
    prob_score = model.predict_proba(pd.DataFrame(input_row))[0][1]
    prediction = "APEC" if prob_score >= 0.5 else "Non-APEC"
    confidence_pct = prob_score * 100 if prediction == "APEC" else (1 - prob_score) * 100

    st.markdown("---")
    
    m_col1, m_col2, m_col3, m_col4 = st.columns(4)
    
    with m_col1:
        if prediction == "APEC":
            st.error(f"### Result: **{prediction}**")
        else:
            st.success(f"### Result: **{prediction}**")
    
    with m_col2:
        st.metric(label="Lineage (MLST)", value=f"ST{mlst_st}", delta=risk_label,
                  delta_color="inverse" if is_high_risk else "normal")
        
    with m_col3:
        st.metric(label="ML Confidence", value=f"{confidence_pct:.1f}%")
    
    with m_col4:
        if serogroup in APEC_SEROGROUPS:
            st.metric(label="Serogroup", value=serogroup, delta=serogroup_risk,
                      delta_color="inverse")
        else:
            st.metric(label="Serogroup", value=serogroup, delta=serogroup_risk)

    with st.expander(" Detailed Serogroup Information", expanded=False):
        col1, col2 = st.columns(2)
        with col1:
            st.write("**Detected Serogroup:**")
            st.write(f"### {serogroup}")
            st.write(f"**Full Serotype:** {serotype}")
        with col2:
            if serogroup in APEC_SEROGROUPS:
                st.success(f"✅ {serogroup} is a known APEC-associated serogroup")
                st.info("This serogroup is commonly found in avian pathogenic E. coli strains")
            elif serogroup not in ["Unknown", "Detection Failed", "Timeout", "Error",
                                   "ECTyper not available", "Detection Error"]:
                st.warning(f"⚠️ {serogroup} is not typically associated with APEC")
            elif serogroup == "ECTyper not available":
                st.info("💡 Install ECTyper for serogroup detection: https://github.com/ssi-dk/ECTyper")
            else:
                st.info("Serogroup could not be determined from the provided sequence")
        
        if ectyper_available and serogroup not in ["Detection Failed", "Timeout", "Error", "Unknown"]:
            st.caption("Serogroup detected using ECTyper")
        elif not ectyper_available:
            st.caption("⚠️ ECTyper is not installed. Install it for accurate serogroup detection")
    
    if resistance_results:
        st.markdown("### Antimicrobial Resistance Profile")
        res_df = pd.DataFrame(resistance_results)
        unique_classes = res_df["Class"].unique()
        st.warning(f"**Predicted Resistance to:** {', '.join(unique_classes)}")
        with st.expander("Show AMR Genes"):
            st.table(res_df)
    else:
        st.success("No AMR markers detected (Susceptible)")

    if detailed_results:
        with st.expander(f"Show {len(detailed_results)} Virulence Factors", expanded=False):
            vir_df = pd.DataFrame(detailed_results)
            csv_data = vir_df.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="Download Results (CSV)",
                data=csv_data,
                file_name=f"apex_ST{mlst_st}_{serogroup}_{prediction}.csv",
                mime="text/csv",
            )
            st.table(vir_df)

    if os.path.exists(temp_path):
        os.remove(temp_path)
    
    for temp_file in ["temp_ectyper_results.txt", "ectyper_output.json"]:
        if os.path.exists(temp_file):
            os.remove(temp_file)
            
