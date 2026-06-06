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
import requests
from datetime import datetime, timezone

# ─── ENVIRONMENT BRIDGE ───────────────────────────────────────────────────────
BIOTOOLS_BIN = os.environ.get(
    "BIOTOOLS_BIN",
    "/home/khera_aycha/miniconda3/envs/apec_a/bin"
)
os.environ["PATH"] = BIOTOOLS_BIN + os.pathsep + os.environ.get("PATH", "")

# ─── SUPABASE CONFIG ──────────────────────────────────────────────────────────
SUPABASE_URL = "https://vfshhvhqlbfwvmzxewcj.supabase.co"
SUPABASE_ANON_KEY = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InZmc2hodmhxbGJmd3Ztenhld2NqIiwicm9sZSI6ImFub24iLCJpYXQiOjE3ODA2NzAwOTgsImV4cCI6MjA5NjI0NjA5OH0.QHcmyhv5siMN4O-Scy_knT_HdxlvsCzDHUCy8opJumI"

def supa_headers(token=None):
    h = {
        "apikey": SUPABASE_ANON_KEY,
        "Content-Type": "application/json"
    }
    if token:
        h["Authorization"] = f"Bearer {token}"
    else:
        h["Authorization"] = f"Bearer {SUPABASE_ANON_KEY}"
    return h

def auth_signup(email, password, full_name, institution):
    r = requests.post(
        f"{SUPABASE_URL}/auth/v1/signup",
        headers=supa_headers(),
        json={"email": email, "password": password,
              "data": {"full_name": full_name, "institution": institution}}
    )
    return r.json()

def auth_login(email, password):
    r = requests.post(
        f"{SUPABASE_URL}/auth/v1/token?grant_type=password",
        headers=supa_headers(),
        json={"email": email, "password": password}
    )
    return r.json()

def get_profile(token, user_id):
    r = requests.get(
        f"{SUPABASE_URL}/rest/v1/profiles?id=eq.{user_id}&select=*",
        headers=supa_headers(token)
    )
    data = r.json()
    return data[0] if data else None

def get_history(token, user_id):
    r = requests.get(
        f"{SUPABASE_URL}/rest/v1/analyses?user_id=eq.{user_id}&select=*&order=created_at.desc",
        headers=supa_headers(token)
    )
    return r.json()

def save_analysis(token, user_id, filename, result, confidence, mlst_st,
                  serogroup, serotype, amr_classes, virulence_gene_count, is_high_risk_st):
    r = requests.post(
        f"{SUPABASE_URL}/rest/v1/analyses",
        headers={**supa_headers(token), "Prefer": "return=minimal"},
        json={
            "user_id": user_id,
            "filename": filename,
            "result": result,
            "confidence": confidence,
            "mlst_st": mlst_st,
            "serogroup": serogroup,
            "serotype": serotype,
            "amr_classes": amr_classes,
            "virulence_gene_count": virulence_gene_count,
            "is_high_risk_st": is_high_risk_st
        }
    )
    return r.status_code in (200, 201)

FREE_TRIAL_LIMIT = 3

def is_subscription_active(profile):
    """Returns True if user has active subscription OR free trial analyses remaining."""
    if not profile:
        return False
    # Paid subscription check
    if profile.get("subscription_active"):
        end = profile.get("subscription_end")
        if end:
            end_dt = datetime.fromisoformat(end.replace("Z", "+00:00"))
            if end_dt >= datetime.now(timezone.utc):
                return True
        else:
            return True  # No expiry = unlimited
    # Free trial check
    used = profile.get("trial_analyses_used", 0) or 0
    if used < FREE_TRIAL_LIMIT:
        return True
    return False

def is_on_trial(profile):
    if not profile:
        return False
    if profile.get("subscription_active"):
        end = profile.get("subscription_end")
        if end:
            return datetime.fromisoformat(end.replace("Z", "+00:00")) < datetime.now(timezone.utc)
        return False
    used = profile.get("trial_analyses_used", 0) or 0
    return used < FREE_TRIAL_LIMIT

def trial_analyses_remaining(profile):
    used = profile.get("trial_analyses_used", 0) or 0
    return max(0, FREE_TRIAL_LIMIT - used)

def increment_trial_usage(token, user_id, current_used):
    requests.patch(
        f"{SUPABASE_URL}/rest/v1/profiles?id=eq.{user_id}",
        headers={**supa_headers(token), "Prefer": "return=minimal"},
        json={"trial_analyses_used": current_used + 1}
    )

# ─── ORIGINAL HELPER FUNCTIONS ────────────────────────────────────────────────
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

# ─── PAGE CONFIG ──────────────────────────────────────────────────────────────
st.set_page_config(page_title=" ApexPredict", page_icon="assets/logo.png", layout="wide")

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

st.markdown(f"""
    <style>
    /* ── Base ── */
    .stApp {{
        {bg_style}
        background-position: bottom right;
        background-repeat: no-repeat;
        background-size: 30%;
        background-color: #f0f2f5 !important;
    }}

    /* ── Header: make transparent, leave collapse arrow alone ── */
    header[data-testid="stHeader"] {{
        background: rgba(240,242,245,0) !important;
        border-bottom: none !important;
        box-shadow: none !important;
    }}

    /* ── Main content: centered with reasonable max-width ── */
    .main .block-container {{
        max-width: 860px !important;
        padding-top: 3rem !important;
        padding-left: 2rem !important;
        padding-right: 2rem !important;
        padding-bottom: 3rem !important;
        margin: 0 auto !important;
    }}
    section.main > div {{ background-color: transparent !important; }}

    /* ── Sidebar ── */
    [data-testid="stSidebar"] {{
        background: #064949 !important;
        border-right: none !important;
    }}
    [data-testid="stSidebar"] > div:first-child {{
        padding-top: 0 !important;
    }}
    [data-testid="stSidebarContent"] {{
        padding: 0 !important;
    }}

    /* ── Sidebar logo area ── */
    .sidebar-logo-block {{
        background: #043333;
        padding: 1.6rem 1.4rem 1.2rem 1.4rem;
        border-bottom: 1px solid rgba(255,255,255,0.08);
        margin-bottom: 0.5rem;
    }}
    .sidebar-logo-block img {{
        height: 38px;
        width: auto;
        display: block;
        margin-bottom: 0.5rem;
    }}
    .sidebar-logo-block .app-name {{
        color: #ffffff !important;
        font-size: 18px;
        font-weight: 700;
        letter-spacing: -0.3px;
        display: block;
    }}
    .sidebar-logo-block .app-tagline {{
        color: rgba(255,255,255,0.5) !important;
        font-size: 11px;
        font-weight: 400;
        display: block;
        margin-top: 2px;
        line-height: 1.4;
    }}

    /* ── Sidebar profile block ── */
    .sidebar-profile {{
        padding: 1rem 1.4rem;
        border-bottom: 1px solid rgba(255,255,255,0.08);
    }}
    .sidebar-profile .profile-name {{
        color: #ffffff !important;
        font-size: 14px;
        font-weight: 600;
        display: block;
        margin-bottom: 2px;
    }}
    .sidebar-profile .profile-institution {{
        color: rgba(255,255,255,0.55) !important;
        font-size: 12px;
        display: block;
        margin-bottom: 2px;
    }}
    .sidebar-profile .profile-email {{
        color: rgba(255,255,255,0.4) !important;
        font-size: 11px;
        display: block;
    }}

    /* ── Sidebar subscription badge ── */
    .sidebar-badge {{
        margin: 0.8rem 1.4rem 0.4rem 1.4rem;
        padding: 0.5rem 0.8rem;
        border-radius: 6px;
        font-size: 12px;
        font-weight: 600;
    }}
    .sidebar-badge.active {{
        background: rgba(22,163,74,0.2);
        color: #4ade80 !important;
        border: 1px solid rgba(74,222,128,0.25);
    }}
    .sidebar-badge.trial {{
        background: rgba(217,119,6,0.2);
        color: #fbbf24 !important;
        border: 1px solid rgba(251,191,36,0.25);
    }}
    .sidebar-badge.inactive {{
        background: rgba(220,38,38,0.2);
        color: #f87171 !important;
        border: 1px solid rgba(248,113,113,0.25);
    }}

    /* ── Sidebar nav buttons ── */
    [data-testid="stSidebar"] .stButton > button {{
        background: transparent !important;
        border: none !important;
        color: rgba(255,255,255,0.75) !important;
        font-size: 13px !important;
        font-weight: 500 !important;
        text-align: left !important;
        padding: 0.55rem 1.4rem !important;
        border-radius: 0 !important;
        width: 100% !important;
        justify-content: flex-start !important;
        transition: all 0.15s ease !important;
        box-shadow: none !important;
    }}
    [data-testid="stSidebar"] .stButton > button:hover {{
        background: rgba(255,255,255,0.08) !important;
        color: #ffffff !important;
    }}
    [data-testid="stSidebar"] .stButton > button[kind="primary"] {{
        background: rgba(255,255,255,0.12) !important;
        color: #ffffff !important;
        border-left: 3px solid #4ade80 !important;
    }}

    /* ── Sidebar alert ── */
    [data-testid="stSidebar"] [data-testid="stAlert"] {{
        background: rgba(255,255,255,0.06) !important;
        border: none !important;
        border-left: 3px solid rgba(251,191,36,0.6) !important;
        color: rgba(255,255,255,0.7) !important;
        font-size: 12px !important;
        margin: 0.5rem 1rem !important;
        border-radius: 4px !important;
    }}

    /* ── Typography ── */
    .main-title, div.main-title {{
        color: #064949 !important;
        font-size: 34px !important;
        font-weight: 700 !important;
        margin-top: 0 !important;
        margin-bottom: 0 !important;
        line-height: 1.15 !important;
        letter-spacing: -0.5px;
    }}
    .sub-title, div.sub-title {{
        color: #4a6060 !important;
        font-size: 14px !important;
        margin-top: 4px !important;
        margin-bottom: 20px !important;
        font-weight: 400;
    }}
    h2, h3 {{ color: #064949 !important; font-weight: 600 !important; }}

    /* ── Section divider ── */
    hr {{
        border: none !important;
        border-top: 1px solid #e2e8f0 !important;
        margin: 1.5rem 0 !important;
    }}

    /* ── Tables ── */
    table {{ width: 100%; border-collapse: collapse; color: #1e3a3a !important; margin-top: 12px; font-size: 13px; }}
    th {{
        background-color: #064949 !important;
        color: white !important;
        text-align: left;
        padding: 10px 14px;
        font-weight: 600;
        font-size: 11px;
        text-transform: uppercase;
        letter-spacing: 0.8px;
    }}
    [data-testid="stTable"] th, [data-testid="stTable"] thead tr th, .stDataFrame thead tr th {{
        background-color: #064949 !important;
        color: white !important;
        font-weight: 600;
    }}
    td {{ padding: 9px 14px; border-bottom: 1px solid #e8f5ee; color: #1e3a3a !important; }}
    [data-testid="stTable"] table td, [data-testid="stTable"] table th,
    .stDataFrame td, .stDataFrame th {{ color: #1e3a3a !important; }}
    [data-testid="stTable"] {{
        background-color: white !important;
        border-radius: 10px;
        padding: 0;
        overflow: hidden;
        box-shadow: 0 1px 6px rgba(0,0,0,0.05);
        border: 1px solid #e8eef0;
    }}

    /* ── Metrics ── */
    [data-testid="stMetricValue"] {{ color: #064949 !important; font-size: 24px; font-weight: 700; }}
    [data-testid="stMetricLabel"] {{ color: #64748b !important; font-size: 11px; text-transform: uppercase; letter-spacing: 0.6px; }}
    [data-testid="stHorizontalBlock"]:has([data-testid="stMetricValue"]) {{
        background-color: white !important;
        border-radius: 10px;
        padding: 1.2rem 1rem;
        box-shadow: 0 1px 6px rgba(6,73,73,0.06);
        border: 1px solid #e8eef0;
    }}

    /* ── Expanders ── */
    .streamlit-expanderContent {{
        background-color: white !important;
        padding: 1rem;
        border-radius: 0 0 8px 8px;
        border: 1px solid #e2e8f0;
        border-top: none;
    }}
    .streamlit-expanderHeader {{
        background-color: #f8fffe !important;
        border-radius: 8px;
        padding: 0.6rem 1rem;
        border: 1px solid #e2e8f0;
        font-weight: 500;
        color: #064949 !important;
        font-size: 14px !important;
    }}

    /* ── Auth form wrapper ── */
    .auth-form {{
        background: white;
        border-radius: 12px;
        padding: 2.5rem 2rem;
        box-shadow: 0 2px 16px rgba(6,73,73,0.08);
        border: 1px solid #e8eef0;
        margin-top: 1rem;
    }}

    /* ── Profile cards ── */
    .profile-card {{
        background: white;
        border-radius: 10px;
        padding: 1.4rem;
        box-shadow: 0 1px 6px rgba(6,73,73,0.06);
        border: 1px solid #e8eef0;
        margin-bottom: 1rem;
    }}
    .profile-card-row {{
        display: flex;
        justify-content: space-between;
        align-items: baseline;
        padding: 0.45rem 0;
        border-bottom: 1px solid #f1f5f5;
        gap: 1rem;
    }}
    .profile-card-row:last-child {{ border-bottom: none; }}
    .profile-card-label {{
        font-size: 12px;
        color: #64748b;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.4px;
        white-space: nowrap;
    }}
    .profile-card-value {{
        font-size: 14px;
        color: #1e3a3a;
        font-weight: 500;
        text-align: right;
    }}

    /* ── Subscription badges ── */
    .subscription-active {{
        background: #f0fdf4;
        border-left: 3px solid #16a34a;
        padding: 0.65rem 1rem;
        border-radius: 6px;
        color: #15803d;
        font-weight: 600;
        font-size: 13px;
        margin-bottom: 0.75rem;
    }}
    .subscription-inactive {{
        background: #fff1f2;
        border-left: 3px solid #dc2626;
        padding: 0.65rem 1rem;
        border-radius: 6px;
        color: #dc2626;
        font-weight: 600;
        font-size: 13px;
        margin-bottom: 0.75rem;
    }}
    .trial-badge {{
        background: #fffbeb;
        border-left: 3px solid #d97706;
        padding: 0.65rem 1rem;
        border-radius: 6px;
        color: #92400e;
        font-weight: 600;
        font-size: 13px;
        margin-bottom: 0.75rem;
    }}

    /* ── Primary buttons ── */
    .stButton > button[kind="primary"] {{
        background-color: #064949 !important;
        border: none !important;
        color: white !important;
        border-radius: 7px !important;
        font-weight: 600 !important;
        font-size: 13px !important;
        transition: background 0.15s ease !important;
    }}
    .stButton > button[kind="primary"]:hover {{
        background-color: #0a5f5f !important;
    }}

    /* ── File uploader ── */
    [data-testid="stFileUploader"] {{
        background: white !important;
        border: 1.5px dashed #c5d8d8 !important;
        border-radius: 10px !important;
        padding: 0.5rem !important;
    }}

    /* ── Alert/info boxes ── */
    [data-testid="stAlert"] {{
        border-radius: 8px !important;
        font-size: 13px !important;
    }}
    </style>
""", unsafe_allow_html=True)

# ─── SESSION STATE ────────────────────────────────────────────────────────────
if "user" not in st.session_state:
    st.session_state.user = None
if "token" not in st.session_state:
    st.session_state.token = None
if "profile" not in st.session_state:
    st.session_state.profile = None
if "page" not in st.session_state:
    st.session_state.page = "login"

# ─── AUTH PAGES ───────────────────────────────────────────────────────────────

def inject_sidebar(show_nav=True):
    """Render the branded sidebar with logo, profile info, and nav."""
    with st.sidebar:
        # ── Logo block ──
        logo_path = "logo.png"
        if os.path.exists(logo_path):
            logo_b64 = get_base64_of_bin_file(logo_path)
            logo_img = f'<img src="data:image/png;base64,{logo_b64}" />'
        else:
            logo_img = ""
        st.markdown(
            f'''<div class="sidebar-logo-block">
                {logo_img}
                <span class="app-name">ApexPredict</span>
                <span class="app-tagline">E. coli Pathotype Intelligence Platform</span>
            </div>''',
            unsafe_allow_html=True
        )

        if show_nav and st.session_state.get("user"):
            profile = st.session_state.get("profile", {}) or {}
            user = st.session_state.get("user", {}) or {}

            # ── Profile info ──
            full_name = profile.get("full_name", "")
            institution = profile.get("institution", "")
            email = user.get("email", "")
            st.markdown(
                f'''<div class="sidebar-profile">
                    <span class="profile-name">{full_name}</span>
                    <span class="profile-institution">{institution}</span>
                    <span class="profile-email">{email}</span>
                </div>''',
                unsafe_allow_html=True
            )

            # ── Subscription badge ──
            if profile.get("subscription_active"):
                end = profile.get("subscription_end", "")
                if end:
                    end_dt = datetime.fromisoformat(end.replace("Z", "+00:00"))
                    if end_dt >= datetime.now(timezone.utc):
                        badge_cls, badge_txt = "active", f"Active — expires {end[:10]}"
                    else:
                        badge_cls, badge_txt = "inactive", "Subscription Expired"
                else:
                    badge_cls, badge_txt = "active", "Active — No expiry"
            elif is_on_trial(profile):
                rem = trial_analyses_remaining(profile)
                badge_cls, badge_txt = "trial", f"Free Trial — {rem}/{FREE_TRIAL_LIMIT} remaining"
            else:
                badge_cls, badge_txt = "inactive", "Trial Ended — No Subscription"
            st.markdown(
                f'<div class="sidebar-badge {badge_cls}">{badge_txt}</div>',
                unsafe_allow_html=True
            )

            st.markdown("<br>", unsafe_allow_html=True)

            # ── Nav buttons ──
            current_page = st.session_state.get("page", "app")
            if st.button("Analysis", use_container_width=True,
                         type="primary" if current_page == "app" else "secondary"):
                st.session_state.page = "app"
                st.rerun()
            if st.button("My Profile", use_container_width=True,
                         type="primary" if current_page == "profile" else "secondary"):
                st.session_state.page = "profile"
                st.rerun()

            st.markdown("<br>" * 8, unsafe_allow_html=True)
            if st.button("Sign Out", use_container_width=True):
                st.session_state.user = None
                st.session_state.token = None
                st.session_state.profile = None
                st.session_state.page = "login"
                st.rerun()

def inject_logo():
    """Legacy shim — now handled by inject_sidebar."""
    pass

def show_login():
    inject_sidebar(show_nav=False)
    st.markdown('<div class="main-title">🦠 ApexPredict</div>', unsafe_allow_html=True)
    st.markdown('<div class="sub-title">Pathotype Classification · MLST · Serogroup Detection · Resistance Profiling</div>', unsafe_allow_html=True)
    st.markdown("---")
    _, center, _ = st.columns([1, 2, 1])
    with center:
        st.markdown("#### Sign In")
        email = st.text_input("Email", key="login_email")
        password = st.text_input("Password", type="password", key="login_password")
        col1, col2 = st.columns(2)
        with col1:
            if st.button("Sign In", use_container_width=True, type="primary"):
                if email and password:
                    with st.spinner("Signing in..."):
                        resp = auth_login(email, password)
                    if "access_token" in resp:
                        st.session_state.token = resp["access_token"]
                        st.session_state.user = resp["user"]
                        st.session_state.profile = get_profile(resp["access_token"], resp["user"]["id"])
                        st.session_state.page = "app"
                        st.rerun()
                    else:
                        msg = resp.get("error_description", resp.get("msg", "Invalid credentials"))
                        st.error(f"{msg}")
                else:
                    st.warning("Please enter your email and password.")
        with col2:
            if st.button("Create Account", use_container_width=True):
                st.session_state.page = "signup"
                st.rerun()

def show_signup():
    inject_sidebar(show_nav=False)
    st.markdown('<div class="main-title">🦠 ApexPredict</div>', unsafe_allow_html=True)
    st.markdown('<div class="sub-title">Create your account</div>', unsafe_allow_html=True)
    st.markdown("---")
    _, center, _ = st.columns([1, 2, 1])
    with center:
        st.markdown("#### Sign Up")
        st.info("New accounts include **3 free trial analyses** — no subscription needed to get started.")
        full_name = st.text_input("Full Name", key="signup_name")
        institution = st.text_input("Institution / Organization", key="signup_inst")
        email = st.text_input("Email", key="signup_email")
        password = st.text_input("Password (min. 6 characters)", type="password", key="signup_pass")
        password2 = st.text_input("Confirm Password", type="password", key="signup_pass2")
        col1, col2 = st.columns(2)
        with col1:
            if st.button("Create Account", use_container_width=True, type="primary"):
                if not all([full_name, email, password, password2]):
                    st.warning("Please fill in all fields.")
                elif password != password2:
                    st.error("Passwords do not match.")
                elif len(password) < 6:
                    st.error("Password must be at least 6 characters.")
                else:
                    with st.spinner("Creating account..."):
                        resp = auth_signup(email, password, full_name, institution)
                    if "id" in resp or ("user" in resp and resp["user"]):
                        st.success("Account created! Please check your email to confirm your account, then sign in.")
                        st.session_state.page = "login"
                        st.rerun()
                    else:
                        msg = resp.get("error_description", resp.get("msg", "Signup failed."))
                        st.error(f"{msg}")
        with col2:
            if st.button("Back to Sign In", use_container_width=True):
                st.session_state.page = "login"
                st.rerun()

def show_profile():
    inject_sidebar(show_nav=True)
    profile = st.session_state.profile
    token = st.session_state.token
    user = st.session_state.user

    st.markdown('<div class="main-title">My Profile</div>', unsafe_allow_html=True)
    st.markdown("---")

    # Profile info — pure HTML so card padding/layout renders correctly
    if profile.get("subscription_active"):
        end = profile.get("subscription_end", "")
        if end:
            end_dt = datetime.fromisoformat(end.replace("Z", "+00:00"))
            if end_dt >= datetime.now(timezone.utc):
                sub_html = f'<div class="subscription-active">Subscription Active<br><small>Expires: {end[:10]}</small></div>'
            else:
                sub_html = '<div class="subscription-inactive">Subscription Expired<br><small>Contact us: apexpredict@contact.com</small></div>'
        else:
            sub_html = '<div class="subscription-active">Subscription Active<br><small>No expiry</small></div>'
    elif is_on_trial(profile):
        remaining = trial_analyses_remaining(profile)
        sub_html = f'<div class="trial-badge">Free Trial — {remaining} of {FREE_TRIAL_LIMIT} remaining<br><small>To subscribe: apexpredict@contact.com</small></div>'
    else:
        sub_html = '<div class="subscription-inactive">Trial Ended<br><small>Contact us: apexpredict@contact.com</small></div>'

    c1, c2 = st.columns(2)
    with c1:
        st.markdown(f'''
        <div class="profile-card">
            <div class="profile-card-row"><span class="profile-card-label">Name</span><span class="profile-card-value">{profile.get("full_name", "N/A")}</span></div>
            <div class="profile-card-row"><span class="profile-card-label">Email</span><span class="profile-card-value">{user.get("email", "N/A")}</span></div>
            <div class="profile-card-row"><span class="profile-card-label">Institution</span><span class="profile-card-value">{profile.get("institution", "N/A")}</span></div>
            <div class="profile-card-row"><span class="profile-card-label">Member since</span><span class="profile-card-value">{profile.get("created_at", "")[:10]}</span></div>
        </div>
        ''', unsafe_allow_html=True)
    with c2:
        st.markdown(f'''
        <div class="profile-card">
            <div style="font-size:12px;text-transform:uppercase;letter-spacing:0.5px;color:#64748b;font-weight:600;margin-bottom:0.75rem;">Subscription</div>
            {sub_html}
        </div>
        ''', unsafe_allow_html=True)

    st.markdown("---")
    st.markdown("## Analysis History")

    history = get_history(token, user["id"])
    if not history:
        st.info("No analyses yet. Run your first analysis to see results here.")
    else:
        st.caption(f"{len(history)} analyses found")
        for rec in history:
            date = rec.get("created_at", "")[:10]
            result = rec.get("result", "N/A")
            color = "#fee2e2" if result == "APEC" else "#dcfce7"
            border = "#dc2626" if result == "APEC" else "#16a34a"
            with st.expander(f"{'[APEC]' if result == 'APEC' else '[Non-APEC]'} {rec.get('filename', 'Unknown')} — {date} — {result}"):
                hc1, hc2, hc3, hc4 = st.columns(4)
                hc1.metric("Result", result)
                hc2.metric("Confidence", f"{rec.get('confidence', 0):.1f}%")
                hc3.metric("ST", f"ST{rec.get('mlst_st', 'N/A')}")
                hc4.metric("Serogroup", rec.get('serogroup', 'N/A'))
                st.caption(f"Serotype: {rec.get('serotype', 'N/A')} | AMR: {rec.get('amr_classes', 'None')} | Virulence genes: {rec.get('virulence_gene_count', 0)}")

def show_subscription_wall():
    st.markdown('<div class="main-title">Subscription Required</div>', unsafe_allow_html=True)
    st.markdown("---")
    st.markdown("""
    <div style='text-align:center; padding: 2rem 1rem;'>
        <h2 style='color:#064949;'>Access Restricted</h2>
        <p style='font-size:16px; color:#4a5568;'>
            Your free trial has ended or your subscription is not active.<br>
            You can still view your analysis history in your profile.
        </p>
        <p style='font-size:14px; color:#064949; margin-top:1.5rem;'>
            <strong>To subscribe, contact us at:<br>
            apexpredict@contact.com</strong>
        </p>
    </div>
    """, unsafe_allow_html=True)
    col1, col2, col3 = st.columns([2, 1, 2])
    with col2:
        if st.button("View My Profile", use_container_width=True, type="primary"):
            st.session_state.page = "profile"
            st.rerun()

# ─── MAIN APP ────────────────────────────────────────────────────────────────
def show_app():
    inject_sidebar(show_nav=True)
    profile = st.session_state.profile
    token = st.session_state.token
    user = st.session_state.user

    # Check subscription / trial
    if not is_subscription_active(profile):
        show_subscription_wall()
        return

    # Show trial banner if user is on trial
    if is_on_trial(profile):
        remaining = trial_analyses_remaining(profile)
        st.markdown(
            f'''<div class="trial-badge">Free Trial — {remaining} analysis remaining. Contact us to subscribe for unlimited access.</div>''',
            unsafe_allow_html=True
        )

    # Header
    st.markdown('<div class="main-title">🦠 ApexPredict</div>', unsafe_allow_html=True)
    st.markdown('<div class="sub-title">Pathotype Classification · MLST · Serogroup Detection · Resistance Profiling</div>', unsafe_allow_html=True)

    @st.cache_resource
    def load_analytical_assets():
        model = joblib.load("apec_random_forest_model.joblib")
        features = joblib.load("model_features.joblib")
        return model, features

    model, model_features = load_analytical_assets()

    ectyper_available = check_ectyper_available()
    if not ectyper_available:
        st.sidebar.warning("ECTyper not detected. Serogroup detection will be limited.")

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
        risk_label = "HIGH RISK" if is_high_risk else "Standard"

        serogroup = "Unknown"
        serotype = "Unknown"
        serogroup_risk = "Unknown"

        if ectyper_available:
            try:
                status_indicator.info("Running Genomic Analysis ...")
                serogroup, serotype = detect_serogroup_with_ectyper(temp_path)
                if serogroup in APEC_SEROGROUPS:
                    serogroup_risk = "APEC-associated"
                elif serogroup not in ["Unknown", "Detection Failed", "Timeout", "Error"]:
                    serogroup_risk = "Other serogroup"
                else:
                    serogroup_risk = "Undetermined"
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

        # ── Save to history ──────────────────────────────────────────────────
        amr_classes_str = ", ".join(set([r["Class"] for r in resistance_results])) if resistance_results else "None"
        save_analysis(
            token=token,
            user_id=user["id"],
            filename=uploaded_file.name,
            result=prediction,
            confidence=round(confidence_pct, 1),
            mlst_st=mlst_st,
            serogroup=serogroup,
            serotype=serotype,
            amr_classes=amr_classes_str,
            virulence_gene_count=len(detailed_results),
            is_high_risk_st=is_high_risk
        )
        # If on free trial, increment usage counter
        if is_on_trial(profile):
            used = profile.get("trial_analyses_used", 0) or 0
            increment_trial_usage(token, user["id"], used)
            st.session_state.profile = get_profile(token, user["id"])

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

        with st.expander("Detailed Serogroup Information", expanded=False):
            col1, col2 = st.columns(2)
            with col1:
                st.write("**Detected Serogroup:**")
                st.write(f"### {serogroup}")
                st.write(f"**Full Serotype:** {serotype}")
            with col2:
                if serogroup in APEC_SEROGROUPS:
                    st.success(f"{serogroup} is a known APEC-associated serogroup")
                    st.info("This serogroup is commonly found in avian pathogenic E. coli strains")
                elif serogroup not in ["Unknown", "Detection Failed", "Timeout", "Error",
                                       "ECTyper not available", "Detection Error"]:
                    st.warning(f"{serogroup} is not typically associated with APEC")
                elif serogroup == "ECTyper not available":
                    st.info("Install ECTyper for serogroup detection: https://github.com/ssi-dk/ECTyper")
                else:
                    st.info("Serogroup could not be determined from the provided sequence")
            if ectyper_available and serogroup not in ["Detection Failed", "Timeout", "Error", "Unknown"]:
                st.caption("Serogroup detected using ECTyper")
            elif not ectyper_available:
                st.caption("ECTyper is not installed. Install it for accurate serogroup detection")

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

# ─── ROUTER ───────────────────────────────────────────────────────────────────
if st.session_state.page == "login":
    show_login()
elif st.session_state.page == "signup":
    show_signup()
elif st.session_state.page == "profile":
    show_profile()
elif st.session_state.page == "app":
    show_app()
