import streamlit as st
import pandas as pd
import joblib
import subprocess
import os
import base64
import shutil
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
FREE_TRIAL_LIMIT = 3

def supa_headers(token=None):
    h = {"apikey": SUPABASE_ANON_KEY, "Content-Type": "application/json"}
    h["Authorization"] = f"Bearer {token}" if token else f"Bearer {SUPABASE_ANON_KEY}"
    return h

def auth_signup(email, password, full_name, institution):
    r = requests.post(f"{SUPABASE_URL}/auth/v1/signup", headers=supa_headers(),
        json={"email": email, "password": password,
              "data": {"full_name": full_name, "institution": institution}})
    return r.json()

def auth_login(email, password):
    r = requests.post(f"{SUPABASE_URL}/auth/v1/token?grant_type=password",
        headers=supa_headers(), json={"email": email, "password": password})
    return r.json()

def get_profile(token, user_id):
    r = requests.get(f"{SUPABASE_URL}/rest/v1/profiles?id=eq.{user_id}&select=*",
        headers=supa_headers(token))
    data = r.json()
    return data[0] if isinstance(data, list) and data else None

def get_history(token, user_id):
    r = requests.get(
        f"{SUPABASE_URL}/rest/v1/analyses?user_id=eq.{user_id}&select=*&order=created_at.desc",
        headers=supa_headers(token))
    return r.json() if isinstance(r.json(), list) else []

def get_all_users(token):
    r = requests.get(f"{SUPABASE_URL}/rest/v1/profiles?select=*&order=created_at.desc",
        headers=supa_headers(token))
    return r.json() if r.status_code == 200 and isinstance(r.json(), list) else []

def get_all_analyses(token):
    r = requests.get(
        f"{SUPABASE_URL}/rest/v1/analyses?select=*,profiles(email,full_name)&order=created_at.desc",
        headers=supa_headers(token))
    return r.json() if r.status_code == 200 and isinstance(r.json(), list) else []

def save_analysis(token, user_id, filename, result, confidence, mlst_st,
                  serogroup, serotype, amr_classes, virulence_gene_count, is_high_risk_st):
    r = requests.post(f"{SUPABASE_URL}/rest/v1/analyses",
        headers={**supa_headers(token), "Prefer": "return=minimal"},
        json={"user_id": user_id, "filename": filename, "result": result,
              "confidence": confidence, "mlst_st": mlst_st, "serogroup": serogroup,
              "serotype": serotype, "amr_classes": amr_classes,
              "virulence_gene_count": virulence_gene_count, "is_high_risk_st": is_high_risk_st})
    return r.status_code in (200, 201)

def is_subscription_active(profile):
    if not profile: return False
    if profile.get("subscription_active"):
        end = profile.get("subscription_end")
        if end:
            if datetime.fromisoformat(end.replace("Z", "+00:00")) >= datetime.now(timezone.utc):
                return True
        else:
            return True
    used = profile.get("trial_analyses_used", 0) or 0
    return used < FREE_TRIAL_LIMIT

def is_on_trial(profile):
    if not profile: return False
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
    requests.patch(f"{SUPABASE_URL}/rest/v1/profiles?id=eq.{user_id}",
        headers={**supa_headers(token), "Prefer": "return=minimal"},
        json={"trial_analyses_used": current_used + 1})

# ─── BIOINFORMATICS HELPERS ───────────────────────────────────────────────────
def get_base64_of_bin_file(bin_file):
    with open(bin_file, 'rb') as f:
        return base64.b64encode(f.read()).decode()

def check_ectyper_available():
    try:
        result = subprocess.run("ectyper --help", shell=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=5)
        return result.returncode in (0, 1)
    except:
        return False

def detect_serogroup_with_ectyper(fasta_path):
    try:
        output_dir = "temp_ectyper_out"
        if os.path.exists(output_dir): shutil.rmtree(output_dir)
        linux_tmp = tempfile.mkdtemp(prefix="/tmp/ectyper_in_")
        linux_fasta = os.path.join(linux_tmp, "sample.fna")
        shutil.copy2(fasta_path, linux_fasta)
        linux_output_dir = tempfile.mkdtemp(prefix="/tmp/ectyper_out_")
        shutil.rmtree(linux_output_dir)
        cmd = f"ectyper --input {linux_fasta} --output {linux_output_dir}"
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, universal_newlines=True, timeout=40)
        output_path = None
        for candidate in [os.path.join(linux_output_dir, "output.csv"),
                          os.path.join(linux_output_dir, "output.tsv")]:
            if os.path.exists(candidate) and os.path.getsize(candidate) > 0:
                output_path = candidate
                break
        if output_path is None:
            sp = (glob.glob(os.path.join(linux_output_dir, "**", "*.csv"), recursive=True) or
                  glob.glob(os.path.join(linux_output_dir, "**", "*.tsv"), recursive=True))
            if sp: output_path = sp[0]
        if output_path and os.path.getsize(output_path) > 0:
            sep = '\t' if output_path.endswith('.tsv') else ','
            df_ec = pd.read_csv(output_path, sep=sep)
            if not df_ec.empty:
                serogroup, serotype = "Unknown", "Unknown"
                for col in ["O_prediction", "O.type", "O_type", "Serogroup"]:
                    if col in df_ec.columns:
                        serogroup = str(df_ec.iloc[0][col]).strip()
                        break
                if "Serotype" in df_ec.columns:
                    serotype = str(df_ec.iloc[0]["Serotype"]).strip()
                elif "O_prediction" in df_ec.columns and "H_prediction" in df_ec.columns:
                    serotype = f"{str(df_ec.iloc[0]['O_prediction']).strip()}:{str(df_ec.iloc[0]['H_prediction']).strip()}"
                if not serogroup or serogroup in ["-", "nan", "NA"]: serogroup = "Unknown"
                if not serotype or serotype in ["-", "nan", "NA"]: serotype = "Unknown"
                for d in [linux_output_dir, linux_tmp, output_dir]:
                    if os.path.exists(d): shutil.rmtree(d)
                return serogroup, serotype
        for d in [linux_output_dir, linux_tmp, output_dir]:
            if os.path.exists(d): shutil.rmtree(d)
        return "Detection Failed", "Detection Failed"
    except subprocess.TimeoutExpired:
        return "Timeout", "Timeout"
    except Exception as e:
        return "Error", f"{str(e)[:50]}"

# ─── PAGE CONFIG ──────────────────────────────────────────────────────────────
st.set_page_config(page_title=" ApexPredict", page_icon="🦠", layout="wide",
                   initial_sidebar_state="expanded")

HIGH_RISK_STS = ["117", "95", "101", "23", "355", "428", "131"]
APEC_SEROGROUPS = ["O1", "O2", "O6", "O7", "O8", "O15", "O18", "O25", "O36",
                   "O45", "O73", "O86", "O88", "O115", "O117", "O119", "O153", "O161", "O166"]

bg_img_path = "bg.png"
bg_style = ""
if os.path.exists(bg_img_path):
    img_base64 = get_base64_of_bin_file(bg_img_path)
    bg_style = f'background-image: url("data:image/png;base64,{img_base64}");'

st.markdown(f"""
<style>
/* ── Main app ── */
.stApp {{
    {bg_style}
    background-position: bottom right;
    background-repeat: no-repeat;
    background-size: 25%;
    background-color: #f4f6f8;
    background-attachment: fixed;
}}
.block-container {{
    padding-top: 4.5rem;
    padding-left: 2.5rem;
    padding-right: 2.5rem;
    padding-bottom: 2rem;
    max-width: 1200px;
}}
section.main > div {{ background-color: transparent !important; }}
    /* Expand main content when sidebar is collapsed */
    [data-testid="collapsedControl"] {{
        display: block !important;
    }}
    section[data-testid="stSidebar"][aria-expanded="false"] {{
        min-width: 0 !important;
        max-width: 0 !important;
        overflow: hidden !important;
    }}
    section[data-testid="stSidebar"][aria-expanded="false"] + section {{
        margin-left: 0 !important;
        width: 100% !important;
    }}

/* ── Sidebar ── */
[data-testid="stSidebar"] {{
    background: linear-gradient(180deg, #064949 0%, #043333 100%) !important;
    min-width: 260px !important;
    max-width: 260px !important;
}}
[data-testid="stSidebar"] * {{ color: white !important; }}
[data-testid="stSidebarContent"] {{ padding: 0 !important; }}
.sidebar-logo {{
    padding: 1.5rem 1.2rem 1rem 1.2rem;
    border-bottom: 1px solid rgba(255,255,255,0.12);
    margin-bottom: 0.5rem;
}}
.sidebar-logo img {{
    height: 52px;
    width: auto;
    margin-bottom: 0.6rem;
}}
.sidebar-appname {{
    font-size: 20px;
    font-weight: 800;
    color: white !important;
    letter-spacing: -0.5px;
    line-height: 1.2;
}}
.sidebar-tagline {{
    font-size: 11px;
    color: rgba(193,244,218,0.8) !important;
    margin-top: 2px;
}}
.sidebar-user {{
    padding: 0.9rem 1.2rem;
    border-bottom: 1px solid rgba(255,255,255,0.10);
    margin-bottom: 0.5rem;
}}
.sidebar-user-label {{
    font-size: 10px;
    text-transform: uppercase;
    letter-spacing: 1px;
    color: rgba(193,244,218,0.6) !important;
    margin-bottom: 3px;
}}
.sidebar-user-name {{
    font-size: 13px;
    font-weight: 600;
    color: white !important;
}}
.sidebar-user-inst {{
    font-size: 11px;
    color: rgba(193,244,218,0.75) !important;
}}
.sidebar-user-email {{
    font-size: 11px;
    color: rgba(193,244,218,0.6) !important;
    margin-top: 1px;
}}
.sidebar-badge-active {{
    background: rgba(22,163,74,0.25);
    border: 1px solid rgba(22,163,74,0.5);
    border-radius: 20px;
    padding: 3px 10px;
    font-size: 11px;
    color: #86efac !important;
    display: inline-block;
    margin-top: 6px;
}}
.sidebar-badge-trial {{
    background: rgba(217,119,6,0.25);
    border: 1px solid rgba(217,119,6,0.5);
    border-radius: 20px;
    padding: 3px 10px;
    font-size: 11px;
    color: #fcd34d !important;
    display: inline-block;
    margin-top: 6px;
}}
.sidebar-badge-inactive {{
    background: rgba(220,38,38,0.25);
    border: 1px solid rgba(220,38,38,0.5);
    border-radius: 20px;
    padding: 3px 10px;
    font-size: 11px;
    color: #fca5a5 !important;
    display: inline-block;
    margin-top: 6px;
}}
/* Sidebar nav buttons */
[data-testid="stSidebar"] .stButton > button {{
    background: transparent !important;
    border: none !important;
    color: rgba(193,244,218,0.85) !important;
    text-align: left !important;
    padding: 0.55rem 1.2rem !important;
    width: 100% !important;
    font-size: 14px !important;
    border-radius: 0 !important;
    font-weight: 400 !important;
    transition: background 0.15s;
}}
[data-testid="stSidebar"] .stButton > button:hover {{
    background: rgba(255,255,255,0.08) !important;
    color: white !important;
}}
.sidebar-nav-active > button {{
    background: rgba(255,255,255,0.12) !important;
    color: white !important;
    font-weight: 700 !important;
    border-left: 3px solid #c1f4da !important;
}}
.sidebar-signout > button {{
    color: rgba(252,165,165,0.8) !important;
    margin-top: 0.5rem;
}}

/* ── Main titles ── */
.main-title {{
    color: #064949 !important;
    font-size: 42px !important;
    font-weight: 800 !important;
    margin-bottom: 0 !important;
    line-height: 1.1 !important;
    letter-spacing: -1px;
}}
.sub-title {{
    color: #064949 !important;
    font-size: 16px !important;
    margin-top: 4px !important;
    margin-bottom: 20px !important;
    opacity: 0.65;
}}

/* ── Auth pages ── */
.auth-header {{
    text-align: center;
    padding: 0.2rem 0 0.6rem 0;
}}
.auth-header .app-name {{
    font-size: 36px;
    font-weight: 800;
    color: #064949;
    letter-spacing: -1px;
}}
.auth-header .app-sub {{
    font-size: 14px;
    color: #64748b;
    margin-top: 4px;
}}

/* ── Tables ── */
table {{ width: 100%; border-collapse: collapse; color: #064949 !important; margin-top: 12px; }}
th {{
    background-color: #064949 !important;
    color: white !important;
    text-align: left;
    padding: 10px 14px;
    font-weight: 600;
    font-size: 12px;
    text-transform: uppercase;
    letter-spacing: 0.5px;
}}
[data-testid="stTable"] th,
    [data-testid="stTable"] thead tr th,
    .stDataFrame thead tr th,
    [data-testid="stTable"] table thead th,
    table thead th {{
    background-color: #064949 !important;
    color: white !important;
    font-weight: 600 !important;
}}
td {{ padding: 9px 14px; border-bottom: 1px solid #e8f5ee; color: #064949 !important; font-size: 13px; }}
[data-testid="stTable"] table td, [data-testid="stTable"] table th,
.stDataFrame td, .stDataFrame th {{ color: #064949 !important; }}
[data-testid="stTable"] {{
    background-color: white !important;
    border-radius: 10px;
    overflow: hidden;
    box-shadow: 0 1px 6px rgba(0,0,0,0.06);
}}

/* ── Metrics ── */
[data-testid="stMetricValue"] {{ color: #064949 !important; font-size: 26px; font-weight: 700; }}
[data-testid="stMetricLabel"] {{ color: #64748b !important; font-size: 11px; text-transform: uppercase; letter-spacing: 0.5px; }}
div[data-testid="stHorizontalBlock"] {{
    background: transparent !important;
    box-shadow: none !important;
    border: none !important;
    padding: 0 !important;
    border-radius: 0 !important;
}}
.streamlit-expanderContent div[data-testid="stHorizontalBlock"] {{
    background-color: white !important;
}}
.metrics-row div[data-testid="stHorizontalBlock"] {{
    background-color: white !important;
    border-radius: 12px !important;
    padding: 1.2rem 1rem !important;
    box-shadow: 0 2px 10px rgba(6,73,73,0.08) !important;
    border: 1px solid #e8f5ee !important;
}}

/* ── Expanders ── */
[data-testid="stExpander"] {{
    background-color: white !important;
    border-radius: 10px !important;
    border: 1px solid #e2e8f0 !important;
    overflow: hidden;
}}
[data-testid="stExpander"] summary {{
    background-color: #f8fffe !important;
    font-weight: 500;
}}
[data-testid="stExpander"] details > div {{
    background-color: white !important;
    padding: 1rem !important;
}}
[data-testid="stExpander"] div[data-testid="stHorizontalBlock"] {{
    background-color: white !important;
}}
.streamlit-expanderContent {{
    background-color: white !important;
    padding: 1rem;
    border-radius: 0 0 10px 10px;
    border: 1px solid #e2e8f0;
    border-top: none;
}}
.streamlit-expanderHeader {{
    background-color: #f8fffe !important;
    border-radius: 10px;
    padding: 0.6rem 1rem;
    border: 1px solid #e2e8f0;
    font-weight: 500;
}}

/* ── Profile cards ── */
.profile-card {{
    background: white;
    border-radius: 12px;
    padding: 1.4rem;
    box-shadow: 0 2px 10px rgba(6,73,73,0.07);
    border: 1px solid #e8f5ee;
    margin-bottom: 1rem;
}}
.subscription-active {{
    background: #f0fdf4; border-left: 4px solid #16a34a;
    padding: 0.75rem 1rem; border-radius: 8px; color: #15803d;
    font-weight: 600; margin-bottom: 0.75rem;
}}
.subscription-inactive {{
    background: #fff1f2; border-left: 4px solid #dc2626;
    padding: 0.75rem 1rem; border-radius: 8px; color: #dc2626;
    font-weight: 600; margin-bottom: 0.75rem;
}}
.trial-badge {{
    background: #fffbeb; border-left: 4px solid #d97706;
    padding: 0.75rem 1rem; border-radius: 8px; color: #92400e;
    font-weight: 600; margin-bottom: 0.75rem;
}}
</style>
""", unsafe_allow_html=True)

# ─── SESSION STATE ────────────────────────────────────────────────────────────
for k, v in [("user", None), ("token", None), ("profile", None), ("page", "login")]:
    if k not in st.session_state:
        st.session_state[k] = v

# ─── SIDEBAR ─────────────────────────────────────────────────────────────────
def render_sidebar():
    profile = st.session_state.profile
    user = st.session_state.user
    if not user or not profile:
        return

    with st.sidebar:
        # Logo + app name
        script_dir = os.path.dirname(os.path.abspath(__file__))
        logo_path = os.path.join(script_dir, "logo.png")
        logo_html = ""
        if os.path.exists(logo_path):
            logo_b64 = get_base64_of_bin_file(logo_path)
            logo_html = f'<img src="data:image/png;base64,{logo_b64}" />'

        st.markdown(f"""
        <div class="sidebar-logo">
            {logo_html}
            <div class="sidebar-appname">ApexPredict</div>
            <div class="sidebar-tagline">E. coli Pathotype Intelligence Platform</div>
        </div>
        """, unsafe_allow_html=True)

        # User info
        full_name = profile.get("full_name") or "User"
        institution = profile.get("institution") or ""
        email = user.get("email") or ""

        if profile.get("subscription_active"):
            end = profile.get("subscription_end", "")
            if end:
                end_dt = datetime.fromisoformat(end.replace("Z", "+00:00"))
                if end_dt >= datetime.now(timezone.utc):
                    badge = f'<span class="sidebar-badge-active">✅ Subscribed</span>'
                else:
                    badge = f'<span class="sidebar-badge-inactive">❌ Expired</span>'
            else:
                badge = f'<span class="sidebar-badge-active">✅ Subscribed</span>'
        elif is_on_trial(profile):
            remaining = trial_analyses_remaining(profile)
            badge = f'<span class="sidebar-badge-trial">🎁 Trial — {remaining} left</span>'
        else:
            badge = f'<span class="sidebar-badge-inactive">❌ Trial Ended</span>'

        st.markdown(f"""
        <div class="sidebar-user">
            <div class="sidebar-user-label">User</div>
            <div class="sidebar-user-name">{full_name}</div>
            <div class="sidebar-user-inst">{institution}</div>
            <div class="sidebar-user-email">{email}</div>
            {badge}
        </div>
        """, unsafe_allow_html=True)

        # Navigation
        st.markdown("<div>", unsafe_allow_html=True)

        page = st.session_state.page
        analysis_class = "sidebar-nav-active" if page == "app" else ""
        profile_class = "sidebar-nav-active" if page == "profile" else ""
        admin_class = "sidebar-nav-active" if page == "admin" else ""

        st.markdown(f'<div class="{analysis_class}">', unsafe_allow_html=True)
        if st.button("  Analysis", key="nav_analysis", use_container_width=True):
            st.session_state.page = "app"
            st.rerun()
        st.markdown("</div>", unsafe_allow_html=True)

        st.markdown(f'<div class="{profile_class}">', unsafe_allow_html=True)
        if st.button("  My Profile", key="nav_profile", use_container_width=True):
            st.session_state.page = "profile"
            st.rerun()
        st.markdown("</div>", unsafe_allow_html=True)

        if profile.get("is_admin"):
            st.markdown(f'<div class="{admin_class}">', unsafe_allow_html=True)
            if st.button("🛠  Admin Dashboard", key="nav_admin", use_container_width=True):
                st.session_state.page = "admin"
                st.rerun()
            st.markdown("</div>", unsafe_allow_html=True)

        st.markdown("</div>", unsafe_allow_html=True)
        st.markdown("<hr style='border-color:rgba(255,255,255,0.1);margin:0.5rem 0;'>", unsafe_allow_html=True)
        st.markdown('<div class="sidebar-signout">', unsafe_allow_html=True)
        if st.button("  Sign Out", key="nav_signout", use_container_width=True):
            for k in ["user", "token", "profile"]:
                st.session_state[k] = None
            st.session_state.page = "login"
            st.rerun()
        st.markdown("</div>", unsafe_allow_html=True)

# ─── AUTH PAGES ───────────────────────────────────────────────────────────────
def show_login():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    logo_path = os.path.join(script_dir, "logo.png")
    _, center, _ = st.columns([1, 1.2, 1])
    with center:
        if os.path.exists(logo_path):
            lc1, lc2, lc3 = st.columns([2,1,2])
            with lc2:
                st.image(logo_path, use_container_width=True)
        st.markdown("""
        <div class="auth-header">
            <div class="app-name">🦠 ApexPredict</div>
            <div class="app-sub">Pathotype Classification, MLST, Serogroup Detection & AMR Profiling</div>
        </div>
        """, unsafe_allow_html=True)
        st.markdown("**Sign In**")
        email = st.text_input("Email", key="login_email")
        password = st.text_input("Password", type="password", key="login_password")
        col1, col2 = st.columns(2)
        with col1:
            if st.button("Sign In", use_container_width=True, type="primary"):
                if email and password:
                    with st.spinner("Signing in..."):
                        resp = auth_login(email, password)
                    if "access_token" in resp:
                        token = resp["access_token"]
                        user = resp["user"]
                        st.session_state.token = token
                        st.session_state.user = user
                        profile = get_profile(token, user["id"])
                        # Sync metadata to profile on first login
                        meta = user.get("user_metadata") or {}
                        if profile and meta:
                            updates = {}
                            if not profile.get("full_name") and meta.get("full_name"):
                                updates["full_name"] = meta["full_name"]
                            if not profile.get("institution") and meta.get("institution"):
                                updates["institution"] = meta["institution"]
                            if updates:
                                requests.patch(
                                    f"{SUPABASE_URL}/rest/v1/profiles?id=eq.{user['id']}",
                                    headers={**supa_headers(token), "Prefer": "return=minimal"},
                                    json=updates)
                                profile = get_profile(token, user["id"])
                        st.session_state.profile = profile
                        st.session_state.page = "admin" if (profile and profile.get("is_admin")) else "app"
                        st.rerun()
                    else:
                        msg = resp.get("error_description", resp.get("msg", "Invalid credentials"))
                        st.error(f"❌ {msg}")
                else:
                    st.warning("Please enter your email and password.")
        with col2:
            if st.button("Create Account", use_container_width=True):
                st.session_state.page = "signup"
                st.rerun()

def show_signup():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    logo_path = os.path.join(script_dir, "logo.png")
    _, center, _ = st.columns([1, 1.2, 1])
    with center:
        if os.path.exists(logo_path):
            st.image(logo_path, width=56)
        st.markdown("""
        <div class="auth-header">
            <div class="app-name">ApexPredict</div>
            <div class="app-sub">Create your account</div>
        </div>
        """, unsafe_allow_html=True)
        st.info(" New accounts include **3 free trial analyses** — no subscription needed to start.")
        st.markdown("#### Sign Up")
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
                    st.error("❌ Passwords do not match.")
                elif len(password) < 6:
                    st.error("❌ Password must be at least 6 characters.")
                else:
                    with st.spinner("Creating account..."):
                        resp = auth_signup(email, password, full_name, institution)
                    if "id" in resp or ("user" in resp and resp.get("user")):
                        st.success("✓ Account created! Check your email to confirm, then sign in.")
                        st.session_state.page = "login"
                        st.rerun()
                    else:
                        msg = resp.get("error_description", resp.get("msg", "Signup failed."))
                        st.error(f"❌ {msg}")
        with col2:
            if st.button("Back to Sign In", use_container_width=True):
                st.session_state.page = "login"
                st.rerun()

# ─── PROFILE PAGE ─────────────────────────────────────────────────────────────
def show_profile():
    render_sidebar()
    token = st.session_state.token
    user = st.session_state.user
    if not st.session_state.profile and token and user:
        st.session_state.profile = get_profile(token, user["id"])
    profile = st.session_state.profile
    if not profile:
        st.error("Could not load profile. Please sign out and sign in again.")
        if st.button("Sign Out"):
            for k in ["user", "token", "profile"]:
                st.session_state[k] = None
            st.session_state.page = "login"
            st.rerun()
        return

    st.markdown('<div class="main-title">My Profile</div>', unsafe_allow_html=True)
    st.markdown("---")

    c1, c2 = st.columns(2)
    with c1:
        st.markdown(f"**Name:** {profile.get('full_name', 'N/A')}")
        st.markdown(f"**Email:** {user.get('email', 'N/A')}")
        st.markdown(f"**Institution:** {profile.get('institution', 'N/A')}")
        st.markdown(f"**Member since:** {profile.get('created_at', '')[:10]}")
    with c2:
        if profile.get("subscription_active"):
            end = profile.get("subscription_end", "")
            if end:
                end_dt = datetime.fromisoformat(end.replace("Z", "+00:00"))
                if end_dt >= datetime.now(timezone.utc):
                    st.markdown(f'<div class="subscription-active">✅ Subscription Active — Expires: {end[:10]}</div>', unsafe_allow_html=True)
                else:
                    st.markdown('<div class="subscription-inactive">❌ Subscription Expired</div>', unsafe_allow_html=True)
                    st.info("📧 **apexpredict@contact.com**")
            else:
                st.markdown('<div class="subscription-active">✅ Subscription Active — No expiry</div>', unsafe_allow_html=True)
        elif is_on_trial(profile):
            remaining = trial_analyses_remaining(profile)
            st.markdown(f'<div class="trial-badge"> Free Trial — {remaining} / {FREE_TRIAL_LIMIT} analyses remaining</div>', unsafe_allow_html=True)
            st.info("📧 To subscribe: **apexpredict@contact.com**")
        else:
            st.markdown('<div class="subscription-inactive">❌ Trial Ended — No Active Subscription</div>', unsafe_allow_html=True)
            st.info("📧 **apexpredict@contact.com**")

    st.markdown("---")
    st.markdown("###  Analysis History")

    history = get_history(token, user["id"])
    if not history or not isinstance(history, list):
        st.info("No analyses yet. Run your first analysis to see results here.")
    else:
        st.caption(f"{len(history)} analyses found")
        for rec in history:
            result = rec.get("result", "N/A")
            date = rec.get("created_at", "")[:10]
            icon = "🔴" if result == "APEC" else "🟢"
            with st.expander(f"{icon} {rec.get('filename', 'Unknown')} — {date} — {result}"):
                hc1, hc2, hc3, hc4 = st.columns(4)
                hc1.metric("Result", result)
                hc2.metric("Confidence", f"{rec.get('confidence', 0):.1f}%")
                hc3.metric("ST", f"ST{rec.get('mlst_st', 'N/A')}")
                hc4.metric("Serogroup", rec.get('serogroup', 'N/A'))
                st.caption(f"Serotype: {rec.get('serotype', 'N/A')} | AMR: {rec.get('amr_classes', 'None')} | Virulence genes: {rec.get('virulence_gene_count', 0)}")

# ─── ADMIN PAGE ───────────────────────────────────────────────────────────────
def show_admin():
    render_sidebar()
    profile = st.session_state.profile
    token = st.session_state.token
    if not profile or not profile.get("is_admin"):
        st.error("❌ Access denied.")
        return

    st.markdown('<div class="main-title">Admin Dashboard</div>', unsafe_allow_html=True)
    st.markdown('<div class="sub-title">Read-only view of all users and analyses</div>', unsafe_allow_html=True)
    st.markdown("---")

    all_users = get_all_users(token)
    all_analyses = get_all_analyses(token)

    total_users = len(all_users)
    active_subs = sum(1 for u in all_users if u.get("subscription_active"))
    trial_users = sum(1 for u in all_users if not u.get("subscription_active") and (u.get("trial_analyses_used") or 0) < FREE_TRIAL_LIMIT)
    total_analyses = len(all_analyses)
    apec_count = sum(1 for a in all_analyses if a.get("result") == "APEC")

    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("Total Users", total_users)
    m2.metric("Active Subscriptions", active_subs)
    m3.metric("On Free Trial", trial_users)
    m4.metric("Total Analyses", total_analyses)
    m5.metric("APEC Results", f"{apec_count}/{total_analyses}" if total_analyses else "0/0")

    st.markdown("---")
    st.markdown("### 👥 All Users")
    if all_users:
        users_df = pd.DataFrame([{
            "Email": u.get("email", ""),
            "Full Name": u.get("full_name", ""),
            "Institution": u.get("institution", ""),
            "Joined": u.get("created_at", "")[:10],
            "Subscription": "✅ Active" if u.get("subscription_active") else "❌ Inactive",
            "Sub Expires": u.get("subscription_end", "")[:10] if u.get("subscription_end") else "—",
            "Trial Used": f"{u.get('trial_analyses_used', 0) or 0}/{FREE_TRIAL_LIMIT}",
            "Admin": "✅" if u.get("is_admin") else "",
        } for u in all_users])
        st.dataframe(users_df, use_container_width=True, hide_index=True)
    else:
        st.info("No users found.")

    st.markdown("---")
    st.markdown("###  All Analyses")
    if all_analyses:
        analyses_df = pd.DataFrame([{
            "Date": a.get("created_at", "")[:10],
            "User": a.get("profiles", {}).get("email", "") if isinstance(a.get("profiles"), dict) else "",
            "File": a.get("filename", ""),
            "Result": a.get("result", ""),
            "Confidence": f"{a.get('confidence', 0):.1f}%",
            "ST": f"ST{a.get('mlst_st', 'N/A')}",
            "Serogroup": a.get("serogroup", ""),
            "AMR": a.get("amr_classes", "None"),
            "Virulence Genes": a.get("virulence_gene_count", 0),
            "High Risk ST": "⚠️" if a.get("is_high_risk_st") else "",
        } for a in all_analyses])
        st.dataframe(analyses_df, use_container_width=True, hide_index=True)
        csv = analyses_df.to_csv(index=False).encode("utf-8")
        st.download_button("⬇️ Download Analyses CSV", csv, "all_analyses.csv", "text/csv")
    else:
        st.info("No analyses yet.")

# ─── SUBSCRIPTION WALL ────────────────────────────────────────────────────────
def show_subscription_wall():
    st.markdown("""
    <div style='text-align:center; padding: 3rem 1rem;'>
        <div style='font-size:48px; margin-bottom:0.5rem;'>🔒</div>
        <h2 style='color:#064949;'>Access Restricted</h2>
        <p style='font-size:16px; color:#4a5568; max-width:480px; margin:0 auto;'>
            Your free trial has ended or your subscription is not active.<br>
            You can still view your analysis history in your profile.
        </p>
        <p style='font-size:15px; color:#064949; margin-top:1.5rem;'>
            <strong>To subscribe, contact us at:<br>apexpredict@contact.com</strong>
        </p>
    </div>
    """, unsafe_allow_html=True)
    _, c, _ = st.columns([2, 1, 2])
    with c:
        if st.button("View My Profile", use_container_width=True, type="primary"):
            st.session_state.page = "profile"
            st.rerun()

# ─── MAIN APP PAGE ────────────────────────────────────────────────────────────
def show_app():
    render_sidebar()
    profile = st.session_state.profile
    token = st.session_state.token
    user = st.session_state.user

    if not is_subscription_active(profile):
        show_subscription_wall()
        return

    st.markdown('<div class="main-title"> 🦠 ApexPredict</div>', unsafe_allow_html=True)
    st.markdown('<div class="sub-title">Pathotype Classification, MLST, Serogroup Detection, and Resistance Profiling</div>', unsafe_allow_html=True)

    if is_on_trial(profile):
        remaining = trial_analyses_remaining(profile)
        st.markdown(
            f'<div class="trial-badge"> Free Trial — {remaining} / {FREE_TRIAL_LIMIT} analyses remaining. '
            f'Contact <strong>apexpredict@contact.com</strong> to subscribe.</div>',
            unsafe_allow_html=True)

    @st.cache_resource
    def load_analytical_assets():
        model = joblib.load("apec_random_forest_model.joblib")
        features = joblib.load("model_features.joblib")
        return model, features

    model, model_features = load_analytical_assets()
    ectyper_available = check_ectyper_available()

    uploaded_file = st.file_uploader("Upload Bacterial FASTA", type=["fna", "fasta"])

    if uploaded_file is not None:
        temp_path = "temp_sample.fna"
        raw_bytes = uploaded_file.getbuffer().tobytes()
        raw_text = raw_bytes.decode("utf-8", errors="replace").replace("\r\n", "\n").replace("\r", "\n")
        clean_lines = [l.strip() for l in raw_text.splitlines() if l.strip()]
        with open(temp_path, "w", newline="\n") as f:
            f.write("\n".join(clean_lines) + "\n")

        status_indicator = st.empty()
        status_indicator.info("Running Genomic Analysis (MLST, Virulence, AMR, Serogroup) ...")

        mlst_st = "Unknown"
        try:
            mlst_proc = subprocess.run(f"mlst {temp_path}", shell=True,
                stdout=subprocess.PIPE, universal_newlines=True, timeout=30)
            if mlst_proc.stdout:
                mlst_data = mlst_proc.stdout.strip().split('\t')
                if len(mlst_data) >= 3:
                    mlst_st = mlst_data[2]
        except Exception as e:
            st.warning(f"MLST failed: {e}")

        is_high_risk = mlst_st in HIGH_RISK_STS
        risk_label = "⚠️ HIGH RISK" if is_high_risk else "Standard"

        serogroup, serotype, serogroup_risk = "Unknown", "Unknown", "Unknown"
        if ectyper_available:
            try:
                serogroup, serotype = detect_serogroup_with_ectyper(temp_path)
                if serogroup in APEC_SEROGROUPS: serogroup_risk = "⚠️ APEC-associated"
                elif serogroup not in ["Unknown", "Detection Failed", "Timeout", "Error"]: serogroup_risk = "ℹ️ Other serogroup"
                else: serogroup_risk = "❓ Undetermined"
            except Exception:
                serogroup, serotype, serogroup_risk = "Detection Error", "Detection Error", "Error"
        else:
            serogroup = "Not available"
            serotype = "Install ECTyper for detection"
            serogroup_risk = "Not available"

        detailed_results, resistance_results = [], []
        current_dir = os.path.dirname(os.path.abspath(__file__))
        db_directory = os.path.join(current_dir, "apec_db")

        for database in ["ecoli_vf", "apec"]:
            cmd = f"abricate --db ecoli_vf {temp_path}" if database == "ecoli_vf" else f"abricate --datadir {db_directory} --db apec {temp_path}"
            try:
                proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                    universal_newlines=True, timeout=60)
                if proc.stdout:
                    lines = proc.stdout.strip().split('\n')
                    for entry in lines[1:]:
                        cols = entry.split('\t')
                        if len(cols) >= 10:
                            clean_name = cols[5].replace('~~~', '|').replace('~', '|').split('|')[0].strip().upper()
                            detailed_results.append({"Gene Name": clean_name, "Contig": cols[1],
                                "Identity (%)": cols[9], "Coverage (%)": cols[6]})
            except Exception as e:
                st.warning(f"Virulence screening for {database} failed: {e}")

        try:
            amr_proc = subprocess.run(f"abricate --db ncbi {temp_path}", shell=True,
                stdout=subprocess.PIPE, universal_newlines=True, timeout=60)
            if amr_proc.stdout:
                for entry in amr_proc.stdout.strip().split('\n')[1:]:
                    cols = entry.split('\t')
                    if len(cols) >= 14:
                        resistance_results.append({"Gene": cols[5], "Class": cols[13], "Identity (%)": cols[9]})
        except Exception as e:
            st.warning(f"AMR screening failed: {e}")

        status_indicator.empty()

        found_genes_list = [r['Gene Name'].strip().upper() for r in detailed_results]
        input_row = {feat: [1 if feat in found_genes_list else 0] for feat in model_features}
        prob_score = model.predict_proba(pd.DataFrame(input_row))[0][1]
        prediction = "APEC" if prob_score >= 0.5 else "Non-APEC"
        confidence_pct = prob_score * 100 if prediction == "APEC" else (1 - prob_score) * 100

        # Save to history
        amr_classes_str = ", ".join(set([r["Class"] for r in resistance_results])) if resistance_results else "None"
        save_analysis(token=token, user_id=user["id"], filename=uploaded_file.name,
            result=prediction, confidence=round(confidence_pct, 1), mlst_st=mlst_st,
            serogroup=serogroup, serotype=serotype, amr_classes=amr_classes_str,
            virulence_gene_count=len(detailed_results), is_high_risk_st=is_high_risk)

        # Increment trial if on trial
        if is_on_trial(profile):
            used = profile.get("trial_analyses_used", 0) or 0
            increment_trial_usage(token, user["id"], used)
            st.session_state.profile = get_profile(token, user["id"])

        st.markdown("---")
        st.markdown('<div class="metrics-row">', unsafe_allow_html=True)
        m_col1, m_col2, m_col3, m_col4 = st.columns(4)
        with m_col1:
            if prediction == "APEC": st.error(f"### Result: **{prediction}**")
            else: st.success(f"### Result: **{prediction}**")
        with m_col2:
            st.metric("Lineage (MLST)", f"ST{mlst_st}", delta=risk_label,
                delta_color="inverse" if is_high_risk else "normal")
        with m_col3:
            st.metric("ML Confidence", f"{confidence_pct:.1f}%")
        with m_col4:
            st.metric("Serogroup", serogroup, delta=serogroup_risk,
                delta_color="inverse" if serogroup in APEC_SEROGROUPS else "normal")

        st.markdown('</div>', unsafe_allow_html=True)
        with st.expander("Detailed Serogroup Information", expanded=False):
            col1, col2 = st.columns(2)
            with col1:
                st.write("**Detected Serogroup:**")
                st.write(f"### {serogroup}")
                st.write(f"**Full Serotype:** {serotype}")
            with col2:
                if serogroup in APEC_SEROGROUPS:
                    st.success(f" {serogroup} is a known APEC-associated serogroup")
                    st.info("This serogroup is commonly found in avian pathogenic E. coli strains")
                elif serogroup not in ["Unknown", "Detection Failed", "Timeout", "Error", "Not available", "Detection Error"]:
                    st.warning(f"⚠️ {serogroup} is not typically associated with APEC")
                else:
                    st.info("Serogroup could not be determined from the provided sequence")

        if resistance_results:
            st.markdown("### Antimicrobial Resistance Profile")
            res_df = pd.DataFrame(resistance_results)
            st.warning(f"**Predicted Resistance to:** {', '.join(res_df['Class'].unique())}")
            with st.expander("Show AMR Genes"):
                st.table(res_df)
        else:
            st.success("No AMR markers detected (Susceptible)")

        if detailed_results:
            with st.expander(f"Show {len(detailed_results)} Virulence Factors", expanded=False):
                vir_df = pd.DataFrame(detailed_results)
                st.download_button("Download Results (CSV)",
                    vir_df.to_csv(index=False).encode('utf-8'),
                    f"apex_ST{mlst_st}_{serogroup}_{prediction}.csv", "text/csv")
                st.table(vir_df)

        if os.path.exists(temp_path): os.remove(temp_path)
        for tf in ["temp_ectyper_results.txt", "ectyper_output.json"]:
            if os.path.exists(tf): os.remove(tf)

# ─── ROUTER ───────────────────────────────────────────────────────────────────
if st.session_state.page == "login":
    show_login()
elif st.session_state.page == "signup":
    show_signup()
elif st.session_state.page == "profile":
    show_profile()
elif st.session_state.page == "admin":
    show_admin()
elif st.session_state.page == "app":
    show_app()
