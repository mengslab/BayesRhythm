import streamlit as st
def render_nature_methods_dashboard(cfg, readiness):
    st.markdown('<style>.nm-hero{background:linear-gradient(135deg,#0f172a 0%,#1d4ed8 48%,#0ea5e9 100%);color:white;border-radius:24px;padding:1.2rem 1.4rem;margin-bottom:1rem;}.block-container{max-width:1500px;padding-top:1.2rem;}</style>', unsafe_allow_html=True)
    st.markdown('<div class="nm-hero"><div style="font-size:2rem;font-weight:800;">BayesRhythm v3.2.0</div><div>Real-data execution hardening for large Chow/TRF single-cell archives.</div></div>', unsafe_allow_html=True)
    cols=st.columns(len(readiness))
    for i,(k,v) in enumerate(readiness.items()): cols[i].metric(k, v.capitalize())
