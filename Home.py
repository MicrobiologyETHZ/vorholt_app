import streamlit as st
import pandas as pd
from check_password import check_password


st.set_page_config(page_title="NCCR Microbiomes ETHZ", layout='wide')


def home_page():
    st.title("RNA-sequencing projects")
    st.table(pd.read_csv("projects.csv"))


if check_password():
    home_page()



