import time

from loguru import logger
import streamlit as st
from pathlib import Path
import tempfile
import os

import grodecoder as gd


class LauncherTime:
    """Class to measure elapsed time.
    
    Modified from:
    https://github.com/Delgan/loguru/issues/867#issuecomment-1537944764
    """
    def __init__(self):
        self._start = time.perf_counter()

    def reset(self):
        self._start = time.perf_counter()

    @property
    def elapsed(self):
        return f"{time.perf_counter() - self._start:.2f}"
    
    
def print_in_log(iterations):
    for i in range(iterations):
        logger.info(f"hello {i}")
        time.sleep(1)


def last_line_log_file():
    with open("test.log", "rb") as file:
        try:
            file.seek(-2, os.SEEK_END)
            while file.read(1) != b'\n':
                file.seek(-2, os.SEEK_CUR)
        except OSError:
            file.seek(0)
        last_line = file.readline().decode()
    return last_line


if __name__ == "__main__":
    timer = LauncherTime()

    logger.remove()
    # logger.add(st.write, format="{time:YYYY-MM-DD HH:mm:ss} | {extra[timer].elapsed} s | {level} | {message}")
    # logger.add("test.log", format="{time:YYYY-MM-DD HH:mm:ss} | {extra[timer].elapsed} s | {level} | {message}")
    logger.add(st.write, format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}")
    logger.add("test.log", format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}")
    # logger = logger.bind(timer=timer)
    st.markdown(
        """
        <style>
        p {
            margin-bottom: 0px;
        }
        </style>
        """, unsafe_allow_html=True)
    st.markdown("# First version of GroDecoder app")
    

    uploaded_file = st.file_uploader("Choose a structural file", type=["pdb", "gro"])

    st.write("Choose some options for the analysis : ")
    draw_graph_option_value = st.checkbox("Draw graph", value=False)
    check_overlapping_residue_value = st.checkbox("Check overlapping residue", value=False)
    check_connectivity_value = st.checkbox("Check connectivity (add degree and number of edge in the fingerprint) ", value=False)
    bond_threshold_value = st.text_input("Enter auto or a threshold value", value="auto")

    if (uploaded_file is not None) and st.button("Run new analysis") : 
        temp_dir = tempfile.mkdtemp()
        path = Path(temp_dir) / uploaded_file.name
        with open(path, "wb") as f:
            f.write(uploaded_file.getvalue())
        st.write(f"Saved to {path}")
        with st.spinner("Running analysis..."):
            gd.main(
                path,
                draw_graph_option=draw_graph_option_value,
                check_overlapping_residue=check_overlapping_residue_value,
                check_connectivity=check_connectivity_value,
                bond_threshold=bond_threshold_value,
                query_pdb=True,
            )

    # Fetch the last line of the log file
    last_line = last_line_log_file()
    last_line_filename = last_line.split()[-1]

    # Test if the last line of the log file is a JSON
    # That also mean that the analysis terminated correctly
    if last_line_filename.split('.')[-1] == "json":
        with open(last_line_filename, 'r') as file:
            file_content = file.read()

        st.download_button(
            label="Download the inventory file for this structural file",
            data=file_content,
            file_name=last_line_filename,
            mime="text/plain"
        )
    
    # if st.button("Run new analysis"):
    #     timer.reset()
    #     with st.spinner("Running analysis..."):
    #         print_in_log(5)
    #     st.write("Done!")