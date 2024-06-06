import time

from loguru import logger
import streamlit as st
from pathlib import Path
import tempfile

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
    st.markdown("# Simple app")
    

    uploaded_file = st.file_uploader("Choose a structural file", type=["pdb", "gro"])
    if (uploaded_file is not None) and st.button("Run new analysis") : 
        temp_dir = tempfile.mkdtemp()
        path = Path(temp_dir) / uploaded_file.name
        with open(path, "wb") as f:
            f.write(uploaded_file.getvalue())
        st.write(f"Saved to {path}")
        with st.spinner("Running analysis..."):
            gd.main(
                path,
                draw_graph_option=False,
                check_overlapping_residue=False,
                check_connectivity=False,
                bond_threshold="auto",
                query_pdb=True,
            )

    
    # if st.button("Run new analysis"):
    #     timer.reset()
    #     with st.spinner("Running analysis..."):
    #         print_in_log(5)
    #     st.write("Done!")