import subprocess
import pandas as pd
import pytest
import os

# Run command to generate results first
@pytest.fixture(scope="module", autouse=True)
def generate_results():
    command = ["nextflow","run","main.nf","-c","test/test.config","-profile","ci"]
    result = subprocess.run(command, capture_output=True, text=True)
    print(result.stdout)
    if result.stderr:
        print(result.stderr)
    return result

# Check if the command was successful
def test_run_completed(generate_results):
    assert generate_results.returncode == 0

def test_ld_block_output(generate_results):
    ld_blocks = pd.read_csv("results/ld_blocks/snps_with_LD_blocks.csv")
    assert len(ld_blocks.RSID) == 3, "Expected 3 SNPs in output LD blocks file."
    assert len(ld_blocks.columns) == 8, "Expected 8 columns in LD blocks file."
    assert ld_blocks.columns.tolist() == ['RSID','CHR','POS','LDblock_lower','LDblock_upper','LDblock_length','lower_bound','upper_bound'], "Column names do not match expected in LD blocks file."

def test_screeplot_exists(generate_results):
    assert os.path.exists("results/pca/ScreePlot.pdf"), "Expected file ScreePlot.pdf from pipeline run."
