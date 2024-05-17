import pandas as pd
from bs4 import BeautifulSoup

URL = "http://argonaute.mit.edu/cgi-bin/genscanw_py.cgi"
ORGANISMS = ["Vertebrate", "Arabidopsis", "Maize"]
EXON_CUTOFFS = [1.00, 0.50, 0.25, 0.10, 0.05, 0.02, 0.01]


def read_seq_from_file(path_to_file: str) -> str:
    seq = ""
    with open(path_to_file) as fna:
        for line in fna.readlines()[1:]:
            line = line.strip()
            seq += line
    return seq


def obtain_genes(found_data: BeautifulSoup) -> str | list:
    genes_raw = found_data.text[
        found_data.text.find("Predicted genes/exons") : found_data.text.find(
            "Suboptimal exons with probability"
        )
    ].split("\n\n")[5:-2]

    exception = "NO EXONS/GENES PREDICTED IN SEQUENCE"
    if exception in genes_raw:
        genes = exception
    else:
        genes = []
        genes_raw = [gene for gene in genes_raw if gene]
        for gene in genes_raw:
            gene_listed = [value for value in gene.split(" ") if value]
            genes.append([gene_listed[0], gene_listed[3], gene_listed[4]])

    return genes


def parse_soup(soup: BeautifulSoup) -> tuple:
    found_data = soup.find("pre")

    predicted_peptide = "".join(
        found_data.text[found_data.text.find("predicted_peptide") :].split("\n\n")[1:-1]
    )
    if predicted_peptide == "":
        predicted_peptide = "NO PEPTIDES PREDICTED"

    genes = obtain_genes(found_data)
    if isinstance(genes, str):
        exons = genes
        introns = genes
    elif isinstance(genes, list):
        df = pd.DataFrame(data=genes, columns=["number", "start", "end"])
        df[["start", "end"]] = df[["start", "end"]].astype(int)

        exons = df.values.tolist()
        introns = []
        for idx in range(df.shape[0] - 1):
            introns.append(
                [
                    df.iloc[idx]["number"],
                    df.iloc[idx]["end"] + 1,
                    df.iloc[idx + 1]["start"] - 1,
                ]
            )
    else:
        raise ValueError("Something went wrong")

    return predicted_peptide, exons, introns
