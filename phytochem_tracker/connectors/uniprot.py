"""UniProt REST: proteins with a given EC number in a given plant taxon.

Complements KEGG for species without a KEGG genome: UniProt includes
TrEMBL entries translated from any sequenced genome or transcript.
"""

from __future__ import annotations

import csv
import io

from .. import config
from ..http import Client

SOURCE = "UniProt REST"
FIELDS = "accession,gene_names,protein_name,organism_name,reviewed,xref_phytozome"


def build_query(ec: str, taxid: str = "", species: str = "") -> str:
    if taxid:
        org = f"taxonomy_id:{taxid}"
    else:
        org = f'organism_name:"{species}"'
    return f"(ec:{ec}) AND ({org})"


def proteins_for_ec(client: Client, ec: str, taxid: str = "", species: str = "", size: int = 10) -> list[dict[str, str]]:
    params = {"query": build_query(ec, taxid, species), "fields": FIELDS, "format": "tsv", "size": size}
    resp = client.request("GET", f"{config.UNIPROT_BASE}/uniprotkb/search", source=SOURCE,
                          params=params, ok_statuses=(200, 400))
    if resp.status_code == 400:  # e.g. unknown field name -> retry with the basic field set
        params["fields"] = "accession,gene_names,protein_name,organism_name,reviewed"
        resp = client.request("GET", f"{config.UNIPROT_BASE}/uniprotkb/search", source=SOURCE, params=params)
    rows = list(csv.DictReader(io.StringIO(resp.text), delimiter="\t"))
    out = []
    for r in rows:
        out.append({
            "accession": r.get("Entry", ""),
            "gene_names": r.get("Gene Names", ""),
            "protein": r.get("Protein names", ""),
            "organism": r.get("Organism", ""),
            "reviewed": r.get("Reviewed", ""),
            "phytozome": r.get("Phytozome", "").strip(";"),
        })
    # Reviewed (Swiss-Prot) entries first
    out.sort(key=lambda r: r["reviewed"] != "reviewed")
    return out
