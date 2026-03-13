#!/usr/bin/env python3
"""
MetaFuncDecoder: Ecologically interpretable functional profiles
from already-annotated (meta)genomes.

Accepts JGI KOfamScan files, JGI/HMMER Pfam outputs, dbCAN overview files,
and pre-merged annotation tables. NCBI-derived annotations can be used via
combined-table mode after harmonizing column names.
Maps KO terms, Pfam IDs, and CAZy families to functional categories defined
in user-editable YAML files (see categories/).

Two input modes:
  Join mode     -- separate annotation files joined on gene_id
                   (--ko-file, --pfam-file, --cazy-file)
  Combined mode -- single pre-merged table with configurable column names
                   (--combined-table, --ko-col, --pfam-col, --cazy-col)
"""

import argparse
import datetime
import hashlib
import json
import re
import sys
from abc import ABC, abstractmethod
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import yaml

VERSION = "1.0.0"

# Script directory — used to resolve default database and category paths
# regardless of the working directory the user invokes the tool from.
_SCRIPT_DIR = Path(__file__).parent


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_patterns(raw_list) -> list:
    """Normalise a raw YAML pattern list to a list of (pattern_str, subcategory) tuples.

    Accepts:
        - plain strings:  "cellulase"            → ("cellulase", "")
        - subcategory dicts: {pattern: "...", subcategory: "..."} → (str, subcat)
        - already-normalised tuples (pass-through)
    """
    result = []
    for item in (raw_list or []):
        if isinstance(item, str):
            result.append((item, ""))
        elif isinstance(item, dict):
            pat = item.get("pattern", "")
            if not pat:
                print(f"  Warning: YAML pattern entry missing 'pattern' key — skipping: {item}")
                continue
            result.append((pat, item.get("subcategory", "")))
        elif isinstance(item, (tuple, list)) and len(item) == 2:
            result.append((item[0], item[1]))
    return result


# ---------------------------------------------------------------------------
# Category loading
# ---------------------------------------------------------------------------

def load_categories(categories_dir: str) -> Dict[str, Any]:
    """Load all YAML category files from directory.

    Returns a dict keyed by category name, each containing:
        kegg_patterns, pfam_patterns, pfam_go_terms, cazy_subcategories
    """
    categories_path = Path(categories_dir)
    if not categories_path.exists():
        print(f"Warning: Categories directory not found: {categories_dir}")
        return {}

    categories = {}
    yaml_files = sorted(categories_path.glob("*.yaml"))

    if not yaml_files:
        print(f"Warning: No .yaml files found in {categories_dir}")
        return {}

    print(f"Loading functional categories from: {categories_dir}")
    for yaml_file in yaml_files:
        try:
            with open(yaml_file, "r") as f:
                cat_def = yaml.safe_load(f)

            cat_name = cat_def.get("name")
            if not cat_name:
                print(f"  Warning: {yaml_file.name} has no 'name' field — skipping")
                continue

            categories[cat_name] = {
                "description": cat_def.get("description", ""),
                "citations": cat_def.get("citations", []),
                "kegg_patterns": cat_def.get("kegg_patterns") or [],
                "pfam_patterns": cat_def.get("pfam_patterns") or [],
                "pfam_go_terms": cat_def.get("pfam_go_terms") or [],
                "cazy_subcategories": cat_def.get("cazy_subcategories") or {},
            }
            n_cazy = len(categories[cat_name]["cazy_subcategories"])
            n_kegg = len(categories[cat_name]["kegg_patterns"])
            print(f"  {cat_name}: {n_kegg} KEGG patterns, {n_cazy} CAZy subcategories")

        except Exception as e:
            print(f"  Error loading {yaml_file.name}: {e}")

    return categories


# ---------------------------------------------------------------------------
# Base mapper
# ---------------------------------------------------------------------------

class FunctionalMapper(ABC):
    """Base class for functional annotation mappers with pattern matching."""

    def __init__(self, patterns: Dict[str, List[str]]):
        self.activity_patterns = patterns
        self.compiled_patterns = self._compile_patterns(patterns)
        self.custom_patterns: Dict[str, List[str]] = {}
        self.compiled_custom_patterns: Dict[str, List[re.Pattern]] = {}

    def _compile_patterns(
        self, pattern_dict: Dict[str, list]
    ) -> Dict[str, List[Tuple[re.Pattern, str]]]:
        """Compile patterns to (regex, subcategory) tuples. Accepts str or (str, subcat)."""
        compiled: Dict[str, List[Tuple[re.Pattern, str]]] = {}
        for cat, pats in pattern_dict.items():
            compiled[cat] = []
            for item in pats:
                if isinstance(item, str):
                    compiled[cat].append((re.compile(item, re.IGNORECASE), ""))
                elif isinstance(item, (tuple, list)) and len(item) == 2:
                    compiled[cat].append((re.compile(item[0], re.IGNORECASE), item[1]))
        return compiled

    def add_custom_patterns(self, custom_patterns_dict: Dict[str, List]) -> None:
        converted = {cat: _parse_patterns(pats) for cat, pats in custom_patterns_dict.items()}
        self.custom_patterns.update(converted)
        self.compiled_custom_patterns.update(self._compile_patterns(converted))
        print(f"Added {len(custom_patterns_dict)} custom pattern categories")

    def _search_patterns_in_text(
        self, text: str, include_custom: bool = True
    ) -> List[Tuple[str, str]]:
        """Return list of (category, subcategory) for all matching patterns."""
        if not text or pd.isna(text):
            return []
        text_lower = str(text).lower()
        results: List[Tuple[str, str]] = []

        for cat, compiled in self.compiled_patterns.items():
            for pattern, subcat in compiled:
                if pattern.search(text_lower):
                    results.append((cat, subcat))
                    break

        if include_custom:
            for cat, compiled in self.compiled_custom_patterns.items():
                for pattern, subcat in compiled:
                    if pattern.search(text_lower):
                        results.append((cat, subcat))
                        break

        return results

    def get_compiled_pattern(self, pattern: str) -> re.Pattern:
        return re.compile(pattern, re.IGNORECASE)

    @abstractmethod
    def search_by_pattern(
        self, search_pattern: str, category_name: str = "custom_search"
    ) -> Dict[str, Any]:
        pass


# ---------------------------------------------------------------------------
# KEGG mapper
# ---------------------------------------------------------------------------

class KEGGMapper(FunctionalMapper):
    """Map KEGG ortholog IDs to functional categories via definition text."""

    def __init__(self, patterns: Dict[str, List[str]]):
        super().__init__(patterns)
        self.kegg_info: Dict[str, Dict] = {}

    def search_by_pattern(
        self, search_pattern: str, category_name: str = "custom_search"
    ) -> Dict[str, Any]:
        compiled = self.get_compiled_pattern(search_pattern)
        matches = {
            k: v for k, v in self.kegg_info.items()
            if compiled.search(v.get("definition", ""))
        }
        print(f"\nKEGG orthologs matching '{search_pattern}': {len(matches)} found")
        for k, v in list(matches.items())[:5]:
            print(f"  {k}: {v['definition'][:80]}")
        return matches

    def load_kegg_info(self, filepath: str = "data/ko_list") -> bool:
        print(f"Loading KEGG ortholog information from: {filepath}")
        if not Path(filepath).exists():
            print(f"  Warning: file not found at {filepath}")
            return False
        try:
            df = pd.read_csv(filepath, sep="\t", low_memory=False)
            for _, row in df.iterrows():
                k_number = row.get("knum", "")
                definition = row.get("definition", "")
                if k_number and definition:
                    ec_match = re.search(
                        r"\[EC:([\d\.-]+(?:\s+[\d\.-]+)*)\]", definition
                    )
                    self.kegg_info[k_number] = {
                        "definition": definition,
                        "ec": ec_match.group(1) if ec_match else "",
                        "threshold": row.get("threshold", ""),
                    }
            print(f"  Loaded {len(self.kegg_info)} KEGG orthologs")
            return True
        except Exception as e:
            print(f"  Error: {e}")
            return False

    def get_kegg_categories(self, kegg_annotation: str) -> List[Tuple[str, str]]:
        """Return list of (category, subcategory) tuples for a KO term."""
        k_match = re.search(r"(K\d+)", str(kegg_annotation))
        if not k_match:
            return [("other", "")]
        info = self.kegg_info.get(k_match.group(1))
        if not info:
            return [("other", "")]
        hits = self._search_patterns_in_text(info["definition"])
        return hits if hits else [("other", "")]

    def get_detailed_kegg_info(self, kegg_annotation: str) -> Optional[Dict]:
        k_match = re.search(r"(K\d+)", str(kegg_annotation))
        if not k_match:
            return None
        return self.kegg_info.get(k_match.group(1))


# ---------------------------------------------------------------------------
# CAZy mapper
# ---------------------------------------------------------------------------

class CAZyMapper(FunctionalMapper):
    """Map CAZy family IDs to functional subcategories and broad categories."""

    def __init__(
        self,
        patterns: Dict[str, List[str]],
        subcategory_to_broad: Dict[str, str],
        family_to_subcategory: Optional[Dict[str, List[str]]] = None,
    ):
        super().__init__(patterns)
        self.family_info: Dict[str, Dict] = {}
        self.subcategory_to_broad = subcategory_to_broad
        # Option C: explicit family ID → subcategory list loaded from YAML.
        # A family may map to multiple subcategories (e.g. heterogeneous families
        # like GH16 that legitimately span two subcategory groups).
        self.family_to_subcategory: Dict[str, List[str]] = family_to_subcategory or {}

    def search_by_pattern(
        self, search_pattern: str, category_name: str = "custom_search"
    ) -> Dict[str, Any]:
        compiled = self.get_compiled_pattern(search_pattern)
        matches = {
            fam: info for fam, info in self.family_info.items()
            if compiled.search(
                f"{info.get('activities', '')} {info.get('note', '')}"
            )
        }
        print(f"\nCAZy families matching '{search_pattern}': {len(matches)} found")
        return matches

    def load_cazy_info(self, filepath: str = "data/FamInfo.txt.04232020.tsv") -> bool:
        print(f"Loading CAZy family information from: {filepath}")
        if not Path(filepath).exists():
            print(f"  Warning: file not found at {filepath}")
            return False
        try:
            df = pd.read_csv(filepath, sep="\t", low_memory=False)
            for _, row in df.iterrows():
                family = row.get("Family", "")
                if family:
                    self.family_info[family] = {
                        "class": row.get("cazy-class", ""),
                        "note": row.get("cazy-note", ""),
                        "activities": row.get("cazy-activities", ""),
                    }
            print(f"  Loaded {len(self.family_info)} CAZy families")
            return True
        except Exception as e:
            print(f"  Error: {e}")
            return False

    def get_family_categories(self, cazy_annotation: str) -> List[str]:
        """Return subcategory names (e.g. 'cellulase'), not broad categories.

        Lookup order:
          1. Explicit family ID list from YAML (families: field) — authoritative.
          2. Text pattern matching against FamInfo activity description — fallback.
        """
        family_match = re.match(r"([A-Z]+\d+(?:_\w+)?)", str(cazy_annotation))
        if not family_match:
            return ["other"]
        family = family_match.group(1)

        # Option C: explicit family ID takes priority; bypasses text matching.
        # Returns all subcategories the family is listed under (may be >1 for
        # heterogeneous families like GH16).
        if family in self.family_to_subcategory:
            return self.family_to_subcategory[family]

        if family not in self.family_info:
            return ["other"]
        hits = self._search_patterns_in_text(self.family_info[family]["activities"])
        cats = [cat for cat, _ in hits]
        return cats if cats else ["other"]

    def get_broad_category(self, subcategory: str) -> str:
        """Map subcategory name to broad category using YAML-defined hierarchy."""
        return self.subcategory_to_broad.get(subcategory, subcategory)

    def get_detailed_family_info(self, cazy_annotation: str) -> Optional[Dict]:
        family_match = re.match(r"([A-Z]+\d+)", str(cazy_annotation))
        if not family_match:
            return None
        return self.family_info.get(family_match.group(1))


# ---------------------------------------------------------------------------
# Pfam mapper
# ---------------------------------------------------------------------------

class Pfam2GOMapper(FunctionalMapper):
    """Map Pfam IDs to functional categories via GO term names and IDs."""

    def __init__(
        self,
        patterns: Dict[str, List[str]],
        go_terms: Dict[str, List[str]],
    ):
        super().__init__(patterns)
        self.pfam2go: Dict[str, List[Dict[str, str]]] = {}
        self.functional_categories = go_terms  # {category: ["GO:XXXXXXX", ...]}

    def search_by_pattern(
        self, search_pattern: str, category_name: str = "custom_search"
    ) -> Dict[str, Any]:
        compiled = self.get_compiled_pattern(search_pattern)
        matches: Dict[str, List] = {}
        for pfam_id, go_terms in self.pfam2go.items():
            for go_term in go_terms:
                if compiled.search(go_term["go_name"]):
                    matches.setdefault(pfam_id, []).append(go_term)
        print(f"\nPfam families matching '{search_pattern}': {len(matches)} found")
        return matches

    def load_pfam2go(self, filepath: str = "data/pfam2go.txt") -> bool:
        print(f"Loading Pfam2GO mapping from: {filepath}")
        if not Path(filepath).exists():
            print(f"  Warning: file not found at {filepath}")
            return False
        try:
            with open(filepath, "r") as f:
                for line in f:
                    if line.startswith("!") or not line.strip():
                        continue
                    m = re.match(
                        r"Pfam:(PF\d+)\s+\S+\s+>\s+GO:(.+?)\s*;\s*(GO:\d+)", line
                    )
                    if m:
                        pfam_id, go_name, go_id = m.groups()
                        self.pfam2go.setdefault(pfam_id, []).append(
                            {"go_id": go_id, "go_name": go_name.strip()}
                        )
            print(f"  Loaded mappings for {len(self.pfam2go)} Pfam families")
            return True
        except Exception as e:
            print(f"  Error: {e}")
            return False

    def map_pfam_to_functional_category(
        self, pfam_id: str
    ) -> Tuple[List[Tuple[str, str]], List[Dict[str, str]]]:
        """Return ([(category, subcategory), ...], go_terms) for a Pfam ID."""
        pfam_clean = pfam_id.split(".")[0] if "." in pfam_id else pfam_id
        if pfam_clean not in self.pfam2go:
            return [("unknown", "")], []
        go_terms = self.pfam2go[pfam_clean]
        found: Dict[str, str] = {}  # cat → best subcategory seen
        for go_term in go_terms:
            for cat, go_list in self.functional_categories.items():
                if go_term["go_id"] in go_list:
                    if cat not in found:
                        found[cat] = ""
            for cat, subcat in self._search_patterns_in_text(go_term["go_name"]):
                if cat not in found or (not found[cat] and subcat):
                    found[cat] = subcat
        result = [(cat, subcat) for cat, subcat in found.items()]
        return (result, go_terms) if result else ([("other", "")], go_terms)


# ---------------------------------------------------------------------------
# Main annotator
# ---------------------------------------------------------------------------

class MetaFuncDecoder:
    """Functional annotation interpreter for annotated (meta)genomes."""

    def __init__(self, categories: Dict[str, Any]):
        self.data: Dict[str, pd.DataFrame] = {}
        self.annotations: Dict[str, pd.DataFrame] = {}
        self.results: Dict[str, Dict] = {}
        self._categories = categories

        # KEGG mapper
        kegg_patterns = {
            cat: _parse_patterns(cfg["kegg_patterns"]) for cat, cfg in categories.items()
        }
        self.kegg_mapper = KEGGMapper(kegg_patterns)

        # Pfam mapper
        pfam_patterns = {
            cat: _parse_patterns(cfg["pfam_patterns"]) for cat, cfg in categories.items()
        }
        go_terms = {
            cat: cfg["pfam_go_terms"] for cat, cfg in categories.items()
        }
        self.pfam_mapper = Pfam2GOMapper(pfam_patterns, go_terms)

        # CAZy mapper: flatten subcategories from all broad categories
        cazy_patterns: Dict[str, list] = {}
        subcategory_to_broad: Dict[str, str] = {}
        family_to_subcategory: Dict[str, List[str]] = {}
        for broad_cat, cfg in categories.items():
            for subcat, subdef in (cfg.get("cazy_subcategories") or {}).items():
                subcategory_to_broad[subcat] = broad_cat
                if subdef and "patterns" in subdef:
                    cazy_patterns[subcat] = _parse_patterns(subdef["patterns"])
                # Option C: explicit family ID list → subcategory list (authoritative).
                # Appending allows one family to appear in multiple subcategories
                # (e.g. GH16 in both beta_glucanase and algal_polysaccharide).
                for family in (subdef.get("families") or [] if subdef else []):
                    family_to_subcategory.setdefault(family, []).append(subcat)
        self.cazy_mapper = CAZyMapper(cazy_patterns, subcategory_to_broad, family_to_subcategory)

    # --- database setup ---------------------------------------------------

    def setup_kegg_info(self, path: str = "data/ko_list") -> bool:
        return self.kegg_mapper.load_kegg_info(path)

    def setup_cazy_info(self, path: str = "data/FamInfo.txt.04232020.tsv") -> bool:
        return self.cazy_mapper.load_cazy_info(path)

    def setup_pfam2go(self, path: str = "data/pfam2go.txt") -> bool:
        return self.pfam_mapper.load_pfam2go(path)

    # --- input: join mode -------------------------------------------------

    def load_join_mode(
        self,
        ko_file: Optional[str] = None,
        pfam_file: Optional[str] = None,
        cazy_file: Optional[str] = None,
        dataset_name: str = "metagenome",
    ) -> None:
        """Load separate annotation files and join on gene_id."""
        print("\nLoading annotation files (join mode)...")
        frames = []

        if ko_file:
            ko_df = self._load_ko_file(ko_file)
            frames.append(ko_df)
            print(f"  KO: {len(ko_df)} genes with threshold-passing hits")

        if pfam_file:
            pfam_df = self._load_pfam_file(pfam_file)
            frames.append(pfam_df)
            print(f"  Pfam: {len(pfam_df)} genes annotated")

        if cazy_file:
            cazy_df = self._load_cazy_file(cazy_file)
            frames.append(cazy_df)
            print(f"  CAZy: {len(cazy_df)} genes annotated")

        if not frames:
            raise ValueError("No annotation files provided.")

        result = frames[0]
        for df in frames[1:]:
            result = pd.merge(result, df, on="gene_id", how="outer")

        print(f"  Combined: {len(result)} unique genes")
        self.data[dataset_name] = result

    def _load_ko_file(self, filepath: str) -> pd.DataFrame:
        """Load JGI KofamScan .assembled.KO format (tab-separated, no header).

        Columns: gene_id, threshold (Yes/No), ko_term, identity,
                 q_start, q_end, s_start, s_end, evalue, bitscore, aln_len
        """
        df = pd.read_csv(
            filepath, sep="\t", header=None, low_memory=False,
            names=[
                "gene_id", "threshold", "kegg_ortholog", "identity",
                "q_start", "q_end", "s_start", "s_end",
                "evalue", "bitscore", "aln_len",
            ],
        )
        df = df[df["threshold"] == "Yes"]
        ko_df = (
            df.groupby("gene_id")["kegg_ortholog"]
            .apply(lambda x: ";".join(x.unique()))
            .reset_index()
        )
        return ko_df

    def _load_pfam_file(self, filepath: str) -> pd.DataFrame:
        """Load a Pfam annotation file.

        Supports two formats, auto-detected from file content:

        1. JGI blout format (tab-separated, no header, no comment lines):
           gene_id  pfam_id  identity  q_start  q_end  s_start  s_end
           evalue  bitscore  aln_len
           Pfam IDs use lowercase prefix (e.g. pfam03891) and are normalised
           to the PF-prefixed format used by Pfam2GO (e.g. PF03891).

        2. HMMER domtblout format (space-separated, comment lines start with #):
           target_acc (col 1) = Pfam accession (e.g. PF00001.23)
           query_name (col 3) = gene ID
        """
        # Detect format from first non-empty line.
        # HMMER domtblout: comment lines start with '#'; data lines have ≥22
        # space-separated fields (target_name acc query_name ...).
        # JGI blout: tab-separated, 2+ fields, no comment lines.
        # Detection uses field count on the first data line so that domtblout
        # files stripped of their header comments are still recognised correctly.
        fmt = "jgi"
        with open(filepath, "r") as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                fmt = "hmmer" if len(stripped.split()) >= 20 else "jgi"
                break

        rows = []
        with open(filepath, "r") as f:
            for line in f:
                if not line.strip():
                    continue
                if fmt == "hmmer":
                    if line.startswith("#"):
                        continue
                    parts = line.split()
                    if len(parts) < 4:
                        continue
                    pfam_id = self._normalise_pfam_id(parts[1])
                    rows.append({"gene_id": parts[3], "pfam_id": pfam_id})
                else:  # jgi blout
                    parts = line.split("\t")
                    if len(parts) < 2:
                        continue
                    pfam_id = self._normalise_pfam_id(parts[1].strip())
                    rows.append({"gene_id": parts[0].strip(), "pfam_id": pfam_id})

        if not rows:
            return pd.DataFrame(columns=["gene_id", "pfam"])

        df = pd.DataFrame(rows)
        pfam_df = (
            df.groupby("gene_id")["pfam_id"]
            .apply(lambda x: ";".join(x.unique()))
            .reset_index()
            .rename(columns={"pfam_id": "pfam"})
        )
        return pfam_df

    @staticmethod
    def _normalise_pfam_id(raw: str) -> str:
        """Normalise Pfam ID to PF-prefixed format used by Pfam2GO.

        Handles:  pfam03891 -> PF03891
                  PF03891.15 -> PF03891   (strip version suffix)
                  PF03891    -> PF03891
        """
        raw = raw.strip()
        # Strip version suffix (e.g. PF03891.15 -> PF03891)
        raw = raw.split(".")[0]
        # Normalise prefix
        if raw.lower().startswith("pfam"):
            return "PF" + raw[4:]
        return raw

    def _load_cazy_file(self, filepath: str) -> pd.DataFrame:
        """Load dbCAN overview.txt format.

        Columns: Gene_ID  EC#  HMMER  dbCAN_sub  DIAMOND  Signalp  #ofTools
        HMMER column contains the CAZy family (e.g. GH18(19-390)) or '-'.
        """
        df = pd.read_csv(filepath, sep="\t")
        df.columns = [c.strip() for c in df.columns]

        if "Gene_ID" in df.columns:
            df = df.rename(columns={"Gene_ID": "gene_id"})

        # Explicit column recognition — prefer HMMER, accept dbCAN, else fail clearly.
        if "HMMER" in df.columns:
            hmmer_col = "HMMER"
        elif "dbCAN" in df.columns:
            hmmer_col = "dbCAN"
        else:
            print(
                f"  Error: CAZy file has no 'HMMER' or 'dbCAN' column.\n"
                f"  Available columns: {', '.join(df.columns.tolist())}"
            )
            sys.exit(1)

        df["cazyme"] = df[hmmer_col].apply(
            lambda x: re.sub(r"\(.*?\)", "", str(x)).replace("+", ";").strip()
            if str(x) not in ("-", "nan") else ""
        )
        df = df[df["cazyme"] != ""].copy()

        df["cazy_ec"] = (
            df["EC#"].apply(lambda x: "" if str(x) in ("-", "nan", "") else str(x).strip())
            if "EC#" in df.columns else ""
        )
        df["cazy_tools"] = (
            df["#ofTools"].astype(str) if "#ofTools" in df.columns else ""
        )

        # Parse dbCAN_sub column — subfamily IDs like GH5_7(19-390) or GH5_7|GH5_40 or -.
        has_sub = "dbCAN_sub" in df.columns

        def _parse_subfamilies(raw) -> list:
            """Return list of valid subfamily IDs from a dbCAN_sub cell value."""
            s = str(raw) if pd.notna(raw) else ""
            if not s or s == "-":
                return []
            ids = []
            for entry in re.split(r"[|+,]", s):
                clean = re.sub(r"\(.*?\)", "", entry).strip()
                if re.match(r"[A-Z]+\d+_\w+", clean):
                    ids.append(clean)
            return ids

        # Expand multi-family rows to one row per family so that EC, tool count,
        # and subfamily stay aligned with their source family after groupby aggregation.
        expanded_rows = []
        for _, row in df.iterrows():
            families = str(row["cazyme"]).split(";")
            ec = row.get("cazy_ec", "")
            tools = row.get("cazy_tools", "")
            row_subfamilies = _parse_subfamilies(row.get("dbCAN_sub", "")) if has_sub else []
            for fam in families:
                fam = fam.strip()
                if fam:
                    # Match subfamilies to this family by prefix (GH5 matches GH5_7, GH5_21)
                    matched_subs = [s for s in row_subfamilies if s.startswith(fam + "_")]
                    expanded_rows.append({
                        "gene_id": row["gene_id"],
                        "cazyme": fam,
                        "cazy_subfamilies": ";".join(matched_subs),
                        "cazy_ec": ec,
                        "cazy_tools": tools,
                    })
        if not expanded_rows:
            return pd.DataFrame(
                columns=["gene_id", "cazyme", "cazy_subfamilies", "cazy_ec", "cazy_tools"]
            )
        df_expanded = pd.DataFrame(expanded_rows)
        df_expanded = df_expanded.drop_duplicates(subset=["gene_id", "cazyme"])
        cazy_df = (
            df_expanded.groupby("gene_id")
            .agg(
                cazyme=("cazyme", lambda x: ";".join(x)),
                cazy_subfamilies=("cazy_subfamilies", lambda x: ";".join(v for v in x if v)),
                cazy_ec=("cazy_ec", lambda x: ";".join(x)),
                cazy_tools=("cazy_tools", lambda x: ";".join(x)),
            )
            .reset_index()
        )
        return cazy_df

    # --- input: combined mode ---------------------------------------------

    def load_combined_mode(
        self,
        filepath: str,
        ko_col: str = "kegg_ortholog",
        pfam_col: str = "pfam",
        cazy_col: str = "cazyme",
        gene_id_col: str = "gene_id",
        allow_missing_cols: bool = False,
        dataset_name: str = "metagenome",
    ) -> None:
        """Load a single pre-merged annotation table with configurable column names."""
        print(f"\nLoading combined annotation table: {filepath}")
        if filepath.endswith(".gz"):
            df = pd.read_csv(filepath, sep="\t", compression="gzip", low_memory=False)
        else:
            df = pd.read_csv(filepath, sep="\t", low_memory=False)
        print(f"  Loaded {len(df)} rows, {len(df.columns)} columns")

        rename = {}
        for src, tgt in [
            (gene_id_col, "gene_id"),
            (ko_col, "kegg_ortholog"),
            (pfam_col, "pfam"),
            (cazy_col, "cazyme"),
        ]:
            if src != tgt and src in df.columns:
                rename[src] = tgt
        if rename:
            df = df.rename(columns=rename)

        for col in ["kegg_ortholog", "pfam", "cazyme"]:
            if col not in df.columns:
                if allow_missing_cols:
                    print(f"  Warning: '{col}' column not found — that database will be skipped")
                else:
                    print(
                        f"  Error: '{col}' column not found in {filepath}.\n"
                        f"  Available columns: {', '.join(df.columns.tolist())}\n"
                        f"  Use --allow-missing-cols to skip missing columns instead of exiting."
                    )
                    sys.exit(1)

        # Validate gene ID column exists after renaming.
        if "gene_id" not in df.columns:
            print(
                f"  Error: gene ID column '{gene_id_col}' not found in {filepath}.\n"
                f"  Available columns: {', '.join(df.columns.tolist())}\n"
                f"  Use --gene-id-col to specify the correct column name."
            )
            sys.exit(1)

        # Aggregate duplicate gene_id rows (warn, then semicolon-join annotation columns).
        if df.duplicated(subset=["gene_id"]).any():
            n_dupes = df.duplicated(subset=["gene_id"]).sum()
            print(
                f"  Warning: {n_dupes} duplicate gene_id rows detected — "
                "aggregating annotation columns by gene_id"
            )
            agg = {}
            for col in ["kegg_ortholog", "pfam", "cazyme"]:
                if col in df.columns:
                    agg[col] = (col, lambda x: ";".join(
                        v for v in x.dropna().astype(str).unique() if v
                    ))
            if agg:
                df = df.groupby("gene_id").agg(**agg).reset_index()
            else:
                df = df.drop_duplicates(subset=["gene_id"])

        self.data[dataset_name] = df
        print(f"  Ready: {len(df)} genes")

    # --- annotation standardization ---------------------------------------

    def standardize_annotations(self, min_cazy_tools: int = 1) -> None:
        """Generate confidence annotations for all loaded datasets."""
        print("\nGenerating confidence annotations...")
        for dataset_name, df in self.data.items():
            print(f"  Processing {dataset_name}...")
            self.annotations[f"{dataset_name}_confidence"] = (
                self._generate_confidence_annotations(df, min_cazy_tools=min_cazy_tools)
            )
            print(
                f"    {len(self.annotations[dataset_name + '_confidence'])} confidence records"
            )

    def _generate_confidence_annotations(
        self, df: pd.DataFrame, min_cazy_tools: int = 1
    ) -> pd.DataFrame:
        """One record per gene per broad category, scored by database agreement."""
        records = []
        for idx, row in df.iterrows():
            base = {"gene_id": str(row.get("gene_id", ""))}

            evidence: Dict[str, Dict] = defaultdict(lambda: {
                "cazy_families": [],
                "cazy_subfamilies": [],
                "subcategories": set(),
                "cazy_subcategories": set(),
                "pfam_subcategories": set(),
                "kegg_subcategories": set(),
                "cazy_description": [],
                "cazy_ec": [],
                "cazy_tool_count": [],
                "pfam_ids": [],
                "pfam_go_names": [],
                "ko_terms": [],
                "ko_definition": [],
            })

            cazyme = row.get("cazyme", "")
            if pd.notna(cazyme) and cazyme:
                cazy_ec_list = str(row.get("cazy_ec", "") or "").split(";")
                cazy_tools_list = str(row.get("cazy_tools", "") or "").split(";")
                cazy_sub_list = str(row.get("cazy_subfamilies", "") or "").split(";")
                for i, ann in enumerate(str(cazyme).split(";")):
                    ann = ann.strip()
                    if not ann:
                        continue
                    ec = cazy_ec_list[i].strip() if i < len(cazy_ec_list) else ""
                    tools = cazy_tools_list[i].strip() if i < len(cazy_tools_list) else ""
                    subfamily = cazy_sub_list[i].strip() if i < len(cazy_sub_list) else ""
                    # Apply --cazy-min-tools filter: skip families below the threshold.
                    # Uses int(float(x)) so both "3" and "3.0" (common in pandas/R
                    # exports) are handled correctly. Genuinely non-numeric values
                    # (empty string, "-", "NA") are passed through with a warning
                    # rather than silently dropped — missing count ≠ failed threshold.
                    if min_cazy_tools > 1 and tools:
                        try:
                            if int(float(tools)) < min_cazy_tools:
                                continue
                        except ValueError:
                            print(
                                f"  Warning: non-numeric cazy_tools value '{tools}' "
                                f"for family '{ann}' — cannot apply --cazy-min-tools "
                                f"filter; hit passed through."
                            )
                    # Use subfamily for lookup when available — more specific than family.
                    # Fall back to family ID if no subfamily was assigned by dbCAN_sub.
                    lookup_id = subfamily if subfamily else ann
                    subcats = self.cazy_mapper.get_family_categories(lookup_id)
                    family_info = self.cazy_mapper.get_detailed_family_info(ann)
                    activities = family_info.get("activities", "") if family_info else ""
                    for subcat in subcats:
                        broad = self.cazy_mapper.get_broad_category(subcat)
                        if broad != "other":
                            evidence[broad]["cazy_families"].append(ann)
                            if subfamily:
                                evidence[broad]["cazy_subfamilies"].append(subfamily)
                            evidence[broad]["subcategories"].add(subcat)
                            evidence[broad]["cazy_subcategories"].add(subcat)
                            if activities:
                                evidence[broad]["cazy_description"].append(activities)
                            if ec:
                                evidence[broad]["cazy_ec"].append(ec)
                            if tools:
                                evidence[broad]["cazy_tool_count"].append(tools)

            pfam = row.get("pfam", "")
            if pd.notna(pfam) and pfam:
                for ann in str(pfam).split(";"):
                    ann = ann.strip()
                    if not ann:
                        continue
                    cats_list, go_terms = self.pfam_mapper.map_pfam_to_functional_category(ann)
                    for cat, subcat in cats_list:
                        if cat not in ("other", "unknown"):
                            evidence[cat]["pfam_ids"].append(ann)
                            if subcat:
                                evidence[cat]["subcategories"].add(subcat)
                                evidence[cat]["pfam_subcategories"].add(subcat)
                            evidence[cat]["pfam_go_names"].extend(
                                g["go_name"] for g in go_terms
                            )

            kegg = row.get("kegg_ortholog", "")
            if pd.notna(kegg) and kegg:
                for ann in str(kegg).split(";"):
                    ann = ann.strip()
                    if not ann:
                        continue
                    hits = self.kegg_mapper.get_kegg_categories(ann)
                    kegg_info = self.kegg_mapper.get_detailed_kegg_info(ann)
                    definition = kegg_info.get("definition", "") if kegg_info else ""
                    for cat, subcat in hits:
                        if cat != "other":
                            evidence[cat]["ko_terms"].append(ann)
                            if subcat:
                                evidence[cat]["subcategories"].add(subcat)
                                evidence[cat]["kegg_subcategories"].add(subcat)
                            if definition:
                                evidence[cat]["ko_definition"].append(definition)

            if not evidence:
                rec = base.copy()
                rec.update({
                    "functional_category": "unannotated",
                    "subcategories": "",
                    "evidence_specificity": "",
                    "cazy_families": "",
                    "cazy_subfamilies": "",
                    "cazy_ec": "",
                    "cazy_tool_count": "",
                    "pfam_ids": "",
                    "ko_terms": "",
                    "cazy_description": "",
                    "pfam_go_names": "",
                    "ko_definition": "",
                    "supporting_databases": "",
                    "n_databases": 0,
                    "confidence_level": "none",
                    "confidence_score": 0.0,
                })
                records.append(rec)
            else:
                for broad_cat, ev in evidence.items():
                    dbs = []
                    if ev["cazy_families"]:
                        dbs.append("CAZy")
                    if ev["pfam_ids"]:
                        dbs.append("Pfam")
                    if ev["ko_terms"]:
                        dbs.append("KEGG")
                    db_count = len(dbs)
                    confidence_level = (
                        "high" if db_count >= 3 else
                        "medium" if db_count == 2 else
                        "low" if db_count == 1 else
                        "none"
                    )
                    # Count databases that agree on each subcategory
                    subcat_db_counts = {}
                    for subcat in ev["subcategories"]:
                        count = sum([
                            subcat in ev["cazy_subcategories"],
                            subcat in ev["pfam_subcategories"],
                            subcat in ev["kegg_subcategories"],
                        ])
                        subcat_db_counts[subcat] = count
                    subcat_n_dbs = max(subcat_db_counts.values()) if subcat_db_counts else 0
                    rec = base.copy()
                    rec.update({
                        "functional_category": broad_cat,
                        "subcategories": ";".join(sorted(ev["subcategories"])),
                        "evidence_specificity": (
                            "specific" if ev["subcategories"] else "broad"
                        ),
                        "subcategory_n_databases": subcat_n_dbs,
                        "cazy_families": ";".join(dict.fromkeys(ev["cazy_families"])),
                        "cazy_subfamilies": ";".join(dict.fromkeys(ev["cazy_subfamilies"])),
                        "cazy_ec": ";".join(ev["cazy_ec"]),
                        "cazy_tool_count": ";".join(ev["cazy_tool_count"]),
                        "pfam_ids": ";".join(ev["pfam_ids"]),
                        "ko_terms": ";".join(ev["ko_terms"]),
                        "cazy_description": " | ".join(ev["cazy_description"]),
                        "pfam_go_names": " | ".join(ev["pfam_go_names"]),
                        "ko_definition": " | ".join(ev["ko_definition"]),
                        "supporting_databases": ";".join(dbs),
                        "n_databases": db_count,
                        "confidence_level": confidence_level,
                        "confidence_score": {3: 1.0, 2: 0.67, 1: 0.33, 0: 0.0}[db_count],
                    })
                    records.append(rec)

        return pd.DataFrame(records)

    # --- analysis ---------------------------------------------------------

    def partition_function(self, functional_category: str) -> Dict:
        print(f"\nPartitioning: {functional_category}")
        partition_results: Dict[str, int] = defaultdict(int)
        total = 0

        for dataset_name, annotations in self.annotations.items():
            if "_confidence" not in dataset_name:
                continue
            base = dataset_name.replace("_confidence", "")
            count = len(
                annotations[
                    (annotations["functional_category"] == functional_category) &
                    (annotations["confidence_level"] != "none")
                ]
            )
            partition_results[base] = count
            total += count
            print(f"  {base}: {count} genes")

        results = {
            "function": functional_category,
            "counts": dict(partition_results),
            "percentages": (
                {k: v / total * 100 for k, v in partition_results.items()}
                if total > 0 else {}
            ),
            "total": total,
        }
        self.results[functional_category] = results
        return results

    def generate_summary_report(self) -> None:
        print("\n" + "=" * 70)
        print("FUNCTIONAL ANNOTATION SUMMARY REPORT")
        print("=" * 70)
        for result in self.results.values():
            print(f"\nFunction: {result['function']}")
            print(f"Total annotations: {result['total']}")
            for dataset, pct in result["percentages"].items():
                count = result["counts"][dataset]
                print(f"  {dataset}: {pct:.1f}% ({count} annotations)")

    def write_run_manifest(
        self,
        output_dir: str,
        args_namespace,
        db_paths: Dict[str, str],
        categories_dir: str,
    ) -> None:
        """Write a JSON sidecar file recording run provenance."""
        def sha256(path: str) -> str:
            try:
                h = hashlib.sha256()
                with open(path, "rb") as f:
                    for chunk in iter(lambda: f.read(65536), b""):
                        h.update(chunk)
                return h.hexdigest()
            except Exception:
                return "unavailable"

        db_hashes = {name: {"path": path, "sha256": sha256(path)}
                     for name, path in db_paths.items()}

        yaml_hashes = {}
        for yaml_file in sorted(Path(categories_dir).glob("*.yaml")):
            yaml_hashes[yaml_file.name] = {
                "path": str(yaml_file),
                "sha256": sha256(str(yaml_file)),
            }

        manifest = {
            "metafuncdecoder_version": VERSION,
            "run_date": datetime.datetime.now().isoformat(timespec="seconds"),
            "command_line": sys.argv,
            "databases": db_hashes,
            "category_files": yaml_hashes,
            "ontology_summary": {
                cat: {
                    "kegg_patterns": len(cfg.get("kegg_patterns", [])),
                    "pfam_patterns": len(cfg.get("pfam_patterns", [])),
                    "cazy_subcategories": len(cfg.get("cazy_subcategories", {})),
                }
                for cat, cfg in self._categories.items()
                if not cfg.get("_custom", False)
            },
        }

        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        manifest_path = out / "run_manifest.json"
        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=2)
        print(f"Run manifest: {manifest_path}")

    def save_results(
        self,
        output_dir: str,
        include_unannotated: bool = False,
        functions: Optional[List[str]] = None,
        exclude_broad: bool = False,
    ) -> None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)

        for name, annotations in self.annotations.items():
            if "_confidence" not in name:
                continue
            base = name.replace("_confidence", "")
            if functions is not None:
                keep = set(functions)
                if include_unannotated:
                    keep.add("unannotated")
                annotations = annotations[
                    annotations["functional_category"].isin(keep)
                ]
            elif not include_unannotated:
                annotations = annotations[annotations["confidence_level"] != "none"]
            if exclude_broad and "evidence_specificity" in annotations.columns:
                before = len(annotations)
                annotations = annotations[annotations["evidence_specificity"] != "broad"]
                dropped = before - len(annotations)
                if dropped:
                    print(f"  --exclude-broad: removed {dropped} broad-evidence rows from {base}")
            outfile = out / f"{base}_confidence_annotations.csv"
            annotations.to_csv(outfile, index=False)
            print(f"{base}: {outfile}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "MetaFuncDecoder: Interpret functional annotations "
            "from annotated (meta)genomes"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:

  # Join mode: separate files from JGI tarball + dbCAN
  python metafuncdecoder.py \\
    --ko-file 80425.assembled.ko \\
    --cazy-file dbcan_overview.txt \\
    --function carbon_cycling

  # Join mode with Pfam and CAZy
  python metafuncdecoder.py \\
    --ko-file 80425.assembled.ko \\
    --pfam-file 80425.assembled.pfam.blout \\
    --cazy-file dbcan_overview.txt \\
    --function carbon_cycling \\
    --output results/

  # Combined mode: single pre-merged table
  python metafuncdecoder.py \\
    --combined-table annotations.tsv \\
    --ko-col ko_id --pfam-col pfam_hits --cazy-col cazy_id \\
    --function nitrogen_cycling

  # Custom categories directory
  python metafuncdecoder.py \\
    --ko-file 80425.assembled.ko \\
    --categories-dir my_categories/
        """,
    )

    # Join mode
    join = parser.add_argument_group("Join mode (separate annotation files)")
    join.add_argument("--ko-file", help="JGI KofamScan .assembled.KO file")
    join.add_argument("--pfam-file", help="HMMER domtblout Pfam annotation file")
    join.add_argument("--cazy-file", help="dbCAN overview.txt CAZy annotation file")

    # Combined mode
    combined = parser.add_argument_group("Combined mode (single pre-merged table)")
    combined.add_argument("--combined-table", help="Pre-merged annotation table (TSV)")
    combined.add_argument(
        "--ko-col", default="kegg_ortholog",
        help="Column name for KO terms (default: kegg_ortholog)",
    )
    combined.add_argument(
        "--pfam-col", default="pfam",
        help="Column name for Pfam IDs (default: pfam)",
    )
    combined.add_argument(
        "--cazy-col", default="cazyme",
        help="Column name for CAZy families (default: cazyme)",
    )
    combined.add_argument(
        "--gene-id-col", default="gene_id",
        help="Column name for gene IDs (default: gene_id)",
    )

    # Reference databases
    dbs = parser.add_argument_group("Reference databases")
    dbs.add_argument(
        "--ko-db", default=str(_SCRIPT_DIR / "data" / "ko_list"),
        help="Path to KOfam ko_list (default: <script_dir>/data/ko_list)",
    )
    dbs.add_argument(
        "--pfam-db", default=str(_SCRIPT_DIR / "data" / "pfam2go.txt"),
        help="Path to pfam2go.txt (default: <script_dir>/data/pfam2go.txt)",
    )
    dbs.add_argument(
        "--cazy-db", default=str(_SCRIPT_DIR / "data" / "FamInfo.txt.04232020.tsv"),
        help="Path to CAZy FamInfo TSV file (default: <script_dir>/data/FamInfo.txt.04232020.tsv)",
    )

    # Categories and analysis
    parser.add_argument(
        "--categories-dir", default=str(_SCRIPT_DIR / "categories"),
        help="Directory with YAML category files (default: <script_dir>/categories/)",
    )
    parser.add_argument(
        "--function", nargs="+", default=["carbon_cycling"],
        metavar="FUNCTION",
        help=(
            "Functional category/categories to analyze. "
            "Pass multiple names to run several at once, or 'all' to run every "
            "category found in --categories-dir. (default: carbon_cycling)"
        ),
    )
    # Custom patterns
    custom = parser.add_argument_group("Custom pattern search")
    custom.add_argument(
        "--custom-ko-pattern", help="Regex to search in KEGG definitions",
    )
    custom.add_argument(
        "--custom-pfam-pattern", help="Regex to search in Pfam GO term names",
    )
    custom.add_argument(
        "--custom-cazy-pattern", help="Regex to search in CAZy activity descriptions",
    )
    custom.add_argument(
        "--custom-category-name", default="custom",
        help="Name for custom category (default: custom)",
    )

    # Output
    parser.add_argument(
        "--dataset-name", default="metagenome",
        help="Dataset name used in output filename: <name>_confidence_annotations.csv (default: metagenome)",
    )
    parser.add_argument(
        "--output", default="results/",
        help="Output directory (default: results/)",
    )
    parser.add_argument(
        "--include-unannotated", action="store_true",
        help="Include genes with no functional category match in the output (default: omit)",
    )
    parser.add_argument(
        "--cazy-min-tools", type=int, default=1, metavar="N",
        help=(
            "Minimum number of dbCAN tools that must agree on a CAZy family assignment "
            "for it to contribute to the confidence model (default: 1, i.e. no filter). "
            "dbCAN recommends 2. Only applies when --cazy-file is used or when the "
            "combined table contains a cazy_tools column."
        ),
    )
    parser.add_argument(
        "--exclude-broad", action="store_true",
        help=(
            "Exclude broad-evidence annotations from the output. Broad hits are genes "
            "matched only by non-specific patterns (e.g. GO:0005975 carbohydrate metabolic "
            "process, or text patterns such as 'antibiotic') that produce no subcategory "
            "label. These correspond to evidence_specificity='broad' in the output. "
            "Equivalent to post-filtering on evidence_specificity != 'broad'."
        ),
    )
    parser.add_argument(
        "--allow-missing-cols", action="store_true",
        help=(
            "Combined mode only: skip missing annotation columns with a warning instead of "
            "exiting with an error (default: error on missing columns)"
        ),
    )

    args = parser.parse_args()

    # Validate input mode
    join_inputs = [args.ko_file, args.pfam_file, args.cazy_file]
    if not any(join_inputs) and not args.combined_table:
        print(
            "Error: Provide join mode inputs (--ko-file / --pfam-file / --cazy-file) "
            "or --combined-table"
        )
        parser.print_help()
        sys.exit(1)
    if any(join_inputs) and args.combined_table:
        print("Error: Use either join mode or --combined-table, not both")
        sys.exit(1)

    # Load categories
    categories = load_categories(args.categories_dir)
    if not categories:
        print("Error: No functional categories loaded. Check --categories-dir")
        sys.exit(1)

    # Initialize
    decoder = MetaFuncDecoder(categories)
    kegg_ok = decoder.setup_kegg_info(args.ko_db)
    cazy_ok = decoder.setup_cazy_info(args.cazy_db)
    pfam_ok = decoder.setup_pfam2go(args.pfam_db)

    if not kegg_ok:
        print(
            f"WARNING: KEGG database not loaded from {args.ko_db!r}. "
            "KEGG annotations will be skipped. Run 'bash download_dbs.sh' to fetch databases."
        )
    if not cazy_ok:
        print(
            f"WARNING: CAZy database not loaded from {args.cazy_db!r}. "
            "CAZy annotations will be skipped."
        )
    if not pfam_ok:
        print(
            f"WARNING: Pfam2GO database not loaded from {args.pfam_db!r}. "
            "Pfam annotations will be skipped."
        )
    if not any([kegg_ok, cazy_ok, pfam_ok]):
        print("ERROR: No reference databases loaded. Cannot produce any annotations.")
        sys.exit(1)

    # Custom patterns
    if args.custom_ko_pattern:
        decoder.kegg_mapper.add_custom_patterns(
            {args.custom_category_name: [args.custom_ko_pattern]}
        )
    if args.custom_pfam_pattern:
        decoder.pfam_mapper.add_custom_patterns(
            {args.custom_category_name: [args.custom_pfam_pattern]}
        )
    if args.custom_cazy_pattern:
        decoder.cazy_mapper.add_custom_patterns(
            {args.custom_category_name: [args.custom_cazy_pattern]}
        )
    if any([args.custom_ko_pattern, args.custom_pfam_pattern, args.custom_cazy_pattern]):
        decoder._categories.setdefault(args.custom_category_name, {})
        if args.custom_category_name not in args.function and args.function != ["all"]:
            args.function.append(args.custom_category_name)

    # Load data
    if args.combined_table:
        decoder.load_combined_mode(
            args.combined_table,
            ko_col=args.ko_col,
            pfam_col=args.pfam_col,
            cazy_col=args.cazy_col,
            gene_id_col=args.gene_id_col,
            allow_missing_cols=args.allow_missing_cols,
            dataset_name=args.dataset_name,
        )
    else:
        decoder.load_join_mode(
            ko_file=args.ko_file,
            pfam_file=args.pfam_file,
            cazy_file=args.cazy_file,
            dataset_name=args.dataset_name,
        )

    # Analyse
    decoder.standardize_annotations(min_cazy_tools=args.cazy_min_tools)

    # Resolve --function all → every loaded category
    if args.function == ["all"]:
        functions = list(decoder._categories.keys())
    else:
        functions = args.function
        unknown = [f for f in functions if f not in decoder._categories]
        if unknown:
            print(f"Warning: unknown function(s): {', '.join(unknown)}")
            print(f"  Available: {', '.join(decoder._categories.keys())}")
            functions = [f for f in functions if f in decoder._categories]

    for func in functions:
        decoder.partition_function(func)

    decoder.write_run_manifest(
        args.output,
        args,
        db_paths={
            "ko_list": args.ko_db,
            "pfam2go": args.pfam_db,
            "cazy_faminfo": args.cazy_db,
        },
        categories_dir=args.categories_dir,
    )
    decoder.generate_summary_report()
    decoder.save_results(
        args.output,
        include_unannotated=args.include_unannotated,
        functions=functions,
        exclude_broad=args.exclude_broad,
    )


if __name__ == "__main__":
    main()
