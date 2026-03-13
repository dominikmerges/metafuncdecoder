"""
Unit tests for MetaFuncDecoder core logic.

Run with:
    pytest tests/test_core.py -v
"""

import sys
import tempfile
import textwrap
from pathlib import Path

import pytest

# Make the parent directory importable regardless of where pytest is run from
sys.path.insert(0, str(Path(__file__).parent.parent))

from metafuncdecoder import (
    MetaFuncDecoder,
    CAZyMapper,
    KEGGMapper,
    Pfam2GOMapper,
    load_categories,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _tmp_file(content: str, suffix: str = ".txt") -> Path:
    """Write content to a named temp file and return its Path."""
    f = tempfile.NamedTemporaryFile(
        mode="w", suffix=suffix, delete=False, encoding="utf-8"
    )
    f.write(textwrap.dedent(content))
    f.flush()
    return Path(f.name)


def _minimal_categories():
    """Minimal categories dict for constructing a MetaFuncDecoder without YAML."""
    return {
        "carbon_cycling": {
            "kegg_patterns": ["cellulase", "chitinase"],
            "pfam_patterns": ["carbohydrate"],
            "pfam_go_terms": ["GO:0005975"],
            "cazy_subcategories": {
                "cellulase": {"patterns": ["cellulase", "cellulose degradation"]},
                "chitinase": {"patterns": ["chitinase", "chitin"]},
            },
        },
        "nitrogen_cycling": {
            "kegg_patterns": ["nitrogenase", "nitrate reductase"],
            "pfam_patterns": ["nitrogen"],
            "pfam_go_terms": [],
            "cazy_subcategories": {},
        },
    }


# ---------------------------------------------------------------------------
# 1. load_categories
# ---------------------------------------------------------------------------

class TestLoadCategories:

    def test_loads_valid_yaml(self, tmp_path):
        yaml_content = """\
            name: test_func
            description: Test functional category
            citations:
              - Author et al. 2020
            kegg_patterns:
              - cellulase
              - chitinase
            pfam_patterns:
              - carbohydrate
            pfam_go_terms:
              - GO:0005975
            cazy_subcategories:
              cellulase:
                patterns:
                  - cellulase
        """
        (tmp_path / "test_func.yaml").write_text(textwrap.dedent(yaml_content))
        cats = load_categories(str(tmp_path))
        assert "test_func" in cats
        assert cats["test_func"]["kegg_patterns"] == ["cellulase", "chitinase"]
        assert cats["test_func"]["pfam_go_terms"] == ["GO:0005975"]
        assert "cellulase" in cats["test_func"]["cazy_subcategories"]

    def test_missing_name_field_skipped(self, tmp_path):
        (tmp_path / "bad.yaml").write_text("description: no name here\n")
        cats = load_categories(str(tmp_path))
        assert cats == {}

    def test_nonexistent_directory_returns_empty(self, tmp_path):
        cats = load_categories(str(tmp_path / "does_not_exist"))
        assert cats == {}

    def test_empty_directory_returns_empty(self, tmp_path):
        cats = load_categories(str(tmp_path))
        assert cats == {}

    def test_optional_fields_default_to_empty(self, tmp_path):
        (tmp_path / "minimal.yaml").write_text("name: minimal\n")
        cats = load_categories(str(tmp_path))
        assert cats["minimal"]["kegg_patterns"] == []
        assert cats["minimal"]["pfam_patterns"] == []
        assert cats["minimal"]["pfam_go_terms"] == []
        assert cats["minimal"]["cazy_subcategories"] == {}

    def test_multiple_yaml_files_all_loaded(self, tmp_path):
        for name in ("alpha", "beta", "gamma"):
            (tmp_path / f"{name}.yaml").write_text(f"name: {name}\n")
        cats = load_categories(str(tmp_path))
        assert set(cats.keys()) == {"alpha", "beta", "gamma"}


# ---------------------------------------------------------------------------
# 2. _normalise_pfam_id
# ---------------------------------------------------------------------------

class TestNormalisePfamId:

    @pytest.fixture
    def decoder(self):
        return MetaFuncDecoder(_minimal_categories())

    def test_lowercase_pfam_prefix(self, decoder):
        assert decoder._normalise_pfam_id("pfam03891") == "PF03891"

    def test_uppercase_pf_prefix_unchanged(self, decoder):
        assert decoder._normalise_pfam_id("PF03891") == "PF03891"

    def test_version_suffix_stripped(self, decoder):
        assert decoder._normalise_pfam_id("PF03891.15") == "PF03891"

    def test_lowercase_pfam_with_version(self, decoder):
        assert decoder._normalise_pfam_id("pfam00001.23") == "PF00001"

    def test_leading_trailing_whitespace(self, decoder):
        assert decoder._normalise_pfam_id("  PF03891  ") == "PF03891"

    def test_already_clean(self, decoder):
        assert decoder._normalise_pfam_id("PF00001") == "PF00001"


# ---------------------------------------------------------------------------
# 3. _load_pfam_file — format auto-detection
# ---------------------------------------------------------------------------

class TestLoadPfamFile:

    @pytest.fixture
    def decoder(self):
        return MetaFuncDecoder(_minimal_categories())

    def test_jgi_blout_format(self, decoder, tmp_path):
        p = tmp_path / "test.blout"
        p.write_text(
            "gene_001\tpfam00001\t95.0\t1\t200\t1\t200\t1e-50\t300\t200\n"
            "gene_001\tpfam00002\t80.0\t1\t150\t1\t150\t1e-30\t200\t150\n"
            "gene_002\tpfam00003\t90.0\t1\t180\t1\t180\t1e-40\t250\t180\n"
        )
        df = decoder._load_pfam_file(str(p))
        assert set(df["gene_id"]) == {"gene_001", "gene_002"}
        # IDs normalised to PF prefix
        gene1_pfams = set(df.loc[df["gene_id"] == "gene_001", "pfam"].values[0].split(";"))
        assert gene1_pfams == {"PF00001", "PF00002"}

    def test_hmmer_domtblout_format(self, decoder, tmp_path):
        # Real domtblout layout (space-separated):
        #   col 0: model name  col 1: Pfam accession  col 2: tlen  col 3: gene_id ...
        p = tmp_path / "test.domtblout"
        p.write_text(
            "# comment line\n"
            "7tm_1  PF00001.23  263  gene_001  -  300  1e-50  300.0  0.0  1  1  1e-50  300  0.0  1  200  1  200  1  200  0.99\n"
            "7tm_2  PF00002.10  200  gene_002  -  250  1e-30  200.0  0.0  1  1  1e-30  200  0.0  1  150  1  150  1  150  0.95\n"
        )
        df = decoder._load_pfam_file(str(p))
        assert set(df["gene_id"]) == {"gene_001", "gene_002"}
        assert df.loc[df["gene_id"] == "gene_001", "pfam"].values[0] == "PF00001"

    def test_jgi_multiple_hits_per_gene_collapsed(self, decoder, tmp_path):
        p = tmp_path / "test.blout"
        p.write_text(
            "gene_001\tpfam00001\t95.0\t1\t200\t1\t200\t1e-50\t300\t200\n"
            "gene_001\tpfam00001\t95.0\t1\t200\t1\t200\t1e-50\t300\t200\n"  # duplicate
        )
        df = decoder._load_pfam_file(str(p))
        assert len(df) == 1
        assert df.iloc[0]["pfam"] == "PF00001"  # not doubled

    def test_empty_file_returns_empty_df(self, decoder, tmp_path):
        p = tmp_path / "empty.blout"
        p.write_text("")
        df = decoder._load_pfam_file(str(p))
        assert len(df) == 0
        assert "gene_id" in df.columns


# ---------------------------------------------------------------------------
# 4. _load_ko_file
# ---------------------------------------------------------------------------

class TestLoadKoFile:

    @pytest.fixture
    def decoder(self):
        return MetaFuncDecoder(_minimal_categories())

    def test_only_threshold_yes_rows_kept(self, decoder, tmp_path):
        p = tmp_path / "test.ko"
        p.write_text(
            "gene_001\tYes\tKO:K01234\t95.0\t1\t200\t1\t200\t1e-50\t300\t200\n"
            "gene_002\tNo\tKO:K05678\t60.0\t1\t150\t1\t150\t1e-10\t100\t150\n"
            "gene_003\tYes\tKO:K09999\t80.0\t1\t180\t1\t180\t1e-40\t250\t180\n"
        )
        df = decoder._load_ko_file(str(p))
        assert set(df["gene_id"]) == {"gene_001", "gene_003"}
        assert "gene_002" not in df["gene_id"].values

    def test_multiple_ko_per_gene_collapsed(self, decoder, tmp_path):
        p = tmp_path / "test.ko"
        p.write_text(
            "gene_001\tYes\tKO:K01234\t95.0\t1\t200\t1\t200\t1e-50\t300\t200\n"
            "gene_001\tYes\tKO:K05678\t80.0\t1\t200\t1\t200\t1e-40\t250\t200\n"
        )
        df = decoder._load_ko_file(str(p))
        assert len(df) == 1
        kos = set(df.iloc[0]["kegg_ortholog"].split(";"))
        assert kos == {"KO:K01234", "KO:K05678"}

    def test_all_no_rows_returns_empty(self, decoder, tmp_path):
        p = tmp_path / "test.ko"
        p.write_text(
            "gene_001\tNo\tKO:K01234\t60.0\t1\t150\t1\t150\t1e-10\t100\t150\n"
        )
        df = decoder._load_ko_file(str(p))
        assert len(df) == 0


# ---------------------------------------------------------------------------
# 5. _load_cazy_file
# ---------------------------------------------------------------------------

class TestLoadCazyFile:

    @pytest.fixture
    def decoder(self):
        return MetaFuncDecoder(_minimal_categories())

    def test_basic_loading(self, decoder, tmp_path):
        p = tmp_path / "overview.txt"
        p.write_text(
            "Gene_ID\tEC#\tHMMER\tdbCAN_sub\tDIAMOND\t#ofTools\n"
            "gene_001\t3.2.1.14\tGH18(1-390)\tGH18_e1\tGH18\t3\n"
            "gene_002\t-\tGT2(1-200)\tGT2_e1\tGT2\t2\n"
        )
        df = decoder._load_cazy_file(str(p))
        assert set(df["gene_id"]) == {"gene_001", "gene_002"}

    def test_hmmer_position_annotation_stripped(self, decoder, tmp_path):
        p = tmp_path / "overview.txt"
        p.write_text(
            "Gene_ID\tEC#\tHMMER\tdbCAN_sub\tDIAMOND\t#ofTools\n"
            "gene_001\t-\tGH18(19-390)\t-\tGH18\t2\n"
        )
        df = decoder._load_cazy_file(str(p))
        assert df.iloc[0]["cazyme"] == "GH18"

    def test_ec_number_captured(self, decoder, tmp_path):
        p = tmp_path / "overview.txt"
        p.write_text(
            "Gene_ID\tEC#\tHMMER\tdbCAN_sub\tDIAMOND\t#ofTools\n"
            "gene_001\t3.2.1.14\tGH18(1-390)\t-\tGH18\t3\n"
        )
        df = decoder._load_cazy_file(str(p))
        assert df.iloc[0]["cazy_ec"] == "3.2.1.14"

    def test_missing_ec_stored_as_empty(self, decoder, tmp_path):
        p = tmp_path / "overview.txt"
        p.write_text(
            "Gene_ID\tEC#\tHMMER\tdbCAN_sub\tDIAMOND\t#ofTools\n"
            "gene_001\t-\tGH18(1-390)\t-\tGH18\t2\n"
        )
        df = decoder._load_cazy_file(str(p))
        assert df.iloc[0]["cazy_ec"] == ""

    def test_tools_count_captured(self, decoder, tmp_path):
        p = tmp_path / "overview.txt"
        p.write_text(
            "Gene_ID\tEC#\tHMMER\tdbCAN_sub\tDIAMOND\t#ofTools\n"
            "gene_001\t-\tGH18(1-390)\t-\tGH18\t3\n"
        )
        df = decoder._load_cazy_file(str(p))
        assert df.iloc[0]["cazy_tools"] == "3"

    def test_dash_hmmer_rows_excluded(self, decoder, tmp_path):
        p = tmp_path / "overview.txt"
        p.write_text(
            "Gene_ID\tEC#\tHMMER\tdbCAN_sub\tDIAMOND\t#ofTools\n"
            "gene_001\t-\t-\t-\tGH18\t1\n"
            "gene_002\t-\tGH18(1-390)\t-\tGH18\t2\n"
        )
        df = decoder._load_cazy_file(str(p))
        assert set(df["gene_id"]) == {"gene_002"}

    def test_multiple_families_per_gene_collapsed(self, decoder, tmp_path):
        p = tmp_path / "overview.txt"
        p.write_text(
            "Gene_ID\tEC#\tHMMER\tdbCAN_sub\tDIAMOND\t#ofTools\n"
            "gene_001\t3.2.1.14\tGH18(1-200)\t-\tGH18\t3\n"
            "gene_001\t2.4.1.-\tGT2(201-400)\t-\tGT2\t2\n"
        )
        df = decoder._load_cazy_file(str(p))
        assert len(df) == 1
        families = df.iloc[0]["cazyme"].split(";")
        assert set(families) == {"GH18", "GT2"}
        ecs = df.iloc[0]["cazy_ec"].split(";")
        assert set(ecs) == {"3.2.1.14", "2.4.1.-"}
        tools = df.iloc[0]["cazy_tools"].split(";")
        assert set(tools) == {"3", "2"}


# ---------------------------------------------------------------------------
# 6. Confidence scoring
# ---------------------------------------------------------------------------

class TestConfidenceScoring:
    """Test _generate_confidence_annotations via a minimal end-to-end run."""

    @pytest.fixture
    def decoder(self):
        d = MetaFuncDecoder(_minimal_categories())
        # Inject minimal KEGG knowledge so pattern matching works
        d.kegg_mapper.kegg_info = {
            "K01234": {"definition": "chitinase [EC:3.2.1.14]", "ec": "3.2.1.14", "threshold": ""},
            "K05678": {"definition": "unrelated enzyme", "ec": "", "threshold": ""},
        }
        # Inject minimal CAZy knowledge
        d.cazy_mapper.family_info = {
            "GH18": {"activities": "chitinase", "class": "GH", "note": ""},
        }
        return d

    def _run(self, decoder, rows):
        import pandas as pd
        df = pd.DataFrame(rows)
        return decoder._generate_confidence_annotations(df)

    def test_single_db_gives_low_confidence(self, decoder):
        result = self._run(decoder, [
            {"gene_id": "g1", "kegg_ortholog": "KO:K01234", "pfam": None, "cazyme": None}
        ])
        hit = result[result["gene_id"] == "g1"]
        assert hit.iloc[0]["confidence_level"] == "low"
        assert hit.iloc[0]["confidence_score"] == 0.33

    def test_two_dbs_give_medium_confidence(self, decoder):
        result = self._run(decoder, [
            {"gene_id": "g1", "kegg_ortholog": "KO:K01234", "pfam": None, "cazyme": "GH18",
             "cazy_ec": "3.2.1.14", "cazy_tools": "3"}
        ])
        hit = result[(result["gene_id"] == "g1") & (result["confidence_level"] != "none")]
        assert hit.iloc[0]["confidence_level"] == "medium"
        assert hit.iloc[0]["confidence_score"] == pytest.approx(0.67)

    def test_unannotated_gene_gets_none_row(self, decoder):
        result = self._run(decoder, [
            {"gene_id": "g_empty", "kegg_ortholog": None, "pfam": None, "cazyme": None}
        ])
        assert result.iloc[0]["confidence_level"] == "none"
        assert result.iloc[0]["functional_category"] == "unannotated"

    def test_supporting_databases_field(self, decoder):
        result = self._run(decoder, [
            {"gene_id": "g1", "kegg_ortholog": "KO:K01234", "pfam": None, "cazyme": "GH18",
             "cazy_ec": "", "cazy_tools": "2"}
        ])
        hit = result[(result["gene_id"] == "g1") & (result["confidence_level"] != "none")]
        dbs = set(hit.iloc[0]["supporting_databases"].split(";"))
        assert "KEGG" in dbs
        assert "CAZy" in dbs

    def test_cazy_ec_and_tools_in_output(self, decoder):
        result = self._run(decoder, [
            {"gene_id": "g1", "kegg_ortholog": None, "pfam": None, "cazyme": "GH18",
             "cazy_ec": "3.2.1.14", "cazy_tools": "3"}
        ])
        hit = result[(result["gene_id"] == "g1") & (result["confidence_level"] != "none")]
        assert hit.iloc[0]["cazy_ec"] == "3.2.1.14"
        assert hit.iloc[0]["cazy_tool_count"] == "3"

    def test_non_matching_kegg_gives_no_annotated_row(self, decoder):
        result = self._run(decoder, [
            {"gene_id": "g1", "kegg_ortholog": "KO:K05678", "pfam": None, "cazyme": None}
        ])
        # K05678 definition ("unrelated enzyme") matches no category pattern
        assert result.iloc[0]["confidence_level"] == "none"

    def test_multiple_genes_processed_independently(self, decoder):
        result = self._run(decoder, [
            {"gene_id": "g1", "kegg_ortholog": "KO:K01234", "pfam": None, "cazyme": None},
            {"gene_id": "g2", "kegg_ortholog": None, "pfam": None, "cazyme": None},
        ])
        assert len(result[result["gene_id"] == "g1"]) >= 1
        assert len(result[result["gene_id"] == "g2"]) == 1

    def test_three_dbs_give_high_confidence(self, decoder):
        """All three databases agreeing on carbon_cycling must produce confidence=high (1.0)."""
        # Inject a Pfam2GO mapping that maps PF00001 → GO:0005975 (carbohydrate metabolic)
        decoder.pfam_mapper.pfam2go = {
            "PF00001": [{"go_id": "GO:0005975", "go_name": "carbohydrate metabolic process"}]
        }
        result = self._run(decoder, [
            {
                "gene_id": "g_high",
                "kegg_ortholog": "KO:K01234",   # definition: "chitinase" → carbon_cycling
                "pfam": "PF00001",              # GO:0005975 → carbon_cycling
                "cazyme": "GH18",              # activities: "chitinase" → carbon_cycling
                "cazy_ec": "3.2.1.14",
                "cazy_tools": "3",
            }
        ])
        hit = result[
            (result["gene_id"] == "g_high") &
            (result["functional_category"] == "carbon_cycling")
        ]
        assert len(hit) == 1, "Expected one carbon_cycling row for three-DB gene"
        assert hit.iloc[0]["confidence_level"] == "high"
        assert hit.iloc[0]["confidence_score"] == pytest.approx(1.0)
        assert hit.iloc[0]["n_databases"] == 3
        dbs = set(hit.iloc[0]["supporting_databases"].split(";"))
        assert dbs == {"CAZy", "Pfam", "KEGG"}


# ---------------------------------------------------------------------------
# 7. Pattern matching
# ---------------------------------------------------------------------------

class TestPatternMatching:

    def test_kegg_pattern_case_insensitive(self):
        mapper = KEGGMapper({"carbon_cycling": ["Cellulase"]})
        assert mapper._search_patterns_in_text("CELLULASE enzyme") == [("carbon_cycling", "")]

    def test_kegg_no_match_returns_empty(self):
        mapper = KEGGMapper({"carbon_cycling": ["cellulase"]})
        assert mapper._search_patterns_in_text("nitrogenase reductase") == []

    def test_kegg_empty_text_returns_empty(self):
        mapper = KEGGMapper({"carbon_cycling": ["cellulase"]})
        assert mapper._search_patterns_in_text("") == []
        assert mapper._search_patterns_in_text(None) == []

    def test_custom_pattern_matched(self):
        mapper = KEGGMapper({})
        mapper.add_custom_patterns({"laccase": ["laccase"]})
        assert mapper._search_patterns_in_text("laccase oxidase") == [("laccase", "")]

    def test_custom_pattern_not_matched_when_excluded(self):
        mapper = KEGGMapper({})
        mapper.add_custom_patterns({"laccase": ["laccase"]})
        assert mapper._search_patterns_in_text("laccase oxidase", include_custom=False) == []

    def test_cazy_broad_category_lookup(self):
        mapper = CAZyMapper(
            {"cellulase": ["cellulase"]},
            {"cellulase": "carbon_cycling"},
        )
        assert mapper.get_broad_category("cellulase") == "carbon_cycling"
        assert mapper.get_broad_category("unknown_subcat") == "unknown_subcat"

    def test_gh5_text_fallback_multi_subcategory(self):
        """GH5 is intentionally absent from families: lists; text-pattern fallback
        must return multiple subcategory labels reflecting its heterogeneous FamInfo
        activities (cellulase AND xylanase AND mannanase AND beta_glucosidase AND
        arabinogalactanase).  The activities string used here is representative of
        the real FamInfo GH5 entry."""
        mapper = CAZyMapper(
            {
                "cellulase":        [("cellulase", ""), ("endoglucanase", "")],
                "xylanase":         [("xylanase", "")],
                "mannanase":        [("mannanase", "")],
                "beta_glucosidase": [("glucosidase", ""), ("beta-glucosidase", "")],
                "arabinogalactanase": [("galactanase", "")],
            },
            {
                "cellulase":          "carbon_cycling",
                "xylanase":           "carbon_cycling",
                "mannanase":          "carbon_cycling",
                "beta_glucosidase":   "carbon_cycling",
                "arabinogalactanase": "carbon_cycling",
            },
        )
        # Activities string representative of the real GH5 FamInfo entry; includes
        # galactanase activity (endo-beta-1,6-galactanase) that triggers arabinogalactanase
        # and beta-glucosidase activity that triggers beta_glucosidase.
        mapper.family_info = {
            "GH5": {
                "activities": (
                    "Cellulase (EC 3.2.1.4); mannan endo-beta-1,4-mannosidase "
                    "(EC 3.2.1.78); endo-beta-1,4-xylanase (EC 3.2.1.8); "
                    "beta-1,3-mannanase (EC 3.2.1.-); beta-glucosidase (EC 3.2.1.21); "
                    "endo-beta-1,6-galactanase (EC 3.2.1.164)"
                ),
                "class": "GH",
                "note": "",
            }
        }
        subcats = mapper.get_family_categories("GH5")
        assert "cellulase" in subcats,          "GH5 must match cellulase via text pattern"
        assert "xylanase" in subcats,           "GH5 must match xylanase via text pattern"
        assert "mannanase" in subcats,          "GH5 must match mannanase via text pattern"
        assert "beta_glucosidase" in subcats,   "GH5 must match beta_glucosidase via text pattern"
        assert "arabinogalactanase" in subcats, "GH5 must match arabinogalactanase via galactanase pattern"
        assert len(subcats) >= 5, "GH5 text-fallback must produce at least 5 subcategory labels"

    def test_ce4_dual_listing(self):
        """CE4 must appear in both chitinase and other_hemicellulase families lists."""
        mapper = CAZyMapper(
            {},
            {"chitinase": "carbon_cycling", "other_hemicellulase": "carbon_cycling"},
            family_to_subcategory={
                "CE4": ["chitinase", "other_hemicellulase"],
            },
        )
        subcats = mapper.get_family_categories("CE4")
        assert "chitinase" in subcats
        assert "other_hemicellulase" in subcats

    def test_ce12_dual_listing(self):
        """CE12 must appear in both pectinase and other_hemicellulase families lists."""
        mapper = CAZyMapper(
            {},
            {"pectinase": "carbon_cycling", "other_hemicellulase": "carbon_cycling"},
            family_to_subcategory={
                "CE12": ["other_hemicellulase", "pectinase"],
            },
        )
        subcats = mapper.get_family_categories("CE12")
        assert "pectinase" in subcats
        assert "other_hemicellulase" in subcats

    def test_gh2_dual_listing(self):
        """GH2 must appear in both mannanase and other_hemicellulase families lists."""
        mapper = CAZyMapper(
            {},
            {"mannanase": "carbon_cycling", "other_hemicellulase": "carbon_cycling"},
            family_to_subcategory={
                "GH2": ["mannanase", "other_hemicellulase"],
            },
        )
        subcats = mapper.get_family_categories("GH2")
        assert "mannanase" in subcats
        assert "other_hemicellulase" in subcats

    def test_denitrification_excludes_assimilatory_nitrate(self):
        """Negative lookbehind must prevent assimilatory nitrate/nitrite reductase
        from being assigned to denitrification."""
        mapper = KEGGMapper(
            {"denitrification": [("(?<!assimilatory )nitrate reductase", "denitrification"),
                                 ("(?<!assimilatory )nitrite reductase", "denitrification")],
             "assimilatory_nitrate_reduction": [("assimilatory nitrate reductase", "assimilatory_nitrate_reduction"),
                                                ("assimilatory nitrite reductase", "assimilatory_nitrate_reduction")]},
        )
        # Dissimilatory (denitrification) — must match
        hits = dict(mapper._search_patterns_in_text("nitrate reductase narG respiratory chain"))
        assert "denitrification" in hits
        # Assimilatory — must NOT match denitrification
        hits2 = dict(mapper._search_patterns_in_text("assimilatory nitrate reductase NADP-dependent"))
        assert "denitrification" not in hits2
        assert "assimilatory_nitrate_reduction" in hits2
        # Assimilatory nitrite reductase — must NOT match denitrification
        hits3 = dict(mapper._search_patterns_in_text("assimilatory nitrite reductase nirA ferredoxin"))
        assert "denitrification" not in hits3

    def test_tetracycline_no_subcategory(self):
        """tetracycline.*resistance pattern must match broad antibiotics category
        but return no subcategory (empty string)."""
        mapper = KEGGMapper(
            {"antibiotics": [("tetracycline.*resistance", "")]},
        )
        hits = mapper._search_patterns_in_text("tetracycline resistance protein tet(A)")
        assert len(hits) == 1
        cat, subcat = hits[0]
        assert cat == "antibiotics"
        assert subcat == ""


# ---------------------------------------------------------------------------
# 8. Parser regression tests
# ---------------------------------------------------------------------------

class TestPfamParserAutodetect:
    """Pfam domtblout files without leading # comment lines must still be
    recognised as HMMER format by field-count heuristic."""

    @pytest.fixture
    def decoder(self):
        return MetaFuncDecoder(_minimal_categories())

    def test_hmmer_no_comment_lines_detected(self, decoder, tmp_path):
        """domtblout data line with ≥20 space-separated fields → HMMER format."""
        p = tmp_path / "no_comments.domtblout"
        p.write_text(
            "7tm_1  PF00001.23  263  gene_001  -  300  1e-50  300.0  0.0  1  1  1e-50  300  0.0  1  200  1  200  1  200  0.99\n"
            "7tm_2  PF00002.10  200  gene_002  -  250  1e-30  200.0  0.0  1  1  1e-30  200  0.0  1  150  1  150  1  150  0.95\n"
        )
        df = decoder._load_pfam_file(str(p))
        assert set(df["gene_id"]) == {"gene_001", "gene_002"}, \
            "domtblout without # header should still be parsed as HMMER format"


class TestIncludeUnannotatedWithFunctionFilter:
    """--include-unannotated must retain unannotated rows even when a function
    filter is applied via save_results(functions=...)."""

    @pytest.fixture
    def decoder(self):
        return MetaFuncDecoder(_minimal_categories())

    def test_unannotated_rows_retained_with_function_filter(self, decoder, tmp_path):
        import pandas as pd
        # Build a tiny annotations table with one annotated + one unannotated row
        rows = [
            {
                "gene_id": "gene_001",
                "functional_category": "carbon_cycling",
                "confidence_level": "low",
                "confidence_score": 0.33,
                "subcategories": "",
                "cazy_families": "",
                "cazy_ec": "",
                "cazy_tool_count": "",
                "pfam_ids": "",
                "ko_terms": "",
                "cazy_description": "",
                "pfam_go_names": "",
                "ko_definition": "",
                "supporting_databases": "KEGG",
                "n_databases": 1,
            },
            {
                "gene_id": "gene_002",
                "functional_category": "unannotated",
                "confidence_level": "none",
                "confidence_score": 0.0,
                "subcategories": "",
                "cazy_families": "",
                "cazy_ec": "",
                "cazy_tool_count": "",
                "pfam_ids": "",
                "ko_terms": "",
                "cazy_description": "",
                "pfam_go_names": "",
                "ko_definition": "",
                "supporting_databases": "",
                "n_databases": 0,
            },
        ]
        decoder.annotations["metagenome_confidence"] = pd.DataFrame(rows)
        decoder.save_results(
            str(tmp_path),
            include_unannotated=True,
            functions=["carbon_cycling"],
        )
        out = pd.read_csv(tmp_path / "metagenome_confidence_annotations.csv")
        assert "gene_001" in out["gene_id"].values, "annotated gene must be present"
        assert "gene_002" in out["gene_id"].values, \
            "unannotated gene must be retained when include_unannotated=True"


class TestCombinedModeStrictColumnValidation:
    """Combined mode must exit with an error on missing annotation columns
    unless --allow-missing-cols is passed."""

    @pytest.fixture
    def decoder(self):
        return MetaFuncDecoder(_minimal_categories())

    def test_missing_column_raises_systemexit(self, decoder, tmp_path):
        import pytest
        p = tmp_path / "missing_col.tsv"
        p.write_text("gene_id\tko_id\n" "gene_001\tK00001\n")
        with pytest.raises(SystemExit):
            decoder.load_combined_mode(
                str(p),
                ko_col="ko_id",
                allow_missing_cols=False,
            )

    def test_missing_column_allowed_with_flag(self, decoder, tmp_path):
        p = tmp_path / "missing_col.tsv"
        p.write_text("gene_id\tko_id\n" "gene_001\tK00001\n")
        # Should not raise
        decoder.load_combined_mode(
            str(p),
            ko_col="ko_id",
            allow_missing_cols=True,
        )
        assert "metagenome" in decoder.data


# ---------------------------------------------------------------------------
# 9. Combined-mode regression tests
# ---------------------------------------------------------------------------

class TestDuplicateGeneAggregation:
    """Combined mode must aggregate duplicate gene_id rows rather than
    emitting them as separate rows, which would inflate downstream counts."""

    @pytest.fixture
    def decoder(self):
        return MetaFuncDecoder(_minimal_categories())

    def test_duplicate_gene_ids_aggregated(self, decoder, tmp_path):
        p = tmp_path / "dupes.tsv"
        p.write_text(
            "gene_id\tkegg_ortholog\tpfam\tcazyme\n"
            "gene_001\tKO:K01234\t\t\n"
            "gene_001\tKO:K05678\t\t\n"   # second row for same gene
            "gene_002\t\tPF00001\t\n"
        )
        decoder.load_combined_mode(str(p), allow_missing_cols=True)
        df = decoder.data["metagenome"]
        # Only one row per gene after aggregation
        assert len(df) == 2, "Duplicate rows must be collapsed to one per gene_id"
        gene1 = df[df["gene_id"] == "gene_001"].iloc[0]
        # Both KO terms must be present in the aggregated cell
        assert "KO:K01234" in gene1["kegg_ortholog"]
        assert "KO:K05678" in gene1["kegg_ortholog"]


class TestSubcategoryNDatabases:
    """subcategory_n_databases must reflect how many databases independently
    support the same subcategory, not just the broad category."""

    @pytest.fixture
    def decoder(self):
        # Use subcategory-labelled kegg_patterns so chitinase hits carry
        # a subcategory label and are tracked in kegg_subcategories.
        cats = {
            "carbon_cycling": {
                "kegg_patterns": [
                    {"pattern": "chitinase", "subcategory": "chitinase"},
                ],
                "pfam_patterns": ["carbohydrate"],
                "pfam_go_terms": ["GO:0005975"],
                "cazy_subcategories": {
                    "chitinase": {"patterns": ["chitinase", "chitin"]},
                },
            },
        }
        d = MetaFuncDecoder(cats)
        d.kegg_mapper.kegg_info = {
            "K01234": {"definition": "chitinase [EC:3.2.1.14]", "ec": "3.2.1.14", "threshold": ""},
        }
        d.cazy_mapper.family_info = {
            "GH18": {"activities": "chitinase", "class": "GH", "note": ""},
        }
        return d

    def _run(self, decoder, rows):
        import pandas as pd
        return decoder._generate_confidence_annotations(pd.DataFrame(rows))

    def test_subcategory_n_databases_two_dbs_agree(self, decoder):
        """KEGG + CAZy both resolve to 'chitinase' subcategory → subcategory_n_databases = 2."""
        result = self._run(decoder, [{
            "gene_id": "g1",
            "kegg_ortholog": "KO:K01234",   # → chitinase
            "pfam": None,
            "cazyme": "GH18",               # → chitinase
            "cazy_ec": "3.2.1.14",
            "cazy_tools": "3",
        }])
        hit = result[
            (result["gene_id"] == "g1") &
            (result["functional_category"] == "carbon_cycling")
        ]
        assert len(hit) == 1
        assert hit.iloc[0]["subcategory_n_databases"] == 2

    def test_subcategory_n_databases_one_db(self, decoder):
        """Single-DB hit → subcategory_n_databases = 1."""
        result = self._run(decoder, [{
            "gene_id": "g2",
            "kegg_ortholog": "KO:K01234",   # → chitinase (KEGG only)
            "pfam": None,
            "cazyme": None,
        }])
        hit = result[
            (result["gene_id"] == "g2") &
            (result["functional_category"] == "carbon_cycling")
        ]
        assert len(hit) == 1
        assert hit.iloc[0]["subcategory_n_databases"] == 1

    def test_subcategory_n_databases_zero_when_no_subcategory(self, decoder):
        """Broad GO-term Pfam hit with no subcategory → subcategory_n_databases = 0."""
        decoder.pfam_mapper.pfam2go = {
            "PF00001": [{"go_id": "GO:0005975", "go_name": "carbohydrate metabolic process"}]
        }
        # carbon_cycling matched via Pfam GO:0005975 — but no subcategory label
        result = self._run(decoder, [{
            "gene_id": "g3",
            "kegg_ortholog": None,
            "pfam": "PF00001",
            "cazyme": None,
        }])
        hit = result[
            (result["gene_id"] == "g3") &
            (result["functional_category"] == "carbon_cycling")
        ]
        assert len(hit) == 1
        assert hit.iloc[0]["subcategory_n_databases"] == 0


class TestRunManifestWritten:
    """write_run_manifest must emit run_manifest.json with required provenance fields."""

    @pytest.fixture
    def decoder(self):
        return MetaFuncDecoder(_minimal_categories())

    def test_run_manifest_written(self, decoder, tmp_path):
        import json
        import argparse

        # Write a dummy DB file so the SHA256 can be computed
        dummy_db = tmp_path / "ko_list"
        dummy_db.write_text("knum\tdefinition\nK00001\tdummy\n")

        decoder.write_run_manifest(
            output_dir=str(tmp_path),
            args_namespace=argparse.Namespace(),
            db_paths={"ko_list": str(dummy_db)},
            categories_dir=str(tmp_path),   # empty dir — no YAML files, yaml_hashes = {}
        )

        manifest_path = tmp_path / "run_manifest.json"
        assert manifest_path.exists(), "run_manifest.json must be created"

        with open(manifest_path) as f:
            manifest = json.load(f)

        assert "metafuncdecoder_version" in manifest
        assert "run_date" in manifest
        assert "command_line" in manifest
        assert "databases" in manifest
        assert "ko_list" in manifest["databases"]
        # SHA256 must be a 64-character hex string (not "unavailable")
        sha = manifest["databases"]["ko_list"]["sha256"]
        assert len(sha) == 64 and all(c in "0123456789abcdef" for c in sha)


# ---------------------------------------------------------------------------
# 10. Evidence specificity and strict evidence mode
# ---------------------------------------------------------------------------

class TestEvidenceSpecificity:
    """evidence_specificity must be 'specific' when a subcategory label is
    present, 'broad' when the category was matched only by a non-specific
    pattern, and empty string for unannotated rows."""

    @pytest.fixture
    def decoder(self):
        cats = {
            "carbon_cycling": {
                "kegg_patterns": [
                    {"pattern": "chitinase", "subcategory": "chitinase"},
                    "carbohydrate.*metabolism",   # broad — no subcategory label
                ],
                "pfam_patterns": ["carbohydrate"],
                "pfam_go_terms": ["GO:0005975"],
                "cazy_subcategories": {
                    "chitinase": {"patterns": ["chitinase"]},
                },
            },
        }
        d = MetaFuncDecoder(cats)
        d.kegg_mapper.kegg_info = {
            "K01234": {"definition": "chitinase [EC:3.2.1.14]", "ec": "3.2.1.14", "threshold": ""},
            "K99999": {"definition": "carbohydrate metabolism general", "ec": "", "threshold": ""},
        }
        return d

    def _run(self, decoder, rows):
        import pandas as pd
        return decoder._generate_confidence_annotations(pd.DataFrame(rows))

    def test_specific_hit_gives_specific(self, decoder):
        """KEGG hit matching a labelled subcategory pattern → evidence_specificity='specific'."""
        result = self._run(decoder, [{
            "gene_id": "g1",
            "kegg_ortholog": "KO:K01234",   # → chitinase subcategory
            "pfam": None,
            "cazyme": None,
        }])
        hit = result[(result["gene_id"] == "g1") & (result["functional_category"] != "unannotated")]
        assert len(hit) == 1
        assert hit.iloc[0]["evidence_specificity"] == "specific"

    def test_broad_hit_gives_broad(self, decoder):
        """KEGG hit matching only a plain-string pattern with no subcategory → 'broad'."""
        result = self._run(decoder, [{
            "gene_id": "g2",
            "kegg_ortholog": "KO:K99999",   # definition: "carbohydrate metabolism general"
            "pfam": None,
            "cazyme": None,
        }])
        hit = result[(result["gene_id"] == "g2") & (result["functional_category"] != "unannotated")]
        assert len(hit) == 1
        assert hit.iloc[0]["evidence_specificity"] == "broad"

    def test_unannotated_row_has_empty_specificity(self, decoder):
        """Unannotated genes must have evidence_specificity=''."""
        result = self._run(decoder, [{
            "gene_id": "g3",
            "kegg_ortholog": None,
            "pfam": None,
            "cazyme": None,
        }])
        assert result.iloc[0]["functional_category"] == "unannotated"
        assert result.iloc[0]["evidence_specificity"] == ""


class TestCazyMinTools:
    """--cazy-min-tools must drop CAZy families whose tool count falls below
    the threshold before they contribute to the confidence model."""

    @pytest.fixture
    def decoder(self):
        d = MetaFuncDecoder(_minimal_categories())
        d.cazy_mapper.family_info = {
            "GH18": {"activities": "chitinase", "class": "GH", "note": ""},
        }
        return d

    def _run(self, decoder, rows, min_cazy_tools=1):
        import pandas as pd
        return decoder._generate_confidence_annotations(
            pd.DataFrame(rows), min_cazy_tools=min_cazy_tools
        )

    def test_cazy_below_threshold_excluded(self, decoder):
        """GH18 with 1 tool must be excluded when min_cazy_tools=2."""
        result = self._run(decoder, [{
            "gene_id": "g1",
            "kegg_ortholog": None,
            "pfam": None,
            "cazyme": "GH18",
            "cazy_ec": "",
            "cazy_tools": "1",
        }], min_cazy_tools=2)
        # CAZy evidence dropped → no annotated row
        assert result.iloc[0]["confidence_level"] == "none"

    def test_cazy_at_threshold_included(self, decoder):
        """GH18 with 2 tools must be included when min_cazy_tools=2."""
        result = self._run(decoder, [{
            "gene_id": "g1",
            "kegg_ortholog": None,
            "pfam": None,
            "cazyme": "GH18",
            "cazy_ec": "",
            "cazy_tools": "2",
        }], min_cazy_tools=2)
        hit = result[result["confidence_level"] != "none"]
        assert len(hit) == 1

    def test_min_tools_1_passes_all(self, decoder):
        """Default min_cazy_tools=1 must pass all CAZy hits regardless of tool count."""
        result = self._run(decoder, [{
            "gene_id": "g1",
            "kegg_ortholog": None,
            "pfam": None,
            "cazyme": "GH18",
            "cazy_ec": "",
            "cazy_tools": "1",
        }], min_cazy_tools=1)
        hit = result[result["confidence_level"] != "none"]
        assert len(hit) == 1

    def test_missing_tool_count_passes_through(self, decoder):
        """When tool count is absent (combined table without that column), hit
        must be allowed through regardless of min_cazy_tools setting."""
        result = self._run(decoder, [{
            "gene_id": "g1",
            "kegg_ortholog": None,
            "pfam": None,
            "cazyme": "GH18",
            "cazy_ec": "",
            "cazy_tools": "",   # empty — no tool count available
        }], min_cazy_tools=3)
        hit = result[result["confidence_level"] != "none"]
        assert len(hit) == 1

    def test_float_formatted_tool_count_filtered(self, decoder):
        """cazy_tools='3.0' (pandas/R float export) must be parsed as 3 and
        filtered correctly — int('3.0') raises ValueError but int(float('3.0'))
        returns 3. A threshold of 4 must exclude it."""
        result = self._run(decoder, [{
            "gene_id": "g1",
            "kegg_ortholog": None,
            "pfam": None,
            "cazyme": "GH18",
            "cazy_ec": "",
            "cazy_tools": "3.0",   # float-formatted, common in pandas/R output
        }], min_cazy_tools=4)
        assert result.iloc[0]["confidence_level"] == "none"

    def test_float_formatted_tool_count_passes_when_above_threshold(self, decoder):
        """cazy_tools='3.0' with min_cazy_tools=2 must pass through (3 >= 2)."""
        result = self._run(decoder, [{
            "gene_id": "g1",
            "kegg_ortholog": None,
            "pfam": None,
            "cazyme": "GH18",
            "cazy_ec": "",
            "cazy_tools": "3.0",
        }], min_cazy_tools=2)
        hit = result[result["confidence_level"] != "none"]
        assert len(hit) == 1


class TestExcludeBroad:
    """--exclude-broad must remove rows with evidence_specificity='broad'
    from the saved CSV while leaving specific and unannotated rows intact."""

    @pytest.fixture
    def decoder(self):
        return MetaFuncDecoder(_minimal_categories())

    def test_exclude_broad_removes_broad_rows(self, decoder, tmp_path):
        import pandas as pd
        rows = [
            {
                "gene_id": "gene_specific",
                "functional_category": "carbon_cycling",
                "evidence_specificity": "specific",
                "confidence_level": "low",
                "confidence_score": 0.33,
                "subcategories": "chitinase",
                "subcategory_n_databases": 1,
                "cazy_families": "", "cazy_subfamilies": "", "cazy_ec": "",
                "cazy_tool_count": "", "pfam_ids": "", "ko_terms": "",
                "cazy_description": "", "pfam_go_names": "", "ko_definition": "",
                "supporting_databases": "KEGG", "n_databases": 1,
            },
            {
                "gene_id": "gene_broad",
                "functional_category": "carbon_cycling",
                "evidence_specificity": "broad",
                "confidence_level": "low",
                "confidence_score": 0.33,
                "subcategories": "",
                "subcategory_n_databases": 0,
                "cazy_families": "", "cazy_subfamilies": "", "cazy_ec": "",
                "cazy_tool_count": "", "pfam_ids": "PF00001", "ko_terms": "",
                "cazy_description": "", "pfam_go_names": "carbohydrate metabolic process",
                "ko_definition": "",
                "supporting_databases": "Pfam", "n_databases": 1,
            },
        ]
        decoder.annotations["metagenome_confidence"] = pd.DataFrame(rows)
        decoder.save_results(str(tmp_path), exclude_broad=True, functions=["carbon_cycling"])
        out = pd.read_csv(tmp_path / "metagenome_confidence_annotations.csv")
        assert "gene_specific" in out["gene_id"].values
        assert "gene_broad" not in out["gene_id"].values

    def test_exclude_broad_false_retains_broad_rows(self, decoder, tmp_path):
        import pandas as pd
        rows = [
            {
                "gene_id": "gene_broad",
                "functional_category": "carbon_cycling",
                "evidence_specificity": "broad",
                "confidence_level": "low",
                "confidence_score": 0.33,
                "subcategories": "",
                "subcategory_n_databases": 0,
                "cazy_families": "", "cazy_subfamilies": "", "cazy_ec": "",
                "cazy_tool_count": "", "pfam_ids": "PF00001", "ko_terms": "",
                "cazy_description": "", "pfam_go_names": "carbohydrate metabolic process",
                "ko_definition": "",
                "supporting_databases": "Pfam", "n_databases": 1,
            },
        ]
        decoder.annotations["metagenome_confidence"] = pd.DataFrame(rows)
        decoder.save_results(str(tmp_path), exclude_broad=False, functions=["carbon_cycling"])
        out = pd.read_csv(tmp_path / "metagenome_confidence_annotations.csv")
        assert "gene_broad" in out["gene_id"].values
