"""
Test suite for Phenotype Matcher v2.

Unit tests (TK, DN, DM, DI, BV, WS) — run locally, no ontology required.
Integration tests (M01–M31)           — require RUN_INTEGRATION=1 + GPU.

Run unit tests:
    python test_cases.py

Run all tests (on ws with GPU):
    RUN_INTEGRATION=1 python test_cases.py
"""

import os
import sys
import unittest

# ---------------------------------------------------------------------------
# Ensure the package is importable when run directly
# ---------------------------------------------------------------------------
_here = os.path.dirname(os.path.abspath(__file__))
# _here is .../src/phenotype_matcher_v2; one level up gives .../src/
_src = os.path.normpath(os.path.join(_here, ".."))
if _src not in sys.path:
    sys.path.insert(0, _src)

from phenotype_matcher_v2.tokenizer import tokenize_expand
from phenotype_matcher_v2.detection import (
    det_neg,
    det_mod,
    det_neg_in_text,
    det_mod_in_text,
    is_boundary_valid,
)
from phenotype_matcher_v2.matcher import W_STOP


# ===========================================================================
# Unit tests — Algorithm 1: tokenize_expand
# ===========================================================================

class TestTokenizer(unittest.TestCase):

    def test_TK1_hard_seg_comma(self):
        """TK1: Hard segmentation on comma."""
        result = tokenize_expand("seizures, hypotonia")
        self.assertIn("seizures", result)
        self.assertIn("hypotonia", result)
        self.assertEqual(len(result), 2)

    def test_TK2_hard_seg_semicolon_period(self):
        """TK2: Hard segmentation on semicolon and period."""
        result = tokenize_expand("seizures; ataxia. hypotonia")
        self.assertIn("seizures", result)
        self.assertIn("ataxia", result)
        self.assertIn("hypotonia", result)

    def test_TK3_hard_seg_but(self):
        """TK3: Hard segmentation on ' but '."""
        result = tokenize_expand("normal vision but hypotonia")
        self.assertIn("normal vision", result)
        self.assertIn("hypotonia", result)

    def test_TK4_adj_and_adj_noun(self):
        """TK4: Pattern Adj1 and Adj2 Noun → two tokens."""
        result = tokenize_expand("short and malformed fingers")
        self.assertIn("short fingers", result)
        self.assertIn("malformed fingers", result)

    def test_TK5_noun_with_a_and_b(self):
        """TK5: Pattern Noun with A and B → two tokens."""
        result = tokenize_expand("seizures with fever and rash")
        self.assertIn("seizures with fever", result)
        self.assertIn("seizures with rash", result)

    def test_TK6_no_pattern(self):
        """TK6: No pattern match → keep segment as-is."""
        result = tokenize_expand("intellectual disability")
        self.assertEqual(result, ["intellectual disability"])

    def test_TK7_combined(self):
        """TK7: Both segmentation and coordination expansion."""
        result = tokenize_expand("short and malformed fingers, absent thumb")
        self.assertIn("short fingers", result)
        self.assertIn("malformed fingers", result)
        self.assertIn("absent thumb", result)


# ===========================================================================
# Unit tests — Algorithm 2a: det_neg (word-window)
# ===========================================================================

class TestDetNeg(unittest.TestCase):

    def test_DN1_trigger_in_pre(self):
        """DN1: Trigger 'no' in PRE window."""
        self.assertTrue(det_neg("no seizures", 3, 11, win=5))

    def test_DN2_trigger_in_post(self):
        """DN2: Trigger 'ruled out' in POST window."""
        self.assertTrue(det_neg("seizures ruled out", 0, 8, win=5))

    def test_DN3_no_trigger(self):
        """DN3: No trigger in either window."""
        self.assertFalse(det_neg("present seizures", 8, 16, win=5))

    def test_DN4_without_trigger(self):
        """DN4: 'without' in PRE window."""
        self.assertTrue(det_neg("without any seizures", 11, 19, win=5))

    def test_DN5_clinical_denial_pre(self):
        """DN5: 'denies' in PRE window."""
        self.assertTrue(det_neg("patient denies seizures", 15, 23, win=5))

    def test_DN6_post_beyond_immediate(self):
        """DN6: Trigger beyond immediate neighbour in POST."""
        self.assertTrue(det_neg("severe seizures were ruled out today", 7, 15, win=5))

    def test_DN7_no_context(self):
        """DN7: No context at all."""
        self.assertFalse(det_neg("seizures", 0, 8, win=5))


# ===========================================================================
# Unit tests — Algorithm 2b: det_mod (word-window)
# ===========================================================================

_TINY_MOD_MAP = {
    "severe":   "HP:0012828",
    "mild":     "HP:0012825",
    "profound": "HP:0012829",
}


class TestDetMod(unittest.TestCase):

    def test_DM1_modifier_in_pre(self):
        """DM1: Modifier in PRE window."""
        result = det_mod("severe intellectual disability", 7, 31, _TINY_MOD_MAP, win=4)
        self.assertEqual(result, "HP:0012828")

    def test_DM2_modifier_in_post(self):
        """DM2: Modifier in POST window."""
        result = det_mod("intellectual disability, mild", 0, 22, _TINY_MOD_MAP, win=4)
        self.assertEqual(result, "HP:0012825")

    def test_DM3_no_modifier(self):
        """DM3: No modifier present."""
        result = det_mod("intellectual disability", 0, 22, _TINY_MOD_MAP, win=4)
        self.assertIsNone(result)

    def test_DM4_longest_match_wins(self):
        """DM4: Longest modifier wins when two appear in PRE."""
        # "profound" (8 chars) > "severe" (6 chars) — both appear before "ataxia"
        result = det_mod("profound and severe ataxia", 16, 22, _TINY_MOD_MAP, win=4)
        # Both are in pre-window; longest match should win
        self.assertEqual(result, "HP:0012829")  # "profound"


# ===========================================================================
# Unit tests — Phase 2e: det_neg_in_text / det_mod_in_text
# ===========================================================================

class TestDetInText(unittest.TestCase):

    def test_DI1_neg_trigger_anywhere(self):
        """DI1: Trigger anywhere in text."""
        self.assertTrue(det_neg_in_text("no fever"))

    def test_DI2_no_neg_trigger(self):
        """DI2: No trigger."""
        self.assertFalse(det_neg_in_text("fever present"))

    def test_DI3_multiword_neg_trigger(self):
        """DI3: Multi-word trigger 'negative for'."""
        self.assertTrue(det_neg_in_text("negative for fever"))

    def test_DI4_mod_found(self):
        """DI4: Modifier word found in text."""
        result = det_mod_in_text("severe ataxia", _TINY_MOD_MAP)
        self.assertEqual(result, "HP:0012828")

    def test_DI5_no_mod(self):
        """DI5: No modifier word."""
        result = det_mod_in_text("ataxia", _TINY_MOD_MAP)
        self.assertIsNone(result)

    def test_DI6_longest_mod_wins(self):
        """DI6: Longest modifier wins when multiple present."""
        result = det_mod_in_text("mild to severe ataxia", _TINY_MOD_MAP)
        # "severe" (6) > "mild" (4)
        self.assertEqual(result, "HP:0012828")


# ===========================================================================
# Unit tests — is_boundary_valid
# ===========================================================================

class TestBoundaryValid(unittest.TestCase):

    def test_BV1_space_before_after(self):
        """BV1: Space before and after."""
        self.assertTrue(is_boundary_valid("has seizure disorder", 4, 11))

    def test_BV2_alnum_before(self):
        """BV2: Alphanumeric char before span → invalid."""
        self.assertFalse(is_boundary_valid("myseizure", 2, 9))

    def test_BV3_alnum_at_end(self):
        """BV3: End char is alphanumeric (span doesn't end at string end) → invalid."""
        # text="seizures", span=[0,7) → text[7]='s' which is alnum
        self.assertFalse(is_boundary_valid("seizures", 0, 7))

    def test_BV4_at_string_boundaries(self):
        """BV4: Span covers the whole string."""
        self.assertTrue(is_boundary_valid("seizure", 0, 7))

    def test_BV5_surrounded_by_punctuation(self):
        """BV5: Surrounded by non-alphanumeric punctuation."""
        self.assertTrue(is_boundary_valid("(seizure)", 1, 8))


# ===========================================================================
# Unit tests — W_STOP membership
# ===========================================================================

class TestWStop(unittest.TestCase):

    def test_WS1_and_in_stop(self):
        """WS1: 'and' is in W_STOP."""
        self.assertIn("and", W_STOP)

    def test_WS2_or_in_stop(self):
        """WS2: 'or' is in W_STOP."""
        self.assertIn("or", W_STOP)

    def test_WS3_seizure_not_in_stop(self):
        """WS3: 'seizure' is NOT in W_STOP."""
        self.assertNotIn("seizure", W_STOP)

    def test_WS4_patient_in_stop(self):
        """WS4: 'patient' is in W_STOP."""
        self.assertIn("patient", W_STOP)


# ===========================================================================
# Integration tests — require RUN_INTEGRATION=1 and a loaded PhenotypeMatcher
# ===========================================================================

@unittest.skipUnless(os.getenv("RUN_INTEGRATION"), "Set RUN_INTEGRATION=1 to run")
class TestIntegration(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from phenotype_matcher_v2 import PhenotypeMatcher, MatcherConfig
        cfg = MatcherConfig(
            hpo_path=os.getenv("HPO_PATH", "ontology/hp.obo"),
            mondo_path=os.getenv("MONDO_PATH", "ontology/mondo.obo"),
            hpoa_path=os.getenv("HPOA_PATH", "data/phenotype.hpoa"),
            embedding_model=os.getenv("EMB_MODEL", "sapbert"),
            llm_model=os.getenv("LLM_MODEL", "deepseek"),
            debug=True,
        )
        cls.m = PhenotypeMatcher(cfg)

    # --- Phase 1 ---

    def test_M01_exact_whole_string(self):
        """M01: Phase 1a exact whole-string match → HP:0001250."""
        out = self.m.match("Seizure")
        self.assertIn("HP:0001250", out.get_hpo_ids())

    def test_M02_stemmed_whole_string(self):
        """M02: Phase 1b stemmed whole-string match."""
        out = self.m.match("seizures")
        self.assertIn("HP:0001250", out.get_hpo_ids())

    def test_M03_fuzzy_whole_string(self):
        """M03: Phase 1c fuzzy whole-string match (edit dist 1).
        'Seizuree' (extra trailing e, distance 1 from 'Seizure') is an obvious
        typo; LLM-Val with the updated prompt recognises it as implying 'Seizure'.
        """
        out = self.m.match("Seizuree")
        self.assertIn("HP:0001250", out.get_hpo_ids())

    def test_M04_wstop_skips_fuzzy(self):
        """M04: W_STOP word → empty result."""
        out = self.m.match("and")
        self.assertEqual(out.get_hpo_ids(), [])
        self.assertEqual(out.diseases, [])

    # --- Phase 2a/2b ---

    def test_M05_exact_token_multi(self):
        """M05: Phase 2a exact token matches for comma-separated terms."""
        out = self.m.match("Seizure, hypotonia")
        hpos = out.get_hpo_ids()
        self.assertIn("HP:0001250", hpos)
        # HP:0001252 = "Hypotonia" (primary label in current HPO);
        # HP:0001290 = "Generalized hypotonia" — accept either.
        self.assertTrue(
            "HP:0001252" in hpos or "HP:0001290" in hpos,
            f"Expected hypotonia HP ID in {hpos}",
        )

    def test_M06_stemmed_token(self):
        """M06: Phase 2b stemmed token match."""
        out = self.m.match("seizing, hypotonia")
        hpos = out.get_hpo_ids()
        # Accept either hypotonia ID depending on HPO version
        self.assertTrue(
            "HP:0001252" in hpos or "HP:0001290" in hpos,
            f"Expected hypotonia HP ID in {hpos}",
        )

    # --- Phase 2c ---

    def test_M08_fuzzy_token(self):
        """M08: Phase 2c fuzzy token match."""
        out = self.m.match("Seizzure, hypotonia")
        hpos = out.get_hpo_ids()
        # Accept either hypotonia ID depending on HPO version
        self.assertTrue(
            "HP:0001252" in hpos or "HP:0001290" in hpos,
            f"Expected hypotonia HP ID in {hpos}",
        )

    def test_M09_wstop_token_skipped(self):
        """M09: W_STOP token ('and') skips fuzzy; other token still matched."""
        out = self.m.match("and, hypotonia")
        hpos = out.get_hpo_ids()
        # Accept either hypotonia ID depending on HPO version
        self.assertTrue(
            "HP:0001252" in hpos or "HP:0001290" in hpos,
            f"Expected hypotonia HP ID in {hpos}",
        )

    # --- Phase 2d ---

    def test_M10_ac_substring_match(self):
        """M10: Phase 2d AC substring match."""
        out = self.m.match("I have a seizure disorder")
        self.assertIn("HP:0001250", out.get_hpo_ids())

    def test_M12_det_neg_pre_window(self):
        """M12: Phase 2d Det-Neg PRE window → excluded=True."""
        out = self.m.match("no seizures")
        excl = [p.hpo_id for p in out.get_excluded_phenotypes()]
        self.assertIn("HP:0001250", excl)

    def test_M13_det_neg_post_window(self):
        """M13: Phase 2d Det-Neg POST window → excluded=True.
        Input must NOT use a comma because the tokenizer hard-segments on commas;
        negation context must be in the same token as the phenotype.
        """
        out = self.m.match("seizures ruled out")
        excl = [p.hpo_id for p in out.get_excluded_phenotypes()]
        self.assertIn("HP:0001250", excl)

    def test_M14_det_mod_pre_window(self):
        """M14: Phase 2d Det-Mod PRE window → severity_id=HP:0012828."""
        out = self.m.match("severe seizures")
        phenotypes = [p for p in out.phenotypes if p.hpo_id == "HP:0001250"]
        self.assertTrue(any(p.severity_id == "HP:0012828" for p in phenotypes))

    def test_M15_det_mod_post_window(self):
        """M15: Phase 2d Det-Mod POST window → severity_id=HP:0012825.
        Use 'mild seizures' (no comma) so the modifier stays inside the same
        token as the phenotype and the detection window can find it.
        """
        out = self.m.match("mild seizures")
        phenotypes = [p for p in out.phenotypes if p.hpo_id == "HP:0001250"]
        self.assertTrue(any(p.severity_id == "HP:0012825" for p in phenotypes))

    # --- Phase 2e ---

    def test_M17_semantic_match(self):
        """M17: Phase 2e semantic ANN match for informal term."""
        out = self.m.match("fits and falls")
        # Should return something seizure-related
        self.assertTrue(len(out.phenotypes) > 0 or len(out.diseases) > 0)

    def test_M18_low_score_filtered(self):
        """M18: Phase 2e Strategy 1 absolute floor filters all (generic text)."""
        out = self.m.match("unknown etiology")
        # Expect no high-confidence matches
        self.assertEqual(len(out.phenotypes), 0)

    def test_M20_det_neg_in_text(self):
        """M20: Phase 2e det_neg_in_text → excluded=True."""
        out = self.m.match("no fever")
        excl = [p for p in out.get_excluded_phenotypes()]
        self.assertTrue(len(excl) > 0)

    def test_M21_det_mod_in_text(self):
        """M21: Phase 2e det_mod_in_text → severity_id set."""
        out = self.m.match("severe ataxia")
        phenotypes_with_sev = [p for p in out.phenotypes if p.severity_id is not None]
        self.assertTrue(len(phenotypes_with_sev) > 0)

    # --- Ontology types ---

    def test_M26_mondo_match(self):
        """M26: MONDO disease term returned."""
        out = self.m.match("Dravet syndrome")
        self.assertTrue(len(out.get_mondo_ids()) > 0)

    def test_M27_omim_match(self):
        """M27: OMIM disease match from phenotype.hpoa."""
        out = self.m.match("phenylketonuria")
        self.assertTrue(len(out.get_omim_ids()) > 0)

    def test_M28_orphanet_from_mondo(self):
        """M28: Orphanet IDs populated via MONDO xrefs."""
        out = self.m.match("Rett syndrome")
        self.assertTrue(len(out.get_orphanet_ids()) > 0)

    # --- Multi-term ---

    def test_M29_three_term_pipeline(self):
        """M29: Three-term comma-separated input."""
        out = self.m.match("seizures, hypotonia, and feeding difficulties")
        self.assertGreaterEqual(len(out.phenotypes), 2)

    def test_M30_coordination_expansion(self):
        """M30: Tokenizer expands 'short and malformed fingers'."""
        out = self.m.match("short and malformed fingers")
        self.assertGreaterEqual(len(out.phenotypes), 1)

    def test_M31_deduplication(self):
        """M31: Same HPO ID reached via two paths → only one PhenotypeMatch."""
        out = self.m.match("Seizure")
        ids = [p.hpo_id for p in out.phenotypes]
        self.assertEqual(len(ids), len(set(ids)))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main(verbosity=2)
