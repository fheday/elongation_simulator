""" Unit-test for the concentraitons generator engine"""
import tempfile
import unittest  # The test framework
import json
import pandas as pd
from concentrations import concentrations_generator


class TestConcentrationsGeneratorPairings(unittest.TestCase):
    """
    Unit-test class for the Concentrations Generator.
    Here we do several condon-anticodon pairing tests.
    """

    def setUp(self):
        self.basepairing_rules = {
            "Watson-Crick": {
                "A": ["U", "&", "3", "1", "~", "N", "S", ")", "{", "V", "}", "P"],
                "C": ["G", "#", "W"],
                "G": ["C", "B"],
                "U": ["A"],
            },
            "Wobble": {
                "A": ["A", "I", "M", "?"],
                "C": ["A", "U", "P", "I", "?", "Q"],
                "G": [
                    "A",
                    "U",
                    "&",
                    "3",
                    "1",
                    "~",
                    "N",
                    "S",
                    ")",
                    "{",
                    "V",
                    "P",
                    "?",
                    "M",
                ],
                "U": ["G", "#", "W", "U", "V", "P", "I", "Q"],
            },
            "Pairing Rules": {"Near-Cognate": {"base-level": [["Wo", "WC", "Wo"]]}},
        }

    def test_default_concentrations_generation(self):

        tRNAs = pd.read_csv("/home/heday/Projects/R_concentrations/data/tRNAs.csv")  # pylint: disable=invalid-name
        codons = pd.read_csv("/home/heday/Projects/R_concentrations/data/codons.csv")
        pairing_rules_string = {
            "Watson-Crick": {
                "A": ["U", "&", "N", "3", "1", "P", "~"],
                "C": ["G", "#", "Q"],
                "G": ["C", "B", "?", "M"],
                "U": ["A", "I"],
            },
            "Wobble": {
                "A": ["A", "U", "&", "N", "3", "1", "P", "~", "I"],
                "C": ["G", "#", "A", "U", "I"],
                "G": [
                    "C",
                    "B",
                    "?",
                    "A",
                    "U",
                    "&",
                    "N",
                    "3",
                    "1",
                    "P",
                    "~",
                    "V",
                    "S",
                    ")",
                    "{",
                ],
                "U": ["A", "G", "U", "I", "#", "Q", "V"],
            },
            "Pairing Rules": {
                "Near-Cognate": {"base-level": [["Wo", "WC", "X"], ["X", "WC", "Wo"]]}
            },
        }
        with tempfile.NamedTemporaryFile(delete=False) as fp:
            with open(fp.name, "w", encoding="utf-8") as json_file:
                json.dump(pairing_rules_string, json_file)
            df = concentrations_generator.make_concentrations(
                concentrations_generator.make_matrix(
                    tRNAs, codons, verbose=False, settings_file_name=fp.name
                ),
                tRNAs,
                codons,
            )
            concentrations_generator.plot_matrix(
                concentrations_generator.make_matrix(tRNAs, codons, verbose=False),
                tRNAs,
                codons,
            )
            df.to_csv("/home/heday/tmp/a.csv", index=False)

    def test_print_codon_anticodon_pairings(self):
        tRNAs = pd.read_csv("/home/heday/Projects/R_concentrations/data/tRNAs.csv") # pylint: disable=invalid-name
        codons = pd.read_csv("/home/heday/Projects/R_concentrations/data/codons.csv")
        concentrations_generator.print_codon_anticodon_pairings(
            concentrations_generator.make_matrix(tRNAs, codons, verbose=False),
            tRNAs,
            codons,
        )

    def test_identify_basepairing(self):
        assert concentrations_generator.identify_3_base_pairing(
            "AAA", "AAA", basepairing_rules=self.basepairing_rules
        ) == ["Wo", "Wo", "Wo"]
        assert concentrations_generator.identify_3_base_pairing(
            "ACG", "CGU", basepairing_rules=self.basepairing_rules
        ) == ["WC", "WC", "WC"]
        assert concentrations_generator.identify_3_base_pairing(
            "ACG", "CG3", basepairing_rules=self.basepairing_rules
        ) == ["WC", "WC", "WC"]

    def test_is_near_cognate_cognate(self):
        assert concentrations_generator.is_near_cognate(
            ["Wo", "WC", "Wo"],
            self.basepairing_rules["Pairing Rules"]["Near-Cognate"]["base-level"],
        )
        assert concentrations_generator.is_near_cognate(
            ["WC", "WC", "WC"],
            self.basepairing_rules["Pairing Rules"]["Near-Cognate"]["base-level"],
        )
        assert concentrations_generator.is_near_cognate(
            ["WC", "WC", "Wo"],
            self.basepairing_rules["Pairing Rules"]["Near-Cognate"]["base-level"],
        )
        assert concentrations_generator.is_near_cognate(
            ["Wo", "WC", "WC"],
            self.basepairing_rules["Pairing Rules"]["Near-Cognate"]["base-level"],
        )
        # new rule
        rule = [["Wo", "WC", "X"]]
        assert concentrations_generator.is_near_cognate(["Wo", "WC", "Wo"], rule)
        assert concentrations_generator.is_near_cognate(["Wo", "WC", "WC"], rule)
        assert concentrations_generator.is_near_cognate(["Wo", "WC", "X"], rule)
        assert concentrations_generator.is_near_cognate(["WC", "WC", "Wo"], rule)
        assert concentrations_generator.is_near_cognate(["WC", "WC", "WC"], rule)
        assert concentrations_generator.is_near_cognate(["WC", "WC", "X"], rule)
        # new rule
        rule = [["X", "WC", "Wo"]]
        assert concentrations_generator.is_near_cognate(["Wo", "WC", "Wo"], rule)
        assert concentrations_generator.is_near_cognate(["WC", "WC", "Wo"], rule)
        assert concentrations_generator.is_near_cognate(["X", "WC", "Wo"], rule)
        assert concentrations_generator.is_near_cognate(["Wo", "WC", "WC"], rule)
        assert concentrations_generator.is_near_cognate(["WC", "WC", "WC"], rule)
        assert concentrations_generator.is_near_cognate(["X", "WC", "WC"], rule)

    def test(self):
        assert concentrations_generator.identify_3_base_pairing(
            "CCA", "AGG", basepairing_rules=self.basepairing_rules
        ) == ["WC", "WC", "Wo"]
        assert concentrations_generator.is_near_cognate(
            ["WC", "WC", "Wo"], [["X", "WC", "Wo"]]
        )
        assert concentrations_generator.is_near_cognate(
            ["WC", "WC", "Wo"], [["Wo", "WC", "X"]]
        )
