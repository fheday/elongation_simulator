from concentrations import concentrations_generator
import pandas as pd
import unittest   # The test framework

class Test_concentrations_generator_pairings(unittest.TestCase):

    def setUp(self):
        self.basepairing_rules ={
                                "Watson-Crick": {
                                    "A": ["U", "&", "3", "1", "~", "N", "S", ")", "{", "V", "}",  "P"],
                                    "C": ["G", "#", "W"],
                                    "G": ["C", "B"],
                                    "U": ["A"]
                                },
                                "Wobble": {
                                    "A": ["A", "I", "M", "?"],
                                    "C": ["A","U", "P", "I", "?", "Q"],
                                    "G": ["A", "U", "&", "3", "1", "~", "N", "S", ")", "{", "V", "P", "?", "M" ],
                                    "U": ["G", "#", "W", "U", "V", "P", "I", "Q"]
                                },
                                "Pairing Rules":{
                                    "Near-Cognate": {
                                        "base-level": [["Wo", "WC", "Wo"]]
                                    }
                                }
                            }

    def test_default_concentrations_generation(self):
        
        tRNAs = pd.read_csv('/home/heday/Projects/R_concentrations/data/tRNAs.csv')
        codons = pd.read_csv('/home/heday/Projects/R_concentrations/data/codons.csv')
        df = concentrations_generator.make_concentrations(concentrations_generator.make_matrix(tRNAs, codons, verbose=False), tRNAs, codons)
        concentrations_generator.plot_matrix(concentrations_generator.make_matrix(tRNAs, codons, verbose=False), tRNAs, codons)

    def test_identify_basepairing(self):
        assert(concentrations_generator.identify_3_base_pairing("AAA", "AAA", basepairing_rules=self.basepairing_rules) == ["Wo", "Wo", "Wo"])
        assert(concentrations_generator.identify_3_base_pairing("ACG", "CGU", basepairing_rules=self.basepairing_rules) == ["WC", "WC", "WC"])
        assert(concentrations_generator.identify_3_base_pairing("ACG", "CG3", basepairing_rules=self.basepairing_rules) == ["WC", "WC", "WC"])
    
    def test_is_near_cognate(self):
        assert(concentrations_generator.is_near_cognate(["Wo", "WC", "Wo"], self.basepairing_rules["Pairing Rules"]["Near-Cognate"]["base-level"]))
        assert(concentrations_generator.is_near_cognate(["WC", "WC", "WC"], self.basepairing_rules["Pairing Rules"]["Near-Cognate"]["base-level"]))
        assert(concentrations_generator.is_near_cognate(["WC", "WC", "Wo"], self.basepairing_rules["Pairing Rules"]["Near-Cognate"]["base-level"]))
        assert(concentrations_generator.is_near_cognate(["Wo", "WC", "WC"], self.basepairing_rules["Pairing Rules"]["Near-Cognate"]["base-level"]))
        # new rule
        rule = [["Wo", "WC", "X"]]
        assert(concentrations_generator.is_near_cognate(["Wo", "WC", "Wo"], rule))
        assert(concentrations_generator.is_near_cognate(["Wo", "WC", "WC"], rule))
        assert(concentrations_generator.is_near_cognate(["Wo", "WC", "X"], rule))
        assert(concentrations_generator.is_near_cognate(["WC", "WC", "Wo"], rule))
        assert(concentrations_generator.is_near_cognate(["WC", "WC", "WC"], rule))
        assert(concentrations_generator.is_near_cognate(["WC", "WC", "X"], rule))
        # new rule
        rule = [["X", "WC", "Wo"]]
        assert(concentrations_generator.is_near_cognate(["Wo", "WC", "Wo"], rule))
        assert(concentrations_generator.is_near_cognate(["WC", "WC", "Wo"], rule))
        assert(concentrations_generator.is_near_cognate(["X", "WC", "Wo"], rule))
        assert(concentrations_generator.is_near_cognate(["Wo", "WC", "WC"], rule))
        assert(concentrations_generator.is_near_cognate(["WC", "WC", "WC"], rule))
        assert(concentrations_generator.is_near_cognate(["X", "WC", "WC"], rule))
