"""Tests the run_metaskat script."""


import os
import sys
import random
import shutil
import logging
import unittest
import collections
from tempfile import mkdtemp

import numpy as np
import pandas as pd

try:
    from unittest.mock import patch, Mock, call
except:
    from mock import patch, Mock, call


# Mocking rpy2 and importing our module to test
mock = Mock()
modules = {
    "rpy2": mock,
    "rpy2.robjects": mock.module,
    "rpy2.robjects.numpy2ri": mock.module,
    "rpy2.robjects.packages": mock.module,
}
with patch.dict("sys.modules", modules):
    import run_metaskat as metaskat


class TestCheckArgs(unittest.TestCase):
    """Tests the 'check_args' function."""
    def setUp(self):
        """Setup the tests."""
        self.tmp_dir = mkdtemp(prefix="metaskat_test_")

        # Creating the configuration file
        self.conf_fn = os.path.join(self.tmp_dir, "conf")
        with open(self.conf_fn, "w") as f:
            pass

        # Creating the gene file
        self.gene_fn = os.path.join(self.tmp_dir, "gene")
        with open(self.gene_fn, "w") as f:
            pass

        # A dummy object (to simulate argparse.Namespace)
        class Dummy(object):
            pass
        self.dummy = Dummy
        self.dummy.yaml_file = self.conf_fn
        self.dummy.gene_list = self.gene_fn
        self.dummy.mac = 4

    def tearDown(self):
        """Finishes the tests."""
        # Cleaning the temporary directory
        shutil.rmtree(self.tmp_dir)

    def _my_compatibility_assertLogs(self, logger=None, level=None):
        """Compatibility 'assertLogs' function for Python < 3.4."""
        if hasattr(self, "assertLogs"):
            return self.assertLogs(logger, level)

        else:
            return AssertLogsContext_Compatibility(self, logger, level)

    def test_normal_functionality(self):
        """Tests when everything is fine."""
        # Calling the function
        metaskat.check_args(self.dummy)

    def test_missing_gene_file(self):
        """Test when the gene file is missing."""
        # Removing the gene file
        os.remove(self.gene_fn)
        with self._my_compatibility_assertLogs(level="CRITICAL") as cm_logs:
            with self.assertRaises(SystemExit) as cm:
                metaskat.check_args(self.dummy)

        # Checking the return code
        self.assertNotEqual(0, cm.exception.code)
        self.assertEqual(
            ["CRITICAL:MetaSKAT Pipeline:{}: no such file".format(
                self.gene_fn
            )],
            cm_logs.output,
        )

    def test_missing_conf_file(self):
        """Test when the conf file is missing."""
        # Removing the gene file
        os.remove(self.conf_fn)
        with self._my_compatibility_assertLogs(level="CRITICAL") as cm_logs:
            with self.assertRaises(SystemExit) as cm:
                metaskat.check_args(self.dummy)

        # Checking the return code
        self.assertNotEqual(0, cm.exception.code)
        self.assertEqual(
            ["CRITICAL:MetaSKAT Pipeline:{}: no such file".format(
                self.conf_fn
            )],
            cm_logs.output,
        )


class TestRead_Conf(unittest.TestCase):
    """Tests the 'read_cohort_configuration' function."""
    def setUp(self):
        """Setup the tests."""
        self.tmp_dir = mkdtemp(prefix="metaskat_test_")

        # Creating the dummy genotype file
        self.prefix = os.path.join(self.tmp_dir, "prefix")
        for ext in (".bed", ".bim", ".fam"):
            with open(self.prefix + ext, "w") as f:
                pass

        # Creating the dummy phenotype file
        self.pheno_fn = os.path.join(self.tmp_dir, "pheno")
        with open(self.pheno_fn, "w") as f:
            pass

        # Writing the configuration
        content = (
            "project_1:\n"
            "    prefix: {prefix}\n"
            "    phenotype_file: {pheno_fn}\n"
            "    phenotype: pheno\n"
            "    covariates:\n"
            "        - covar_1\n"
            "        - covar_2\n"
        ).format(prefix=self.prefix, pheno_fn=self.pheno_fn)
        self.conf_fn = os.path.join(self.tmp_dir, "conf.yaml")
        with open(self.conf_fn, "w") as f:
            f.write(content)

    def tearDown(self):
        """Finishes the tests."""
        # Cleaning the temporary directory
        shutil.rmtree(self.tmp_dir)

    def _my_compatibility_assertLogs(self, logger=None, level=None):
        """Compatibility 'assertLogs' function for Python < 3.4."""
        if hasattr(self, "assertLogs"):
            return self.assertLogs(logger, level)

        else:
            return AssertLogsContext_Compatibility(self, logger, level)

    def test_normal_functionality(self):
        """Tests when everything is fine."""
        # Testing the function
        metaskat.read_cohort_configuration(self.conf_fn)

    def test_missing_pheno_file(self):
        """Tests when the phenotype file is missing."""
        # Deleting the phenotype file
        os.remove(self.pheno_fn)

        # Testing the function
        with self._my_compatibility_assertLogs(level="CRITICAL") as cm_logs:
            with self.assertRaises(SystemExit) as cm:
                metaskat.read_cohort_configuration(self.conf_fn)

        # Checking the return code
        self.assertNotEqual(0, cm.exception.code)
        self.assertEqual(
            ["CRITICAL:MetaSKAT Pipeline:{}: no such file".format(
                self.pheno_fn
            )],
            cm_logs.output,
        )

    def test_missing_geno_file(self):
        """Tests when the genotype file is missing."""
        # Deleting the phenotype file
        ext_to_remove = random.choice((".bed", ".bim", ".fam"))
        os.remove(self.prefix + ext_to_remove)

        # Testing the function
        with self._my_compatibility_assertLogs(level="CRITICAL") as cm_logs:
            with self.assertRaises(SystemExit) as cm:
                metaskat.read_cohort_configuration(self.conf_fn)

        # Checking the return code
        self.assertNotEqual(0, cm.exception.code)
        self.assertEqual(
            ["CRITICAL:MetaSKAT Pipeline:{}: no such file".format(
                self.prefix + ext_to_remove
            )],
            cm_logs.output,
        )

    def test_missing_kw(self):
        """Tests when an important keyword is missing."""
        # Writing a new configuration
        content = (
            "project_1:\n"
            "    phenotype_file: {pheno_fn}\n"
            "    phenotype: pheno\n"
            "    covariates:\n"
            "        - covar_1\n"
            "        - covar_2\n"
        ).format(pheno_fn=self.pheno_fn)
        self.conf_fn = os.path.join(self.tmp_dir, "conf.yaml")
        with open(self.conf_fn, "w") as f:
            f.write(content)

        # Testing the function
        with self._my_compatibility_assertLogs(level="CRITICAL") as cm_logs:
            with self.assertRaises(SystemExit) as cm:
                metaskat.read_cohort_configuration(self.conf_fn)

        # Checking the return code
        self.assertNotEqual(0, cm.exception.code)
        self.assertEqual(
            ["CRITICAL:MetaSKAT Pipeline:{}: missing keyword '{}'".format(
                "project_1",
                "prefix",
            )],
            cm_logs.output,
        )


class TestExtractPlink(unittest.TestCase):
    """Tests the 'extract_plink' function."""
    def setUp(self):
        """Setup the tests."""
        self.tmp_dir = mkdtemp(prefix="metaskat_test_")

    def tearDown(self):
        """Finishes the tests."""
        # Cleaning the temporary directory
        shutil.rmtree(self.tmp_dir)

    def _my_compatibility_assertLogs(self, logger=None, level=None):
        """Compatibility 'assertLogs' function for Python < 3.4."""
        if hasattr(self, "assertLogs"):
            return self.assertLogs(logger, level)

        else:
            return AssertLogsContext_Compatibility(self, logger, level)

    @patch.object(metaskat, "execute_command")
    def test_normal_functionality(self, mocked):
        """Tests when everything is fine."""
        # The required data
        samples = (
            ("f1", "s1"),
            ("f1", "s2"),
            ("f2", "s2"),
            ("f3", "s3"),
            ("f4", "s4"),
            ("f5", "s5"),
            ("f6", "s6"),
            ("f7", "s7"),
        )
        i_prefix = os.path.join(self.tmp_dir, "i_prefix")
        o_prefix = os.path.join(self.tmp_dir, "o_prefix")

        # Executing the function
        metaskat.extract_plink(i_prefix, o_prefix, samples)

        # Checking the file sample to extracts is the same
        self.assertTrue(os.path.isfile(o_prefix + ".to_extract"))
        with open(o_prefix + ".to_extract", "r") as f:
            expected = f.read()
        self.assertEqual(
            "f1 s1\nf1 s2\nf2 s2\nf3 s3\nf4 s4\nf5 s5\nf6 s6\nf7 s7\n",
            expected,
        )

        # Testing the command that should have been executed is the right
        self.assertTrue(mocked.called)
        mocked.assert_called_once_with(
            ["plink", "--noweb", "--bfile", i_prefix, "--keep",
             o_prefix + ".to_extract", "--make-bed", "--out", o_prefix]
        )


class TestGetAnalysisData(unittest.TestCase):
    """Tests the 'get_analysis_data' function."""
    def setUp(self):
        """Setup the tests."""
        self.tmp_dir = mkdtemp(prefix="metaskat_test_")

        # create a fam
        self.fam = pd.DataFrame(
            [("f1", "i1", 0, 0, 1, -9),
             ("f2", "i2", 0, 0, 1, -9),
             ("f3", "i3", 0, 0, 1, -9),
             ("f4", "i4", 0, 0, 1, -9),
             ("f5", "i5", 0, 0, 1, -9),
             ("f6", "i6", 0, 0, 1, -9),
             ("f7", "i7", 0, 0, 1, -9),
             ("f8", "i8", 0, 0, 1, -9),
             ("f9", "i9", 0, 0, 1, -9),
             ("f10", "i10", 0, 0, 1, -9)],
            columns=["FID", "IID", "FAT", "MOT", "SEX", "PHENO"]
        )

        # create a pheno-covariates file
        self.cov = pd.DataFrame(
            [("f1", "i1", 25, 1, 0.1),
             ("f2", "i2", 25, 1, 0.2),
             ("f3", "i3", 25, 1, 0.3),
             ("f4", "i4", 25, 1, 0.4),
             ("f5", "i5", 25, 1, 0.5),
             ("f6", "i6", 25, 1, 0.5),
             ("f7", "i7", 25, 1, 0.4),
             ("f8", "i8", 25, 1, 0.3),
             ("f9", "i9", 25, 1, 0.2),
             ("f10", "i10", 25, 1, 0.1)],
            columns=["FID", "IID", "AGE", "SEX", "PHENO"]
        )

    def tearDown(self):
        """Finishes the tests."""
        # Cleaning the temporary directory
        shutil.rmtree(self.tmp_dir)

    def _my_compatibility_assertLogs(self, logger=None, level=None):
        """Compatibility 'assertLogs' function for Python < 3.4."""
        if hasattr(self, "assertLogs"):
            return self.assertLogs(logger, level)

        else:
            return AssertLogsContext_Compatibility(self, logger, level)

    def test_normal_functionality(self):
        """Tests when everything is fine."""
        prefix = os.path.join(self.tmp_dir, "test")
        fam = self.fam.reindex(np.random.permutation(self.fam.index))
        fam.to_csv(prefix + ".fam", sep=" ", index=False, header=False)

        phenofile = os.path.join(self.tmp_dir, "pheno.txt")
        cov = self.cov.reindex(np.random.permutation(self.cov.index))
        cov.to_csv(phenofile, sep="\t", index=False)

        obs_prefix, obs_pheno = metaskat.get_analysis_data(
            prefix, "PHENO",
            ["SEX", "AGE"],
            "FID", "IID",
            phenofile, os.path.join(self.tmp_dir, "foo"))

        self.assertEqual(prefix, obs_prefix)
        self.assertTrue(isinstance(obs_pheno, pd.DataFrame))
        self.assertEqual((10, 3), obs_pheno.shape)
        self.assertEqual(["PHENO", "SEX", "AGE"], list(obs_pheno.columns))
        self.assertEqual(fam.set_index(["FID", "IID"]).index.tolist(),
                         obs_pheno.index.tolist())

    def test_subset_pheno(self):
        """Tests when there are more sample in phenofile than famfile."""
        prefix = os.path.join(self.tmp_dir, "test")
        fam = self.fam.reindex(np.random.permutation(self.fam.index)).head(n=8)
        fam.to_csv(prefix + ".fam", sep=" ", index=False, header=False)

        phenofile = os.path.join(self.tmp_dir, "pheno.txt")
        cov = self.cov.reindex(np.random.permutation(self.cov.index))
        cov.to_csv(phenofile, sep="\t", index=False)

        obs_prefix, obs_pheno = metaskat.get_analysis_data(
            prefix, "PHENO",
            ["SEX", "AGE"],
            "FID", "IID",
            phenofile, os.path.join(self.tmp_dir, "foo"))

        self.assertEqual(prefix, obs_prefix)
        self.assertTrue(isinstance(obs_pheno, pd.DataFrame))
        self.assertEqual((8, 3), obs_pheno.shape)
        self.assertEqual(["PHENO", "SEX", "AGE"], list(obs_pheno.columns))
        self.assertEqual(fam.set_index(["FID", "IID"]).index.tolist(),
                         obs_pheno.index.tolist())

    def test_subset_fam(self):
        """Tests when there are more sample in famfile than phenofile."""
        prefix = os.path.join(self.tmp_dir, "test")
        fam = self.fam.reindex(np.random.permutation(self.fam.index))
        fam.to_csv(prefix + ".fam", sep=" ", index=False, header=False)

        phenofile = os.path.join(self.tmp_dir, "pheno.txt")
        cov = self.cov.reindex(np.random.permutation(self.cov.index)).head(n=8)
        cov.to_csv(phenofile, sep="\t", index=False)

        samples = fam.set_index(["FID", "IID"]).index.intersection(
            cov.set_index(["FID", "IID"]).index
        )

        new_prefix = prefix + "_new"
        fam = fam.loc[fam.IID.isin(cov.IID), ]
        fam = fam.reindex(np.random.permutation(fam.index))
        fam.to_csv(new_prefix + ".fam", sep=" ", index=False, header=False)

        with patch.object(metaskat, "extract_plink",
                          return_value=new_prefix) as mock:
            obs_prefix, obs_pheno = metaskat.get_analysis_data(
                prefix, "PHENO",
                ["SEX", "AGE"],
                "FID", "IID",
                phenofile, os.path.join(self.tmp_dir, "foo"))

        self.assertEqual(new_prefix, obs_prefix)
        self.assertTrue(isinstance(obs_pheno, pd.DataFrame))
        self.assertEqual((8, 3), obs_pheno.shape)
        self.assertEqual(["PHENO", "SEX", "AGE"], list(obs_pheno.columns))
        self.assertEqual(fam.set_index(["FID", "IID"]).index.tolist(),
                         obs_pheno.index.tolist())

        self.assertTrue(mock.called)
        mock.assert_called_once_with(
            prefix, os.path.join(os.path.join(self.tmp_dir, "foo"), "subset"),
            samples.tolist()
        )

    def test_subset_fam_pheno(self):
        """Tests when there are different sample in famfile and phenofile."""
        fam_sample = {"i{}".format(i) for i in range(1, 5)}
        cov_sample = {"i{}".format(i) for i in range(3, 7)}

        prefix = os.path.join(self.tmp_dir, "test")
        fam = self.fam.reindex(np.random.permutation(self.fam.index))
        fam = fam[fam.IID.isin(fam_sample)]
        fam.to_csv(prefix + ".fam", sep=" ", index=False, header=False)

        phenofile = os.path.join(self.tmp_dir, "pheno.txt")
        cov = self.cov.reindex(np.random.permutation(self.cov.index))
        cov = cov[cov.IID.isin(cov_sample)]
        cov.to_csv(phenofile, sep="\t", index=False)

        samples = fam.set_index(["FID", "IID"]).index.intersection(
            cov.set_index(["FID", "IID"]).index
        )

        new_prefix = prefix + "_new"
        fam = fam.loc[fam.IID.isin({"i3", "i4"}), ]
        fam = fam.reindex(np.random.permutation(fam.index))
        fam.to_csv(new_prefix + ".fam", sep=" ", index=False, header=False)

        with patch.object(metaskat, "extract_plink",
                          return_value=new_prefix) as mock:
            obs_prefix, obs_pheno = metaskat.get_analysis_data(
                prefix, "PHENO",
                ["SEX", "AGE"],
                "FID", "IID",
                phenofile, os.path.join(self.tmp_dir, "foo"))

        self.assertEqual(new_prefix, obs_prefix)
        self.assertTrue(isinstance(obs_pheno, pd.DataFrame))
        self.assertEqual((len(samples.tolist()), 3), obs_pheno.shape)
        self.assertEqual(["PHENO", "SEX", "AGE"], list(obs_pheno.columns))
        self.assertEqual(fam.set_index(["FID", "IID"]).index.tolist(),
                         obs_pheno.index.tolist())

        self.assertTrue(mock.called)
        mock.assert_called_once_with(
            prefix, os.path.join(os.path.join(self.tmp_dir, "foo"), "subset"),
            samples.tolist()
        )

    def test_no_sample_fam_pheno(self):
        """Tests when there is no sample in common in famfile and phenofile."""
        fam_sample = {"i{}".format(i) for i in range(1, 4)}
        cov_sample = {"i{}".format(i) for i in range(5, 7)}

        prefix = os.path.join(self.tmp_dir, "test")
        fam = self.fam.reindex(np.random.permutation(self.fam.index))
        fam = fam[fam.IID.isin(fam_sample)]
        fam.to_csv(prefix + ".fam", sep=" ", index=False, header=False)

        phenofile = os.path.join(self.tmp_dir, "pheno.txt")
        cov = self.cov.reindex(np.random.permutation(self.cov.index))
        cov = cov[cov.IID.isin(cov_sample)]
        cov.to_csv(phenofile, sep="\t", index=False)

        # Testing the function
        with self._my_compatibility_assertLogs(level="CRITICAL") as cm_logs:
            with self.assertRaises(SystemExit) as cm:
                metaskat.get_analysis_data(prefix, "PHENO", ["SEX", "AGE"],
                                           "FID", "IID", phenofile,
                                           os.path.join(self.tmp_dir, "foo"))

        # Checking the return code
        self.assertNotEqual(0, cm.exception.code)
        self.assertEqual(
            ["CRITICAL:MetaSKAT Pipeline:{} vs {}: no sample in common".format(
                prefix + ".fam",
                phenofile,
            )],
            cm_logs.output,
        )

    def test_missing_column(self):
        """Tests when there is a missing column in the phenotype file."""
        prefix = os.path.join(self.tmp_dir, "test")
        fam = self.fam.reindex(np.random.permutation(self.fam.index))
        fam.to_csv(prefix + ".fam", sep=" ", index=False, header=False)

        phenofile = os.path.join(self.tmp_dir, "pheno.txt")
        cov = self.cov.reindex(np.random.permutation(self.cov.index))
        cov.to_csv(phenofile, sep="\t", index=False)

        # Testing the function
        with self._my_compatibility_assertLogs(level="CRITICAL") as cm_logs:
            with self.assertRaises(SystemExit) as cm:
                metaskat.get_analysis_data(prefix, "PHENO", ["SEX", "A"],
                                           "FID", "IID", phenofile,
                                           os.path.join(self.tmp_dir, "foo"))

        # Checking the return code
        self.assertNotEqual(0, cm.exception.code)
        self.assertEqual(
            ["CRITICAL:MetaSKAT Pipeline:{}: missing column '{}'".format(
                phenofile,
                "A",
            )],
            cm_logs.output,
        )

    def test_na_pheno(self):
        """Tests when there are na values in phenofile."""
        prefix = os.path.join(self.tmp_dir, "test")
        fam = self.fam.reindex(np.random.permutation(self.fam.index))
        fam.to_csv(prefix + ".fam", sep=" ", index=False, header=False)

        phenofile = os.path.join(self.tmp_dir, "pheno.txt")
        cov = self.cov.reindex(np.random.permutation(self.cov.index))
        cov.loc[cov.IID.isin({"i3", "i5"}), "PHENO"] = np.nan
        cov.to_csv(phenofile, sep="\t", index=False)
        cov = cov.dropna()

        samples = fam.set_index(["FID", "IID"]).index.intersection(
            cov.set_index(["FID", "IID"]).index
        )

        new_prefix = prefix + "_new"
        fam = fam.loc[fam.IID.isin(cov.IID), ]
        fam = fam.reindex(np.random.permutation(fam.index))
        fam.to_csv(new_prefix + ".fam", sep=" ", index=False, header=False)

        with patch.object(metaskat, "extract_plink",
                          return_value=new_prefix) as mock:
            obs_prefix, obs_pheno = metaskat.get_analysis_data(
                prefix, "PHENO",
                ["SEX", "AGE"],
                "FID", "IID",
                phenofile, os.path.join(self.tmp_dir, "foo"))

        self.assertEqual(new_prefix, obs_prefix)
        self.assertTrue(isinstance(obs_pheno, pd.DataFrame))
        self.assertEqual((8, 3), obs_pheno.shape)
        self.assertEqual(["PHENO", "SEX", "AGE"], list(obs_pheno.columns))
        self.assertEqual(fam.set_index(["FID", "IID"]).index.tolist(),
                         obs_pheno.index.tolist())

        self.assertTrue(mock.called)
        mock.assert_called_once_with(
            prefix, os.path.join(os.path.join(self.tmp_dir, "foo"), "subset"),
            samples.tolist()
        )


class TestExecuteSKAT(unittest.TestCase):
    """Tests the 'execute_skat' function."""
    def setUp(self):
        """Setup the tests."""
        self.tmp_dir = mkdtemp(prefix="metaskat_test_")

        # create a pheno-covariates file
        self.cov = pd.DataFrame(
            [("f1", "i1", 25, 1, 33, 0.1),
             ("f2", "i2", 25, 1, 34, 0.2),
             ("f3", "i3", 25, 1, 35, 0.3),
             ("f4", "i4", 25, 1, 30, 0.4),
             ("f5", "i5", 25, 1, 36, 0.5),
             ("f6", "i6", 25, 1, 38, 0.5),
             ("f7", "i7", 25, 1, 23, 0.4),
             ("f8", "i8", 25, 1, 13, 0.3),
             ("f9", "i9", 25, 1, 31, 0.2),
             ("f10", "i10", 25, 1, 29, 0.1)],
            columns=["FID", "IID", "AGE", "SEX", "FOO", "PHENO"]
        ).set_index(["FID", "IID"], verify_integrity=True)

    def tearDown(self):
        """Finishes the tests."""
        # Cleaning the temporary directory
        shutil.rmtree(self.tmp_dir)

    def _my_compatibility_assertLogs(self, logger=None, level=None):
        """Compatibility 'assertLogs' function for Python < 3.4."""
        if hasattr(self, "assertLogs"):
            return self.assertLogs(logger, level)

        else:
            return AssertLogsContext_Compatibility(self, logger, level)

    def test_normal_functionality(self):
        """Tests the normal functionality of the 'execute_skat' function."""
        # The cohort information
        cohort_info = collections.OrderedDict()
        cohort_info["cohort_1"] = dict(
                prefix=os.path.join(self.tmp_dir, "prefix_1"),
                phenotype="PHENO",
                covariates=["SEX", "AGE"],
                phenotype_file=os.path.join(self.tmp_dir, "pheno"),
        )
        cohort_info["cohort_2"] = dict(
                prefix=os.path.join(self.tmp_dir, "prefix_2"),
                phenotype="PHENO",
                covariates=["SEX", "FOO"],
                phenotype_file=os.path.join(self.tmp_dir, "pheno_2"),
                family_id="fam_id",
                individual_id="ind_id",
        )

        # Creating the mocks and executing the function
        prefix = os.path.join(self.tmp_dir, "new_prefix")

        class FormulaEnv(object):
            environment = {}

        formula_env = FormulaEnv()

        # What the get analysis will return
        analysis_data_return = [
            (prefix + "_1", self.cov[["SEX", "AGE", "PHENO"]]),
            (prefix + "_2", self.cov[["SEX", "FOO", "PHENO"]]),
        ]

        # Executing the function
        with patch.object(metaskat, "get_analysis_data",
                          side_effect=analysis_data_return) as mock_get_data, \
             patch.object(metaskat, "write_valid_segments",
                          return_value=None) as mock_write_segments, \
             patch.object(metaskat, "rformula",
                          return_value=formula_env) as mock_rformula, \
             patch.object(metaskat.skat, "SKAT_Null_Model",
                          return_value="dummy_model") as mock_skat, \
             patch.object(metaskat.metaskat, "Generate_Meta_Files",
                          return_value=None) as mock_metaskat, \
             patch.object(metaskat.robjects, "r",
                          return_value="NULL") as mock_rob:
            metaskat.execute_skat(
                cohort_info,
                os.path.join(self.tmp_dir, "genes"),
                os.path.join(self.tmp_dir, "o_prefix"),
                self.tmp_dir,
                mac=6,
            )

        # Testing the first mock was called with the right argument
        self.assertEqual(2, mock_get_data.call_count)
        mock_get_data.assert_has_calls([
            call(
                plink_prefix=os.path.join(self.tmp_dir, "prefix_1"),
                pheno="PHENO",
                covariates=["SEX", "AGE"],
                pheno_fn=os.path.join(self.tmp_dir, "pheno"),
                fid="FID",
                iid="IID",
                tmp_dir=os.path.join(self.tmp_dir, "cohort_1"),
            ),
            call(
                plink_prefix=os.path.join(self.tmp_dir, "prefix_2"),
                pheno="PHENO",
                covariates=["SEX", "FOO"],
                pheno_fn=os.path.join(self.tmp_dir, "pheno_2"),
                fid="fam_id",
                iid="ind_id",
                tmp_dir=os.path.join(self.tmp_dir, "cohort_2"),
            ),
        ])

        # Testing the third mock was called with the right argument
        self.assertEqual(2, mock_rformula.call_count)
        mock_rformula.assert_any_call(
            "cohort_1.PHENO ~ cohort_1.SEX + cohort_1.AGE"
        )
        mock_rformula.assert_any_call(
            "cohort_2.PHENO ~ cohort_2.SEX + cohort_2.FOO"
        )

        # Testing what we have in the formula environment is OK
        self.assertEqual({"cohort_1.PHENO", "cohort_1.AGE", "cohort_1.SEX",
                          "cohort_2.PHENO", "cohort_2.FOO", "cohort_2.SEX"},
                         set(formula_env.environment.keys()))
        self.assertEqual(list(self.cov.PHENO),
                         list(formula_env.environment["cohort_1.PHENO"]))
        self.assertEqual(list(self.cov.AGE),
                         list(formula_env.environment["cohort_1.AGE"]))
        self.assertEqual(list(self.cov.SEX),
                         list(formula_env.environment["cohort_1.SEX"]))
        self.assertEqual(list(self.cov.PHENO),
                         list(formula_env.environment["cohort_2.PHENO"]))
        self.assertEqual(list(self.cov.FOO),
                         list(formula_env.environment["cohort_2.FOO"]))
        self.assertEqual(list(self.cov.SEX),
                         list(formula_env.environment["cohort_2.SEX"]))

        # Checking what SKAT was called correctly
        self.assertEqual(2, mock_skat.call_count)
        mock_skat.assert_called_with(
            formula_env,
            type="C",
        )

        # Checking that MetaSKAT was called correctly
        self.assertEqual(2, mock_metaskat.call_count)
        mock_metaskat.assert_has_calls([
            call(
                "dummy_model",
                prefix + "_1.bed",
                prefix + "_1.bim",
                os.path.join(self.tmp_dir, "o_prefix", "valid_segments.txt"),
                os.path.join(self.tmp_dir, "o_prefix", "cohort_1.MSSD"),
                os.path.join(self.tmp_dir, "o_prefix", "cohort_1.MInfo"),
                10,
                File_Permu="NULL",
            ),
            call(
                "dummy_model",
                prefix + "_2.bed",
                prefix + "_2.bim",
                os.path.join(self.tmp_dir, "o_prefix", "valid_segments.txt"),
                os.path.join(self.tmp_dir, "o_prefix", "cohort_2.MSSD"),
                os.path.join(self.tmp_dir, "o_prefix", "cohort_2.MInfo"),
                10,
                File_Permu="NULL",
            ),
        ])

    def test_with_metaskat_error(self):
        """Tests when there is a problem with MetaSKAT."""
        # The cohort information
        cohort_info = dict(
            cohort_1=dict(
                prefix=os.path.join(self.tmp_dir, "prefix_1"),
                phenotype="PHENO",
                covariates=["SEX", "AGE"],
                phenotype_file=os.path.join(self.tmp_dir, "pheno"),
            ),
            cohort_2=dict(
                prefix=os.path.join(self.tmp_dir, "prefix_2"),
                phenotype="PHENO",
                covariates=["SEX", "FOO"],
                phenotype_file=os.path.join(self.tmp_dir, "pheno_2"),
                family_id="fam_id",
                individual_id="ind_id",
            ),
        )

        # Creating the mocks and executing the function
        prefix = os.path.join(self.tmp_dir, "new_prefix")

        class FormulaEnv(object):
            environment = {}

        formula_env = FormulaEnv()

        with patch.object(metaskat, "get_analysis_data",
                          return_value=(prefix, self.cov)) as mock_get_data, \
             patch.object(metaskat, "write_valid_segments",
                          return_value=None) as mock_write_segments, \
             patch.object(metaskat, "rformula",
                          return_value=formula_env) as mock_rformula, \
             patch.object(metaskat.skat, "SKAT_Null_Model",
                          return_value="dummy_model") as mock_skat, \
             patch.object(metaskat.metaskat, "Generate_Meta_Files",
                          return_value=None) as mock_metaskat, \
             patch.object(metaskat.robjects, "r",
                          side_effect=Exception("An error")) as mock_robj, \
             self._my_compatibility_assertLogs(level="CRITICAL") as cm_logs, \
             self.assertRaises(SystemExit) as cm:
            metaskat.execute_skat(
                cohort_info,
                os.path.join(self.tmp_dir, "genes"),
                os.path.join(self.tmp_dir, "o_prefix"),
                self.tmp_dir,
            )

        # Checking the log
        self.assertNotEqual(0, cm.exception.code)
        self.assertEqual(
            ["CRITICAL:MetaSKAT Pipeline:problem with SKAT\nAn error"],
            cm_logs.output,
        )


class TestExecuteMetaAnalysis(unittest.TestCase):
    """Tests the 'execute_meta_analysis' function."""
    def setUp(self):
        """Setup the tests."""
        self.tmp_dir = mkdtemp(prefix="metaskat_test_")

    def tearDown(self):
        """Finishes the tests."""
        # Cleaning the temporary directory
        shutil.rmtree(self.tmp_dir)

    def _my_compatibility_assertLogs(self, logger=None, level=None):
        """Compatibility 'assertLogs' function for Python < 3.4."""
        if hasattr(self, "assertLogs"):
            return self.assertLogs(logger, level)

        else:
            return AssertLogsContext_Compatibility(self, logger, level)

    def test_normal_functionality(self):
        """Tests the normal functionality of the function."""
        # The cohort information
        cohort_info = dict(
            cohort_1=dict(
                prefix=os.path.join(self.tmp_dir, "prefix_1"),
                phenotype="PHENO",
                covariates=["SEX", "AGE"],
                phenotype_file=os.path.join(self.tmp_dir, "pheno"),
            ),
            cohort_2=dict(
                prefix=os.path.join(self.tmp_dir, "prefix_2"),
                phenotype="PHENO",
                covariates=["SEX", "AGE"],
                phenotype_file=os.path.join(self.tmp_dir, "pheno"),
            ),
            cohort_3=dict(
                prefix=os.path.join(self.tmp_dir, "prefix_3"),
                phenotype="PHENO",
                covariates=["SEX", "AGE"],
                phenotype_file=os.path.join(self.tmp_dir, "pheno"),
            ),
        )

        # Creating a dummy cohort information
        class DummyMetaCohortInfo(object):
            def __init__(self, values, names):
                self.values = values
                self.names = names

            def __getitem__(self, k):
                return self.values[k]

        # A mock result
        meta_result = Mock()

        meta_cohort = DummyMetaCohortInfo(
            values=[[
                DummyMetaCohortInfo(values=[meta_result], names=["Info"]),
                DummyMetaCohortInfo(values=[meta_result], names=["Info"]),
                DummyMetaCohortInfo(values=[meta_result], names=["Info"]),
            ]],
            names=["EachInfo"],
        )

        # Creating a dummy value for the final cohort information
        with patch.object(metaskat.metaskat, "MetaSKAT_MSSD_ALL",
                          return_value=meta_result) as mock_metaskat, \
             patch.object(metaskat.metaskat, "Open_MSSD_File_2Read",
                          return_value=meta_cohort) as mock_file2read, \
             patch.object(metaskat.np, "array",
                          side_effect=list) as mock_np_array, \
             patch.object(metaskat.robjects, "r",
                          return_value="NULL") as mock_rob:
            metaskat.execute_meta_analysis(
                cohort_info,
                os.path.join(self.tmp_dir, "genes"),
                self.tmp_dir,
            )

        # Checking what have been called for 'Open_MSSD_File_2Read'
        self.assertEqual(1, mock_file2read.call_count)
        mock_file2read.assert_called_with(
            [os.path.join(self.tmp_dir,
                          "cohort_{}.MSSD".format(i+1)) for i in range(3)],
            [os.path.join(self.tmp_dir,
                          "cohort_{}.MInfo".format(i+1)) for i in range(3)],
        )

        # Checking what have been called for 'MetaSKAT_MSSD_ALL'
        self.assertEqual(2, mock_metaskat.call_count)
        mock_metaskat.assert_has_calls([
            call(**{
                "Cohort.Info": meta_cohort,
                "combined.weight": True,
                "weights.beta": [1, 25],
                "method": "davies",
                "r.corr": 0,
                "is.separate": False,
                "Group_Idx": "NULL",
                "MAF.cutoff": 5,
            }),
            call(**{
                "Cohort.Info": meta_cohort,
                "combined.weight": True,
                "weights.beta": [1, 25],
                "method": "davies",
                "r.corr": 0,
                "is.separate": True,
                "Group_Idx": "NULL",
                "MAF.cutoff": 5,
            }),
        ])

        # Checking what has been called for the MetaSKAT results
        self.assertEqual(5, meta_result.to_csvfile.call_count)
        meta_result.to_csvfile.assert_has_calls([
            call(os.path.join(self.tmp_dir, "metaSKAT.homo.txt"),
                 quote=False, sep="\t", row_names=False, col_names=True),
            call(os.path.join(self.tmp_dir, "metaSKAT.hetero.txt"),
                 quote=False, sep="\t", row_names=False, col_names=True),
            call(os.path.join(self.tmp_dir, "cohort_1.MetaInfo.txt"),
                 sep="\t", quote=False, row_names=False),
            call(os.path.join(self.tmp_dir, "cohort_2.MetaInfo.txt"),
                 sep="\t", quote=False, row_names=False),
            call(os.path.join(self.tmp_dir, "cohort_3.MetaInfo.txt"),
                 sep="\t", quote=False, row_names=False),
        ])


class BaseTestCaseContext_Compatibility:

    def __init__(self, test_case):
        self.test_case = test_case

    def _raiseFailure(self, standardMsg):
        msg = self.test_case._formatMessage(self.msg, standardMsg)
        raise self.test_case.failureException(msg)


_LoggingWatcher = collections.namedtuple("_LoggingWatcher",
                                         ["records", "output"])


class CapturingHandler_Compatibility(logging.Handler):
    """
    A logging handler capturing all (raw and formatted) logging output.
    """

    def __init__(self):
        logging.Handler.__init__(self)
        self.watcher = _LoggingWatcher([], [])

    def flush(self):
        pass

    def emit(self, record):
        self.watcher.records.append(record)
        msg = self.format(record)
        self.watcher.output.append(msg)


class AssertLogsContext_Compatibility(BaseTestCaseContext_Compatibility):
    """A context manager used to implement TestCase.assertLogs()."""

    LOGGING_FORMAT = "%(levelname)s:%(name)s:%(message)s"

    def __init__(self, test_case, logger_name, level):
        BaseTestCaseContext_Compatibility.__init__(self, test_case)
        self.logger_name = logger_name
        if level:
            # Python < 3.4 logging doesn't have a _nameToLevel dictionary
            nameToLevel = {
                'CRITICAL': logging.CRITICAL,
                'ERROR': logging.ERROR,
                'WARN': logging.WARNING,
                'WARNING': logging.WARNING,
                'INFO': logging.INFO,
                'DEBUG': logging.DEBUG,
                'NOTSET': logging.NOTSET,
            }
            self.level = nameToLevel.get(level, level)
        else:
            self.level = logging.INFO
        self.msg = None

    def __enter__(self):
        if isinstance(self.logger_name, logging.Logger):
            logger = self.logger = self.logger_name
        else:
            logger = self.logger = logging.getLogger(self.logger_name)
        formatter = logging.Formatter(self.LOGGING_FORMAT)
        handler = CapturingHandler_Compatibility()
        handler.setFormatter(formatter)
        self.watcher = handler.watcher
        self.old_handlers = logger.handlers[:]
        self.old_level = logger.level
        self.old_propagate = logger.propagate
        logger.handlers = [handler]
        logger.setLevel(self.level)
        logger.propagate = False
        return handler.watcher

    def __exit__(self, exc_type, exc_value, tb):
        self.logger.handlers = self.old_handlers
        self.logger.propagate = self.old_propagate
        self.logger.setLevel(self.old_level)
        if exc_type is not None:
            # let unexpected exceptions pass through
            return False
        if len(self.watcher.records) == 0:
            self._raiseFailure(
                "no logs of level {} or higher triggered on {}"
                .format(logging.getLevelName(self.level), self.logger.name))


if __name__ == "__main__":
    unittest.main()
