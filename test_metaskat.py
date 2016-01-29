"""Tests the run_metaskat script."""


import os
import random
import shutil
import logging
import unittest
import collections
from tempfile import mkdtemp

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
