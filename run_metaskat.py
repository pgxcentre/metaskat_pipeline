#!/usr/bin/env python3

from __future__ import print_function

import os
import sys
import logging
import argparse

import yaml
import numpy as np
import pandas as pd

import rpy2.robjects as robjects
from rpy2.robjects import Formula as rformula
import rpy2.robjects.numpy2ri as numpy2ri
from rpy2.robjects.packages import importr as rimport


# Activating numpy array to rpy2 objects automatic conversion
numpy2ri.activate()

# Loading required R libraries
metaskat = rimport("MetaSKAT")
skat = rimport("SKAT")


__copyright__ = "Copyright 2016, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "MIT"
__version__ = "0.1"


# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s %(name)s %(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("MetaSKAT Pipeline")


def main():
    """The main function."""
    args = parse_args()
    check_args(args)

    # Reading the YAML configuration
    cohort_information = read_cohort_configuration(args.yaml_file)

    # Creating the null model
    generate_meta_files(cohort_information, args.gene_list,
                        out_dir=args.output)


def generate_meta_files(cohorts, genes, out_dir):
    """Generating the SKAT meta files."""
    logger.info("Generating the meta files")

    for cohort, cohort_info in cohorts.items():
        logger.info("  - " + cohort)

        # Creating the required variables
        prefix = cohort_info["prefix"]
        pheno_name = cohort_info["phenotype"]
        pheno_filename = cohort_info["phenotype_file"]
        covariates = cohort_info["covariates"]

        # Plink binary files
        bed = prefix + ".bed"
        bim = prefix + ".bim"
        fam = prefix + ".fam"

        # Reading the FAM file
        fam = pd.read_csv(fam, sep=" ",
                          names=["fid", "iid", "father", "mother", "gender",
                                 "status"])
        fam = fam.set_index("iid", verify_integrity=True)

        # Reading the phenotype and covariates
        pheno = pd.read_csv(pheno_filename, sep="\t")
        pheno = pheno.set_index("IID", drop=False, verify_integrity=True)

        # Keeping only samples with genotypes
        pheno = pheno.loc[fam.index, ]

        # Checking we have all required information
        if (len(pheno) != len(fam)) or (not all(pheno.index == fam.index)):
            logger.critical("{}: missing sample information".format(cohort))
            sys.exit(1)

        # Checking the covariates
        for covariate in covariates:
            if covariate not in pheno.columns:
                logger.critical("{}: missing covariate '{}'".format(
                    pheno_filename,
                    covariate,
                ))
                sys.exit(1)

        # Checking the phenotype
        if pheno_name not in pheno.columns:
            logger.critical("{}: missing phenotype '{}'".format(
                pheno_filename,
                pheno_name,
            ))
            sys.exit(1)

        # Creating the formula
        formula = rformula("{} ~ {}".format(
            cohort + "." + pheno_name,
            " + ".join(["{}.{}".format(cohort, c) for c in covariates]),
        ))

        # Adding vectors to the formula environment
        env = formula.environment
        for covariate in covariates:
            env[cohort + "." + covariate] = pheno[covariate].values
        env[cohort + "." + pheno_name] = pheno[pheno_name].values

        # The null model
        model = skat.SKAT_Null_Model(formula, type="C")

        # The output files
        mssd = os.path.join(out_dir, cohort + ".MSSD")
        minfo = os.path.join(out_dir, cohort + ".MInfo")

        # Generating the meta files
        metaskat.Generate_Meta_Files(model, bed, bim, genes, mssd, minfo,
                                     fam.shape[0],
                                     File_Permu=robjects.r("NULL"))


def read_cohort_configuration(fn):
    """Reads the cohort information using YAML.

    Args:
        fn (str): the name of the YAML file.

    Returns:
        dict: a dictionary containing the cohort information.

    Some keywords are necessary for cohort description.

    +----------------+---------------------------------------------------+
    | Keyword        | Description                                       |
    +----------------+---------------------------------------------------+
    | prefix         | the prefix of the genotype dataset (Plink binary  |
    |                | format).                                          |
    +----------------+---------------------------------------------------+
    | phenotype_file | The name of the file containing the phenotype and |
    |                | covariates.                                       |
    +----------------+---------------------------------------------------+
    | phenotype      | The outcome of the analysis (phenotype).          |
    +----------------+---------------------------------------------------+
    | covariates     | A list of covariates.                             |
    +----------------+---------------------------------------------------+

    """
    # The configuration
    conf = None
    with open(fn, "r") as f:
        conf = yaml.load(f)

    req_kw = ("prefix", "phenotype_file", "phenotype", "covariates")
    for cohort, cohort_info in conf.items():
        # Checking the keywords
        for kw in req_kw:
            if kw not in conf[cohort]:
                logger.critical("{}: missing keyword '{}'".format(cohort, kw))
                sys.exit(1)

        # Checking the BED/BIM/FAM files
        for suffix in (".bed", ".bim", ".fam"):
            fn = cohort_info["prefix"] + suffix
            if not os.path.isfile(fn):
                logger.critical("{}: no such file".format(fn))
                sys.exit(1)

        # Checking the phenotype file
        pheno_file = cohort_info["phenotype_file"]
        if not os.path.isfile(pheno_file):
            logger.critical("{}: no such file".format(pheno_file))
            sys.exit(1)

    return conf


def check_args(args):
    """Checks the arguments and options."""
    if not os.path.isfile(args.yaml_file):
        logger.critical("{}: no such file".format(args.yaml_file))
        sys.exit(1)

    if not os.path.isfile(args.gene_list):
        logger.critical("{}: no such file".format(args.gene_list))
        sys.exit(1)

    if not os.path.isdir(args.output):
        os.mkdir(args.output)


def parse_args():
    """Parses the arguments and options."""
    parser = argparse.ArgumentParser(
        description="Automatically execute MetaSKAT on multiple cohorts "
                    "(using rpy2).",
    )

    # The version of the script
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s version " + __version__)

    # The input options
    group = parser.add_argument_group("INPUT OPTIONS")
    group.add_argument("-c", "--conf", required=True, metavar="YAML",
                       dest="yaml_file", help="The YAML file containing the "
                                              "information about the cohorts.")
    group.add_argument("-g", "--gene-list", required=True, metavar="FILE",
                       help="The gene list (with markers) required by "
                            "MetaSKAT")

    # The output options
    group = parser.add_argument_group("OUTPUT OPTIONS")
    group.add_argument("-o", "--output", default="metaskat_output",
                       metavar="DIR",
                       help="The output directory [%(default)s]")

    return parser.parse_args()


if __name__ == "__main__":
    main()
