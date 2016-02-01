#!/usr/bin/env python3

from __future__ import print_function

import os
import sys
import shutil
import logging
import argparse
import subprocess
from tempfile import mkdtemp

import yaml
import numpy as np
import pandas as pd

from rpy2 import rinterface
import rpy2.robjects as robjects
from rpy2.robjects import Formula as rformula
import rpy2.robjects.numpy2ri as numpy2ri
from rpy2.robjects.packages import importr as rimport


# Activating numpy array to rpy2 objects automatic conversion
numpy2ri.activate()

# Removing R output
rinterface.set_writeconsole(None)

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

    # Creating the output directory, if required
    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    # Creating a temporary directory
    tmp_dir = mkdtemp(prefix="metaskat_tmp_", dir=args.output)

    # Reading the YAML configuration
    cohort_information = read_cohort_configuration(args.yaml_file)

    # Running SKAT on individual cohorts
    execute_skat(cohort_information, args.gene_list, o_prefix=args.output,
                 tmp_dir=tmp_dir)

    # Deleting temporary directory
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)


def execute_skat(cohorts, genes, o_prefix, tmp_dir):
    """Generating the SKAT meta files.

    Args:
        cohorts (dict): information about individual cohorts.
        genes (str): name of the file containing the SNP set (by gene).
        o_prefix (str): the output prefix.
        tmp_dir (str): the name of the temporary directory.

    """
    for cohort, cohort_info in cohorts.items():
        # Getting the plink file prefix (str) and the phenotype (DataFrame)
        plink_prefix, pheno = get_analysis_data(
            plink_prefix=cohort_info["prefix"],
            pheno=cohort_info["phenotype"],
            covariates=cohort_info["covariates"],
            pheno_fn=cohort_info["phenotype_file"],
            fid=cohort_info.get("family_id", "FID"),
            iid=cohort_info.get("individual_id", "IID"),
            tmp_dir=os.path.join(tmp_dir, cohort),
        )

        # Getting the required information
        covariates = cohort_info["covariates"]
        pheno_name = cohort_info["phenotype"]

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
        mssd = os.path.join(o_prefix, cohort + ".MSSD")
        minfo = os.path.join(o_prefix, cohort + ".MInfo")

        # Generating the meta files
        try:
            metaskat.Generate_Meta_Files(model, plink_prefix + ".bed",
                                         plink_prefix + ".bim", genes, mssd,
                                         minfo, pheno.shape[0],
                                         File_Permu=robjects.r("NULL"))
        except Exception as e:
            logger.critical("problem with SKAT\n{}".format(e))
            sys.exit(1)


def get_analysis_data(plink_prefix, pheno, covariates, fid, iid, pheno_fn,
                      tmp_dir):
    """Extract phenotypes and genomic information about the sample.

    Args:
        plink_prefix (str): the prefix of the Plink files.
        pheno (str): the name of the phenotype (in the phenotype file).
        covariates (list): the list of covariate names (in the pheno file).
        fid (str): the name of the family ID column (in the pheno file).
        iid (str): the name of the individual ID column (in the pheno file).
        pheno_fn (str): the name of the phenotype file.
        tmp_dir (str): the name of the temporary file.

    Returns:
        tuple: first element is the Plink prefix for the genotype data, the
               second element is the DataFrame containing the phenotype and
               covariate data.

    """
    # Reading the FAM file
    fam = pd.read_csv(plink_prefix + ".fam", sep=" ",
                      names=["fid", "iid", "father", "mother", "gender",
                             "status"])
    fam = fam.set_index(["fid", "iid"], verify_integrity=True)

    # Reading the phenotype and covariates
    phenotypes = pd.read_csv(pheno_fn, sep="\t")

    # Checking that we have all the required column in the phenotype file
    req_col = [fid, iid, pheno] + covariates
    for col in req_col:
        if col not in phenotypes.columns:
            logger.critical("{}: missing column '{}'".format(pheno_fn, col))
            sys.exit(1)

    # Keeping only the required column (and setting the index)
    phenotypes = phenotypes[req_col].set_index(
        [fid, iid],
        verify_integrity=True,
    ).dropna().sortlevel()

    # Finding the intersection between phenotype and sample file
    same_samples = fam.index.intersection(phenotypes.index)
    nb_same = len(same_samples.tolist())

    if nb_same == 0:
        logger.critical("{} vs {}: no sample in common".format(
            plink_prefix + ".fam",
            pheno_fn,
        ))
        sys.exit(1)

    # Extracting the samples from the phenotype file
    phenotypes = phenotypes.loc[same_samples, ]

    if (nb_same == fam.shape[0]) and (nb_same == phenotypes.shape[0]):
        # We have the same samples, hence, we only return the ordered phenotype
        # data (same order as FAM file)
        phenotypes = phenotypes.loc[fam.index, ]
        assert fam.index.tolist() == phenotypes.index.tolist()
        return plink_prefix, phenotypes

    # We need extraction, so we create the directory, if required
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)

    # Extracting required samples using plink
    new_prefix = extract_plink(plink_prefix, os.path.join(tmp_dir, "subset"),
                               same_samples.tolist())

    # Rerunning the function
    return get_analysis_data(new_prefix, pheno, covariates, fid, iid, pheno_fn,
                             tmp_dir)


def extract_plink(i_prefix, o_prefix, samples):
    """Extract a list of samples from a Plink binary file.

    Args:
        i_prefix (str): the prefix of the input Plink file.
        o_prefix (str): the prefix of the output Plink file.
        samples (list): the list of samples to extract

    """
    # Creating the list of samples to extract
    with open(o_prefix + ".to_extract", "w") as f:
        for fid, iid in samples:
            print(fid, iid, file=f)

    # Executing Plink (for extraction)
    plink_command = [
        "plink",
        "--noweb",
        "--bfile", i_prefix,
        "--keep", o_prefix + ".to_extract",
        "--make-bed",
        "--out", o_prefix,
    ]

    # Executing the command
    execute_command(plink_command)

    # Returning the prefix
    return o_prefix


def execute_command(command):   # pragma: no cover
    """Execute a command.

    Args:
        command (list): the command to execute.

    """
    try:
        proc = subprocess.Popen(command, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    except OSError:
        logger.critical("plink: missing binary")
        sys.exit(0)

    try:
        stdout, stderr = proc.communicate()

    except:
        logger.critical("problem with plink extraction")
        sys.exit(0)

    if proc.returncode != 0:
        logger.critical("there was a problem with plink\n" + stderr.decode())
        sys.exit(1)


def execute_meta_analysis(cohorts, genes, out_dir):
    # Combining the results
    mssd_files = np.array(list_mssd)
    minfo_files = np.array(list_minfo)
    meta_cohort_info = metaskat.Open_MSSD_File_2Read(mssd_files, minfo_files)

    # Computing the stats
    meta_hom = metaskat.MetaSKAT_MSSD_ALL(**{
        "Cohort.Info": meta_cohort_info,
        "combined.weight": True,
        "weights.beta": np.array([1, 25]),
        "method": "davies",
        "r.corr": 0,
        "is.separate": False,
        "Group_Idx": robjects.r("NULL"),
        "MAF.cutoff": 5,
    })

    meta_hom.to_csvfile(os.path.join(out_dir, "metaSKAT.homo.txt"),
                        quote=False, sep="\t", row_names=False, col_names=True)

    meta_het = metaskat.MetaSKAT_MSSD_ALL(**{
        "Cohort.Info": meta_cohort_info,
        "combined.weight": True,
        "weights.beta": np.array([1, 25]),
        "method": "davies",
        "r.corr": 0,
        "is.separate": True,
        "Group_Idx": robjects.r("NULL"),
        "MAF.cutoff": 5,
    })

    meta_het.to_csvfile(os.path.join(out_dir, "metaSKAT.hetero.txt"),
                        quote=False, sep="\t", row_names=False, col_names=True)

    # Printing final results per cohort
    each_cohort_info = meta_cohort_info[
        meta_cohort_info.names.index("EachInfo")
    ]
    for i, cohort in enumerate(cohort_names):
        final_info = each_cohort_info[i]
        final_info = final_info[final_info.names.index("Info")]
        final_info.to_csvfile(os.path.join(out_dir, cohort + ".MetaInfo.txt"),
                              sep="\t", quote=False, row_names=False)


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
