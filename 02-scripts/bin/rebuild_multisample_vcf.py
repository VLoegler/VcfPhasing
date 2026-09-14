#!/usr/bin/env python3

"""
Author: Victor Loegler
Date: 2026-09-11

Description
-----------
Rebuild a multisample VCF by transferring phased genotype information from
single-sample WhatsHap VCFs back into the original multisample VCF.

Input:
    - Original multisample VCF
    - One or more phased single-sample VCFs

Output:
    - Multisample VCF containing phased genotypes for samples that were phased

Notes
-----
All phased VCFs are assumed to derive from the same original multisample VCF
and therefore contain identical variant coordinates.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from cyvcf2 import VCF, Writer

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Rebuild multisample VCF from phased single-sample VCFs."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Original multisample VCF (.vcf.gz)"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output phased multisample VCF (.vcf.gz)"
    )

    parser.add_argument(
        "--phased-vcfs",
        required=True,
        nargs="+",
        help="Single-sample phased VCFs"
    )

    return parser.parse_args()


def load_phased_genotypes(
    phased_vcfs: list[str]
) -> dict[str, dict[tuple[str, int], list[int]]]:
    """
    Load phased genotype calls from Whatshap output VCFs.

    Returns
    -------
    dict
        sample -> {(chrom, pos): genotype}
    """

    phased_data: dict[str, dict[tuple[str, int], list[int]]] = {}

    for vcf_path in phased_vcfs:

        logger.info("Loading %s", vcf_path)

        vcf = VCF(vcf_path)

        if len(vcf.samples) != 1:
            raise ValueError(
                f"{vcf_path} contains {len(vcf.samples)} samples, expected 1"
            )

        sample = vcf.samples[0]

        phased_data[sample] = {}

        n_records = 0

        for record in vcf:

            key = (record.CHROM, record.POS)

            # cyvcf2 genotype structure:
            # [allele1, allele2, phased]
            phased_data[sample][key] = record.genotypes[0]

            n_records += 1

        logger.info(
            "Loaded %d variants for sample %s",
            n_records,
            sample
        )

    return phased_data


def rebuild_vcf(
    input_vcf: str,
    output_vcf: str,
    phased_data: dict[str, dict[tuple[str, int], list[int]]]
) -> None:
    """
    Inject phased GT fields into the original multisample VCF.
    """

    vcf = VCF(input_vcf)

    sample_to_index = {
        sample: idx
        for idx, sample in enumerate(vcf.samples)
    }

    missing_samples = (
        set(phased_data.keys()) - set(sample_to_index.keys())
    )

    if missing_samples:
        raise ValueError(
            f"Samples absent from original VCF: "
            f"{', '.join(sorted(missing_samples))}"
        )

    writer = Writer(output_vcf, vcf)

    updated_variants = 0

    for record in vcf:

        key = (record.CHROM, record.POS)

        genotypes = record.genotypes

        modified = False

        for sample, sample_variants in phased_data.items():

            gt = sample_variants.get(key)

            if gt is None:
                continue

            sample_idx = sample_to_index[sample]

            genotypes[sample_idx] = gt

            modified = True

        if modified:
            record.set_genotypes(genotypes)

            updated_variants += 1

        writer.write_record(record)

    writer.close()
    vcf.close()

    logger.info(
        "Updated phase information for %d variants",
        updated_variants
    )


def main() -> None:

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s"
    )

    args = parse_args()

    phased_data = load_phased_genotypes(
        args.phased_vcfs
    )

    rebuild_vcf(
        input_vcf=args.input,
        output_vcf=args.output,
        phased_data=phased_data
    )

    logger.info("Finished")


if __name__ == "__main__":
    main()