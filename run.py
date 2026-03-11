import polars as pl
import numpy as np
import argparse



def parse_args():
    parser = argparse.ArgumentParser(
        prog='ctnna2_mutation_rate',
        description='',
    )

    parser.add_argument('--exons', required=True)
    parser.add_argument('--dnms', required=True)
    parser.add_argument('--roulette', required=True)
    parser.add_argument('--seed', default=0)
    parser.add_argument('--out', required=True)

    return parser.parse_args()


def multinomial_resampling(N_exons: int, N_observations: int, normalized_mutation_rates, N_iterations: int, seed: int = 0):
    rng = np.random.default_rng(seed = seed)
    events = rng.choice(a=N_exons, size=(N_iterations, N_observations), p=normalized_mutation_rates)
    return np.eye(N_exons)[events].sum(axis=1).astype(int) # Convert from sparse to dense matrix (counts per exon)


def main():
    args = parse_args()

    df_exons = pl.read_csv(args.exons, separator="\t") \
        .with_columns((pl.col("chromEnd") - pl.col("chromStart")).alias("exon_length"))

    df_dnms = pl.read_csv(args.dnms, separator="\t")
    df_roulette = pl.read_csv(args.roulette, separator="\t")

    df_dnms_per_exon = df_exons \
        .join_where(df_dnms, ((pl.col("chromStart") <= pl.col("POS")) & (pl.col("chromEnd") > pl.col("POS"))) ) \
        .group_by("name") \
        .len() \
        .rename({"len": "observed"})

    # Missense mutation rates per exon
    df_roulette_per_exon = df_exons \
        .join_where(df_roulette, ((pl.col("chromStart") <= pl.col("POS")) & (pl.col("chromEnd") > pl.col("POS"))) ) \
        .select(["name", "MR"]) \
        .group_by("name") \
        .sum()

    # Observed missense mutations per exon
    df_per_exon = df_dnms_per_exon \
        .join(df_roulette_per_exon, how="full", on="name", coalesce=True) \
        .with_columns(pl.col("observed").fill_null(0), pl.col("MR").fill_null(0)) \
        .sort(by="name")

    # Extract values from dataframe
    N_exons = df_per_exon.height
    observed = df_per_exon["observed"].to_numpy()
    N_observed = observed.sum()
    mutation_rates = df_per_exon["MR"].to_numpy()
    mutation_rates_norm = mutation_rates / mutation_rates.sum()

    # Multinomial resampling
    N_samplings = 100000
    counts = multinomial_resampling(N_exons, N_observed, mutation_rates_norm, N_samplings, seed=args.seed)

    # Calculate p values
    observed_deviation = np.abs(observed - mutation_rates_norm)
    expected_deviations = np.abs(counts - mutation_rates_norm)

    p_values = (1 + (counts > observed).sum(axis=0)) / (N_samplings + 1)

    # Add p values back to dataframe and adjust with bonferroni
    df = df_per_exon \
        .with_columns(pl.Series(p_values).alias("p_value")) \
        .with_columns((df_per_exon.height * pl.col("p_value")).alias("adj_p_value")) \
        .with_columns(pl.when(pl.col("adj_p_value") > 1).then(pl.lit(1)).otherwise(pl.col("adj_p_value")).alias("adj_p_value"))

    df.write_csv(args.out, separator="\t")


if __name__ == "__main__":
    main()
