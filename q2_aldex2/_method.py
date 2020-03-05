import os
import qiime2
import pandas as pd
import tempfile
import subprocess

from q2_aldex2._visualizer import _effect_statistic_functions


def run_commands(cmds, verbose=True):
    if verbose:
        print("Running external command line application(s). This may print "
              "messages to stdout and/or stderr.")
        print("The command(s) being run are below. These commands cannot "
              "be manually re-run as they will depend on temporary files that "
              "no longer exist.")
    for cmd in cmds:
        if verbose:
            print("\nCommand:", end=' ')
            print(" ".join(cmd), end='\n\n')
        subprocess.run(cmd, check=True)


def aldex2(table: pd.DataFrame,
           metadata: qiime2.CategoricalMetadataColumn,
           mc_samples: int = 128,
           test: str = 't',
           denom: str = 'all') -> pd.DataFrame:

    # create series from the metadata column
    meta = metadata.to_series()

    # The condition is just the only column in the passed metadata column
    condition = metadata.name

    # filter the metadata so only the samples present in the table are used
    # this also reorders it for the correct condition selection
    # it has to be re ordered for aldex to correctly input the conditions
    meta = meta.loc[list(table.index)]

    # force reorder based on the data to ensure conds are selected correctly

    with tempfile.TemporaryDirectory() as temp_dir_name:
        biom_fp = os.path.join(temp_dir_name, 'input.tsv.biom')
        map_fp = os.path.join(temp_dir_name, 'input.map.txt')
        summary_fp = os.path.join(temp_dir_name, 'output.summary.txt')

        # Need to manually specify header=True for Series (i.e. "meta"). It's
        # already the default for DataFrames (i.e. "table"), but we manually
        # specify it here anyway to alleviate any potential confusion.
        table.to_csv(biom_fp, sep='\t', header=True)
        meta.to_csv(map_fp, sep='\t', header=True)

        cmd = ['run_aldex2.R', biom_fp, map_fp, condition, mc_samples,
               test, denom, summary_fp]
        cmd = list(map(str, cmd))

        try:
            run_commands([cmd])
        except subprocess.CalledProcessError as e:
            raise Exception("An error was encountered while running ALDEx2"
                            " in R (return code %d), please inspect stdout"
                            " and stderr to learn more." % e.returncode)

        summary = pd.read_csv(summary_fp, index_col=0)
        #differentials = summary[['effect']]
	# hack to fix column name for features because aldex removes
	#it in R because of row.names = 1

        summary.index.name = "featureid"
        summary.rename(index=str, inplace=True)
        return summary

def extract_differences(table: pd.DataFrame, sig_threshold: float = 0.1, effect_threshold: float = 1, difference_threshold: float = 1, test: str = 'welch') -> pd.DataFrame:

    # checks to make sure there is no error
    # ensure max or min, depending on case

    effect_statistic_function = _effect_statistic_functions[test]

    if sig_threshold < table[effect_statistic_function].min():
        raise ValueError("You have selected a significance threshold that "
        "is lower than minimum Q score (-p--sig-threshold). Select a "
        "higher threshold.")

    # absolute values needed for effect or difference to see change in either
    # condition
    if effect_threshold > abs(table['effect']).max():
        raise ValueError("You have selected an effect threshold that exceeds "
        "maximum effect size (-p--effect-threshold). Choose a lower "
        "threshold, or be aware that there there will be no features "
        "in the output.")

    if difference_threshold > abs(table['diff.btw']).max():
        raise ValueError("You have selected a difference threshold that "
        "exceeds maximum difference (-p--difference-threshold). Choose a "
        "lower threshold, or be aware that there will be no features in "
        "the output.")

    # subset the table if it psases all the threshold
    differentials_sig = table[(table[effect_statistic_function]
    <= sig_threshold) & (abs(table['effect']) > effect_threshold) &
    (abs(table['diff.btw']) > difference_threshold)]

    return differentials_sig
