############################################
# imports
############################################

import os
import re
import pandas as pd
import numpy as np

from pathlib import Path
from joblib import Parallel, delayed
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

from . import SpecificityFilterBase
from ..utils import get_sequence_from_annotation

############################################
# Oligo Blast Filter Classes
############################################


class Blastn(SpecificityFilterBase):
    """This class filters oligos based on the blast alignment tool. All the oligos which have a match with a percentage identity higher
    than the one given in input are filtered out.Moreover, an additional filter with AI models
     can be used to finely select which off-target match is relevant.

    :param dir_specificity: directory where alignement temporary files can be written
    :type dir_specificity: str
    :param word_size: word size for the blastn seed (exact match to target)
    :type word_size: int
    :param percent_identity: maximum similarity between oligos and target sequences, ranging from 0 to 100% (no missmatch)
    :type percent_identity: int
    :param strand: strand of the query sequence to search
    :type strand: str
    :param coverage: minimum coverage between oligos and target sequence, ranging from 0 to 100% (full coverage)
    :type coverage: float
    :param min_alignment_length: minimum alignment length between oligos and target sequence
    :type min_alignment_length: int
    :param ai_filter: machine learning model used to filter the oligos {None, 'hybridization_probability'}, defaults to None
    :type ai_filter: str
    :param ai_filter_thershold: threshold above (>=) which the oligos are filtered, defaults to 0.
    :type ai_filter_thershold: float, optional
    :param ai_filter_path: path to the machine learning model used to filter the oligos, if None the pretrained model provided will be used, defaults to None
    :type ai_filter_path: str, optional
    """

    def __init__(
        self,
        dir_specificity: str,
        word_size: int,
        percent_identity: float,
        strand: str,
        coverage: float = None,
        min_alignment_length: int = None,
        ai_filter: str = None,
        ai_filter_threshold: float = 0.5,
        ai_filter_path : str = None,
    ):
        """Constructor."""
        super().__init__(dir_specificity)

        self.word_size = word_size
        self.percent_identity = percent_identity
        self.strand = strand
        self.ai_filter = ai_filter
        self.ai_filter_threshold = ai_filter_threshold
        self.ai_filter_path = ai_filter_path

        if coverage is None and min_alignment_length is not None:
            self.min_alignment_length = min_alignment_length
        elif min_alignment_length is None and coverage is not None:
            self.coverage = coverage

        else:
            raise Exception(
                "Please provide either coverage or a minimum alignment length"
            )

        self.dir_blast = os.path.join(self.dir_specificity, "blast")
        Path(self.dir_blast).mkdir(parents=True, exist_ok=True)

        self.dir_fasta = os.path.join(self.dir_specificity, "fasta")
        Path(self.dir_fasta).mkdir(parents=True, exist_ok=True)
        self.filtered = pd.DataFrame(columns = ["Vanilla", "AI_filter", "Difference"])
        self.ai_time = pd.Series()

    def apply(self, database: dict, file_reference: str, n_jobs: int):
        """Apply the blastn filter in parallel on the given ``database``. Each jobs filters a single region, and  at the same time are generated at most ``n_job`` jobs.
        The filtered database is returned.

        :param database: database containing the oligos and their features
        :type database: dict
        :param file_reference: path to the file that will be used as reference for the alignement
        :type file_reference: str
        :param n_jobs: number of simultaneous parallel computations
        :type n_jobs: int
        :return: oligo info of user-specified regions
        :rtype: dict
        """

        # create blast database
        database_exists = False
        database_name = os.path.basename(file_reference)
        # Check if blast database exists
        for file in os.listdir(self.dir_blast):
            if re.search(f"^{database_name}.*", file):
                database_exists = True
                break

        if not database_exists:
            cmd = NcbimakeblastdbCommandline(
                input_file=file_reference,
                dbtype="nucl",
                out=os.path.join(self.dir_blast, database_name),
            )
            out, err = cmd()

        # run the balst search
        regions = list(database.keys())
        filtered_database_regions = Parallel(n_jobs=n_jobs)(
            delayed(self._run_blast)(database[region], region, database_name, file_reference)
            for region in regions
        )

        # reconstruct the oligos db and return it
        for region, filtered_database_region in zip(regions, filtered_database_regions):
            database[region] = filtered_database_region

        return database

    def _run_blast(self, database_region, region, database_name, file_reference):
        """Run BlastN alignment tool to find regions of local similarity between sequences, where sequences are oligos and background sequences (e.g. transcript, genome, etc.).
        BlastN identifies the transcript regions where oligos match with a certain coverage and similarity.

        :param database_region: database containing the oligos form one region
        :type database_region: dict
        :param region: id of the region processed
        :type region: str
        :param database_name: path to the blatn database
        :type database_name: str
        """
        # run the blast search and write the results
        file_oligo_fasta_gene = self._create_fasta_file(
            database_region, self.dir_fasta, region
        )
        file_blast_gene = os.path.join(self.dir_blast, f"blast_{region}.txt")
        cmd = NcbiblastnCommandline(
            query=file_oligo_fasta_gene,
            db=os.path.join(self.dir_blast, database_name),
            outfmt="6 qseqid sseqid length qlen qstart qend sstart send qseq sseq sstrand",
            out=file_blast_gene,
            strand=self.strand,
            word_size=self.word_size,
            perc_identity=self.percent_identity,
            num_threads=1,  # ????
        )
        out, err = cmd()

        # read the reuslts of the blast seatch
        blast_results = self._read_blast_output(file_blast_gene)
        # filter the DB based on the blast results
        matching_oligos = self._find_matching_oligos(blast_results)
        if self.ai_filter is not None:
            matching_oligos = self._ai_filter_matching_oligos(matching_oligos, file_reference, database_region, region)
        filtered_database_region = self._filter_matching_oligos(
            database_region, matching_oligos
        )
        # remove temporary files
        os.remove(file_blast_gene)
        os.remove(file_oligo_fasta_gene)
        return filtered_database_region

    def _read_blast_output(self, file_blast_gene):
        """Load the output of the BlastN alignment search into a DataFrame and process the results."""

        blast_results = pd.read_csv(
            file_blast_gene,
            header=None,
            sep="\t",
            low_memory=False,
            names=[
                "query",
                "target",
                "alignment_length",
                "query_length",
                "query_start",
                "query_end",
                "target_start",
                "target_end",
                "query_sequence",
                "target_sequence",
                "target_strand"
            ],
            engine="c",
            dtype={
                "query": str,
                "target": str,
                "alignment_length": int,
                "query_length": int,
                "query_start": int,
                "query_end": int,
                "target_start": int,
                "target_end": int,
                "query_sequence": str,
                "target_sequence": str,
                "target_strand": str,
            },
        )
        # return the real matches, that is the ones not belonging to the same region of the query oligo
        blast_results["query_gene_id"] = blast_results["query"].str.split("_").str[0]
        blast_results["target_gene_id"] = blast_results["target"].str.split("::").str[0]
        return blast_results

    def _find_matching_oligos(self, blast_results):
        """Use the results of the BlastN alignement search to remove oligos with high similarity,
        oligo coverage and ligation site coverage based on user defined thresholds.

        :param blast_results: DataFrame with processed blast alignment search results.
        :type blast_results: pandas.DataFrame
        """

        blast_matches = blast_results[
            blast_results["query_gene_id"] != blast_results["target_gene_id"]
        ]
        if self.coverage is not None:
            values = blast_matches["query_length"] * self.coverage / 100
            blast_matches.insert(
                len(blast_matches.columns), "min_alignment_length", values
            )
        else:
            blast_matches["min_alignment_length"] = self.min_alignment_length

        blast_matches_filtered = blast_matches.loc[
            blast_matches.alignment_length > blast_matches.min_alignment_length
        ]

        return blast_matches_filtered
    
    def _ai_filter_matching_oligos(self, matching_oligos: pd.DataFrame, file_reference: str, database_region: dict, region: str) -> pd.DataFrame:
        
        from odt_ai_filters.api import APIHybridizationProbability #rewrite the import with the new

        # check if there are any oligos to filter
        if len(matching_oligos) == 0:
            return matching_oligos
        
        # filter the oligos based on the ai filter outcomes
        targets = self._get_targets_fasta(matching_oligos, file_reference, region)

        # retrive original sequences of query and target adn encode the insertion and deletions
        queries = [database_region[query_id]["sequence"] for query_id in matching_oligos["query"]]
        assert len(queries) == len(targets), "The targets haven't been correctly retrieved."
        gapped_queries, gapped_targets = self._add_alignement_gaps(matching_oligos=matching_oligos, queries=queries, targets=targets)

        # create the dataset
        dataset = self._create_ai_filter_dataset(queries=queries, gapped_queries=gapped_queries, targets=targets, gapped_targets=gapped_targets)
        
        # define the AI filter
        if self.ai_filter == "hybridization_probability":
            model = APIHybridizationProbability(ai_filter_path = self.ai_filter_path)
        else:
            raise ValueError(f"The AI filter {self.ai_filter} is not supported.")

        # filter the database, keep only the oligos above the threshold
        predictions = model.predict(data = dataset)
        matching_oligos.reset_index(drop=True, inplace=True)
        ids_vanilla = len(matching_oligos["query"].unique())
        below_threshold = np.where(predictions < self.ai_filter_threshold)[0]
        matching_oligos.drop(index=below_threshold, inplace=True)
        ids_ai_filter = len(matching_oligos["query"].unique())
        self.filtered.loc[region] = [ids_vanilla, ids_ai_filter, ids_vanilla - ids_ai_filter]
        return matching_oligos


    def _get_targets_fasta(self, matching_oligos: pd.DataFrame, file_reference: str, region: str):
        """Create a FASTA file containing the sequences of the targets of the oligos that match the blast search.

        :param matching_oligos: DataFrame with the oligos that match the blast search.
        :type matching_oligos: pd.DataFrame
        :param file_reference: _description_
        :type file_reference: str
        :param region: _description_
        :type region: str
        :return: _description_
        :rtype: _type_
        """

        # set the positions to a 0-based index
        targets_plus = matching_oligos[matching_oligos["target_strand"] == "plus"]
        targets_plus["query_start"] = targets_plus["query_start"] - 1
        targets_plus["target_start"] = targets_plus["target_start"] - 1
        targets_minus = matching_oligos[matching_oligos["target_strand"] == "minus"]
        targets_minus["query_start"] = targets_minus["query_start"] - 1
        targets_minus["target_end"] = targets_minus["target_end"] - 1
        matching_oligos = pd.concat([targets_plus, targets_minus])
        matching_oligos.sort_index(inplace=True) # restore the intial indexes

        # create bed file for the plus strand
        bed_plus = pd.DataFrame(columns=["chr", "start", "end", "name", "score", "strand"], index=targets_plus.index)
        bed_plus["chr"] = targets_plus["target"]
        bed_plus["start"] = targets_plus["target_start"] - targets_plus["query_start"]
        bed_plus["end"] = targets_plus["target_end"] + (targets_plus["query_length"] - targets_plus["query_end"])
        bed_plus["name"] = targets_plus["query"]
        bed_plus["score"] = 0
        bed_plus["strand"] = '+'
        
        # create bed file for the minus strand
        bed_minus = pd.DataFrame(columns=["chr", "start", "end", "name", "score", "strand"], index=targets_minus.index)
        bed_minus["chr"] = targets_minus["target"]
        bed_minus["start"] = targets_minus["target_end"] - (targets_minus["query_length"] - targets_minus["query_end"])
        bed_minus["end"] = targets_minus["target_start"] + targets_minus["query_start"]
        bed_minus["name"] = targets_minus["query"]
        bed_minus["score"] = 0
        bed_minus["strand"] = '-'

        # concatenate the bed files
        bed = pd.concat([bed_plus, bed_minus])
        # adjust for possible overflows (e.g. new coordinates are not included in the gene boundaries)
        # additionally we store how muchpadding we have to do to have two seqeunces of the same length
        bed["overflow_start"] = bed["start"].apply(lambda x : -x if x < 0  else 0)
        bed["start"] = bed["start"].apply(lambda x : x if x >= 0  else 0)
        # what about having the leght of a region as a method of the reference database?
        records = SeqIO.index(file_reference, "fasta")
        bed["len_region"] = bed["chr"].apply(lambda x : len(records[x].seq))
        bed["overflow_end"] = bed[["end", "len_region"]].apply(lambda x : x["end"] - x["len_region"] if x["end"] > x["len_region"]  else 0, axis=1)
        bed["end"] = bed[["end", "len_region"]].apply(lambda x : x["end"] if x["end"] <= x["len_region"]  else x["len_region"], axis=1)
        bed.sort_index(inplace=True)
        file_bed = os.path.join(self.dir_specificity, f"targets_{region}.bed")
        bed.to_csv(file_bed, sep="\t", index=False, header=False, columns=["chr", "start", "end", "name", "score", "strand"])

        # generate the fasta file
        targets_fasta_file = os.path.join(self.dir_specificity, f"targets_{region}.fasta")
        get_sequence_from_annotation(file_bed, file_reference,targets_fasta_file, strand=True, nameOnly=True)
        
        targets = [off_target.seq for off_target in SeqIO.parse(targets_fasta_file, "fasta")]
        targets_padded = []
        for target, overflow_start, overflow_end in zip(targets, bed["overflow_start"], bed["overflow_end"]):
            if overflow_start != 0:
                # add padding in front
                target = '-' * overflow_start + target
            if overflow_end != 0:
                # add padding at the end
                target = target + '-' * overflow_end
            targets_padded.append(target)
        return targets_padded
    

    def _add_alignement_gaps(self, matching_oligos: pd.DataFrame, queries: list, targets: list):
        """Adjust the sequences of the oligos to remove the gaps introduced by the alignement search.

        :param matching_oligos: DataFrame with the oligos that have a match with the query oligo.
        :type matching_oligos: pandas.DataFrame
        :param database_region: Dictionary with the information of the oligo sequences.
        :type database_region: dict
        """

        def add_gaps(seq, gaps):
            for i, gap in enumerate(gaps):
                seq = seq[:gap+i] + '-' + seq[gap+i:] 
            return seq

        matching_oligos["query_gaps"] = matching_oligos["query_sequence"].apply(lambda x: np.where(np.array(list(x)) == '-')[0]) + matching_oligos["query_start"]
        matching_oligos["target_gaps"] = matching_oligos["target_sequence"].apply(lambda x: np.where(np.array(list(x)) == '-')[0]) + matching_oligos["query_start"]
        gapped_queries = [add_gaps(query, gaps) for query, gaps in zip(queries, matching_oligos["query_gaps"])]
        gapped_targets = [add_gaps(target, gaps) for target, gaps in zip(targets, matching_oligos["target_gaps"])]
        return gapped_queries, gapped_targets
    
    def _create_ai_filter_dataset(self, queries, gapped_queries, targets, gapped_targets):
        """Create a database with the information of the oligos that match the blast search.

        :param queries: List with the sequences of the query oligos.
        :type queries: list
        :param gapped_queries: List with the sequences of the query oligos with gaps.
        :type gapped_queries: list
        :param targets: List with the sequences of the target oligos.
        :type targets: list
        :param gapped_targets: List with the sequences of the target oligos with gaps.
        :type gapped_targets: list
        :return: dataset
        :rtype: pd.DataFrame
        """


        dataset = pd.DataFrame(columns=["query_sequence", "query_length", "query_GC_content", "off_target_sequence", "off_target_length", "off_target_GC_content", "number_mismatches"])
        dataset["query_sequence"] = gapped_queries
        dataset["query_length"] = [len(query) for query in queries]
        dataset["query_GC_content"] = [gc_fraction(query) for query in queries]
        dataset["off_target_sequence"] = gapped_targets
        dataset["off_target_length"] = [len(target) for target in targets]
        dataset["off_target_GC_content"] = [gc_fraction(target) for target in targets]
        dataset["number_mismatches"] = [sum(query != target for query, target in zip(query, target)) for query, target in zip(gapped_queries, gapped_targets)]
        return dataset
        