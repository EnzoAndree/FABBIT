import subprocess
import shlex
import pandas as pd
from io import StringIO

class DIAMOND:
    def __init__(self, diamond_path="diamond"):
        self.diamond_path = diamond_path
        self.sensitivity_mapping = {
            'faster': '--faster',
            'fast': '--fast',
            'mid-sensitive': '--mid-sensitive',
            'sensitive': '--sensitive',
            'more-sensitive': '--more-sensitive',
            'very-sensitive': '--very-sensitive',
            'ultra-sensitive': '--ultra-sensitive'
        }

    def _run_command(self, command):
        if isinstance(command, str):
            command = [self.diamond_path] + shlex.split(command)
        elif isinstance(command, list):
            command = [self.diamond_path] + command
        else:
            raise ValueError("Command must be a string or a list")
        process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if process.returncode != 0:
            raise Exception(f"Command '{' '.join(command)}' failed with return code {process.returncode}:\n{process.stderr}")
        return process.stdout

    def _format_command_argument(self, key, value):
        key = key.replace('_', '-')
        
        if key == 'sensitivity':
            return self._handle_sensitivity(value)
        elif isinstance(value, bool):
            return f" --{key}" if value else ""
        elif key == 'threads':
            return f" --{key} 1"  # DIAMOND runs with 1 core
        else:
            return f" --{key} {shlex.quote(str(value))}"

    def _handle_sensitivity(self, value):
        if value == 'default':
            return ""
        elif value in self.sensitivity_mapping:
            return f" {self.sensitivity_mapping[value]}"
        else:
            raise ValueError(f"Invalid sensitivity value: {value}")

    def makedb(self, in_file, db_file, **kwargs):
        command = ["makedb", "--in", in_file, "--db", db_file]
        
        for key, value in kwargs.items():
            command += shlex.split(self._format_command_argument(key, value))

        return self._run_command(command)

    def blastp(self, db, query, outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore", **kwargs):
        command = f"blastp --ignore-warnings --db {shlex.quote(db)} --query {shlex.quote(query)} --outfmt '{shlex.quote(outfmt)}'"
        
        for key, value in kwargs.items():
            command += self._format_command_argument(key, value)
        
        results_string = self._run_command(command)
        return pd.read_csv(StringIO(results_string), sep='\t', names=outfmt.split()[1:])

    def bidirectional_best_hit(self, fasta_file1, fasta_file2, **kwargs):
        results_1_to_2 = self.blastp(db=fasta_file2, query=fasta_file1, **kwargs)
        results_2_to_1 = self.blastp(db=fasta_file1, query=fasta_file2, **kwargs)

        best_hits_1_to_2 = self._get_best_hits(results_1_to_2)
        best_hits_2_to_1 = self._get_best_hits(results_2_to_1)

        bbh_df = pd.merge(
            best_hits_1_to_2,
            best_hits_2_to_1,
            left_on=['qseqid', 'sseqid'],
            right_on=['sseqid', 'qseqid'],
            suffixes=('', '_reverse')
        )

        relevant_reverse_cols = ['qseqid_reverse', 'sseqid_reverse', 'pident_reverse', 'evalue_reverse', 'bitscore_reverse']
        return bbh_df[results_1_to_2.columns.tolist() + relevant_reverse_cols]

    @staticmethod
    def _get_best_hits(results):
        idx = results.groupby('qseqid')['evalue'].idxmin()
        return results.loc[idx].reset_index(drop=True)