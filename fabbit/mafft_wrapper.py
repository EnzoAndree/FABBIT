import subprocess
import shlex
import os

class MAFFT:
    def __init__(self, mafft_path="mafft"):
        self.mafft_path = mafft_path
        self.algorithm_presets = {
            'auto': '--auto',
            'L-INS-i': '--localpair --maxiterate 1000',
            'G-INS-i': '--globalpair --maxiterate 1000',
            'E-INS-i': '--ep 0 --genafpair --maxiterate 1000',
            'FFT-NS-i': '--retree 2 --maxiterate 1000',
            'FFT-NS-2': '--retree 2 --maxiterate 0',
            'FFT-NS-1': '--retree 1 --maxiterate 0',
            'NW-NS-i': '--retree 2 --maxiterate 2 --nofft',
            'NW-NS-2': '--retree 2 --maxiterate 0 --nofft',
            'NW-NS-PartTree-1': '--retree 1 --maxiterate 0 --nofft --parttree'
        }

    def _run_command(self, command, input_data=None):
        if isinstance(command, str):
            command = [self.mafft_path] + shlex.split(command)
        elif isinstance(command, list):
            command = [self.mafft_path] + command
        else:
            raise ValueError("Command must be a string or a list")
        
        process = subprocess.Popen(
            command,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        stdout, stderr = process.communicate(input=input_data)
        
        if process.returncode != 0:
            raise Exception(f"Command '{' '.join(command)}' failed with return code {process.returncode}:\n{stderr}")
        
        return stdout

    def _format_command_argument(self, key, value):
        key = key.replace('_', '-')
        
        if isinstance(value, bool):
            return f"--{key}" if value else ""
        else:
            return f"--{key} {shlex.quote(str(value))}"

    def align(self, input_sequences, output_file=None, algorithm='auto', **kwargs):
        command = self.algorithm_presets.get(algorithm, algorithm)
        command = shlex.split(command)
        
        for key, value in kwargs.items():
            command.append(self._format_command_argument(key, value))
        
        command.append("-")  # Tell MAFFT to read from stdin
        
        alignment_result = self._run_command(command, input_data=input_sequences)
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write(alignment_result)
            return f"Alignment saved to {output_file}"
        else:
            return alignment_result

    def align_from_file(self, input_file, output_file=None, algorithm='auto', **kwargs):
        with open(input_file, 'r') as f:
            input_sequences = f.read()
        
        return self.align(input_sequences, output_file, algorithm, **kwargs)

    def align_and_save(self, input_sequences, output_file, algorithm='auto', **kwargs):
        return self.align(input_sequences, output_file, algorithm, **kwargs)

    def profile_align(self, group1, group2, output_file=None, **kwargs):
        command = ["--maxiterate", "1000", "--seed", group1, "--seed", group2, "/dev/null"]
        
        for key, value in kwargs.items():
            command.append(self._format_command_argument(key, value))
        
        alignment_result = self._run_command(command)
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write(alignment_result)
            return f"Profile alignment saved to {output_file}"
        else:
            return alignment_result

# Example usage:
# mafft = MAFFT()
# result = mafft.align(sequences, algorithm="L-INS-i", op=1.53, ep=0.123)
# print(result)
#
# mafft.align_from_file("input.fasta", "output.aln", algorithm="G-INS-i", bl=45)
#
# mafft.align_and_save(sequences, "output.aln", algorithm="E-INS-i", lop=-2.00, lep=0.1)
#
# mafft.profile_align("group1.fasta", "group2.fasta", "profile_output.aln", op=1.53, ep=0.123)