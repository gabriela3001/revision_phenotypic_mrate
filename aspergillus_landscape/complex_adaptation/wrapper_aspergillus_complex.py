import subprocess

script_path = 'aspergillus_complex_adaptation.py'
cmd = 'bsub -q new-short -R "rusage[mem=1000]" -u /dev/null python ' + script_path

for param in range(10000):
  paramstr = ' '.join(map(str,[param]))
  _cmd = ' '.join([cmd, paramstr])
  subprocess.run(_cmd, shell=True)