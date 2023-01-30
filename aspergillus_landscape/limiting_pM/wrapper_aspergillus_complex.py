import subprocess

script_path = 'script_limiting_pM_aspergillus.py'
cmd = 'bsub -q new-short -R "rusage[mem=1024]" -u /dev/null python ' + script_path

for param in range(10000):
  paramstr = ' '.join(map(str,[param]))
  _cmd = ' '.join([cmd, paramstr])
  subprocess.run(_cmd, shell=True)