import subprocess

script_path = 'script_large_U.py'
cmd = 'bsub -q new-medium -R "rusage[mem=600]" -u /dev/null python ' + script_path

for p1 in range(12):
  for p in range(100):
      _cmd = ' '.join([cmd, str(p1), str(p)])
      subprocess.run(_cmd, shell=True)