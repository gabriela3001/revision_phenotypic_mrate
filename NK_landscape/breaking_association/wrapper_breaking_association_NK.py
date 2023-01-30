import subprocess
import pickle

script_path = 'script_breaking_association_NK.py'
cmd = 'bsub -q new-short -R "rusage[mem=200]" -u /dev/null python ' + script_path

for NKind in [1,3,5]:
  for param in range(10000):
#for param, NKind in missing:
      paramstr = ' '.join(map(str,[param, NKind]))
      _cmd = ' '.join([cmd, paramstr])
      subprocess.run(_cmd, shell=True)