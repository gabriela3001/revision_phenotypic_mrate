import subprocess

script_path = 'script_complex_adaptation_NK.py'
cmd = 'bsub -q new-short -R "rusage[mem=10000]" -u /dev/null python ' + script_path

for NKind in [1]:
  for param in range(10000):
      paramstr = ' '.join(map(str,[param, NKind]))
      _cmd = ' '.join([cmd, paramstr])
      subprocess.run(_cmd, shell=True)