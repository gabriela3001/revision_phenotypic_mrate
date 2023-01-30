import subprocess

script_path = 'script_NK_msb.py'
cmd = 'bsub -q new-short -R "rusage[mem=1000]" -u /dev/null python ' + script_path

for NKind in [1,3,5]:
  for param in range(400):
      paramstr = ' '.join(map(str,[param, NKind]))
      _cmd = ' '.join([cmd, paramstr])
      subprocess.run(_cmd, shell=True)