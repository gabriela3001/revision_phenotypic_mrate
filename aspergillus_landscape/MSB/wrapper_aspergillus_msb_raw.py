import subprocess

script_path = 'script_aspergillus_msb.py'
cmd = 'bsub -q new-short -R "rusage[mem=1000]" -u /dev/null python ' + script_path

for param in range(400):
#for param in {20, 104, 113, 186, 188, 190, 196, 199, 379, 389, 390, 395}:
    paramstr = ' '.join(map(str,[param]))
    _cmd = ' '.join([cmd, paramstr])
    subprocess.run(_cmd, shell=True)