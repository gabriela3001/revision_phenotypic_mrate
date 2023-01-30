import subprocess

script_path = 'script_stepwise_MSB.py'
cmd = 'bsub -q new-medium -R "rusage[mem=600]" -u /dev/null python ' + script_path

for p1 in range(36):
    for p in range(100):
        _cmd = ' '.join([cmd, str(p1), str(p)])
        subprocess.run(_cmd, shell=True)