import shlex, subprocess
args = shlex.split("sesansfit.sh sesans_parameters_sphere.py --edit")
subprocess.Popen(args, shell="True")