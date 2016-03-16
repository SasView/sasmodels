import shlex, subprocess
args = shlex.split("sesansfit.sh sesans_parameters_css-hs.py --edit")
subprocess.Popen(args, shell="True")