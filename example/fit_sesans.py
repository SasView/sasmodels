import shlex, subprocess, sys
args = shlex.split("sesansfit.sh blah --edit")
args[1] = sys.argv[1]
subprocess.Popen(args, shell="True")