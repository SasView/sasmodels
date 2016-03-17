import shlex, subprocess, sys
import platform
args = shlex.split("sesansfit.sh blah --edit")
args[1] = sys.argv[1]
if platform.system() != "Windows":
    args = ["sh"] + args
    shellFlag = False
else:
    shellFlag = True
subprocess.Popen(args, shell=shellFlag)