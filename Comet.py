#! usr/bin/python
import subprocess
import sys


sshProcess = subprocess.Popen(['ssh', 'seroot@comet.sdsc.edu'],stdin =subprocess.PIPE, stdout=subprocess.PIPE)
sshProcess.stdin.write("go")
