import os
import subprocess

cmd = ['sphinx-build', os.getcwd(), '_build']

subprocess.run(cmd, check=True)
