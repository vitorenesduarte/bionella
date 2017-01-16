import os

site_path = os.path.join(os.getcwd(), "site")
os.chdir(site_path)
files = [f for f in os.listdir() if os.path.isfile(f)]

for f in files:
    if f.endswith(".ipynb"):
        # para os ficheiros ipython notebook

        if os.name == "nt":
            cmd = "ipython"
        elif os.name == "posix":
            cmd = "jupyter"

        # converter para html
        os.system(cmd + " nbconvert --to html --template full " + f)
