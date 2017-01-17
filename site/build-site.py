import os

[
    os.system("jupyter nbconvert --to html --template full " + f)
    for f in os.listdir() 
    if os.path.isfile(f) and f.endswith(".ipynb")
]
